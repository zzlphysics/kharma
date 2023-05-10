/* 
 *  File: kharma.cpp
 *  
 *  BSD 3-Clause License
 *  
 *  Copyright (c) 2020, AFD Group at UIUC
 *  All rights reserved.
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1. Redistributions of source code must retain the above copyright notice, this
 *     list of conditions and the following disclaimer.
 *  
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "kharma.hpp"

#include <iostream>

#include <parthenon/parthenon.hpp>

#include "decs.hpp"

// Packages
#include "b_flux_ct.hpp"
#include "b_cd.hpp"
#include "b_cleanup.hpp"
#include "current.hpp"
#include "kharma_driver.hpp"
#include "electrons.hpp"
#include "implicit.hpp"
#include "floors.hpp"
#include "grmhd.hpp"
#include "reductions.hpp"
#include "emhd.hpp"
#include "wind.hpp"

#include "bondi.hpp"
#include "boundaries.hpp"
#include "resize_restart.hpp"
#include "resize_restart_kharma.hpp"

std::shared_ptr<KHARMAPackage> KHARMA::InitializeGlobals(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    // All truly global state.  Mostly mutable state in order to avoid scope creep
    auto pkg = std::make_shared<KHARMAPackage>("Globals");
    Params &params = pkg->AllParams();
    // Current time in the simulation.  For ramping things up, ramping things down,
    // or preventing bad outcomes at known times
    params.Add("time", 0.0, true);
    // Last step's dt (Parthenon SimTime tm.dt), which must be preserved to output jcon
    params.Add("dt_last", 0.0, true);
    // Whether we are computing initial outputs/timestep, or versions in the execution loop
    params.Add("in_loop", false, true);

    // Log levels, the other acceptable global
    // Made mutable in case we want to bump global log level on certain events
    // TODO allow a "go_verbose" file watch
    int verbose = pin->GetOrAddInteger("debug", "verbose", 0);
    params.Add("verbose", verbose, true);
    int flag_verbose = pin->GetOrAddInteger("debug", "flag_verbose", 0);
    params.Add("flag_verbose", flag_verbose, true);
    int extra_checks = pin->GetOrAddInteger("debug", "extra_checks", 0);
    params.Add("extra_checks", extra_checks, true);

    // Record the problem name, just in case we need to special-case for different problems.
    // Please favor packages & options before using this, and modify problem-specific code
    // to be more general as it matures.
    std::string problem_name = pin->GetString("parthenon/job", "problem_id");
    params.Add("problem", problem_name);

    // Update the times with callbacks
    pkg->MeshPreStepUserWorkInLoop = KHARMA::MeshPreStepUserWorkInLoop;
    pkg->MeshPostStepUserWorkInLoop = KHARMA::MeshPostStepUserWorkInLoop;

    return pkg;
}
void KHARMA::ResetGlobals(ParameterInput *pin, Mesh *pmesh)
{
    // The globals package was loaded & exists, retrieve it
    auto pkg = pmesh->packages.Get("Globals");
    Params &params = pkg->AllParams();
    // This needs to be reset to guarantee that EstimateTimestep doesn't try to
    // calculate a new dt from a blank 'ctop' variable,
    // just uses whatever the next step was going to be at reset
    params.Update("in_loop", false);

    // Everything else is a per-step variable, not per-run, so they're fine
    // to be restored by Parthenon
}

void KHARMA::MeshPreStepUserWorkInLoop(Mesh *pmesh, ParameterInput *pin, const SimTime &tm)
{
    Flag("KHARMA Pre-step");
    auto& globals = pmesh->packages.Get("Globals")->AllParams();
    if (!globals.Get<bool>("in_loop")) {
        globals.Update<bool>("in_loop", true);
    }
    globals.Update<double>("dt_last", tm.dt);
    globals.Update<double>("time", tm.time);
}

void KHARMA::MeshPostStepUserWorkInLoop(Mesh *pmesh, ParameterInput *pin, const SimTime &tm)
{
    Flag("KHARMA Post-step");
    // Knowing this works took a little digging into Parthenon's EvolutionDriver.
    // The order of operations after calling Step() is:
    // 1. Call PostStepUserWorkInLoop and PostStepDiagnostics (this function and following)
    // 2. Set the timestep tm.dt to the minimum from the EstimateTimestep calls
    // 3. Generate any outputs, e.g. jcon
    // Thus we preserve tm.dt (which has not yet been reset) as dt_last for Current::FillOutput
    auto& globals = pmesh->packages.Get("Globals")->AllParams();
    globals.Update<double>("dt_last", tm.dt);
    globals.Update<double>("time", tm.time);
}

void KHARMA::FixParameters(std::unique_ptr<ParameterInput>& pin)
{
    Flag("Fixing parameters");
    // Parthenon sets 2 ghost zones as a default.
    // We can't override that default while allowing a file-specified value.
    // Fine for now because we crash with 2. (Flux CT)
    // TODO add under different name?  Better precedence/origin code?
    pin->SetInteger("parthenon/mesh", "nghost", 4);
    Globals::nghost = pin->GetInteger("parthenon/mesh", "nghost");
    // Warn if using less than 4 ghost zones in any circumstances, it's still not tested well
    // if (Globals::nghost < 4) {
    //     std::cerr << "WARNING: Using less than 4 ghost zones is untested!" << std::endl;
    // }

    // If we're restarting (not via Parthenon), read the restart file to get most parameters
    std::string prob = pin->GetString("parthenon/job", "problem_id");
    if (prob == "resize_restart") {
        ReadIharmRestartHeader(pin->GetString("resize_restart", "fname"), pin);
    }
    if (prob == "resize_restart_kharma") {
        ReadKharmaRestartHeader(pin->GetString("resize_restart", "fname"), pin);
    }

    // Construct a CoordinateEmbedding object.  See coordinate_embedding.hpp for supported systems/tags
    CoordinateEmbedding tmp_coords(pin.get());
    // Record whether we're in spherical as we'll need that
    pin->SetBoolean("coordinates", "spherical", tmp_coords.is_spherical());

    // Do a bunch of autodetection/setting in spherical coordinates
    // Note frequent use of "GetOrAddX": this sets a default if not present but allows overriding
    if (tmp_coords.is_spherical()) {
        // Spherical systems can specify r_out and optionally r_in,
        // instead of xNmin/max.
        if (!pin->DoesParameterExist("parthenon/mesh", "x1min") ||
            !pin->DoesParameterExist("parthenon/mesh", "x1max")) {
            // Outer radius is always specified
            GReal Rout = pin->GetReal("coordinates", "r_out");
            GReal x1max = tmp_coords.r_to_native(Rout);
            pin->GetOrAddReal("parthenon/mesh", "x1max", x1max);

            if (mpark::holds_alternative<SphMinkowskiCoords>(tmp_coords.base)) {
                // In Minkowski coordinates, require Rin so the singularity is at user option
                GReal Rin = pin->GetReal("coordinates", "r_in");
                GReal x1min = tmp_coords.r_to_native(Rin);
                pin->GetOrAddReal("parthenon/mesh", "x1min", x1min);
            } else { // Any spherical BH metric: KS, BL, and derivatives
                // Set inner radius if not specified
                if (pin->DoesParameterExist("coordinates", "r_in")) {
                    GReal Rin = pin->GetReal("coordinates", "r_in");
                    GReal x1min = tmp_coords.r_to_native(Rin);
                    pin->GetOrAddReal("parthenon/mesh", "x1min", x1min);
                    if (Rin < 2.0){ // warn if there are fewer than 5 zones inside the event horizon
                        GReal dx = (x1max - x1min) / pin->GetInteger("parthenon/mesh", "nx1");
                        if (tmp_coords.X1_to_embed(x1min + 5*dx) > tmp_coords.get_horizon()) {
                            std::cerr << "WARNING: inner radius is near/in the EH, but does not allow 5 zones inside!" << std::endl;
                        }
                    }
                } else {
                    int nx1 = pin->GetInteger("parthenon/mesh", "nx1");
                    // Allow overriding Rhor for bondi_viscous problem
                    const GReal Rhor = pin->GetOrAddReal("coordinates", "Rhor", tmp_coords.get_horizon());
                    const GReal x1hor = tmp_coords.r_to_native(Rhor);

                    // Set Rin such that we have 5 zones completely inside the event horizon
                    // If xeh = log(Rhor), xin = log(Rin), and xout = log(Rout),
                    // then we want xeh = xin + 5.5 * (xout - xin) / N1TOT:
                    const GReal x1min = (nx1 * x1hor / 5.5 - x1max) / (-1. + nx1 / 5.5);
                    if (x1min < 0.0) {
                        throw std::invalid_argument("Not enough radial zones were specified to put 5 zones inside EH!");
                    }
                    pin->GetOrAddReal("parthenon/mesh", "x1min", x1min);
                    pin->GetOrAddReal("coordinates", "r_in", tmp_coords.X1_to_embed(Rhor));
                }
            }
        }

        // If the simulation domain extends inside the EH, we change some boundary options
        pin->SetBoolean("coordinates", "domain_intersects_eh", pin->GetReal("coordinates", "r_in") < tmp_coords.get_horizon());

        // Spherical systems will also want KHARMA's spherical boundary conditions.
        // Note boundaries are now exclusively set by KBoundaries package
        pin->GetOrAddString("boundaries", "inner_x1", "outflow");
        pin->GetOrAddString("boundaries", "outer_x1", "outflow");
        pin->GetOrAddString("boundaries", "inner_x2", "reflecting");
        pin->GetOrAddString("boundaries", "outer_x2", "reflecting");
        pin->GetOrAddString("boundaries", "inner_x3", "periodic");
        pin->GetOrAddString("boundaries", "outer_x3", "periodic");
    } else {
        // We can set reasonable default boundary conditions for Cartesian sims,
        // but not default domain bounds
        pin->GetOrAddString("boundaries", "inner_x1", "periodic");
        pin->GetOrAddString("boundaries", "outer_x1", "periodic");
        pin->GetOrAddString("boundaries", "inner_x2", "periodic");
        pin->GetOrAddString("boundaries", "outer_x2", "periodic");
        pin->GetOrAddString("boundaries", "inner_x3", "periodic");
        pin->GetOrAddString("boundaries", "outer_x3", "periodic");
    }

    // Default boundaries are to cover the domain of our native coordinate system
    // std::cout << "Coordinate transform has boundaries: "
    //             << tmp_coords.startx(1) << " "
    //             << tmp_coords.startx(2) << " "
    //             << tmp_coords.startx(3) << " to "
    //             << tmp_coords.stopx(1) << " "
    //             << tmp_coords.stopx(2) << " "
    //             << tmp_coords.stopx(3) << std::endl;
    // TODO(BSP) is this worth looping?  I say probably no.
    if (tmp_coords.startx(1) >= 0)
        pin->GetOrAddReal("parthenon/mesh", "x1min", tmp_coords.startx(1));
    if (tmp_coords.stopx(1) >= 0)
        pin->GetOrAddReal("parthenon/mesh", "x1max", tmp_coords.stopx(1));
    if (tmp_coords.startx(2) >= 0)
        pin->GetOrAddReal("parthenon/mesh", "x2min", tmp_coords.startx(2));
    if (tmp_coords.stopx(2) >= 0)
        pin->GetOrAddReal("parthenon/mesh", "x2max", tmp_coords.stopx(2));
    if (tmp_coords.startx(3) >= 0)
        pin->GetOrAddReal("parthenon/mesh", "x3min", tmp_coords.startx(3));
    if (tmp_coords.stopx(3) >= 0)
        pin->GetOrAddReal("parthenon/mesh", "x3max", tmp_coords.stopx(3));

    Flag("Fixed");
}

TaskStatus KHARMA::AddPackage(std::shared_ptr<Packages_t>& packages,
                              std::function<std::shared_ptr<KHARMAPackage>(ParameterInput*, std::shared_ptr<Packages_t>&)> package_init,
                              ParameterInput *pin)
{
    Flag("AddPackage");
    const auto& pkg = package_init(pin, packages);
    packages->Add(pkg);
    EndFlag("AddPackage "+pkg->label());
    return TaskStatus::complete;
}

Packages_t KHARMA::ProcessPackages(std::unique_ptr<ParameterInput> &pin)
{
    // See above.  Only run if 
    //if ()
    FixParameters(pin);

    Flag("ProcessPackages");

    // Allocate the packages list as a shared pointer, to be updated in various tasks
    auto packages = std::make_shared<Packages_t>();

    TaskCollection tc;
    auto& tr = tc.AddRegion(1);
    auto& tl = tr[0];
    TaskID t_none(0);
    // The globals package will never have dependencies
    auto t_globals = tl.AddTask(t_none, KHARMA::AddPackage, packages, KHARMA::InitializeGlobals, pin.get());
    // Driver package is the foundation
    auto t_driver = tl.AddTask(t_none, KHARMA::AddPackage, packages, KHARMADriver::Initialize, pin.get());
    // Floors package has no dependencies
    if (!pin->GetOrAddBoolean("floors", "disable_floors", false)) {
        auto t_floors = tl.AddTask(t_none, KHARMA::AddPackage, packages, Floors::Initialize, pin.get());
    }
    // GRMHD needs globals to mark packages
    auto t_grmhd = tl.AddTask(t_globals | t_driver, KHARMA::AddPackage, packages, GRMHD::Initialize, pin.get());
    // Inverter (TODO: split out fixups, then don't load this when GRMHD isn't loaded)
    auto t_inverter = tl.AddTask(t_grmhd, KHARMA::AddPackage, packages, Inverter::Initialize, pin.get());

    // B field solvers, to ensure divB ~= 0.
    // Bunch of logic here: basically we want to load <=1 solver with an encoded order of preference
    auto t_b_field = t_none;
    std::string b_field_solver = pin->GetOrAddString("b_field", "solver", "flux_ct");
    if (b_field_solver == "none" || b_field_solver == "b_cleanup") {
        // Don't add a B field
    } else if (b_field_solver == "constraint_damping" || b_field_solver == "b_cd") {
        // Constraint damping, probably only useful for non-GR MHD systems
        t_b_field = tl.AddTask(t_grmhd, KHARMA::AddPackage, packages, B_CD::Initialize, pin.get());
    } else {
        // Don't even error on bad values.  This is probably what you want
        t_b_field = tl.AddTask(t_grmhd, KHARMA::AddPackage, packages, B_FluxCT::Initialize, pin.get());
    }
    // Cleanup for the B field, using an elliptic solve for eliminating divB
    // Almost always loaded explicitly in addition to another transport, just for cleaning at simulation start
    // Enable b_cleanup package if we want it explicitly
    bool b_cleanup_package = pin->GetOrAddBoolean("b_cleanup", "on", (b_field_solver == "b_cleanup"));
    // OR if we need it for resizing a dump
    bool is_resize = pin->GetString("parthenon/job", "problem_id") == "resize_restart" &&
                     !pin->GetOrAddBoolean("resize_restart", "skip_b_cleanup", false);
    // OR if we ordered an initial cleanup pass for some other reason
    bool initial_cleanup = pin->GetOrAddBoolean("b_field", "initial_cleanup", false);
    bool use_b_cleanup = b_cleanup_package || is_resize || initial_cleanup;
    pin->SetBoolean("b_cleanup", "on", use_b_cleanup);
    auto t_b_cleanup = t_none;
    if (use_b_cleanup) {
        t_b_cleanup = tl.AddTask(t_grmhd, KHARMA::AddPackage, packages, B_Cleanup::Initialize, pin.get());
        if (t_b_field == t_none) t_b_field = t_b_cleanup;
    }

    // Enable calculating jcon iff it is in any list of outputs (and there's even B to calculate it)
    // Since it is never required to restart, this is the only time we'd write (hence, need) it
    // TODO use GetVector & == when available
    if (FieldIsOutput(pin.get(), "jcon") && t_b_field != t_none) {
        auto t_current = tl.AddTask(t_b_field, KHARMA::AddPackage, packages, Current::Initialize, pin.get());
    }
    // Electrons are usually boring but not impossible without a B field (TODO add a test?)
    if (pin->GetOrAddBoolean("electrons", "on", false)) {
        auto t_electrons = tl.AddTask(t_grmhd, KHARMA::AddPackage, packages, Electrons::Initialize, pin.get());
    }
    if (pin->GetOrAddBoolean("emhd", "on", false)) {
        auto t_electrons = tl.AddTask(t_grmhd, KHARMA::AddPackage, packages, EMHD::Initialize, pin.get());
    }
    if (pin->GetOrAddBoolean("wind", "on", false)) {
        auto t_electrons = tl.AddTask(t_grmhd, KHARMA::AddPackage, packages, Wind::Initialize, pin.get());
    }

    // Execute the whole collection (just in case we do something fancy?)
    while (!tr.Execute()); // TODO this will inf-loop on error

    // The boundaries package may need to know variable counts for allocating memory,
    // so we initialize it after the main dependency tree
    // TODO only init if at least one boundary is "user"
    KHARMA::AddPackage(packages, KBoundaries::Initialize, pin.get());

    // Load the implicit package *last*, if there are any variables which need implicit evolution
    // TODO print what we're doing here & do some sanity checks, if verbose
    int n_implicit = packages->Get("Driver")->Param<int>("n_implicit_vars");
    if (n_implicit > 0) {
        KHARMA::AddPackage(packages, Implicit::Initialize, pin.get());
        // Implicit evolution must use predictor-corrector i.e. "vl2" integrator
        pin->SetString("parthenon/time", "integrator", "vl2");
    }

    EndFlag("ProcessPackages"); // TODO print full package list way up here?
    return std::move(*packages);
}
