// Fishbone-Moncrief torus initialization functions
#pragma once

#include "decs.hpp"
#include "types.hpp"

/**
 * Initialize a wide variety of different fishbone-moncrief torii.
 *
 * @param rin is the torus innermost radius, in r_g
 * @param rmax is the radius of maximum density of the F-M torus in r_g
 */
TaskStatus InitializeKZTorus(std::shared_ptr<MeshBlockData<Real>>& rc, ParameterInput *pin);

/**
 * Torus solution for ln h, See Fishbone and Moncrief eqn. 3.6. 
 */
KOKKOS_INLINE_FUNCTION Real lnh_calc(const GReal a, const Real kzeta, const Real l, const GReal rin, const GReal r, const GReal th)
{
    // TODO this isn't faster than splitting into two evaluations of a sub-function,
    // and it doesn't matter anyway.  Make it clearer
    Real sth = m::sin(th);
    Real cth = m::cos(th);

    // Real r2 = r*r;
    // Real a2 = a*a;
    // // Metric 
    // Real DD = r2 - 2. * r + a2;
    // Real AA = m::pow(r2 + a2, 2) - DD * a2 * sth * sth;
    // Real SS = r2 + a2 * cth * cth;

    // Real thin = M_PI / 2.;
    // Real sthin = m::sin(thin);
    // Real cthin = m::cos(thin);

    // Real rin2 = m::pow(rin, 2);
    // Real DDin = rin2 - 2. * rin + a2;
    // Real AAin = m::pow(rin2 + a2, 2) - DDin * a2 * sthin * sthin;
    // Real SSin = rin2 + a2 * cthin * cthin;

    if (r >= rin) {
        Real sth = m::sin(th);
        Real cth = m::cos(th);
        Real sth2 = sth*sth;
        Real cth2 = 1.0 - sth2;
        Real s2th = 2.0*sth*cth;
        Real c2th = m::cos(2.0*th);
        Real cscth = 1.0/sth;
        Real cotth = cth/sth;
        Real res = (-m::sqrt(1.0 + (l*l*r*(a*a*r - 2.0*r*r + m::pow(r,3) - kzeta)*
                m::pow(a*a + 2.0*r*r + a*a*c2th,2)*m::pow(cscth,6))/
                m::pow(a*a*(2.0*r*r + kzeta) + a*a*r*(a*a + r*r)*m::pow(cotth,2) + 
                m::pow(r,3)*(a*a + r*r)*m::pow(cscth,2),2)) + 
            m::log(1.0 + m::sqrt(1.0 + (l*l*r*(a*a*r - 2.0*r*r + m::pow(r,3) - kzeta)*
                    m::pow(a*a + 2.0*r*r + a*a*c2th,2)*m::pow(cscth,6))/
                m::pow(a*a*(2.0*r*r + kzeta) + a*a*r*(a*a + r*r)*m::pow(cotth,2) + 
                    m::pow(r,3)*(a*a + r*r)*m::pow(cscth,2),2))) - 
            m::log(((a*a*r - 2.0*r*r + m::pow(r,3) - kzeta)*(a*a + 2.0*r*r + a*a*c2th))/
            (2.*(m::pow(r,3)*(a*a + r*r) + a*a*r*(a*a + r*r)*cth2 + 
                a*a*(2.0*r*r + kzeta)*sth2))) - 
            (2.0*a*l*(2.0*r*r + kzeta))/
            (m::pow(r,3)*(a*a + r*r) + a*a*r*(a*a + r*r)*cth2 + 
                a*a*(2.0*r*r + kzeta)*sth2))/2.;

        Real th_in = M_PI / 2.;
        Real sth_in = m::sin(th_in);
        Real cth_in = m::cos(th_in);
        Real sth2_in = sth_in*sth_in;
        Real cth2_in = 1.0 - sth2_in;
        Real s2th_in = 2.0*sth_in*cth_in;
        Real c2th_in = m::cos(2.0*th_in);
        Real cscth_in = 1.0/sth_in;
        Real cotth_in = cth_in/sth_in;
        Real res_in = (-m::sqrt(1.0 + (l*l*rin*(a*a*rin - 2.0*rin*rin + m::pow(rin,3) - kzeta)*
                m::pow(a*a + 2.0*rin*rin + a*a*c2th_in,2)*m::pow(cscth_in,6))/
                m::pow(a*a*(2.0*rin*rin + kzeta) + a*a*rin*(a*a + rin*rin)*m::pow(cotth_in,2) + 
                m::pow(rin,3)*(a*a + rin*rin)*m::pow(cscth_in,2),2)) + 
            m::log(1.0 + m::sqrt(1.0 + (l*l*rin*(a*a*rin - 2.0*rin*rin + m::pow(rin,3) - kzeta)*
                    m::pow(a*a + 2.0*rin*rin + a*a*c2th_in,2)*m::pow(cscth_in,6))/
                m::pow(a*a*(2.0*rin*rin + kzeta) + a*a*rin*(a*a + rin*rin)*m::pow(cotth_in,2) + 
                    m::pow(rin,3)*(a*a + rin*rin)*m::pow(cscth_in,2),2))) - 
            m::log(((a*a*rin - 2.0*rin*rin + m::pow(rin,3) - kzeta)*(a*a + 2.0*rin*rin + a*a*c2th_in))/
            (2.*(m::pow(rin,3)*(a*a + rin*rin) + a*a*rin*(a*a + rin*rin)*cth2_in + 
                a*a*(2.0*rin*rin + kzeta)*sth2_in))) - 
            (2.0*a*l*(2.0*rin*rin + kzeta))/
            (m::pow(rin,3)*(a*a + rin*rin) + a*a*rin*(a*a + rin*rin)*cth2_in + 
                a*a*(2.0*rin*rin + kzeta)*sth2_in))/2.;
        return res - res_in;   
        // return
        //     0.5 *
        //         m::log((1. +
        //                 m::sqrt(1. +
        //                     4. * (l * l * SS * SS) * DD / (AA * AA * sth * sth))) /
        //             (SS * DD / AA)) -
        //     0.5 * m::sqrt(1. +
        //                 4. * (l * l * SS * SS) * DD /
        //                     (AA * AA * sth * sth)) -
        //     2. * a * r * l / AA -
        //         (0.5 *
        //             m::log((1. +
        //                 m::sqrt(1. +
        //                     4. * (l * l * SSin * SSin) * DDin /
        //                         (AAin * AAin * sthin * sthin))) /
        //                 (SSin * DDin / AAin)) -
        //         0.5 * m::sqrt(1. +
        //                 4. * (l * l * SSin * SSin) * DDin / (AAin * AAin * sthin * sthin)) -
        //         2. * a * rin * l / AAin);
    } else {
        return 1.;
    }
}

/**
 * This function calculates specific the angular momentum of the
 * Fishbone-Moncrief solution in the midplane, as a function of radius.
 * (see Fishbone & Moncrief eqn. 3.8)
 * It improves on (3.8) by requiring no sign changes for
 * co-rotating (a > 0) vs counter-rotating (a < 0) disks.
 */
// KOKKOS_INLINE_FUNCTION Real lfish_calc(const GReal a, const GReal r)
// {
//     GReal sqtr = m::sqrt(r);
//     return ((a*a - 2. * a * sqtr + r*r) *
//              ((-2. * a * r * (a*a - 2. * a * sqtr + r*r)) /
//                   m::sqrt(2. * a * sqtr + (-3. + r) * r) +
//               ((a + (-2. + r) * sqtr) * (r*r*r + a*a * (2. + r))) /
//                   m::sqrt(1 + (2. * a) / m::pow(r, 1.5) - 3. / r))) /
//             (r*r*r * m::sqrt(2. * a * sqtr + (-3. + r) * r) *
//              (a*a + (-2. + r) * r));
// }

/**
 * Torus solution for density at a given location.
 * 
 * This function is *not* used for the actual initialization (where rho is calculated
 * alongside the other primitive variables).  Rather, it is for:
 * 1. Normalization, in which the max of this function over the domain is calculated.
 * 2. B field initialization, which requires density of the untilted disk for simplicity
 */
KOKKOS_INLINE_FUNCTION Real kz_torus_rho(const GReal a, const GReal kzeta, const GReal rin, const GReal rmax, const GReal l, const Real gam,
                                         const Real kappa, const GReal r, const GReal th)
{
    // Real l = lfish_calc(a, rmax);
    Real lnh = lnh_calc(a, kzeta, l, rin, r, th);
    if (lnh >= 0. && r >= rin) {
        // Calculate rho
        Real hm1 = m::exp(lnh) - 1.;
        return m::pow(hm1 * (gam - 1.) / (kappa * gam),
                            1. / (gam - 1.));
    } else {
        return 0;
    }
}
