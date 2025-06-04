if [[ $HOST == *".gravity" ]]; then
    MPI_EXE=mpirun
    NPROC=16
    
    if [[ "$ARGS" == *"hip"* ]]; then
        module purge
        module load cmake/3.29.2 ucc/1.3-rocm ucx/1.19-rocm openmpi/5.0.5-rocm rocm/6.1.3 hdf5/1.14.3-openmpi-5.0.5-rocm
        
        # OpenMPI-GPU支持
        export UCX_MEMTYPE_CACHE=n
        export UCX_IB_GPU_DIRECT_RDMA=n  # 添加这行
        export UCX_TLS=sm,self,rocm_copy,rocm_ipc
        export UCX_RNDV_SCHEME=put_zcopy
        export UCX_RNDV_THRESH=16384
        export UCX_MAX_RNDV_RAILS=1
        export UCX_MEMTYPE_REG_WHOLE_ALLOC_TYPES=rocm  # 添加这行
        # MPI_EXTRA_ARGS="--gpu-bind=closest"
        
        # MPI运行参数
        MPI_NUM_PROCS=${MPI_NUM_PROCS:-4}
        # MPI_EXTRA_ARGS="--mca pml ucx --mca btl ^vader,tcp,openib --mca opal_hip_support 1 --bind-to none"
        
        # 架构设置
        HOST_ARCH=Kokkos_ARCH_SKX
        DEVICE_ARCH=AMD_GFX906
        C_NATIVE=hipcc
        CXX_NATIVE=hipcc
        
        # 额外编译标志
        EXTRA_FLAGS="-DPARTHENON_ENABLE_GPU_MPI_CHECKS=OFF $EXTRA_FLAGS"
    fi
fi