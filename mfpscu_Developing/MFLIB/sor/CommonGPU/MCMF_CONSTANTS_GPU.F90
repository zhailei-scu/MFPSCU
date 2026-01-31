!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MCMF_CONSTANTS_GPU
    use cudafor
    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
    use MCMF_CONSTANTS

    implicit none

    !***Box Info

    real(kind=KMCSF),constant::dm_BOXBOUNDARY(3,2)                                ! siumulation box boundary, in unit of atomic radiua
    integer,constant::dm_BDCTYPE(3,2)
    real(kind=KMCSF),constant::dm_BOXSIZE(3)                                      ! simulation box size
    real(kind=KMCSF),constant::dm_HBOXSIZE(3)                                     ! half box size
    integer,constant::dm_PERIOD(3)
    real(kind=KMCDF),constant::dm_TKB


    !********************GPU Parameters******************************
    integer,parameter::p_BLOCKSIZE = 128
    integer,parameter::p_BLOCKDIMX = 256
    integer,parameter::p_NUMINSEG = 10000

    !integer,parameter::p_ReduceAllStatu_BLOCKSIZE = 256
    integer,parameter::p_Reduce_BLOCKSIZE = 512

    contains


    subroutine copyInBoxParamsConstant(BOXBOUNDARY,BOXSIZE,HBOXSIZE)
        implicit none
        !---dummy Vars----
        real(kind=KMCDF),intent(in)::BOXBOUNDARY(3,2)
        real(kind=KMCDF),intent(in)::BOXSIZE(3)
        real(kind=KMCDF),intent(in)::HBOXSIZE(3)
        !---Body----
        !*** copy to device constant memory
        dm_BOXBOUNDARY = BOXBOUNDARY
        dm_BOXSIZE  = BOXSIZE
        dm_HBOXSIZE = HBOXSIZE

        return
    end subroutine

    subroutine copyInPhyParamsConstant(Host_SimuCtrlParam)
        implicit none
        !---dummy Vars----
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Body----
        !*** copy to device constant memory
        dm_PERIOD  = Host_SimuCtrlParam%PERIOD
        dm_BDCTYPE = Host_SimuCtrlParam%BDCTYPE
        dm_TKB = Host_SimuCtrlParam%TKB

        return
    end subroutine

end module MCMF_CONSTANTS_GPU
