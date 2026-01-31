!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MCMF_CONSTANTS

    implicit none

     !Date conversion
    character(len=1),parameter::KEYWORD_HEAD = "&"
    ! Data File type
    character(len=10),parameter::OKMC_OUTCFG_FORMAT18 = "&BOXOKMC18"
    character(len=8),parameter::MF_OUTCFG_FORMAT18 = "&BOXMF18"
    character(len=10),parameter::SPMF_OUTCFG_FORMAT18 = "&BOXSPMF18"

    ! Data type
    integer,parameter::KMCDF  = 8
    integer,parameter::KMCSF   = 4
    integer,parameter::KMCLINT = 8
    integer,parameter::KMCINT  = 4
    integer,parameter::KMCINT2 = 2

    integer,parameter::KMCDF_SIZE_PGI = 8             ! 8 bit for a KMCDF kind
    integer,parameter::DEFAULT_SIZE_PGI = 4


    integer,parameter::mp_CalcNeighborList_NNEAREST = 1
    integer,parameter::mp_CalcNeighborList_RCUT = 2

    integer,parameter::mp_NEIGHBORUPDATEBYSTEP    =  0
    integer,parameter::mp_NEIGHBORUPDATEBYNCREMIND = 1

    integer,parameter::mp_TermTimeFlag_ByStep = 0
    integer,parameter::mp_TermTimeFlag_ByRealTime = 1

    integer,parameter::mp_SelfAdjustlStep_NearestSep = 0   ! the time step is determined by average distance of nearest cluster
    integer,parameter::mp_FixedTimeStep = 1
    integer,parameter::mp_SelfAdjustlStep_AveSep = 2       ! the time step is determined by volume average distance and suppose the clusters distribute uniform in the box

    integer,parameter::mp_UpdateStatisFlag_ByIntervalSteps = 0
    integer,parameter::mp_UpdateStatisFlag_ByIntervalRealTime = 1

    integer,parameter::mp_OutTimeFlag_ByIntervalSteps = 0
    integer,parameter::mp_OutTimeFlag_ByIntervalRealTime = 1
    integer,parameter::mp_OutTimeFlag_ByIntervalTimeMagnification = 2


    integer, parameter::p_DEPT_DIS_Layer = 0                            ! uniform distribution in the box
    integer, parameter::p_DEPT_DIS_BOX = 1                              ! boxed uniform distribution. The box is a subbox of simulated box
    integer, parameter::p_DEPT_DIS_GAS = 2                              ! Gauss distribution in Z-depth


    !*** data structure defining an cluster
    #ifdef ATOMSGROUP
    integer, parameter::p_ATOMS_GROUPS_NUMBER = ATOMSGROUP
    #else
    integer, parameter::p_ATOMS_GROUPS_NUMBER = 3
    #endif

    !---Boundary condition
    integer, parameter::p_Neumann_BDC = 0
    integer, parameter::p_Dirichlet_BDC = 1


    character*9, parameter::p_CDirichlet_BDC = "Dirichlet"
    character*7, parameter::p_CNeumann_BDC = "Neumann"

    character(len=1),parameter::p_ElementsTypeSpe = "@"
    character(len=1),parameter::p_ElementsNumSpe = "#"
    character(len=1),parameter::p_NumRangeSpe = "-"
    character(len=3),parameter::p_InfStr = "INF"


    integer,parameter::p_ReactionCoefficientTypesNum = 2
    integer,parameter::p_ReactionCoefficient_ByValue = 1
    integer,parameter::p_ReactionCoefficient_ByArrhenius = 2

    integer,parameter::p_DiffuseCoefficientTypesNum = 3
    integer,parameter::p_DiffuseCoefficient_ByValue = 1
    integer,parameter::p_DiffuseCoefficient_ByArrhenius = 2
    integer,parameter::p_DiffuseCoefficient_ByBCluster = 3

    integer,parameter::p_ECRTypesNum = 2
    integer,parameter::p_ECR_ByValue = 1
    integer,parameter::p_ECR_ByBCluster = 2

    !*** The clusters type
    integer,parameter::p_NUMBER_OF_STATU = 6
    integer,parameter::p_ACTIVEFREE_STATU = 1
    integer,parameter::p_ACTIVEINGB_STATU = 2
    integer,parameter::p_OUT_DESTROY_STATU = 3
    integer,parameter::p_EXP_DESTROY_STATU = 4
    integer,parameter::p_MIS_DESTROY_STATU = 5
    integer,parameter::p_ABSORBED_STATU = 6
    integer,parameter::p_Empty = 0
    character*20,parameter::p_CStatu(p_NUMBER_OF_STATU) = (/"ACTIVEFREE","ACTIVEINGB","OUT_DESTROY","EXP_DESTRO","MIS_DESTROY","ABSORBED"/)


    real(kind=KMCDF),parameter::p_GAMMA = 4.D0                          ! the parameter for cluster diffusion

    !--- numbers
    real(kind=KMCDF), parameter::ZERO=0
    real(kind=KMCDF), parameter::ONE=1
    real(kind=KMCDF), parameter::TWO=2
    real(kind=KMCDF), parameter::THREE=3
    real(kind=KMCDF), parameter::TEN=10
    real(kind=KMCDF), parameter::HUNDRED=100
    integer(kind=KMCINT), parameter::IHUNDRED=100
    integer(kind=KMCLINT),parameter::TENPOWTHREE = 1*10**3
    integer(kind=KMCLINT),parameter::TENPOWFOUR = 1*10**4
    integer(kind=KMCLINT),parameter::TENPOWFIVE = 1*10**5
    integer(kind=KMCLINT),parameter::TENPOWSIX = 1*10**6
    integer(kind=KMCLINT),parameter::TENPOWSEVEN = 1*10**7
    integer(kind=KMCLINT),parameter::TENPOWEIGHT = 1*10**8

    real(kind=KMCDF),parameter::ZERO_PROXIMITY = 1.D-16


    !*** Math. and  Phys. constants used
    real(kind=KMCDF), parameter::PI=3.1415926535897932
    real(kind=KMCDF), parameter::HALFPI=0.5D0*PI
    real(kind=KMCDF), parameter::TWOPI=2.0D0*PI
    real(kind=KMCDF), parameter::FOURPI=4.0D0*PI
    real(kind=KMCDF), parameter::FRADEG=180.0D0/PI
    real(kind=KMCDF), parameter::A0B=5.29177249D-9, &      !BOHR RADIU
                              AVOG=6.0221367D23            !Avigado constants

    real(kind=KMCDF), parameter::C_FOURBYTHREE = 4.D0/3.D0
    real(kind=KMCDF), parameter::C_4PI_3 = (4.D0*PI)/3.D0
    real(kind=KMCDF), parameter::C_3_4PI = 3.D0/(4.D0*PI)
    real(kind=KMCDF), parameter::C_ONEBYTHREE = 1.D0/3.D0
    real(kind=KMCDF), parameter::C_UM2CM = 1.D-4
    real(kind=KMCDF), parameter::C_CM2UM = 1.D4
    real(kind=KMCDF), parameter::C_NM2CM = 1.D-7
    real(kind=KMCDF), parameter::C_CM2NM = 1.D7
    real(kind=KMCDF), parameter::C_AM2CM = 1.D-8
    real(kind=KMCDF), parameter::C_CM2AM = 1.D8
    real(kind=KMCDF), parameter::C_JPERM2_TO_ERGPERCM2 = 1.D3
    real(kind=KMCDF), parameter::C_KB      =  1.38054D-16              !Boltzmann constant, in ERG/K
    real(kind=KMCDF), parameter::C_EV2ERG   = 1.60219D-12

    !*** Memory management *****************
    integer, parameter::C_BYTE = 8    ! (8 bits)
    integer, parameter::C_KBYTES = 1024*C_BYTE
    integer, parameter::C_MBYTES = 1024*C_KBYTES
    integer, parameter::C_GBYTES = 1024*C_MBYTES

    !***The directory associated parameters*************
    ! "/"  In fact, the "/" and "\" is same in windows and in cygwin bash or command windows, but in linux, should be "/",
    ! and in cygwin environment, if we call the system(mkdir ) command in fortran program, we cannot use the "/", it should be "\"
    ! so, in fortran program, it is necessary to use "/" in linux and use "\" in windows(cygwin)
    #ifdef CYGWIN
    character(len=1), parameter::FolderSpe = achar(92)   ! "\"
    #else
    character(len=1), parameter::FolderSpe = achar(47)   ! "/"
    #endif

    character(len=1), parameter::RelativeHead = achar(46) ! "."


end module MCMF_CONSTANTS
