module MFLIB_CONSTANTS
    implicit none

    integer,parameter::KMCSF = 4
    integer,parameter::KMCDF = 8



    real(kind=KMCDF), parameter::PI=3.1415926535897932


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
    real(kind=KMCDF),parameter::CP_KB = 1.38054D-16                  ! ERG/K
    real(kind=KMCDF), parameter::C_EV2ERG   = 1.60219D-12

    !*** Memory management *****************
    integer, parameter::C_BYTE = 8    ! (8 bits)
    integer, parameter::C_KBYTES = 1024*C_BYTE
    integer, parameter::C_MBYTES = 1024*C_KBYTES
    integer, parameter::C_GBYTES = 1024*C_MBYTES



    #ifdef WIN_FILE
    character(len=1),parameter::FolderSpe = achar(94)
    #else
    character(len=1),parameter::FolderSpe = achar(47)
    #endif


end module MFLIB_CONSTANTS
