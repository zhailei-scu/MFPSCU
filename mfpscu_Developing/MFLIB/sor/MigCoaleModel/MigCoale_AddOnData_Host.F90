!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MIGCOALE_ADDONDATA_HOST
    use MSM_TYPEDEF_InputPaser
    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
    implicit none


    real(kind=KMCDF)::m_FREEDIFCOES(3) = 0.D0                       ! the surface diffusion coeffciency
    real(kind=KMCDF)::m_FREEDIFCOESPRE(3) = 0.D0                    ! the prefactor for the surface diffusion coeffciency
    real(kind=KMCDF)::m_FREEDIFCOESES(3) = 0.D0                     ! the surface ative energy for surface diffusion coeffciency

    real(kind=KMCDF)::m_FREESURDIFPRE = 0.D0

    real(kind=KMCDF)::m_GBDIFCOES(3) = 0.D0                     ! the surface diffusion coeffciency in GB
    real(kind=KMCDF)::m_GBDIFCOESPRE(3) = 0.D0                  ! the prefactor for the surface diffusion coeffciency in GB
    real(kind=KMCDF)::m_GBDIFCOESES(3) = 0.D0                   ! the surface ative energy for surface diffusion coeffciency in GB

    real(kind=KMCDF)::m_GBSURDIFPRE = 0.D0

    real(kind=KMCDF)::m_SURFE= C_JPERM2_TO_ERGPERCM2         ! the surface energy. In ERG/cm**2
    real(kind=KMCDF)::m_RNFACTOR                             ! the factor in relation between radiius and number of atom
                                                             ! 8*PI*SURFE/(3*KB*TEMP)
    logical::m_DumplicateBox = .true.

    contains

    subroutine resolveAddOnData(Host_Boxes,Host_SimuCtrlParam)
        !---Dummy Vars---
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Local Vars---
        integer::LINE
        character*256::STR
        character*32::KEYWORD
        character*20::STRTEMP(10)
        integer::N
        integer::I
        !---Body---

        KEYWORD = "&DUMPLICATEBOX"
        call Get_StatementList(KEYWORD(1:LENTRIM(KEYWORD)), Host_SimuCtrlParam%AddOnData, STR, LINE)
        call RemoveComments(STR,"!")
        call EXTRACT_NUMB(STR,1,N,STRTEMP)
        if(N .LT. 1) then
            write(*,*) "MFPSCUERROR: Too few parameters for dumplicate box strategy at line: ",LINE
            write(*,*) STR
            write(*,*) "You should special: &DUMPLICATEBOX  If use the dumplicate box strategy = "
            pause
            stop
        end if
        if(ISTR(STRTEMP(1)) .eq. 0) then
            m_DumplicateBox = .false.
        else
            m_DumplicateBox = .true.
        end if


        KEYWORD = "&SURENG"
        call Get_StatementList(KEYWORD(1:LENTRIM(KEYWORD)), Host_SimuCtrlParam%AddOnData, STR, LINE)
        call EXTRACT_NUMB(STR,1,N,STRTEMP)
        if(N .LT. 1) then
            write(*,*) "MFPSCUERROR: Too few parameters for surface energy at line: ",LINE
            write(*,*) STR
            write(*,*) "You should special: &SURENG THE SURFACE ENERGY OF A BUBBLE = ! (ERG/CM^2))"
            pause
            stop
        end if
        m_SURFE = DRSTR(STRTEMP(1))
        m_RNFACTOR = 8.D0*PI*m_SURFE/(3.D0*Host_SimuCtrlParam%TKB)


        KEYWORD = "&SURDIF"
        call Get_StatementList(KEYWORD(1:LENTRIM(KEYWORD)), Host_SimuCtrlParam%AddOnData, STR, LINE)
        call EXTRACT_NUMB(STR,6,N,STRTEMP)
        if(N .LT. 6) then
            write(*,*) "MFPSCUERROR: Too few parameters for surface diffusion oarameters at line: ",LINE
            write(*,*) STR
            write(*,*) "You should special: &SURDIF  THE Surface Diffusion coefficiens, prefactor (cm^2/s) and ES(ev): 0.0012, 1.0, 0.0012, 1.0, 0.0012, 1.0"
            pause
            stop
        end if
        DO I=1,3
            m_FREEDIFCOESPRE(I) = DRSTR(STRTEMP(2*I-1))
            m_FREEDIFCOESES(I) =  DRSTR(STRTEMP(2*I))
            m_FREEDIFCOES(I) = m_FREEDIFCOESPRE(I)*DEXP(-m_FREEDIFCOESES(I)*C_EV2ERG/Host_SimuCtrlParam%TKB)
        END DO
        m_FREESURDIFPRE = (3.D0/(2.D0*PI))*(Host_Boxes%MatrixAtom%m_Volum**C_FOURBYTHREE)*m_FREEDIFCOES(1)

        KEYWORD = "&GBSURDIF"
        call Get_StatementList(KEYWORD(1:LENTRIM(KEYWORD)), Host_SimuCtrlParam%AddOnData, STR, LINE)
        call EXTRACT_NUMB(STR,6,N,STRTEMP)
        if(N .LT. 6) then
            write(*,*) "MFPSCUERROR: Too few parameters for surface diffusion parameters in GB at line: ",LINE
            write(*,*) STR
            write(*,*) "You should special: &GBSURDIF THE Surface Diffusion coefficiens in GB, prefactor (cm^2/s) and ES(ev): 0.0012, 1.0, 0.0012, 1.0, 0.0012, 1.0"
            pause
            stop
        end if
        DO I=1,3
            m_GBDIFCOESPRE(I) = DRSTR(STRTEMP(2*I-1))
            m_GBDIFCOESES(I) =  DRSTR(STRTEMP(2*I))
            m_GBDIFCOES(I) = m_GBDIFCOESPRE(I)*DEXP(-m_GBDIFCOESES(I)*C_EV2ERG/Host_SimuCtrlParam%TKB)
        END DO

        m_GBSURDIFPRE = (3.D0/(2.D0*PI))*(Host_Boxes%MatrixAtom%m_Volum**C_FOURBYTHREE)*m_GBDIFCOES(1)

        return
    end subroutine resolveAddOnData

end module MIGCOALE_ADDONDATA_HOST
