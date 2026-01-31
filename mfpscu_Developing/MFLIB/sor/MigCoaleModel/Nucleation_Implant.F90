!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module NUCLEATION_IMPLANT
!    use MFLIB_GLOBAL
!    use MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
!    implicit none
!    !---Local Vars---
!    character*256::fileName
!    !---Body---
!
!    real(kind=KMCDF),dimension(:),allocatable::NBPV
!
!    real(kind=KMCDF),dimension(:),allocatable::NATOMS
!
!    integer::Dumplicate = 1
!
!    integer::m_StatisticFile
!    integer::LastOutPutConifgIndex = 0
!    contains
!
!    subroutine InitSimu_IMPLANT(Host_SimuCtrlParam)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        !---Local Vars---
!        integer::I
!        !---Body---
!
!        call AvailableIOUnit(m_StatisticFile)
!
!        fileName = Host_SimuCtrlParam%OutFilePath(1:len_trim(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"Statistic.dat"
!
!        pause
!
!        open(Unit=m_StatisticFile,file=fileName(1:len_trim(fileName)))
!
!        write(m_StatisticFile,FMT="(15(A15,1x))") "Step","Time","TStep","ImplantedNum","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"
!
!        return
!    end subroutine
!
!    !***************************************************
!    subroutine NucleationSimu_IMPLANT(Host_SimuCtrlParam)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        !---Local Vars---
!        integer::CKind
!        integer::IKind
!        integer::JKind
!        integer::ITIME
!        real(kind=KMCDF)::TSTEP
!        real(kind=KMCDF)::TTIME
!        real(kind=KMCDF)::deta
!        real(kind=KMCDF),dimension(:),allocatable::tempNBPV
!        real(kind=KMCDF),dimension(:),allocatable::tempNBPVChangeRate
!        real(kind=KMCDF)::Pre
!        real(kind=KMCDF)::NPOWER0Ave
!        real(kind=KMCDF)::NPOWER1DIV2Ave
!        real(kind=KMCDF)::NPOWER1Ave
!        real(kind=KMCDF)::NPOWER3DIV2Ave
!        real(kind=KMCDF)::MeanR
!        real(kind=KMCDF)::MeanNA
!        real(kind=KMCDF)::MeanVolum
!        real(kind=KMCDF)::NewAddedAtoms
!        real(kind=KMCDF)::ImplantedNum
!        logical::startAnnealing
!        integer::I
!        real(kind=KMCDF)::Factor
!        real(kind=KMCDF)::tempTimeStep
!        real(kind=KMCDF)::minTimeStep
!        logical::adjustlTime
!        !---Body---
!        allocate(tempNBPV(CKind))
!
!        allocate(tempNBPVChangeRate(CKind))
!
!        ITIME = 0
!
!        TTIME = 0.D0
!
!        TSTEP = 0.01
!
!        ImplantedNum = 0.D0
!
!        Pre = 4*PI*m_SURDIFPRE*(m_RNFACTOR)**(1.5D0)
!
!        pause
!
!        startAnnealing = .false.
!
!        call Cal_Statistic_IMPLANT(NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
!
!        DO While(.true.)
!
!            adjustlTime = .false.
!
!            tempNBPVChangeRate = 0.D0
!
!            minTimeStep = 1.D32
!
!            NewAddedAtoms = 0.D0
!
!            if((NPOWER1Ave .GE. m_TargetConCentrat .AND. startAnnealing .eq. .false.) .or. (m_Flux .LE. 0.D0 .AND. startAnnealing .eq. .false.)) then
!                startAnnealing = .true.
!                write(*,*) "Start Annealing..."
!            end if
!
!            DO IKind = 1,CKind
!
!                if(IKind .eq. 1 .AND. startAnnealing .eq. .false.) then
!                    NewAddedAtoms = TSTEP*m_Flux
!                    NBPV(IKind) = NBPV(IKind) + NewAddedAtoms
!                end if
!
!                DO JKind = IKind,CKind
!
!                        deta = Dumplicate*NBPV(IKind)*NBPV(JKind)*(DSQRT(NATOMS(IKind)) + DSQRT(NATOMS(JKind)))*(NATOMS(IKind)**(-2) + NATOMS(JKind)**(-2))*Pre
!
!                        if(IKind .eq. JKind) then
!
!                            Factor = 0.5D0
!
!                            tempNBPVChangeRate(IKind) =  tempNBPVChangeRate(IKind) - deta
!
!                        else
!
!                            Factor = 1.D0
!
!                            tempNBPVChangeRate(IKind) =  tempNBPVChangeRate(IKind) - deta
!
!                            tempNBPVChangeRate(JKind) =  tempNBPVChangeRate(JKind) - deta
!
!                        end if
!
!                        if((IKind + JKind) .LE. CKind) then
!                            tempNBPVChangeRate(IKind + JKind) = tempNBPVChangeRate(IKind + JKind) + Factor*deta*(NATOMS(IKind) + NATOMS(JKind))/NATOMS(IKind+JKind)
!                        end if
!
!                END DO
!
!            END DO
!
!            TSTEP = maxval(NBPV)*m_MaxReactChangeRate/maxval(dabs(tempNBPVChangeRate))
!
!            tempNBPV = NBPV
!
!            NBPV = NBPV + tempNBPVChangeRate*TSTEP
!
!            DO I = 1,CKind
!                if(NBPV(I) .LT. 0) then
!                    tempTimeStep = tempNBPV(I)/dabs(tempNBPVChangeRate(I))
!                    if(tempTimeStep .LE. minTimeStep) then
!                        minTimeStep = tempTimeStep
!                    end if
!                    adjustlTime = .true.
!                end if
!            END DO
!
!            if(adjustlTime .eq. .true.) then
!
!                NBPV = tempNBPV
!                TSTEP = minTimeStep
!                NBPV = NBPV + tempNBPVChangeRate*TSTEP
!            end if
!
!            TTIME = TTIME + TSTEP
!            ITIME = ITIME + 1
!
!            ImplantedNum = ImplantedNum + NewAddedAtoms
!
!            if(mod(ITIME,1) .eq. 0) then
!
!                call Cal_Statistic_IMPLANT(NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
!
!                call Put_Out_IMPLANT(Host_SimuCtrlParam,ITIME,TTIME,TSTEP,ImplantedNum,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
!            end if
!
!            !if(NBPV(CKind) .GT. 1.D-10) then
!            if(DSQRT(NATOMS(CKind))*NBPV(CKind) .GT. NPOWER1DIV2Ave*m_DumplicateFactor) then
!
!                DO I = 1,(CKind -1)/2 + 1
!                    if(2*I .LE. CKind) then
!                        NBPV(I) = (NBPV(2*I - 1)*NATOMS(2*I - 1) + NBPV(2*I)*NATOMS(2*I))/(2.D0*NATOMS(2*I))
!                    else
!                        NBPV(I) = NBPV(2*I - 1)
!                    end if
!                END DO
!
!                NBPV((CKind -1)/2+2:CKind) = 0.D0
!
!                DO I = 1,(CKind -1)/2 + 1
!                    if(2*I .LE. CKind) then
!                        NATOMS(I) = NATOMS(2*I)
!                    else
!                        NATOMS(I) = NATOMS(2*I - 1)
!                    end if
!                END DO
!
!                DO I = (CKind -1)/2 + 2,CKind
!                    NATOMS(I) = 2*NATOMS(I)
!                END DO
!
!                Dumplicate = Dumplicate*2
!                write(*,*) "Dumplicate",Dumplicate
!
!            end if
!
!
!            if(TTIME .GT. m_TermTime) then
!                exit
!            end if
!
!        END DO
!
!    end subroutine
!
!    !---------------------------------------------------------------------
!    subroutine Cal_Statistic_IMPLANT(NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
!        implicit none
!        !---Dummy Vars---
!        real(kind=KMCDF)::NPOWER0Ave
!        real(kind=KMCDF)::NPOWER1DIV2Ave
!        real(kind=KMCDF)::NPOWER1Ave
!        real(kind=KMCDF)::NPOWER3DIV2Ave
!        !---Local Vars---
!        real(kind=KMCDF)::Temp
!        integer::I
!        !---Body---
!        NPOWER0Ave = Dumplicate*sum(NBPV)
!
!        Temp = 0.D0
!        DO I = 1,CKind
!            Temp = Temp + (NATOMS(I)**(0.5D0))*NBPV(I)
!        END DO
!        NPOWER1DIV2Ave = Dumplicate*Temp
!
!        Temp = 0.D0
!        DO I = 1,CKind
!            Temp = Temp + NATOMS(I)*NBPV(I)
!        END DO
!        NPOWER1Ave = Dumplicate*Temp
!
!        Temp = 0.D0
!        DO I = 1,CKind
!            Temp = Temp + (NATOMS(I)**(1.5D0))*NBPV(I)
!        END DO
!        NPOWER3DIV2Ave = Dumplicate*Temp
!
!        return
!    end subroutine
!
!    !---------------------------------------------
!    subroutine Put_Out_IMPLANT(Host_SimuCtrlParam,Step,TTime,TStep,ImplantedNum,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        integer,intent(in)::Step
!        real(kind=KMCDF),intent(in)::TTime
!        real(kind=KMCDF),intent(in)::TStep
!        real(kind=KMCDF),intent(in)::ImplantedNum
!        real(kind=KMCDF)::NPOWER0Ave
!        real(kind=KMCDF)::NPOWER1DIV2Ave
!        real(kind=KMCDF)::NPOWER1Ave
!        real(kind=KMCDF)::NPOWER3DIV2Ave
!        !---Local Vars---
!        integer::hFile
!        real(kind=KMCDF)::N1
!        real(kind=KMCDF)::N2
!        real(kind=KMCDF)::N3
!        character*256::C_Index
!        character*256::fileName
!        integer::I
!        real(kind=KMCDF)::Rave
!        !---Body---
!
!        N1 = (NPOWER1DIV2Ave/NPOWER0Ave)**2
!
!        N2 = NPOWER1Ave/NPOWER0Ave
!
!        N3 = (NPOWER3DIV2Ave/NPOWER0Ave)**(dble(2)/dble(3))
!
!        Rave = DSQRT(N1/m_RNFACTOR)
!
!        write(6,FMT="(15(A15,1x))") "Step","Time","TStep","ImplantedNum","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"
!
!        write(6,FMT="(I15,1x,15(E15.5,1x))") Step,TTime,TStep,ImplantedNum,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave*C_CM2NM
!
!        write(m_StatisticFile,FMT="(I15,1x,15(E15.5,1x))") Step,TTime,TStep,ImplantedNum,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave*C_CM2NM
!
!        Flush(m_StatisticFile)
!
!
!        LastOutPutConifgIndex = LastOutPutConifgIndex + 1
!
!
!            call AvailableIOUnit(hFile)
!
!            write(C_Index,*) LastOutPutConifgIndex
!
!            C_Index = adjustl(C_Index)
!
!            fileName = Host_SimuCtrlParam%OutFilePath(1:len_trim(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"Config"//C_Index(1:len_trim(C_Index))//".dat"
!
!            open(Unit=hFile,file=trim(fileName))
!
!            write(hFile,*) "Time",TTime
!
!            write(hFile,FMT="(4(A15))") "n","NBPV","u","Z(u)"
!
!            DO I = 1,CKind
!                write(hFile,FMT="(15(E15.5))") NATOMS(I),NBPV(I),NATOMS(I)*TTime**(-dble(2)/dble(5)),NBPV(I)/(TTime**(-dble(4)/dble(5)))
!            END DO
!
!            close(hFile)
!
!
!        return
!    end subroutine
!
end module NUCLEATION_IMPLANT
