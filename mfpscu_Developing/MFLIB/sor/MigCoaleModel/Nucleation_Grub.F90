module NUCLEATION_GRUB
!    implicit none
!
!
!    integer,parameter::KMCDF = 8
!    integer,parameter::p_BKind = 200
!    real(kind=KMCDF),parameter::CP_KB = 1.38054D-16                  ! ERG/K
!
!    real(kind=KMCDF),dimension(:),allocatable::NBPAPV
!    real(kind=KMCDF),dimension(:),allocatable::NATOMS
!
!    real(kind=KMCDF)::MaxReactChangeRate = 0.05
!
!    real(kind=KMCDF)::Temperature = 1000                ! (K)
!
!    real(kind=KMCDF)::SurfEnergy = 1700                 ! (ERG)/cm^2
!
!    real(kind=KMCDF)::JumpLength = 2.55D-8              ! (cm)
!
!    real(kind=KMCDF)::SurfDiffCoffe = 1.D-5             ! (cm^2)/s
!
!    real(kind=KMCDF)::Concentration = 1.D20             ! (cm^(-2))
!
!    integer::Dumplicate = 1
!
!    real(kind=KMCDF)::m_DumplicateFactor = 1.D-6
!
!    integer::m_StatisticFile
!    integer::LastOutPutConifgIndex = 0
!    integer::LastOutPutConfigTime = 1
!
!    contains
!
!    subroutine InitSimu()
!        implicit none
!        !---Local Vars---
!        logical::hasOpened
!        integer::I
!        integer::J
!        !---Body---
!        hasOpened = .true.
!
!        if(allocated(NBPAPV)) then
!            deallocate(NBPAPV)
!        end if
!        allocate(NBPAPV(p_BKind))
!
!        NBPAPV = 0.D0
!
!        NBPAPV(1) = 1.D0
!
!        if(allocated(NATOMS)) then
!            deallocate(NATOMS)
!        end if
!        allocate(NATOMS(p_BKind))
!
!        DO I = 1,p_BKind
!            NATOMS(I) = I
!        END DO
!
!
!        DO m_StatisticFile = 10,99
!            INQUIRE(unit=m_StatisticFile,opened=hasOpened)
!
!            if(.not. hasOpened) then
!                exit
!            end if
!        END DO
!
!        open(Unit=m_StatisticFile,file="/home/zhail/MeanFiledResult/Statistic.dat")
!
!        write(m_StatisticFile,FMT="(15(A15,1x))") "Step","Time","ReduceTime","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3"
!
!        return
!    end subroutine InitSimu
!
!
!    subroutine NucleationSimu()
!        implicit none
!        !---Local Vars---
!        integer::CKind
!        integer::IKind
!        integer::JKind
!        integer::ITIME
!        real(kind=KMCDF)::ReduceTime
!        real(kind=KMCDF)::TSTEP
!        real(kind=KMCDF)::TTIME
!        real(kind=KMCDF)::deta
!        real(kind=KMCDF)::ReduceTimeStep
!        real(kind=KMCDF)::tempReduceTimeStep
!        real(kind=KMCDF),dimension(:),allocatable::tempNBPAPVChangeRate
!        real(kind=KMCDF),dimension(:),allocatable::tempNBPAPV
!        real(kind=KMCDF)::Pre
!        real(kind=KMCDF)::NPOWER0Ave
!        real(kind=KMCDF)::NPOWER1DIV2Ave
!        real(kind=KMCDF)::NPOWER1Ave
!        real(kind=KMCDF)::NPOWER3DIV2Ave
!        real(kind=KMCDF)::MeanR
!        real(kind=KMCDF)::MeanNA
!        real(kind=KMCDF)::MeanVolum
!        logical::adjustlTime
!        real(kind=KMCDF)::minReduceTimeStep
!        real(kind=KMCDF)::Factor
!        integer::I
!        !---Body---
!
!        allocate(tempNBPAPVChangeRate(p_BKind))
!
!        allocate(tempNBPAPV(p_BKind))
!
!        ITIME = 0
!
!        TTIME = 0.D0
!
!        ReduceTime = 0.D0
!
!        Pre = 91.7*(JumpLength**4)*SurfDiffCoffe*Concentration*(SurfEnergy/(CP_KB*Temperature))**(1.5D0)
!
!        DO While(.true.)
!
!            tempNBPAPVChangeRate = 0.D0
!
!            adjustlTime = .false.
!
!            minReduceTimeStep = 1.D32
!
!            DO IKind = 1,p_BKind
!
!                DO JKind = IKind,p_BKind
!
!                        deta = Dumplicate*NBPAPV(IKind)*NBPAPV(JKind)*(DSQRT(NATOMS(IKind)) + DSQRT(NATOMS(JKind)))*(NATOMS(IKind)**(-2) + NATOMS(JKind)**(-2))
!
!                        if(IKind .eq. JKind) then
!
!                            Factor = 0.5D0
!
!                            tempNBPAPVChangeRate(IKind) =  tempNBPAPVChangeRate(IKind) - deta
!
!                        else
!
!                            Factor = 1.D0
!
!                            tempNBPAPVChangeRate(IKind) =  tempNBPAPVChangeRate(IKind) - deta
!
!                            tempNBPAPVChangeRate(JKind) =  tempNBPAPVChangeRate(JKind) - deta
!
!                        end if
!
!                        if((IKind + JKind) .LE. p_BKind) then
!                            tempNBPAPVChangeRate(IKind + JKind) = tempNBPAPVChangeRate(IKind + JKind) + Factor*deta*(NATOMS(IKind) + NATOMS(JKind))/NATOMS(IKind+JKind)
!                        end if
!
!
!                END DO
!
!            END DO
!
!            ReduceTimeStep = maxval(NBPAPV)*MaxReactChangeRate/maxval(dabs(tempNBPAPVChangeRate))
!
!            tempNBPAPV = NBPAPV
!
!            NBPAPV = NBPAPV + tempNBPAPVChangeRate*ReduceTimeStep
!
!            DO I = 1,p_BKind
!                if(NBPAPV(I) .LT. 0) then
!                    tempReduceTimeStep = tempNBPAPV(I)/dabs(tempNBPAPVChangeRate(I))
!                    if(tempReduceTimeStep .LE. minReduceTimeStep) then
!                        minReduceTimeStep = tempReduceTimeStep
!                    end if
!                    adjustlTime = .true.
!                end if
!            END DO
!
!            if(adjustlTime .eq. .true.) then
!
!                NBPAPV = tempNBPAPV
!                ReduceTimeStep = minReduceTimeStep
!                NBPAPV = NBPAPV + tempNBPAPVChangeRate*ReduceTimeStep
!            end if
!
!            ReduceTime = ReduceTime + ReduceTimeStep
!            TSTEP = ReduceTimeStep/Pre
!
!            TTIME = TTIME + TSTEP
!            ITIME = ITIME + 1
!
!            if(mod(ITIME,1) .eq. 0) then
!
!                call Cal_Statistic(NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
!
!                call Put_Out(ITIME,TTIME,ReduceTime,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
!            end if
!
!            !if(NBPAPV(p_BKind) .GT. 1.D-10) then
!            if(DSQRT(NATOMS(p_BKind))*NBPAPV(p_BKind)*Concentration .GT. NPOWER1DIV2Ave*m_DumplicateFactor) then
!
!                DO I = 1,(p_BKind -1)/2 + 1
!                    if(2*I .LE. p_BKind) then
!                        NBPAPV(I) = (NBPAPV(2*I - 1)*NATOMS(2*I - 1) + NBPAPV(2*I)*NATOMS(2*I))/(2.D0*NATOMS(2*I))
!                    else
!                        NBPAPV(I) = NBPAPV(2*I - 1)
!                    end if
!                END DO
!
!                NBPAPV((p_BKind -1)/2+2:p_BKind) = 0.D0
!
!                DO I = 1,(p_BKind -1)/2 + 1
!                    if(2*I .LE. p_BKind) then
!                        NATOMS(I) = NATOMS(2*I)
!                    else
!                        NATOMS(I) = NATOMS(2*I - 1)
!                    end if
!                END DO
!
!                DO I = (p_BKind -1)/2 + 2,p_BKind
!                    NATOMS(I) = 2*NATOMS(I)
!                END DO
!
!                Dumplicate = Dumplicate*2
!                write(*,*) "Dumplicate",Dumplicate
!
!            end if
!
!            if(ReduceTime .GT. 1.D6) then
!                exit
!            end if
!
!        END DO
!
!    end subroutine
!
!    !---------------------------------------------------------------------
!    subroutine Cal_Statistic(NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
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
!        NPOWER0Ave = Dumplicate*Concentration*sum(NBPAPV)
!
!        Temp = 0.D0
!        DO I = 1,p_BKind
!            Temp = Temp + (NATOMS(I)**(0.5D0))*NBPAPV(I)
!        END DO
!        NPOWER1DIV2Ave = Dumplicate*Concentration*Temp
!
!        Temp = 0.D0
!        DO I = 1,p_BKind
!            Temp = Temp + NATOMS(I)*NBPAPV(I)
!        END DO
!        NPOWER1Ave = Dumplicate*Concentration*Temp
!
!        Temp = 0.D0
!        DO I = 1,p_BKind
!            Temp = Temp + (NATOMS(I)**(1.5D0))*NBPAPV(I)
!        END DO
!        NPOWER3DIV2Ave = Dumplicate*Concentration*Temp
!
!        return
!    end subroutine
!
!    !---------------------------------------------
!    subroutine Put_Out(Step,TTime,ReduceTime,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
!        implicit none
!        !---Dummy Vars---
!        integer,intent(in)::Step
!        real(kind=KMCDF),intent(in)::TTime
!        real(kind=KMCDF),intent(in)::ReduceTime
!        real(kind=KMCDF)::NPOWER0Ave
!        real(kind=KMCDF)::NPOWER1DIV2Ave
!        real(kind=KMCDF)::NPOWER1Ave
!        real(kind=KMCDF)::NPOWER3DIV2Ave
!        !---Local Vars---
!        integer::hFile
!        real(kind=KMCDF)::N1
!        real(kind=KMCDF)::N2
!        real(kind=KMCDF)::N3
!        logical::hasOpened
!        character*256::C_Index
!        character*256::fileName
!        integer::I
!        !---Body---
!
!        hasOpened = .true.
!
!        N1 = (NPOWER1DIV2Ave/NPOWER0Ave)**2
!
!        N2 = NPOWER1Ave/NPOWER0Ave
!
!        N3 = (NPOWER3DIV2Ave/NPOWER0Ave)**(dble(2)/dble(3))
!
!        write(6,FMT="(15(A15,1x))") "Step","Time","ReduceTime","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3"
!
!        write(6,FMT="(I15,1x,15(E15.8,1x))") Step,TTime,ReduceTime,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3
!
!
!        write(m_StatisticFile,FMT="(I15,1x,15(E15.10,1x))") Step,TTime,ReduceTime,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3
!
!        Flush(m_StatisticFile)
!
!            LastOutPutConifgIndex = LastOutPutConifgIndex + 1
!
!            DO hFile = 10,99
!                INQUIRE(unit=hFile,opened=hasOpened)
!
!                if(.not. hasOpened) then
!                    exit
!                end if
!            END DO
!
!            write(C_Index,*) LastOutPutConifgIndex
!
!            C_Index = adjustl(C_Index)
!
!            fileName = "/home/zhail/MeanFiledResult/Config"//C_Index(1:len_trim(C_Index))//".dat"
!
!            open(Unit=hFile,file=trim(fileName))
!
!            write(hFile,*) "Time",TTime
!
!            write(hFile,*) "ReduceTime",ReduceTime
!
!            write(hFile,FMT="(4(A15))") "n","NBPAPV","u","Z(u)"
!
!            DO I = 1,p_BKind
!                write(hFile,FMT="(15(E15.8,1x))") NATOMS(I),NBPAPV(I),NATOMS(I)*ReduceTime**(-dble(2)/dble(5)),NBPAPV(I)/(ReduceTime**(-dble(4)/dble(5)))
!            END DO
!
!            close(hFile)
!
!
!        return
!    end subroutine
!
end module NUCLEATION_GRUB
