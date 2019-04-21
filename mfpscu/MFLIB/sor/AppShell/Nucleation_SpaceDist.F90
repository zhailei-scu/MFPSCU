module NUCLEATION_SPACEDIST
    use MFLIB_GLOBAL
    implicit none
    !---Local Vars---
    character*256::fileName
    !---Body---

    real(kind=KMCDF),dimension(:,:),allocatable::NBPV

    integer::m_ImplantLayerNum = 10
    real(kind=KMCDF),dimension(:),allocatable::m_ImplantSpaceDist
    real(kind=KMCDF),dimension(:),allocatable::m_ImplantDist

    real(kind=KMCDF),dimension(:),allocatable::NATOMS

    real(kind=KMCDF),dimension(:),allocatable::MatA
    real(kind=KMCDF),dimension(:),allocatable::MatB
    real(kind=KMCDF),dimension(:),allocatable::MatC
    real(kind=KMCDF),dimension(:),allocatable::MatD
    real(kind=KMCDF),dimension(:),allocatable::MatW
    real(kind=KMCDF),dimension(:),allocatable::MatH

    integer::Dumplicate = 1

    integer::m_StatisticFile
    integer::LastOutPutConifgIndex = 0
    contains

    subroutine InitSimu_SpaceDist()
        implicit none
        !---Dummy Vars---
        integer::I
        !---Body---

        call Load_SimulationParams()

        if(allocated(NBPV)) then
            deallocate(NBPV)
        end if
        allocate(NBPV(m_BKind,m_NNodes))

        NBPV = 0.D0

        !NBPV(1,1:m_NNodes) = m_ConCentrat0
        NBPV(1,1) = m_ConCentrat0

        if(allocated(NATOMS)) then
            deallocate(NATOMS)
        end if
        allocate(NATOMS(m_BKind))

        DO I = 1,m_BKind
            NATOMS(I) = I
        END DO

        if(allocated(MatA)) then
            deallocate(MatA)
        end if
        allocate(MatA(m_NNodes))

        if(allocated(MatB)) then
            deallocate(MatB)
        end if
        allocate(MatB(m_NNodes))

        if(allocated(MatC)) then
            deallocate(MatC)
        end if
        allocate(MatC(m_NNodes))

        if(allocated(MatD)) then
            deallocate(MatD)
        end if
        allocate(MatD(m_NNodes))

        if(allocated(MatW)) then
            deallocate(MatW)
        end if
        allocate(MatW(m_NNodes))

        if(allocated(MatH)) then
            deallocate(MatH)
        end if
        allocate(MatH(m_NNodes))

        if(allocated(m_ImplantSpaceDist)) then
            deallocate(m_ImplantSpaceDist)
        end if
        allocate(m_ImplantSpaceDist(m_ImplantLayerNum))
        m_ImplantSpaceDist = m_NodeSpace

        if(allocated(m_ImplantDist)) then
            deallocate(m_ImplantDist)
        end if
        allocate(m_ImplantDist(m_ImplantLayerNum))
        m_ImplantDist = 1.D0/m_ImplantLayerNum


        m_StatisticFile = AvailableIOUnit()

        fileName = OutFilePath(1:len_trim(OutFilePath))//FolderSpe//"Statistic.dat"

        write(*,*) "OutFilePath",OutFilePath

        write(*,*) "fileName",fileName

        open(Unit=m_StatisticFile,file=fileName(1:len_trim(fileName)))

        write(m_StatisticFile,FMT="(15(A15,1x))") "Step","Time","TStep","ImplantedNum","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"

        return
    end subroutine InitSimu_SpaceDist

    !***************************************************
    subroutine NucleationSimu_SpaceDist()
        implicit none
        !---Local Vars---
        integer::CKind
        integer::IKind
        integer::JKind
        integer::ITIME
        real(kind=KMCDF)::TSTEP
        real(kind=KMCDF)::TTIME
        real(kind=KMCDF)::deta
        real(kind=KMCDF),dimension(:,:),allocatable::tempNBPV
        real(kind=KMCDF),dimension(:,:),allocatable::tempNBPVChangeRate
        real(kind=KMCDF)::Pre
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::N1
        real(kind=KMCDF),dimension(:),allocatable::N2
        real(kind=KMCDF),dimension(:),allocatable::N3
        real(kind=KMCDF),dimension(:),allocatable::Rave
        real(kind=KMCDF)::NewAddedAtoms
        real(kind=KMCDF)::ImplantedNum
        logical::startAnnealing
        integer::INode
        real(kind=KMCDF)::Factor
        real(kind=KMCDF)::tempTimeStep
        real(kind=KMCDF)::minTimeStep
        logical::adjustlTime
        real(kind=KMCDF)::DiffGradient1
        real(kind=KMCDF)::DiffGradient2
        integer::IImplantLayer
        !---Body---
        allocate(tempNBPV(m_BKind,m_NNodes))

        allocate(tempNBPVChangeRate(m_BKind,m_NNodes))

        allocate(NPOWER0Ave(m_NNodes),NPOWER1DIV2Ave(m_NNodes),NPOWER1Ave(m_NNodes),NPOWER3DIV2Ave(m_NNodes))

        allocate(N1(m_NNodes),N2(m_NNodes),N3(m_NNodes),Rave(m_NNodes))

        ITIME = 0

        TTIME = 0.D0

        TSTEP = 0.01

        ImplantedNum = 0.D0

        Pre = 4*PI*m_SURDIFPRE*(m_RNFACTOR)**(1.5D0)

        startAnnealing = .false.

        call Cal_Statistic_IMPLANT(NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

        DO While(.true.)

            adjustlTime = .false.

            tempNBPVChangeRate = 0.D0

            minTimeStep = 1.D32

            NewAddedAtoms = 0.D0

            if((sum(NPOWER1Ave)/m_NNodes .GE. m_TargetConCentrat .AND. startAnnealing .eq. .false.) .or. (m_Flux .LE. 0.D0 .AND. startAnnealing .eq. .false.)) then
                startAnnealing = .true.
                write(*,*) "Start Annealing..."
            end if

            DO INode = 1,m_NNodes

                DO IKind = 1,m_BKind

                    !if(IKind .eq. 1 .AND. startAnnealing .eq. .false.) then
                    !    NewAddedAtoms = TSTEP*m_Flux
                    !    NBPV(IKind,INode) = NBPV(IKind,INode) + NewAddedAtoms
                    !end if

                    DO JKind = IKind,m_BKind

                        deta = Dumplicate*NBPV(IKind,INode)*NBPV(JKind,INode)*(DSQRT(NATOMS(IKind)) + DSQRT(NATOMS(JKind)))*(NATOMS(IKind)**(-2) + NATOMS(JKind)**(-2))*Pre

                        if(IKind .eq. JKind) then

                            Factor = 0.5D0

                            tempNBPVChangeRate(IKind,INode) =  tempNBPVChangeRate(IKind,INode) - deta

                        else

                            Factor = 1.D0

                            tempNBPVChangeRate(IKind,INode) =  tempNBPVChangeRate(IKind,INode) - deta

                            tempNBPVChangeRate(JKind,INode) =  tempNBPVChangeRate(JKind,INode) - deta

                        end if

                        if((IKind + JKind) .LE. m_BKind) then
                            tempNBPVChangeRate(IKind + JKind,INode) = tempNBPVChangeRate(IKind + JKind,INode) + Factor*deta*(NATOMS(IKind) + NATOMS(JKind))/NATOMS(IKind+JKind)
                        end if

                    END DO

                END DO
            END DO

            TSTEP = maxval(NBPV)*m_MaxChangeRate/maxval(dabs(tempNBPVChangeRate))

            tempNBPV = NBPV

            tempNBPV = tempNBPV + tempNBPVChangeRate*TSTEP

            DO IKind = 1,m_BKind
                DO INode = 1,m_NNodes
                    if(tempNBPV(IKind,INode) .LT. 0) then
                        tempTimeStep = NBPV(IKind,INode)/dabs(tempNBPVChangeRate(IKind,INode))
                        if(tempTimeStep .LE. minTimeStep) then
                            minTimeStep = tempTimeStep
                        end if
                        adjustlTime = .true.
                    end if
                END DO
            END DO

            if(adjustlTime .eq. .true.) then
                TSTEP = minTimeStep
            end if

            DO IKind = 1,m_BKind

                MatA = 0.D0
                MatB = 0.D0
                MatC = 0.D0
                MatD = 0.D0

                DO INode = 1,m_NNodes
                    if(INode .eq. 1) then
                        !For surface, the Dirichlet boundary condition is applied
                        DiffGradient1 = m_SURDIFPRE*(NATOMS(IKind)**(-2))*(m_RNFACTOR**2)/m_NodeSpace
                        DiffGradient2 = m_SURDIFPRE*(NATOMS(IKind)**(-2) + NATOMS(IKind)**(-2))*(m_RNFACTOR**2)/(m_NodeSpace + m_NodeSpace)

!                        write(*,*) "m_SURDIFPRE",m_SURDIFPRE
!                        write(*,*) "m_NodeSpace",m_NodeSpace
!                        write(*,*) "m_RNFACTOR",m_RNFACTOR
!                        write(*,*) "NATOMS(IKind)",NATOMS(IKind)
!                        write(*,*) "DiffGradient1",DiffGradient1
!                        write(*,*) "DiffGradient2",DiffGradient2
!
!                        pause

                        MatA(INode) = 0.D0
                        MatB(INode) = m_NodeSpace/TSTEP + (DiffGradient1 + DiffGradient2)
                        MatC(INode) = -DiffGradient2
                        MatD(INode) = NBPV(IKind,INode)*m_NodeSpace/TSTEP + tempNBPVChangeRate(IKind,INode)*TSTEP
                    else
                        !For surface, the Dirichlet boundary condition is applied
                        DiffGradient1 = m_SURDIFPRE*(NATOMS(IKind)**(-2) + NATOMS(IKind)**(-2))*(m_RNFACTOR**2)/(m_NodeSpace + m_NodeSpace)
                        DiffGradient2 = m_SURDIFPRE*(NATOMS(IKind)**(-2) + NATOMS(IKind)**(-2))*(m_RNFACTOR**2)/(m_NodeSpace + m_NodeSpace)
                        MatA(INode) = -DiffGradient1
                        MatB(INode) = m_NodeSpace/TSTEP + (DiffGradient1 + DiffGradient2)
                        MatC(INode) = -DiffGradient2
                        MatD(INode) = NBPV(IKind,INode)*m_NodeSpace/TSTEP + tempNBPVChangeRate(IKind,INode)*TSTEP
                    end if

                    if(IKind .eq. 1 .AND. startAnnealing .eq. .false.) then
                        DO IImplantLayer = 1,m_ImplantLayerNum
                            if(sum(m_ImplantSpaceDist(1:IImplantLayer)) .GT. (INode-1)*m_NodeSpace .AND.  &
                                sum(m_ImplantSpaceDist(1:IImplantLayer)) .LE. INode*m_NodeSpace) then
                                    MatD(INode) = MatD(INode) + TSTEP*m_Flux*m_ImplantDist(IImplantLayer)
                                    NewAddedAtoms = TSTEP*m_Flux*m_ImplantDist(IImplantLayer)
                            end if
                        END DO
                    end if

                END DO


                call ReloveTridag(MatA,MatB,MatC,MatD,NBPV(IKind,1:m_NNodes),m_NNodes,MatW,MatH)

                if(IKind .eq. 1) then

                    DO INode = 1,m_NNodes
                        write(*,*) "INode",INode,"NBPV(1,INode)",NBPV(1,INode)
                    END DO
                end if

!                write(*,*) "************",IKind,"********************"
!                DO INode = 1,m_NNodes
!                    write(*,*) "NBPV(I,INode)",NBPV(IKind,INode)
!                    write(*,*) "MatA(INode)",MatA(INode)
!                    write(*,*) "MatB(INode)",MatB(INode)
!                    write(*,*) "MatC(INode)",MatC(INode)
!                    write(*,*) "MatD(INode)",MatD(INode)
!                END DO


            END DO

            TTIME = TTIME + TSTEP
            ITIME = ITIME + 1

            ImplantedNum = ImplantedNum + NewAddedAtoms

            if(mod(ITIME,1) .eq. 0) then

                call Cal_Statistic_IMPLANT(NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

                call Put_Out_IMPLANT(ITIME,TTIME,TSTEP,ImplantedNum,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
            end if

            DO INode = 1,m_NNodes
                !if(NBPV(m_BKind) .GT. 1.D10) then
                if(DSQRT(NATOMS(m_BKind))*NBPV(m_BKind,INode) .GT. NPOWER1DIV2Ave(INode)*m_DumplicateFactor) then
                    write(*,*) "INode",INode
                    write(*,*) "NBPV(m_BKind,INode) ",NBPV(m_BKind,INode)
                    write(*,*) "NPOWER1DIV2Ave(INode)",NPOWER1DIV2Ave(INode)


                    DO IKind = 1,(m_BKind -1)/2 + 1
                        if(2*IKind .LE. m_BKind) then
                            if(NBPV(2*IKind - 1,INode) .eq. 0) then
                                NATOMS(IKind) = NATOMS(2*IKind)
                            else if(NBPV(2*IKind,INode) .eq. 0) then
                                NATOMS(IKind) = NATOMS(2*IKind - 1)
                            else
                                NATOMS(IKind) = (NATOMS(2*IKind - 1)*NBPV(2*IKind - 1,INode) + NATOMS(2*IKind)*NBPV(2*IKind,INode))/(NBPV(2*IKind - 1,INode) + NBPV(2*IKind,INode))
                            end if
                        else
                            NATOMS(IKind) = NATOMS(2*IKind - 1)
                        end if
                    END DO

                    DO IKind = (m_BKind -1)/2 + 2,m_BKind
                        NATOMS(IKind) = 2*NATOMS(IKind)
                    END DO

                    DO IKind = 1,(m_BKind -1)/2 + 1
                        if(2*IKind .LE. m_BKind) then
                            NBPV(IKind,INode) = (NBPV(2*IKind - 1,INode) + NBPV(2*IKind,INode))/2.D0
                        else
                            NBPV(IKind,INode) = NBPV(2*IKind - 1,INode)
                        end if
                    END DO

                    NBPV((m_BKind -1)/2+2:m_BKind,INode) = 0.D0

                    Dumplicate = Dumplicate*2
                    write(*,*) "Dumplicate",Dumplicate
                end if
            END DO


            if(TTIME .GT. m_TermTime) then
                exit
            end if

        END DO

    end subroutine NucleationSimu_SpaceDist

    !----------------------------------------------------------------------
    subroutine ReloveTridag(MatrixA,MatrixB,MatrixC,MatrixD,Solver,MatrixSize,w,h)
        implicit none
        !---Dummy Vars---
        real(kind=KMCDF),dimension(:),allocatable::MatrixA
        real(kind=KMCDF),dimension(:),allocatable::MatrixB
        real(kind=KMCDF),dimension(:),allocatable::MatrixC
        real(kind=KMCDF),dimension(:),allocatable::MatrixD
        real(kind=KMCDF)::Solver(:)
        integer,intent(in)::MatrixSize
        real(kind=KMCDF),dimension(:),allocatable::w
        real(kind=KMCDF),dimension(:),allocatable::h
        !---Local Vars---
        integer::I
        !---Body---
        w = 0.D0
        h = 0.D0

        w(1) = MatrixB(1)
        h(1) = MatrixD(1)/w(1)
        DO I = 2,MatrixSize
            w(I) = MatrixB(I) - MatrixC(I)*MatrixA(I)/w(I-1)
            h(I) = (MatrixD(I) - MatrixA(I-1)*h(I-1))/w(I)
        END DO
        Solver(MatrixSize) = h(MatrixSize)

        DO I = MatrixSize-1,1
            Solver(I) = h(I) - MatrixD(I)*Solver(I+1)/w(I)
        END DO


        return
    end subroutine

    !---------------------------------------------------------------------
    subroutine Cal_Statistic_IMPLANT(NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
        implicit none
        !---Dummy Vars---
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        !---Local Vars---
        real(kind=KMCDF)::Temp
        integer::I
        integer::INode
        !---Body---
        DO INode = 1,m_NNodes
            NPOWER0Ave(INode) = Dumplicate*sum(NBPV(1:m_BKind,INode))

            Temp = 0.D0
            DO I = 1,m_BKind
                Temp = Temp + (NATOMS(I)**(0.5D0))*NBPV(I,INode)
            END DO
            NPOWER1DIV2Ave(INode) = Dumplicate*Temp

            Temp = 0.D0
            DO I = 1,m_BKind
                Temp = Temp + NATOMS(I)*NBPV(I,INode)
            END DO
            NPOWER1Ave(INode) = Dumplicate*Temp

            Temp = 0.D0
            DO I = 1,m_BKind
                Temp = Temp + (NATOMS(I)**(1.5D0))*NBPV(I,INode)
            END DO
            NPOWER3DIV2Ave(INode) = Dumplicate*Temp
        END DO

        return
    end subroutine

    !---------------------------------------------
    subroutine Put_Out_IMPLANT(Step,TTime,TStep,ImplantedNum,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::Step
        real(kind=KMCDF),intent(in)::TTime
        real(kind=KMCDF),intent(in)::TStep
        real(kind=KMCDF),intent(in)::ImplantedNum
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::N1
        real(kind=KMCDF),dimension(:),allocatable::N2
        real(kind=KMCDF),dimension(:),allocatable::N3
        real(kind=KMCDF),dimension(:),allocatable::Rave
        !---Local Vars---
        integer::hFile
        character*256::C_Index
        character*256::fileName
        integer::I
        integer::INode
        !---Body---

        N1 = 0.D0
        N2 = 0.D0
        N3 = 0.D0
        Rave = 0.D0

        DO INode = 1,m_NNodes

            if(NPOWER0Ave(INode) .GT. 0) then
                N1(INode) = (NPOWER1DIV2Ave((INode))/NPOWER0Ave(INode))**2

                N2(INode) = NPOWER1Ave(INode)/NPOWER0Ave(INode)

                N3(INode) = (NPOWER3DIV2Ave(INode)/NPOWER0Ave(INode))**(dble(2)/dble(3))

                Rave(INode) = DSQRT(N1(INode)/m_RNFACTOR)
            end if
        END DO


        write(6,FMT="(15(A15,1x))") "Step","Time","TStep","ImplantedNum","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"

        write(6,FMT="(I15,1x,15(E15.5,1x))") Step,TTime,TStep,ImplantedNum,sum(NPOWER0Ave)/m_NNodes,sum(NPOWER1DIV2Ave)/m_NNodes,sum(NPOWER1Ave)/m_NNodes,sum(NPOWER3DIV2Ave)/m_NNodes, &
                                             sum(N1)/m_NNodes,sum(N2)/m_NNodes,sum(N3)/m_NNodes,sum(Rave)*C_CM2NM/m_NNodes

        write(m_StatisticFile,FMT="(I15,1x,15(E15.5,1x))") Step,TTime,TStep,ImplantedNum,sum(NPOWER0Ave)/m_NNodes,sum(NPOWER1DIV2Ave)/m_NNodes,sum(NPOWER1Ave)/m_NNodes,sum(NPOWER3DIV2Ave)/m_NNodes, &
                                             sum(N1)/m_NNodes,sum(N2)/m_NNodes,sum(N3)/m_NNodes,sum(Rave)*C_CM2NM/m_NNodes

        Flush(m_StatisticFile)


        LastOutPutConifgIndex = LastOutPutConifgIndex + 1


            hFile = AvailableIOUnit()

            write(C_Index,*) LastOutPutConifgIndex

            C_Index = adjustl(C_Index)

            fileName = OutFilePath(1:len_trim(OutFilePath))//FolderSpe//"Config"//C_Index(1:len_trim(C_Index))//".dat"

            open(Unit=hFile,file=trim(fileName))

            write(hFile,*) "Time",TTime

            write(hFile,FMT="(4(A15))") "n","NBPV","u","Z(u)"

            DO I = 1,m_BKind
                write(hFile,FMT="(15(E15.5))") NATOMS(I),sum(NBPV(I,1:m_NNodes))/m_NNodes,NATOMS(I)*TTime**(-dble(2)/dble(5)),(sum(NBPV(I,1:m_NNodes))/m_NNodes)/(TTime**(-dble(4)/dble(5)))
            END DO

            close(hFile)


        return
    end subroutine

end module NUCLEATION_SPACEDIST
