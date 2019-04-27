module NUCLEATION_SPACEDIST
    use MFLIB_GLOBAL
    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
    use MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
    use MIGCOALE_ADDONDATA_HOST
    implicit none
    !---Local Vars---
    character*256::fileName
    !---Body---
    real(kind=KMCDF),dimension(:),allocatable::MatA
    real(kind=KMCDF),dimension(:),allocatable::MatB
    real(kind=KMCDF),dimension(:),allocatable::MatC
    real(kind=KMCDF),dimension(:),allocatable::MatD
    real(kind=KMCDF),dimension(:),allocatable::MatW
    real(kind=KMCDF),dimension(:),allocatable::MatH

    real(kind=KMCDF),dimension(:),allocatable::m_ImplantDist
    integer::m_ImplantLayerNum = 1

    integer::Dumplicate = 1

    integer::m_StatisticFile
    integer::LastOutPutConifgIndex = 0
    contains

    subroutine InitSimu_SpaceDist(Host_Boxes,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer::I
        integer::NNodes
        integer::CKind
        !---Body---

        NNodes = Host_Boxes%NNodes

        CKind = Host_Boxes%CKind

        if(allocated(MatA)) then
            deallocate(MatA)
        end if
        allocate(MatA(NNodes))

        if(allocated(MatB)) then
            deallocate(MatB)
        end if
        allocate(MatB(NNodes))

        if(allocated(MatC)) then
            deallocate(MatC)
        end if
        allocate(MatC(NNodes))

        if(allocated(MatD)) then
            deallocate(MatD)
        end if
        allocate(MatD(NNodes))

        if(allocated(MatW)) then
            deallocate(MatW)
        end if
        allocate(MatW(NNodes))

        if(allocated(MatH)) then
            deallocate(MatH)
        end if
        allocate(MatH(NNodes))

        if(allocated(m_ImplantDist)) then
            deallocate(m_ImplantDist)
        end if
        allocate(m_ImplantDist(m_ImplantLayerNum))
        m_ImplantDist = 1.D0/m_ImplantLayerNum


        call AvailableIOUnit(m_StatisticFile)

        fileName = Host_SimuCtrlParam%OutFilePath(1:len_trim(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"Statistic.dat"

        open(Unit=m_StatisticFile,file=fileName(1:len_trim(fileName)))

        write(m_StatisticFile,FMT="(15(A15,1x))") "Step","Time","TStep","ImplantedNum","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"

        return
    end subroutine InitSimu_SpaceDist

    !***************************************************
    subroutine NucleationSimu_SpaceDist_Old(Host_SimBoxes,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
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
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::N1
        real(kind=KMCDF),dimension(:),allocatable::N2
        real(kind=KMCDF),dimension(:),allocatable::N3
        real(kind=KMCDF),dimension(:),allocatable::Rave
        real(kind=KMCDF),dimension(:),allocatable::NewAddedAtoms
        real(kind=KMCDF),dimension(:),allocatable::ImplantedNum
        real(kind=KMCDF),dimension(:),allocatable::FSurfAccum
        real(kind=KMCDF),dimension(:),allocatable::FOutAccum
        real(kind=KMCDF),dimension(:),allocatable::CSurfAccum
        real(kind=KMCDF),dimension(:),allocatable::COutAccum
        real(kind=KMCDF),dimension(:),allocatable::FSurfEachStep
        real(kind=KMCDF),dimension(:),allocatable::FOutEachStep
        real(kind=KMCDF),dimension(:),allocatable::CSurfEachStep
        real(kind=KMCDF),dimension(:),allocatable::COutEachStep
        logical::startAnnealing
        integer::I
        integer::INode
        integer::NNodes
        real(kind=KMCDF)::Factor
        real(kind=KMCDF)::tempTimeStep
        real(kind=KMCDF)::minTimeStep
        logical::adjustlTime
        real(kind=KMCDF)::DiffGradient1
        real(kind=KMCDF)::DiffGradient2
        integer::IImplantLayer
        type(DiffusorValue)::TheDiffusorValue
        integer::AtomuNumbSubject
        integer::AtomuNumbObject
        integer::AtomuNumbProductor
        real(kind=KMCDF)::ConCentrat0
        !---Body---
        CKind = Host_SimBoxes%CKind
        NNodes = Host_SimBoxes%NNodes

        allocate(tempNBPV(CKind,NNodes))

        allocate(tempNBPVChangeRate(CKind,NNodes))

        allocate(NPOWER0Ave(NNodes),NPOWER1DIV2Ave(NNodes),NPOWER1Ave(NNodes),NPOWER3DIV2Ave(NNodes))

        allocate(N1(NNodes),N2(NNodes),N3(NNodes),Rave(NNodes))

        allocate(NewAddedAtoms(NNodes))

        allocate(ImplantedNum(NNodes))

        allocate(FSurfAccum(NNodes))

        allocate(FOutAccum(NNodes))

        allocate(CSurfAccum(NNodes))

        allocate(COutAccum(NNodes))

        allocate(FSurfEachStep(NNodes))

        allocate(FOutEachStep(NNodes))

        allocate(CSurfEachStep(NNodes))

        allocate(COutEachStep(NNodes))

        ITIME = 0

        TTIME = 0.D0

        TSTEP = 0.01

        ImplantedNum = 0.D0

        FSurfAccum = 0.D0

        FOutAccum = 0.D0

        CSurfAccum = 0.D0

        COutAccum = 0.D0

        startAnnealing = .false.

        call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

        Associate(ClustersKind=>Host_SimBoxes%ClustersKindArray,Concent=>Host_SimBoxes%Concentrate,NodeSpace=>Host_SimBoxes%NodeSpace)

          ConCentrat0 = sum(Concent)

          DO While(.true.)

            adjustlTime = .false.

            tempNBPVChangeRate = 0.D0

            minTimeStep = 1.D32

            NewAddedAtoms = 0.D0

            FSurfEachStep = 0.D0

            FOutEachStep = 0.D0

            CSurfEachStep = 0.D0

            COutEachStep = 0.D0

            DO INode = 1,NNodes

                DO IKind = 1,CKind

                    !if(IKind .eq. 1 .AND. startAnnealing .eq. .false.) then
                    !    NewAddedAtoms = TSTEP*m_Flux
                    !    Concent(IKind,INode)(IKind,INode) = Concent(IKind,INode) + NewAddedAtoms
                    !end if

                    DO JKind = IKind,CKind

                        deta = Dumplicate*Concent(IKind,INode)*Concent(JKind,INode)* &
                               (ClustersKind(IKind)%m_RAD + ClustersKind(JKind)%m_RAD)* &
                               (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(JKind)%m_DiffCoeff)

                        if(IKind .eq. JKind) then

                            Factor = 0.5D0

                            tempNBPVChangeRate(IKind,INode) =  tempNBPVChangeRate(IKind,INode) - deta

                        else

                            Factor = 1.D0

                            tempNBPVChangeRate(IKind,INode) =  tempNBPVChangeRate(IKind,INode) - deta

                            tempNBPVChangeRate(JKind,INode) =  tempNBPVChangeRate(JKind,INode) - deta

                        end if

                        if((IKind + JKind) .LE. CKind) then

                            AtomuNumbSubject = sum(ClustersKind(IKind)%m_Atoms(:)%m_NA)
                            AtomuNumbObject = sum(ClustersKind(JKind)%m_Atoms(:)%m_NA)
                            AtomuNumbProductor = sum(ClustersKind(IKind+JKind)%m_Atoms(:)%m_NA)

                            tempNBPVChangeRate(IKind + JKind,INode) = tempNBPVChangeRate(IKind + JKind,INode) + &
                                                                      Factor*deta*(AtomuNumbSubject+AtomuNumbObject)/AtomuNumbProductor
                        end if

                    END DO

                END DO
            END DO

            TSTEP = maxval(Concent)*Host_SimuCtrlParam%MaxChangeRate/maxval(dabs(tempNBPVChangeRate))

            TSTEP = min(TSTEP,1D-8)

            tempNBPV = Concent

            tempNBPV = tempNBPV + tempNBPVChangeRate*TSTEP

            DO IKind = 1,CKind
                DO INode = 1,NNodes
                    if(tempNBPV(IKind,INode) .LT. 0) then
                        tempTimeStep = Concent(IKind,INode)/dabs(tempNBPVChangeRate(IKind,INode))
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

            DO IKind = 1,CKind

                MatA = 0.D0
                MatB = 0.D0
                MatC = 0.D0
                MatD = 0.D0

                DO INode = 1,NNodes
                    if(INode .eq. 1) then  ! upper surface
                        !For surface, the Dirichlet boundary condition is applied
                        DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                        MatA(INode) = 0.D0
                        MatB(INode) = NodeSpace(INode)/TSTEP + (DiffGradient1 + DiffGradient2)
                        MatC(INode) = -DiffGradient2
                        MatD(INode) = Concent(IKind,INode)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(IKind,INode)*TSTEP
                    else if(INode .eq. NNodes) then  ! Low surface
                        !For surface, the Dirichlet boundary condition is applied
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)
                        MatA(INode) = -DiffGradient1
                        MatB(INode) = NodeSpace(INode)/TSTEP + (DiffGradient1 + DiffGradient2)
                        MatC(INode) = 0.D0
                        MatD(INode) = Concent(IKind,INode)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(IKind,INode)*TSTEP
                    else
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))
                        MatA(INode) = -DiffGradient1
                        MatB(INode) = NodeSpace(INode)/TSTEP + (DiffGradient1 + DiffGradient2)
                        MatC(INode) = -DiffGradient2
                        MatD(INode) = Concent(IKind,INode)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(IKind,INode)*TSTEP
                    end if

!                    if(IKind .eq. 1) then
!                        DO IImplantLayer = 1,m_ImplantLayerNum
!                            if(sum(m_ImplantSpaceDist(1:IImplantLayer)) .GT. (INode-1)*m_NodeSpace .AND.  &
!                                sum(m_ImplantSpaceDist(1:IImplantLayer)) .LE. INode*m_NodeSpace) then
!                                    MatD(INode) = MatD(INode) + TSTEP*m_Flux*m_ImplantDist(IImplantLayer)
!                                    NewAddedAtoms(IImplantLayer) = TSTEP*m_Flux*m_ImplantDist(IImplantLayer)
!                            end if
!                        END DO
!                    end if

                END DO

                call SolveTridag(IKind,MatA,MatB,MatC,MatD,Host_SimBoxes%Concentrate,NNodes,MatW,MatH)

                DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(NNodes)
                FOutEachStep(IKind) = DiffGradient2*Concent(IKind,NNodes)
                COutEachStep(IKind) = DiffGradient2*Concent(IKind,NNodes)*TSTEP/NodeSpace(NNodes)
                FOutAccum(IKind) = FOutAccum(IKind) + FOutEachStep(IKind)
                COutAccum(IKind) = COutAccum(IKind) + COutEachStep(IKind)


                DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(1)
                FSurfEachStep(IKind) = DiffGradient1*Concent(IKind,1)
                CSurfEachStep(IKind) = DiffGradient1*Concent(IKind,1)*TSTEP/NodeSpace(1)
                FSurfAccum(IKind) = FSurfAccum(IKind) + FSurfEachStep(IKind)
                CSurfAccum(IKind) = CSurfAccum(IKind) + CSurfEachStep(IKind)

                if(IKind .eq. 1) then

                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out flux from up surface",FSurfAccum(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out concentrate from up surface",CSurfAccum(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Out flux from up surface in current step",FSurfEachStep(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Out concentrate from up surface in current step",CSurfEachStep(IKind)
                    DO INode = 1,NNodes
                        write(*,fmt="(A10,1x,I10,1x,A15,1x,1ES14.6)") "INode",INode,"Concent(1,INode)",Concent(1,INode)
                    END DO
                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out flux from lower surface",FOutAccum(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out concentrate from lower surface",COutAccum(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Out flux from lower surface in current step",FOutEachStep(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Out concentrate from lower surface in current step",COutEachStep(IKind)

                    write(*,*) "----------------------------------------"
                    write(*,*) "sum(Concent(1,:))",sum(Concent(1,:))
                    write(*,*) "sum(Concent(1,:)) + CSurfAccum(1) + COutAccum(1)",sum(Concent(1,:)) + CSurfAccum(1) + COutAccum(1)
                    write(*,*) "(sum(Concent(1,:)) + CSurfAccum(1) + COutAccum(1))/ConCentrat0",(ConCentrat0 - (sum(Concent(1,:)) + CSurfAccum(1) + COutAccum(1)))/ConCentrat0
                    write(*,*) "----------------------------------------"
                end if


            END DO

            TTIME = TTIME + TSTEP
            ITIME = ITIME + 1

            ImplantedNum = ImplantedNum + NewAddedAtoms

            if(mod(ITIME,1) .eq. 0) then

                call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

                call Put_Out_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,ITIME,TTIME,TSTEP,sum(ImplantedNum),NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
            end if

            !if(Concent(CKind) .GT. 1.D-10) then
            if(DSQRT(dble(sum(ClustersKind(CKind)%m_Atoms(:)%m_NA)))*sum(Concent(CKind,1:NNodes)) .GT. &
               sum(NPOWER1DIV2Ave)*Host_SimuCtrlParam%DumplicateFactor) then

                DO INode = 1,NNodes

                    DO I = 1,(CKind -1)/2 + 1
                        if(2*I .LE. CKind) then
                            Concent(I,INode) = (Concent(2*I-1,INode)*sum(ClustersKind(2*I-1)%m_Atoms(:)%m_NA)+Concent(2*I,INode)*sum(ClustersKind(2*I)%m_Atoms(:)%m_NA)) &
                                               /(2.D0*sum(ClustersKind(2*I)%m_Atoms(:)%m_NA))
                        else
                            Concent(I,INode) = Concent(2*I - 1,INode)
                        end if
                    END DO

                    Concent((CKind -1)/2+2:CKind,INode) = 0.D0

                    DO I = 1,(CKind -1)/2 + 1
                        if(2*I .LE. CKind) then
                            ClustersKind(I) = ClustersKind(2*I)
                        else
                            ClustersKind(I) = ClustersKind(2*I - 1)
                        end if
                    END DO

                    DO I = (CKind -1)/2 + 2,CKind
                        ClustersKind(I)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA = 2*ClustersKind(I)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA
                        TheDiffusorValue = Host_SimBoxes%m_DiffusorTypesMap%get(ClustersKind(I))

                        select case(TheDiffusorValue%ECRValueType_Free)
                            case(p_ECR_ByValue)
                                ClustersKind(I)%m_RAD = TheDiffusorValue%ECR_Free
                            case(p_ECR_ByBCluster)
                                ClustersKind(I)%m_RAD = DSQRT(sum(ClustersKind(I)%m_Atoms(:)%m_NA)/m_RNFACTOR)
                        end select

                        select case(TheDiffusorValue%DiffusorValueType_Free)
                            case(p_DiffuseCoefficient_ByValue)
                                ClustersKind(I)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                            case(p_DiffuseCoefficient_ByArrhenius)
                                ClustersKind(I)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
                            case(p_DiffuseCoefficient_ByBCluster)
                                ! Here we adopt a model that D=D0*(1/R)**Gama
                                ClustersKind(I)%m_DiffCoeff = m_FREESURDIFPRE*(ClustersKind(I)%m_RAD**(-p_GAMMA))
                        end select
                    END DO

                END DO

                Dumplicate = Dumplicate*2
                write(*,*) "Dumplicate",Dumplicate
            end if


            if(TTIME .GT. Host_SimuCtrlParam%TermTValue) then
                exit
            end if

          END DO
        END Associate

    end subroutine NucleationSimu_SpaceDist_Old

    !***************************************************
    subroutine NucleationSimu_SpaceDist(Host_SimBoxes,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
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
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::N1
        real(kind=KMCDF),dimension(:),allocatable::N2
        real(kind=KMCDF),dimension(:),allocatable::N3
        real(kind=KMCDF),dimension(:),allocatable::Rave
        real(kind=KMCDF),dimension(:),allocatable::NewAddedAtoms
        real(kind=KMCDF),dimension(:),allocatable::ImplantedNum
        real(kind=KMCDF),dimension(:),allocatable::FSurfAccum
        real(kind=KMCDF),dimension(:),allocatable::FOutAccum
        real(kind=KMCDF),dimension(:),allocatable::CSurfAccum
        real(kind=KMCDF),dimension(:),allocatable::COutAccum
        real(kind=KMCDF),dimension(:),allocatable::FSurfEachStep
        real(kind=KMCDF),dimension(:),allocatable::FOutEachStep
        real(kind=KMCDF),dimension(:),allocatable::CSurfEachStep
        real(kind=KMCDF),dimension(:),allocatable::COutEachStep
        logical::startAnnealing
        integer::NNodes
        integer::INode
        real(kind=KMCDF)::Factor
        real(kind=KMCDF)::tempTimeStep
        real(kind=KMCDF)::minTimeStep
        logical::adjustlTime
        real(kind=KMCDF)::DiffGradient1
        real(kind=KMCDF)::DiffGradient2
        integer::IImplantLayer
        integer::I
        type(DiffusorValue)::TheDiffusorValue
        integer::AtomuNumbSubject
        integer::AtomuNumbObject
        integer::AtomuNumbProductor
        real(kind=KMCDF)::ConCentrat0
        !---Body---
        CKind = Host_SimBoxes%CKind

        NNodes = Host_SimBoxes%NNodes

        allocate(tempNBPV(CKind,NNodes))

        allocate(tempNBPVChangeRate(CKind,NNodes))

        allocate(NPOWER0Ave(NNodes),NPOWER1DIV2Ave(NNodes),NPOWER1Ave(NNodes),NPOWER3DIV2Ave(NNodes))

        allocate(N1(NNodes),N2(NNodes),N3(NNodes),Rave(NNodes))

        allocate(NewAddedAtoms(NNodes))

        allocate(ImplantedNum(NNodes))

        allocate(FSurfAccum(NNodes))

        allocate(FOutAccum(NNodes))

        allocate(CSurfAccum(NNodes))

        allocate(COutAccum(NNodes))

        allocate(FSurfEachStep(NNodes))

        allocate(FOutEachStep(NNodes))

        allocate(CSurfEachStep(NNodes))

        allocate(COutEachStep(NNodes))

        ITIME = 0

        TTIME = 0.D0

        TSTEP = 0.01

        ImplantedNum = 0.D0

        FSurfAccum = 0.D0

        FOutAccum = 0.D0

        CSurfAccum = 0.D0

        COutAccum = 0.D0

        startAnnealing = .false.

        call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

        Associate(ClustersKind=>Host_SimBoxes%ClustersKindArray,Concent=>Host_SimBoxes%Concentrate,NodeSpace=>Host_SimBoxes%NodeSpace)

          ConCentrat0 = sum(Concent)

          DO While(.true.)

            adjustlTime = .false.

            tempNBPVChangeRate = 0.D0

            minTimeStep = 1.D32

            NewAddedAtoms = 0.D0

            FSurfEachStep = 0.D0

            FOutEachStep = 0.D0

            CSurfEachStep = 0.D0

            COutEachStep = 0.D0

            DO INode = 1,NNodes

                DO IKind = 1,CKind

                    !if(IKind .eq. 1 .AND. startAnnealing .eq. .false.) then
                    !    NewAddedAtoms = TSTEP*m_Flux
                    !    Concent(IKind,INode) = Concent(IKind,INode) + NewAddedAtoms
                    !end if

                    DO JKind = IKind,CKind

                        deta = Dumplicate*Concent(IKind,INode)*Concent(JKind,INode)* &
                               (ClustersKind(IKind)%m_RAD + ClustersKind(JKind)%m_RAD)* &
                               (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(JKind)%m_DiffCoeff)

                        if(IKind .eq. JKind) then

                            Factor = 0.5D0

                            tempNBPVChangeRate(IKind,INode) =  tempNBPVChangeRate(IKind,INode) - deta

                        else

                            Factor = 1.D0

                            tempNBPVChangeRate(IKind,INode) =  tempNBPVChangeRate(IKind,INode) - deta

                            tempNBPVChangeRate(JKind,INode) =  tempNBPVChangeRate(JKind,INode) - deta

                        end if

                        if((IKind + JKind) .LE. CKind) then
                            AtomuNumbSubject = sum(ClustersKind(IKind)%m_Atoms(:)%m_NA)
                            AtomuNumbObject = sum(ClustersKind(JKind)%m_Atoms(:)%m_NA)
                            AtomuNumbProductor = sum(ClustersKind(IKind+JKind)%m_Atoms(:)%m_NA)

                            tempNBPVChangeRate(IKind + JKind,INode) = tempNBPVChangeRate(IKind + JKind,INode) + &
                                                                      Factor*deta*(AtomuNumbSubject+AtomuNumbObject)/AtomuNumbProductor
                        end if

                    END DO

                END DO
            END DO

            TSTEP = maxval(Concent)*Host_SimuCtrlParam%MaxChangeRate/maxval(dabs(tempNBPVChangeRate))

            TSTEP = min(TSTEP,1D-8)

            tempNBPV = Concent

            tempNBPV = tempNBPV + tempNBPVChangeRate*TSTEP

            DO IKind = 1,CKind
                DO INode = 1,NNodes
                    if(tempNBPV(IKind,INode) .LT. 0) then
                        tempTimeStep = Concent(IKind,INode)/dabs(tempNBPVChangeRate(IKind,INode))
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

            DO IKind = 1,CKind

                MatA = 0.D0
                MatB = 0.D0
                MatC = 0.D0
                MatD = 0.D0

                DO INode = 1,NNodes
                    if(INode .eq. 1) then  ! upper surface
                        !For surface, the Dirichlet boundary condition is applied
                        DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                        MatA(INode) = 0.D0
                        MatB(INode) = (NodeSpace(INode)/TSTEP - (DiffGradient1 + DiffGradient2))*Concent(IKind,INode)
                        MatC(INode) = DiffGradient2*Concent(IKind,INode+1)
                        MatD(INode) = tempNBPVChangeRate(IKind,INode)*TSTEP
                    else if(INode .eq. NNodes) then  ! Low surface
                        !For surface, the Dirichlet boundary condition is applied
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)
                        MatA(INode) = DiffGradient1*Concent(IKind,INode-1)
                        MatB(INode) = (NodeSpace(INode)/TSTEP - (DiffGradient1 + DiffGradient2))*Concent(IKind,INode)
                        MatC(INode) = 0.D0
                        MatD(INode) = tempNBPVChangeRate(IKind,INode)*TSTEP
                    else
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))
                        MatA(INode) = DiffGradient1*Concent(IKind,INode-1)
                        MatB(INode) = (NodeSpace(INode)/TSTEP - (DiffGradient1 + DiffGradient2))*Concent(IKind,INode)
                        MatC(INode) = DiffGradient2*Concent(IKind,INode+1)
                        MatD(INode) = tempNBPVChangeRate(IKind,INode)*TSTEP
                    end if

!                    if(IKind .eq. 1) then
!                        DO IImplantLayer = 1,m_ImplantLayerNum
!                            if(sum(m_ImplantSpaceDist(1:IImplantLayer)) .GT. (INode-1)*m_NodeSpace .AND.  &
!                                sum(m_ImplantSpaceDist(1:IImplantLayer)) .LE. INode*m_NodeSpace) then
!                                    MatD(INode) = MatD(INode) + m_NodeSpace*TSTEP*m_Flux*m_ImplantDist(IImplantLayer)
!                                    NewAddedAtoms(IImplantLayer) = TSTEP*m_Flux*m_ImplantDist(IImplantLayer)
!                            end if
!                        END DO
!                    end if

                END DO

!                call SolveTridag(IKind,MatA,MatB,MatC,MatD,Concent,NNodes,MatW,MatH)

                DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(NNodes)
                FOutEachStep(IKind) = DiffGradient2*Concent(IKind,NNodes)
                COutEachStep(IKind) = DiffGradient2*Concent(IKind,NNodes)*TSTEP/NodeSpace(NNodes)
                FOutAccum(IKind) = FOutAccum(IKind) + FOutEachStep(IKind)
                COutAccum(IKind) = COutAccum(IKind) + COutEachStep(IKind)


                DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(1)
                FSurfEachStep(IKind) = DiffGradient1*Concent(IKind,1)
                CSurfEachStep(IKind) = DiffGradient1*Concent(IKind,1)*TSTEP/NodeSpace(1)
                FSurfAccum(IKind) = FSurfAccum(IKind) + FSurfEachStep(IKind)
                CSurfAccum(IKind) = CSurfAccum(IKind) + CSurfEachStep(IKind)


                DO INode = 1,NNodes
                    Concent(IKind,INode) = (MatA(INode) + MatB(INode) + MatC(INode) + MatD(INode))*TSTEP/NodeSpace(INode)
                END DO


                if(IKind .eq. 1) then

                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out flux from up surface",FSurfAccum(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out concentrate from up surface",CSurfAccum(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Out flux from up surface in current step",FSurfEachStep(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Out concentrate from up surface in current step",CSurfEachStep(IKind)
                    DO INode = 1,NNodes
                        write(*,fmt="(A10,1x,I10,1x,A15,1x,1ES14.6)") "INode",INode,"Concent(1,INode)",Concent(1,INode)
                    END DO
                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out flux from lower surface",FOutAccum(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out concentrate from lower surface",COutAccum(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Out flux from lower surface in current step",FOutEachStep(IKind)
                    write(*,fmt="(A50,1x,1ES14.6)") "Out concentrate from lower surface in current step",COutEachStep(IKind)

                    write(*,*) "----------------------------------------"
                    write(*,*) "sum(Concent(1,:))",sum(Concent(1,:))
                    write(*,*) "sum(Concent(1,:)) + CSurfAccum(1) + COutAccum(1)",sum(Concent(1,:)) + CSurfAccum(1) + COutAccum(1)
                    write(*,*) "(sum(Concent(1,:)) + CSurfAccum(1) + COutAccum(1))/ConCentrat0",(ConCentrat0 - (sum(Concent(1,:)) + CSurfAccum(1) + COutAccum(1)))/ConCentrat0
                    write(*,*) "----------------------------------------"
                end if


            END DO

            TTIME = TTIME + TSTEP
            ITIME = ITIME + 1

            ImplantedNum = ImplantedNum + NewAddedAtoms

            if(mod(ITIME,1) .eq. 0) then

                call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

                call Put_Out_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,ITIME,TTIME,TSTEP,sum(ImplantedNum),NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
            end if

            !if(Concent(CKind) .GT. 1.D-10) then
            if(DSQRT(dble(sum(ClustersKind(CKind)%m_Atoms(:)%m_NA)))*sum(Concent(CKind,1:NNodes)) .GT. &
               sum(NPOWER1DIV2Ave)*Host_SimuCtrlParam%DumplicateFactor) then

                DO INode = 1,NNodes

                    DO I = 1,(CKind -1)/2 + 1
                        if(2*I .LE. CKind) then
                            Concent(I,INode) = (Concent(2*I-1,INode)*sum(ClustersKind(2*I-1)%m_Atoms(:)%m_NA)+Concent(2*I,INode)*sum(ClustersKind(2*I)%m_Atoms(:)%m_NA)) &
                                               /(2.D0*sum(ClustersKind(2*I)%m_Atoms(:)%m_NA))
                        else
                            Concent(I,INode) = Concent(2*I - 1,INode)
                        end if
                    END DO

                    Concent((CKind -1)/2+2:CKind,INode) = 0.D0

                    DO I = 1,(CKind -1)/2 + 1
                        if(2*I .LE. CKind) then
                            ClustersKind(I) = ClustersKind(2*I)
                        else
                            ClustersKind(I) = ClustersKind(2*I - 1)
                        end if
                    END DO

                    DO I = (CKind -1)/2 + 2,CKind
                        ClustersKind(I)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA = 2*ClustersKind(I)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA
                        TheDiffusorValue = Host_SimBoxes%m_DiffusorTypesMap%get(ClustersKind(I))

                        select case(TheDiffusorValue%ECRValueType_Free)
                            case(p_ECR_ByValue)
                                ClustersKind(I)%m_RAD = TheDiffusorValue%ECR_Free
                            case(p_ECR_ByBCluster)
                                ClustersKind(I)%m_RAD = DSQRT(sum(ClustersKind(I)%m_Atoms(:)%m_NA)/m_RNFACTOR)
                        end select

                        select case(TheDiffusorValue%DiffusorValueType_Free)
                            case(p_DiffuseCoefficient_ByValue)
                                ClustersKind(I)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                            case(p_DiffuseCoefficient_ByArrhenius)
                                ClustersKind(I)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
                            case(p_DiffuseCoefficient_ByBCluster)
                                ! Here we adopt a model that D=D0*(1/R)**Gama
                                ClustersKind(I)%m_DiffCoeff = m_FREESURDIFPRE*(ClustersKind(I)%m_RAD**(-p_GAMMA))
                        end select
                    END DO

                END DO

                Dumplicate = Dumplicate*2
                write(*,*) "Dumplicate",Dumplicate
            end if


            if(TTIME .GT. Host_SimuCtrlParam%TermTValue) then
                exit
            end if

          END DO
        END Associate

    end subroutine NucleationSimu_SpaceDist

    !----------------------------------------------------------------------
    subroutine SolveTridag(IKind,MatrixA,MatrixB,MatrixC,MatrixD,Solver,MatrixSize,w,h)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::IKind
        real(kind=KMCDF),dimension(:),allocatable::MatrixA
        real(kind=KMCDF),dimension(:),allocatable::MatrixB
        real(kind=KMCDF),dimension(:),allocatable::MatrixC
        real(kind=KMCDF),dimension(:),allocatable::MatrixD
        real(kind=KMCDF),dimension(:,:),allocatable::Solver
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
            h(I) = (MatrixD(I) - MatrixA(I)*h(I-1))/w(I)
        END DO
        Solver(IKind,MatrixSize) = h(MatrixSize)

        DO I = MatrixSize-1,1,-1
            Solver(IKind,I) = h(I) - MatrixC(I)*Solver(IKind,I+1)/w(I)
        END DO


        return
    end subroutine SolveTridag

    !---------------------------------------------------------------------
    subroutine Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        !---Local Vars---
        integer::CKind
        real(kind=KMCDF)::Temp
        integer::I
        integer::NNodes
        integer::INode
        !---Body---
        NNodes = Host_SimBoxes%NNodes

        CKind = Host_SimBoxes%CKind

        Associate(ClustersKind=>Host_SimBoxes%ClustersKindArray,Concent=>Host_SimBoxes%Concentrate)

          DO INode = 1,NNodes
            NPOWER0Ave(INode) = Dumplicate*sum(Concent(1:CKind,INode))

            Temp = 0.D0
            DO I = 1,CKind
                Temp = Temp + (sum(ClustersKind(I)%m_Atoms(:)%m_NA)**(0.5D0))*Concent(I,INode)
            END DO
            NPOWER1DIV2Ave(INode) = Dumplicate*Temp

            Temp = 0.D0
            DO I = 1,CKind
                Temp = Temp + sum(ClustersKind(I)%m_Atoms(:)%m_NA)*Concent(I,INode)
            END DO
            NPOWER1Ave(INode) = Dumplicate*Temp

            Temp = 0.D0
            DO I = 1,CKind
                Temp = Temp + (sum(ClustersKind(I)%m_Atoms(:)%m_NA)**(1.5D0))*Concent(I,INode)
            END DO
            NPOWER3DIV2Ave(INode) = Dumplicate*Temp
          END DO

        END Associate
        return
    end subroutine

    !---------------------------------------------
    subroutine Put_Out_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,Step,TTime,TStep,ImplantedNum,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
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
        integer::CKind
        character*256::C_Index
        character*256::fileName
        integer::I
        integer::NNodes
        integer::INode
        !---Body---



        NNodes = Host_SimBoxes%NNodes

        CKind = Host_SimBoxes%CKind

        N1 = 0.D0
        N2 = 0.D0
        N3 = 0.D0
        Rave = 0.D0

        Associate(ClustersKind=>Host_SimBoxes%ClustersKindArray,Concent=>Host_SimBoxes%Concentrate)

          DO INode = 1,NNodes

            if(NPOWER0Ave(INode) .GT. 0) then
                N1(INode) = (NPOWER1DIV2Ave((INode))/NPOWER0Ave(INode))**2

                N2(INode) = NPOWER1Ave(INode)/NPOWER0Ave(INode)

                N3(INode) = (NPOWER3DIV2Ave(INode)/NPOWER0Ave(INode))**(dble(2)/dble(3))

                Rave(INode) = DSQRT(N1(INode)/m_RNFACTOR)
            end if
          END DO


          write(6,FMT="(15(A15,1x))") "Step","Time","TStep","ImplantedNum","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"

          write(6,FMT="(I15,1x,15(E15.5,1x))") Step,TTime,TStep,ImplantedNum,sum(NPOWER0Ave)/NNodes,sum(NPOWER1DIV2Ave)/NNodes,sum(NPOWER1Ave)/NNodes,sum(NPOWER3DIV2Ave)/NNodes, &
                                             sum(N1)/NNodes,sum(N2)/NNodes,sum(N3)/NNodes,sum(Rave)*C_CM2NM/NNodes

          write(m_StatisticFile,FMT="(I15,1x,15(E15.5,1x))") Step,TTime,TStep,ImplantedNum,sum(NPOWER0Ave)/NNodes,sum(NPOWER1DIV2Ave)/NNodes,sum(NPOWER1Ave)/NNodes,sum(NPOWER3DIV2Ave)/NNodes, &
                                             sum(N1)/NNodes,sum(N2)/NNodes,sum(N3)/NNodes,sum(Rave)*C_CM2NM/NNodes

          Flush(m_StatisticFile)


          LastOutPutConifgIndex = LastOutPutConifgIndex + 1


          call AvailableIOUnit(hFile)

          write(C_Index,*) LastOutPutConifgIndex

          C_Index = adjustl(C_Index)

          fileName = Host_SimuCtrlParam%OutFilePath(1:len_trim(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"Config"//C_Index(1:len_trim(C_Index))//".dat"

          open(Unit=hFile,file=trim(fileName))

          write(hFile,*) "Time",TTime

          write(hFile,FMT="(4(A15))") "n","NBPV","u","Z(u)"

          DO I = 1,CKind
              write(hFile,FMT="(15(E15.5))") sum(ClustersKind(I)%m_Atoms(:)%m_NA),                              &
                                             sum(Concent(I,1:NNodes))/NNodes,                                   &
                                             sum(ClustersKind(I)%m_Atoms(:)%m_NA)*TTime**(-dble(2)/dble(5)),    &
                                             (sum(Concent(I,1:NNodes))/NNodes)/(TTime**(-dble(4)/dble(5)))
          END DO

          close(hFile)
         END associate

        return
    end subroutine

end module NUCLEATION_SPACEDIST
