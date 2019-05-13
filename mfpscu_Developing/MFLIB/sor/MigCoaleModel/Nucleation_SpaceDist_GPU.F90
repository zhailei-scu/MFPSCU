module NUCLEATION_SPACEDIST_GPU
    use cudafor
    use MFLIB_GLOBAL
    use MCMF_CONSTANTS_GPU
    use MIGCOALE_IMPLANTATION
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

    integer::Dumplicate = 1

    integer::LastOutPutConifgIndex = 0

    integer::m_StatisticFile

    !--------
    type(ACluster),device,dimension(:),allocatable::dm_ClustersKindArray
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_Concentrate
    real(kind=KMCDF),device,dimension(:),allocatable::dm_NodeSpace

    contains

    subroutine InitSimu_SpaceDist_GPU(Host_Boxes,Host_SimuCtrlParam)
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

        if(allocated(dm_ClustersKindArray)) then
            deallocate(dm_ClustersKindArray)
        end if
        allocate(dm_ClustersKindArray(CKind))
        dm_ClustersKindArray = Host_Boxes%m_ClustersInfo_CPU%ClustersKindArray

        if(allocated(dm_Concentrate)) then
            deallocate(dm_Concentrate)
        end if
        allocate(dm_Concentrate(NNodes,CKind))
        dm_Concentrate = Host_Boxes%m_ClustersInfo_CPU%Concentrate

        if(allocated(dm_NodeSpace)) then
            deallocate(dm_NodeSpace)
        end if
        allocate(dm_NodeSpace(NNodes))
        dm_NodeSpace = Host_Boxes%NodeSpace

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


        call AvailableIOUnit(m_StatisticFile)

        fileName = Host_SimuCtrlParam%OutFilePath(1:len_trim(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"Statistic.dat"

        open(Unit=m_StatisticFile,file=fileName(1:len_trim(fileName)))

        write(m_StatisticFile,FMT="(15(A15,1x))") "Step","Time","TStep","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"

        return
    end subroutine InitSimu_SpaceDist_GPU

    !***************************************************
    subroutine NucleationSimu_SpaceDist_Balance_GrubDumplicate_GPU(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,TheImplantSection)
        use RAND32_MODULE
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        type(MigCoaleStatInfoWrap)::TheMigCoaleStatInfoWrap
        type(MigCoalClusterRecord)::Record
        type(ImplantSection)::TheImplantSection
        !---Local Vars---
        integer::CKind
        integer::IKind
        integer::JKind
        real(kind=KMCDF)::TSTEP
        real(kind=KMCDF)::deta
        real(kind=KMCDF),dimension(:,:),allocatable::tempNBPVChangeRate
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::N1
        real(kind=KMCDF),dimension(:),allocatable::N2
        real(kind=KMCDF),dimension(:),allocatable::N3
        real(kind=KMCDF),dimension(:),allocatable::Rave
        real(kind=KMCDF),dimension(:,:),allocatable::ImplantedRate
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
        real(kind=KMCDF)::DiffGradient1
        real(kind=KMCDF)::DiffGradient2
        integer::IImplantLayer
        type(DiffusorValue)::TheDiffusorValue
        integer::AtomuNumbSubject
        integer::AtomuNumbObject
        integer::AtomuNumbProductor
        real(kind=KMCDF)::ConCentrat0
        type(ReactionValue)::TheReactionValue
        real(kind=KMCDF)::ReactionCoeff
        real(kind=KMCDF)::SFlux
        !---Body---
        CKind = Host_SimBoxes%CKind
        NNodes = Host_SimBoxes%NNodes

        allocate(tempNBPVChangeRate(NNodes,CKind))

        allocate(NPOWER0Ave(NNodes),NPOWER1DIV2Ave(NNodes),NPOWER1Ave(NNodes),NPOWER3DIV2Ave(NNodes))

        allocate(N1(NNodes),N2(NNodes),N3(NNodes),Rave(NNodes))

        allocate(ImplantedRate(NNodes,CKind))

        allocate(FSurfAccum(NNodes))

        allocate(FOutAccum(NNodes))

        allocate(CSurfAccum(NNodes))

        allocate(COutAccum(NNodes))

        allocate(FSurfEachStep(NNodes))

        allocate(FOutEachStep(NNodes))

        allocate(CSurfEachStep(NNodes))

        allocate(COutEachStep(NNodes))

        TSTEP = 0.01

        FSurfAccum = 0.D0

        FOutAccum = 0.D0

        CSurfAccum = 0.D0

        COutAccum = 0.D0

        startAnnealing = .false.

        call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

        Associate(ClustersKind=>Host_SimBoxes%m_ClustersInfo_CPU%ClustersKindArray,Concent=>Host_SimBoxes%m_ClustersInfo_CPU%Concentrate,NodeSpace=>Host_SimBoxes%NodeSpace)

          ConCentrat0 = sum(Concent)

          DO While(.true.)

            call Record%IncreaseOneSimuStep()

            tempNBPVChangeRate = 0.D0

            FSurfEachStep = 0.D0

            FOutEachStep = 0.D0

            CSurfEachStep = 0.D0

            COutEachStep = 0.D0

            ImplantedRate = 0.D0

            call TheImplantSection%Cal_ImplantClustersRate(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,ImplantedRate)

            DO INode = 1,NNodes

                DO IKind = 1,CKind

                    DO JKind = IKind,CKind

                        TheReactionValue = Host_SimBoxes%m_ReactionsMap%get(ClustersKind(IKind),ClustersKind(JKind))

                        ReactionCoeff = 0.D0
                        select case(TheReactionValue%ReactionCoefficientType)
                            case(p_ReactionCoefficient_ByValue)
                                ReactionCoeff = TheReactionValue%ReactionCoefficient_Value
                            case(p_ReactionCoefficient_ByArrhenius)
                                ReactionCoeff = TheReactionValue%PreFactor*exp(-C_EV2ERG*TheReactionValue%ActEnergy/Host_SimuCtrlParam%TKB)
                        end select

                        if(ReactionCoeff .GE. DRAND32()) then

                            deta = Dumplicate*4*PI*Concent(INode,IKind)*Concent(INode,JKind)* &
                                    (ClustersKind(IKind)%m_RAD + ClustersKind(JKind)%m_RAD)* &
                                    (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(JKind)%m_DiffCoeff)

                            if(IKind .eq. JKind) then

                                Factor = 0.5D0

                                tempNBPVChangeRate(INode,IKind) =  tempNBPVChangeRate(INode,IKind) - deta
                            else
                                Factor = 1.D0

                                tempNBPVChangeRate(INode,IKind) =  tempNBPVChangeRate(INode,IKind) - deta

                                tempNBPVChangeRate(INode,JKind) =  tempNBPVChangeRate(INode,JKind) - deta
                            end if

                            if((IKind + JKind) .LE. CKind) then

                                AtomuNumbSubject = sum(ClustersKind(IKind)%m_Atoms(:)%m_NA)
                                AtomuNumbObject = sum(ClustersKind(JKind)%m_Atoms(:)%m_NA)
                                AtomuNumbProductor = sum(ClustersKind(IKind+JKind)%m_Atoms(:)%m_NA)

                                tempNBPVChangeRate(INode,IKind + JKind) = tempNBPVChangeRate(INode,IKind + JKind) + &
                                                                            Factor*deta*(AtomuNumbSubject+AtomuNumbObject)/AtomuNumbProductor
                            end if

                        end if

                    END DO

                END DO
            END DO

            TSTEP = maxval(Concent)*Host_SimuCtrlParam%MaxReactChangeRate/maxval(dabs(tempNBPVChangeRate))

            DO IKind = 1,CKind
                DO INode = 1,NNodes
                    tempTimeStep = Concent(INode,IKind)/dabs(tempNBPVChangeRate(INode,IKind))
                    if(tempNBPVChangeRate(INode,IKind) .LT. 0.D0 .AND. Concent(INode,IKind) .GT. 0.D0) then
                        TSTEP = min(TSTEP,tempTimeStep)
                    end if

                    if(INode .eq. 1) then  ! upper surface

                        DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        if(NNodes .LE. 1) then
                            DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            SFlux = DiffGradient1*Concent(INode,IKind) + DiffGradient2*Concent(INode,IKind)
                            if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                                TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                            end if
                        else
                            DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                            SFlux = DiffGradient1*Concent(INode,IKind) - DiffGradient2*(Concent(INode+1,IKind) - Concent(INode,IKind))
                            if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                                TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                            end if

                        end if

                    else if(INode .eq. NNodes) then  ! Low surface

                        DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        if(NNodes .LE. 1) then
                            DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            SFlux = DiffGradient1*Concent(INode,IKind) + DiffGradient2*Concent(INode,IKind)
                            if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                                TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                            end if
                        else
                            DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))

                            SFlux = DiffGradient1*(Concent(INode,IKind) - Concent(INode-1,IKind)) + DiffGradient2*Concent(INode,IKind)
                            if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                                TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                            end if
                        end if

                    else
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                        SFlux = DiffGradient1*(Concent(INode,IKind) - Concent(INode-1,IKind)) - DiffGradient2*(Concent(INode+1,IKind) - Concent(INode,IKind))
                        if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                            TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                        end if
                    end if

                END DO
            END DO

            DO IKind = 1,CKind

                MatA = 0.D0
                MatB = 0.D0
                MatC = 0.D0
                MatD = 0.D0

                DO INode = 1,NNodes
                    if(INode .eq. 1) then  ! upper surface

                        DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        MatA(INode) = 0.D0
                        if(NNodes .LE. 1) then
                            DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            select case(Host_SimuCtrlParam%BDCTYPE(3,1))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                                case(p_Neumann_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,1)
                                    pause
                                    stop
                            end select

                            MatC(INode) = 0.D0
                        else
                            DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                            select case(Host_SimuCtrlParam%BDCTYPE(3,1))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                                case(p_Neumann_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient2
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,1)
                                    pause
                                    stop
                            end select

                            MatC(INode) = -DiffGradient2
                        end if

                        MatD(INode) = Concent(INode,IKind)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    else if(INode .eq. NNodes) then  ! Low surface

                        DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        if(NNodes .LE. 1) then
                            DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            MatA(INode) = 0.D0

                            select case(Host_SimuCtrlParam%BDCTYPE(3,2))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                                case(p_Neumann_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,2)
                                    pause
                                    stop
                            end select
                        else
                            DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))

                            MatA(INode) = -DiffGradient1

                            select case(Host_SimuCtrlParam%BDCTYPE(3,2))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                                case(p_Neumann_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,2)
                                    pause
                                    stop
                            end select

                        end if

                        MatC(INode) = 0.D0
                        MatD(INode) = Concent(INode,IKind)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    else
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))
                        MatA(INode) = -DiffGradient1
                        MatB(INode) = NodeSpace(INode)/TSTEP + (DiffGradient1 + DiffGradient2)
                        MatC(INode) = -DiffGradient2
                        MatD(INode) = Concent(INode,IKind)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    end if

                END DO

                call SolveTridag(IKind,MatA,MatB,MatC,MatD,Concent,NNodes,MatW,MatH)

!                DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(NNodes)
!                FOutEachStep(IKind) = DiffGradient2*Concent(NNodes,IKind)
!                COutEachStep(IKind) = DiffGradient2*Concent(NNodes,IKind)*TSTEP/NodeSpace(NNodes)
!                FOutAccum(IKind) = FOutAccum(IKind) + FOutEachStep(IKind)
!                COutAccum(IKind) = COutAccum(IKind) + COutEachStep(IKind)
!
!                DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(1)
!                FSurfEachStep(IKind) = DiffGradient1*Concent(1,IKind)
!                CSurfEachStep(IKind) = DiffGradient1*Concent(1,IKind)*TSTEP/NodeSpace(1)
!                FSurfAccum(IKind) = FSurfAccum(IKind) + FSurfEachStep(IKind)
!                CSurfAccum(IKind) = CSurfAccum(IKind) + CSurfEachStep(IKind)

                if(IKind .eq. 1) then

                    write(*,*) "Accumulated out flux from up surface",FSurfAccum(IKind)
                    write(*,*) "Accumulated out concentrate from up surface",CSurfAccum(IKind)
                    write(*,*) "Out flux from up surface in current step",FSurfEachStep(IKind)
                    write(*,*) "Out concentrate from up surface in current step",CSurfEachStep(IKind)
                    DO INode = 1,NNodes
                        write(*,*) "INode",INode,"Concent(INode,1)",Concent(INode,1)
                    END DO
                    write(*,*) "Accumulated out flux from lower surface",FOutAccum(IKind)
                    write(*,*) "Accumulated out concentrate from lower surface",COutAccum(IKind)
                    write(*,*) "Out flux from lower surface in current step",FOutEachStep(IKind)
                    write(*,*) "Out concentrate from lower surface in current step",COutEachStep(IKind)

                    write(*,*) "----------------------------------------"
                    write(*,*) "sum(Concent(:,1))",sum(Concent(:,1))
                    write(*,*) "sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)",sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)
                    write(*,*) "(sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1))/ConCentrat0",(ConCentrat0 - (sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)))/ConCentrat0
                    write(*,*) "----------------------------------------"
                end if


            END DO

            call Record%AddSimuTimes(TSTEP)

            call OutPutCurrent(Host_SimBoxes,Host_SimuCtrlParam,Record)

            if(mod(Record%GetSimuSteps(),1) .eq. 0) then

                call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

                call Put_Out_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,Record%GetSimuSteps(),Record%GetSimuTimes(),TSTEP,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
            end if

            !if(Concent(CKind) .GT. 1.D-10) then
            if(DSQRT(dble(sum(ClustersKind(CKind)%m_Atoms(:)%m_NA)))*sum(Concent(1:NNodes,CKind)) .GT. &
               sum(NPOWER1DIV2Ave)*Host_SimuCtrlParam%DumplicateFactor) then

                DO INode = 1,NNodes

                    DO I = 1,(CKind -1)/2 + 1
                        if(2*I .LE. CKind) then
                            Concent(INode,I) = (Concent(INode,2*I-1)*sum(ClustersKind(2*I-1)%m_Atoms(:)%m_NA)+Concent(INode,2*I)*sum(ClustersKind(2*I)%m_Atoms(:)%m_NA)) &
                                               /(2.D0*sum(ClustersKind(2*I)%m_Atoms(:)%m_NA))
                        else
                            Concent(INode,I) = Concent(INode,2*I - 1)
                        end if
                    END DO

                    Concent(INode,(CKind -1)/2+2:CKind) = 0.D0

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

            if(Record%GetSimuTimes() .GT. Host_SimuCtrlParam%TermTValue) then
                exit
            end if

            TSTEP = TSTEP*1.5

          END DO
        END Associate

    end subroutine NucleationSimu_SpaceDist_Balance_GrubDumplicate_GPU

    !***************************************************
    attributes(global) subroutine Kernel_CalReaction_Vanish(NNode,CKind,Dev_Concent,Dev_ClusterKindArray,Dev_tempNBPVChangeRate)
        implicit none
        !---Dummy Vars---
        integer,value::NNode
        integer,value::CKind
        real(kind=KMCDF),device::Dev_Concent(:,:)
        type(ACluster),device::Dev_ClusterKindArray(:)
        real(kind=KMCDF),device::Dev_tempNBPVChangeRate(:,:)
        !---Local Vars---
        integer::tid
        integer::bid
        integer::cid
        integer::IKind
        integer::JKInd
        integer::INode
        real(kind=KMCDF)::deta
        real(kind=KMCDF)::Factor
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x +threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*p_BlockSize + tid

        IKind = (cid - 1)/NNode + 1
        INode = cid - (IKind-1)*NNode

        DO JKind = 1,CKind

!            TheReactionValue = Host_SimBoxes%m_ReactionsMap%get(ClustersKind(IKind),ClustersKind(JKind))
!
!            ReactionCoeff = 0.D0
!            select case(TheReactionValue%ReactionCoefficientType)
!                case(p_ReactionCoefficient_ByValue)
!                    ReactionCoeff = TheReactionValue%ReactionCoefficient_Value
!                case(p_ReactionCoefficient_ByArrhenius)
!                    ReactionCoeff = TheReactionValue%PreFactor*exp(-C_EV2ERG*TheReactionValue%ActEnergy/Host_SimuCtrlParam%TKB)
!            end select

!            if(ReactionCoeff .GE. DRAND32()) then


                deta =  4*PI*Dev_Concent(INode,IKind)*Dev_Concent(INode,JKind)* &
                                    (Dev_ClusterKindArray(IKind)%m_RAD + Dev_ClusterKindArray(JKind)%m_RAD)* &
                                    (Dev_ClusterKindArray(IKind)%m_DiffCoeff + Dev_ClusterKindArray(JKind)%m_DiffCoeff)

                Dev_tempNBPVChangeRate(INode,IKind) =  Dev_tempNBPVChangeRate(INode,IKind) - deta

!            end if

        END DO

        return
    end subroutine Kernel_CalReaction_Vanish

    !***************************************************
    attributes(global) subroutine Kernel_CalReaction_Generate(NNode,CKind,Dev_Concent,Dev_ClusterKindArray,Dev_tempNBPVChangeRate)
        implicit none
        !---Dummy Vars---
        integer,value::NNode
        integer,value::CKind
        real(kind=KMCDF),device::Dev_Concent(:,:)
        type(ACluster),device::Dev_ClusterKindArray(:)
        real(kind=KMCDF),device::Dev_tempNBPVChangeRate(:,:)
        !---Local Vars---
        integer::tid
        integer::bid
        integer::cid
        integer::IKind
        integer::JKInd
        integer::LKind
        integer::INode
        real(kind=KMCDF)::deta
        real(kind=KMCDF)::Factor
        integer::AtomuNumbSubject
        integer::AtomuNumbObject
        integer::AtomuNumbProductor
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x +threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*p_BlockSize + tid

        LKind = (cid - 1)/NNode + 1
        INode = cid - (LKind-1)*NNode

        DO Ikind = 1,LKind - 1

            JKind = LKind - IKind

!            TheReactionValue = Host_SimBoxes%m_ReactionsMap%get(ClustersKind(IKind),ClustersKind(JKind))
!
!            ReactionCoeff = 0.D0
!            select case(TheReactionValue%ReactionCoefficientType)
!                case(p_ReactionCoefficient_ByValue)
!                    ReactionCoeff = TheReactionValue%ReactionCoefficient_Value
!                case(p_ReactionCoefficient_ByArrhenius)
!                    ReactionCoeff = TheReactionValue%PreFactor*exp(-C_EV2ERG*TheReactionValue%ActEnergy/Host_SimuCtrlParam%TKB)
!            end select

!            if(ReactionCoeff .GE. DRAND32()) then

                Factor = 1.D0
                if(IKind .eq. JKInd) then
                    Factor = 0.5D0
                end if

                deta =  4*PI*Dev_Concent(INode,IKind)*Dev_Concent(INode,JKind)* &
                                    (Dev_ClusterKindArray(IKind)%m_RAD + Dev_ClusterKindArray(JKind)%m_RAD)* &
                                    (Dev_ClusterKindArray(IKind)%m_DiffCoeff + Dev_ClusterKindArray(JKind)%m_DiffCoeff)

                AtomuNumbSubject = sum(Dev_ClusterKindArray(IKind)%m_Atoms(:)%m_NA,dim=1)
                AtomuNumbObject = sum(Dev_ClusterKindArray(JKind)%m_Atoms(:)%m_NA,dim=1)
                AtomuNumbProductor = sum(Dev_ClusterKindArray(LKind)%m_Atoms(:)%m_NA,dim=1)

                Dev_tempNBPVChangeRate(INode,LKind) =  Dev_tempNBPVChangeRate(INode,LKind) + &
                                                       Factor*deta*(AtomuNumbSubject+AtomuNumbObject)/AtomuNumbProductor

!            end if

        END DO

        return
    end subroutine Kernel_CalReaction_Generate

    !***************************************************
    attributes(global) subroutine Kernel_DoReaction_And_NodeDiffusion_Balance(NNode,CKind,Dev_Concent,Dev_ClusterKindArray,Dev_tempNBPVChangeRate)
        implicit none
        !---Dummy Vars---
        integer,value::NNode
        integer,value::CKind
        real(kind=KMCDF),device::Dev_Concent(:,:)
        type(ACluster),device::Dev_ClusterKindArray(:)
        real(kind=KMCDF),device::Dev_tempNBPVChangeRate(:,:)
        !---Local Vars---
        integer::tid
        integer::bid
        integer::cid
        integer::IKind
        integer::JKInd
        integer::LKind
        integer::INode
        real(kind=KMCDF)::deta
        real(kind=KMCDF)::Factor
        integer::AtomuNumbSubject
        integer::AtomuNumbObject
        integer::AtomuNumbProductor
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x +threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*p_BlockSize + tid

        LKind = (cid - 1)/NNode + 1
        INode = cid - (LKind-1)*NNode

        MatA = 0.D0
        MatB = 0.D0
        MatC = 0.D0
        MatD = 0.D0

        DO INode = 1,NNodes
                    if(INode .eq. 1) then  ! upper surface

                        DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        MatA(INode) = 0.D0
                        if(NNodes .LE. 1) then
                            DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            select case(Host_SimuCtrlParam%BDCTYPE(3,1))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                                case(p_Neumann_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,1)
                                    pause
                                    stop
                            end select

                            MatC(INode) = 0.D0
                        else
                            DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                            select case(Host_SimuCtrlParam%BDCTYPE(3,1))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                                case(p_Neumann_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient2
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,1)
                                    pause
                                    stop
                            end select

                            MatC(INode) = -DiffGradient2
                        end if

                        MatD(INode) = Concent(INode,IKind)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    else if(INode .eq. NNodes) then  ! Low surface

                        DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        if(NNodes .LE. 1) then
                            DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            MatA(INode) = 0.D0

                            select case(Host_SimuCtrlParam%BDCTYPE(3,2))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                                case(p_Neumann_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,2)
                                    pause
                                    stop
                            end select
                        else
                            DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))

                            MatA(INode) = -DiffGradient1

                            select case(Host_SimuCtrlParam%BDCTYPE(3,2))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                                case(p_Neumann_BDC)
                                    MatB(INode) = NodeSpace(INode)/TSTEP + DiffGradient1
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,2)
                                    pause
                                    stop
                            end select

                        end if

                        MatC(INode) = 0.D0
                        MatD(INode) = Concent(INode,IKind)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    else
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))
                        MatA(INode) = -DiffGradient1
                        MatB(INode) = NodeSpace(INode)/TSTEP + (DiffGradient1 + DiffGradient2)
                        MatC(INode) = -DiffGradient2
                        MatD(INode) = Concent(INode,IKind)*NodeSpace(INode)/TSTEP + tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    end if

                END DO

        call SolveTridag(IKind,MatA,MatB,MatC,MatD,Concent,NNodes,MatW,MatH)


        return
    end subroutine Kernel_DoReaction_And_NodeDiffusion_Balance


    !***************************************************
    subroutine NucleationSimu_SpaceDist_Transient_GrubDumplicate_GPU(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,TheImplantSection)
        use RAND32_MODULE
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        type(MigCoaleStatInfoWrap)::TheMigCoaleStatInfoWrap
        type(MigCoalClusterRecord)::Record
        type(ImplantSection)::TheImplantSection
        !---Local Vars---
        integer::CKind
        integer::IKind
        integer::JKind
        real(kind=KMCDF)::TSTEP
        real(kind=KMCDF)::deta
        real(kind=KMCDF),dimension(:,:),allocatable::tempNBPVChangeRate
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::N1
        real(kind=KMCDF),dimension(:),allocatable::N2
        real(kind=KMCDF),dimension(:),allocatable::N3
        real(kind=KMCDF),dimension(:),allocatable::Rave
        real(kind=KMCDF),dimension(:,:),allocatable::ImplantedRate
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
        real(kind=KMCDF)::DiffGradient1
        real(kind=KMCDF)::DiffGradient2
        integer::IImplantLayer
        integer::I
        type(DiffusorValue)::TheDiffusorValue
        integer::AtomuNumbSubject
        integer::AtomuNumbObject
        integer::AtomuNumbProductor
        real(kind=KMCDF)::ConCentrat0
        real(kind=KMCDF)::ReactionCoeff
        type(ReactionValue)::TheReactionValue
        real(kind=KMCDF)::SFlux
        !---Body---
        CKind = Host_SimBoxes%CKind

        NNodes = Host_SimBoxes%NNodes

        allocate(tempNBPVChangeRate(NNodes,CKind))

        allocate(ImplantedRate(NNodes,CKind))

        allocate(NPOWER0Ave(NNodes),NPOWER1DIV2Ave(NNodes),NPOWER1Ave(NNodes),NPOWER3DIV2Ave(NNodes))

        allocate(N1(NNodes),N2(NNodes),N3(NNodes),Rave(NNodes))

        allocate(FSurfAccum(NNodes))

        allocate(FOutAccum(NNodes))

        allocate(CSurfAccum(NNodes))

        allocate(COutAccum(NNodes))

        allocate(FSurfEachStep(NNodes))

        allocate(FOutEachStep(NNodes))

        allocate(CSurfEachStep(NNodes))

        allocate(COutEachStep(NNodes))

        TSTEP = 0.01

        FSurfAccum = 0.D0

        FOutAccum = 0.D0

        CSurfAccum = 0.D0

        COutAccum = 0.D0

        startAnnealing = .false.

        call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

        Associate(ClustersKind=>Host_SimBoxes%m_ClustersInfo_CPU%ClustersKindArray,Concent=>Host_SimBoxes%m_ClustersInfo_CPU%Concentrate,NodeSpace=>Host_SimBoxes%NodeSpace)

          ConCentrat0 = sum(Concent)

          DO While(.true.)

            call Record%IncreaseOneSimuStep()

            tempNBPVChangeRate = 0.D0

            FSurfEachStep = 0.D0

            FOutEachStep = 0.D0

            CSurfEachStep = 0.D0

            COutEachStep = 0.D0

            ImplantedRate = 0.D0

            call TheImplantSection%Cal_ImplantClustersRate(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,ImplantedRate)

            DO INode = 1,NNodes

                DO IKind = 1,CKind

                    DO JKind = IKind,CKind

                        TheReactionValue = Host_SimBoxes%m_ReactionsMap%get(ClustersKind(IKind),ClustersKind(JKind))

                        ReactionCoeff = 0.D0
                        select case(TheReactionValue%ReactionCoefficientType)
                            case(p_ReactionCoefficient_ByValue)
                                ReactionCoeff = TheReactionValue%ReactionCoefficient_Value
                            case(p_ReactionCoefficient_ByArrhenius)
                                ReactionCoeff = TheReactionValue%PreFactor*exp(-C_EV2ERG*TheReactionValue%ActEnergy/Host_SimuCtrlParam%TKB)
                        end select

                        if(ReactionCoeff .GE. DRAND32()) then


                            deta = Dumplicate*4*PI*Concent(INode,IKind)*Concent(INode,JKind)* &
                                    (ClustersKind(IKind)%m_RAD + ClustersKind(JKind)%m_RAD)* &
                                    (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(JKind)%m_DiffCoeff)

                            if(IKind .eq. JKind) then

                                Factor = 0.5D0

                                tempNBPVChangeRate(INode,IKind) =  tempNBPVChangeRate(INode,IKind) - deta

                            else

                                Factor = 1.D0

                                tempNBPVChangeRate(INode,IKind) =  tempNBPVChangeRate(INode,IKind) - deta

                                tempNBPVChangeRate(INode,JKind) =  tempNBPVChangeRate(INode,JKind) - deta

                            end if

                            if((IKind + JKind) .LE. CKind) then
                                AtomuNumbSubject = sum(ClustersKind(IKind)%m_Atoms(:)%m_NA)
                                AtomuNumbObject = sum(ClustersKind(JKind)%m_Atoms(:)%m_NA)
                                AtomuNumbProductor = sum(ClustersKind(IKind+JKind)%m_Atoms(:)%m_NA)

                                tempNBPVChangeRate(INode,IKind + JKind) = tempNBPVChangeRate(INode,IKind + JKind) + &
                                                                            Factor*deta*(AtomuNumbSubject+AtomuNumbObject)/AtomuNumbProductor
                            end if

                        end if

                    END DO

                END DO
            END DO

            TSTEP = maxval(Concent)*Host_SimuCtrlParam%MaxReactChangeRate/maxval(dabs(tempNBPVChangeRate))

            DO IKind = 1,CKind
                DO INode = 1,NNodes
                    tempTimeStep = Concent(INode,IKind)/dabs(tempNBPVChangeRate(INode,IKind))
                    if(tempNBPVChangeRate(INode,IKind) .LT. 0.D0 .AND. Concent(INode,IKind) .GT. 0.D0) then
                        TSTEP = min(TSTEP,tempTimeStep)
                    end if

                    if(INode .eq. 1) then  ! upper surface

                        DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        if(NNodes .LE. 1) then
                            DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            SFlux = DiffGradient1*Concent(INode,IKind) + DiffGradient2*Concent(INode,IKind)
                            if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                                TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                            end if
                        else
                            DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                            SFlux = DiffGradient1*Concent(INode,IKind) - DiffGradient2*(Concent(INode+1,IKind) - Concent(INode,IKind))
                            if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                                TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                            end if
                        end if

                    else if(INode .eq. NNodes) then  ! Low surface

                        DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        if(NNodes .LE. 1) then
                            DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            SFlux = DiffGradient1*Concent(INode,IKind) + DiffGradient2*Concent(INode,IKind)
                            if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                                TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                            end if
                        else
                            DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))

                            SFlux = DiffGradient1*(Concent(INode,IKind) - Concent(INode-1,IKind)) + DiffGradient2*Concent(INode,IKind)
                            if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                                TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                            end if

                        end if

                    else
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                        SFlux = DiffGradient1*(Concent(INode,IKind) - Concent(INode-1,IKind)) - DiffGradient2*(Concent(INode+1,IKind) - Concent(INode,IKind))
                        if(SFlux .GT. 0.D0 .AND. Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind) .GT. 0) then
                            TSTEP = min(TSTEP,Host_SimuCtrlParam%MaxDiffuseChangeRate*Concent(INode,IKind)/(dabs(SFlux)/NodeSpace(INode)))
                        end if
                    end if

                END DO
            END DO

            DO IKind = 1,CKind

                MatA = 0.D0
                MatB = 0.D0
                MatC = 0.D0
                MatD = 0.D0

                DO INode = 1,NNodes
                    if(INode .eq. 1) then  ! upper surface

                        DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        MatA(INode) = 0.D0
                        if(NNodes .LE. 1) then
                            DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            select case(Host_SimuCtrlParam%BDCTYPE(3,1))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = (NodeSpace(INode)/TSTEP - DiffGradient1 - DiffGradient2)*Concent(INode,IKind)
                                case(p_Neumann_BDC)
                                    MatB(INode) = (NodeSpace(INode)/TSTEP)*Concent(INode,IKind)
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,1)
                                    pause
                                    stop
                            end select

                            MatC(INode) = 0.D0
                        else
                            DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))

                            select case(Host_SimuCtrlParam%BDCTYPE(3,1))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = (NodeSpace(INode)/TSTEP - DiffGradient1 - DiffGradient2)*Concent(INode,IKind)
                                case(p_Neumann_BDC)
                                    MatB(INode) = (NodeSpace(INode)/TSTEP - DiffGradient2)*Concent(INode,IKind)
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,1)
                                    pause
                                    stop
                            end select

                            MatC(INode) = DiffGradient2*Concent(INode+1,IKind)
                        end if

                        MatD(INode) = tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    else if(INode .eq. NNodes) then  ! Low surface

                        DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                        if(NNodes .LE. 1) then
                            DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(INode)

                            MatA(INode) = 0.D0
                            select case(Host_SimuCtrlParam%BDCTYPE(3,2))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = (NodeSpace(INode)/TSTEP - DiffGradient1 - DiffGradient2)*Concent(INode,IKind)
                                case(p_Neumann_BDC)
                                    MatB(INode) = (NodeSpace(INode)/TSTEP)*Concent(INode,IKind)
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,2)
                                    pause
                                    stop
                            end select
                        else
                            DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))

                            MatA(INode) = DiffGradient1*Concent(INode-1,IKind)
                            select case(Host_SimuCtrlParam%BDCTYPE(3,2))
                                case(p_Dirichlet_BDC)
                                    MatB(INode) = (NodeSpace(INode)/TSTEP - DiffGradient1 - DiffGradient2)*Concent(INode,IKind)
                                case(p_Neumann_BDC)
                                    MatB(INode) = (NodeSpace(INode)/TSTEP - DiffGradient1)*Concent(INode,IKind)
                                case default
                                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",Host_SimuCtrlParam%BDCTYPE(3,2)
                                    pause
                                    stop
                            end select

                        end if

                        MatC(INode) = 0.D0
                        MatD(INode) = tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    else
                        DiffGradient1 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode-1) + NodeSpace(INode))
                        DiffGradient2 = (ClustersKind(IKind)%m_DiffCoeff + ClustersKind(IKind)%m_DiffCoeff)/(NodeSpace(INode) + NodeSpace(INode+1))
                        MatA(INode) = DiffGradient1*Concent(INode-1,IKind)
                        MatB(INode) = (NodeSpace(INode)/TSTEP - (DiffGradient1 + DiffGradient2))*Concent(INode,IKind)
                        MatC(INode) = DiffGradient2*Concent(INode+1,IKind)
                        MatD(INode) = tempNBPVChangeRate(INode,IKind)*NodeSpace(INode) + ImplantedRate(INode,IKind)*NodeSpace(INode)
                    end if

                END DO

!                DiffGradient2 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(NNodes)
!                FOutEachStep(IKind) = DiffGradient2*Concent(NNodes,IKind)
!                COutEachStep(IKind) = DiffGradient2*Concent(NNodes,IKind)*TSTEP/NodeSpace(NNodes)
!                FOutAccum(IKind) = FOutAccum(IKind) + FOutEachStep(IKind)
!                COutAccum(IKind) = COutAccum(IKind) + COutEachStep(IKind)
!
!
!                DiffGradient1 = ClustersKind(IKind)%m_DiffCoeff/NodeSpace(1)
!                FSurfEachStep(IKind) = DiffGradient1*Concent(1,IKind)
!                CSurfEachStep(IKind) = DiffGradient1*Concent(1,IKind)*TSTEP/NodeSpace(1)
!                FSurfAccum(IKind) = FSurfAccum(IKind) + FSurfEachStep(IKind)
!                CSurfAccum(IKind) = CSurfAccum(IKind) + CSurfEachStep(IKind)


                DO INode = 1,NNodes
                    Concent(INode,IKind) = (MatA(INode) + MatB(INode) + MatC(INode) + MatD(INode))*TSTEP/NodeSpace(INode)
                END DO


                if(IKind .eq. 1) then

!                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out flux from up surface",FSurfAccum(IKind)
!                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out concentrate from up surface",CSurfAccum(IKind)
!                    write(*,fmt="(A50,1x,1ES14.6)") "Out flux from up surface in current step",FSurfEachStep(IKind)
!                    write(*,fmt="(A50,1x,1ES14.6)") "Out concentrate from up surface in current step",CSurfEachStep(IKind)
!                    DO INode = 1,NNodes
!                        write(*,fmt="(A10,1x,I10,1x,A15,1x,1ES14.6)") "INode",INode,"Concent(INode,1)",Concent(INode,1)
!                    END DO
!                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out flux from lower surface",FOutAccum(IKind)
!                    write(*,fmt="(A50,1x,1ES14.6)") "Accumulated out concentrate from lower surface",COutAccum(IKind)
!                    write(*,fmt="(A50,1x,1ES14.6)") "Out flux from lower surface in current step",FOutEachStep(IKind)
!                    write(*,fmt="(A50,1x,1ES14.6)") "Out concentrate from lower surface in current step",COutEachStep(IKind)
!
!                    write(*,*) "----------------------------------------"
!                    write(*,*) "sum(Concent(:,1))",sum(Concent(:,1))
!                    write(*,*) "sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)",sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)
!                    write(*,*) "(sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1))/ConCentrat0",(ConCentrat0 - (sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)))/ConCentrat0
!                    write(*,*) "----------------------------------------"
                end if


            END DO

            call Record%AddSimuTimes(TSTEP)

            call OutPutCurrent(Host_SimBoxes,Host_SimuCtrlParam,Record)

            if(mod(Record%GetSimuSteps(),1) .eq. 0) then

                call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

                call Put_Out_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,Record%GetSimuSteps(),Record%GetSimuTimes(),TSTEP,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
            end if

            !if(Concent(CKind) .GT. 1.D-10) then
            if(DSQRT(dble(sum(ClustersKind(CKind)%m_Atoms(:)%m_NA)))*sum(Concent(1:NNodes,CKind)) .GT. &
               sum(NPOWER1DIV2Ave)*Host_SimuCtrlParam%DumplicateFactor) then

                DO INode = 1,NNodes

                    DO I = 1,(CKind -1)/2 + 1
                        if(2*I .LE. CKind) then
                            Concent(INode,I) = (Concent(INode,2*I-1)*sum(ClustersKind(2*I-1)%m_Atoms(:)%m_NA)+Concent(INode,2*I)*sum(ClustersKind(2*I)%m_Atoms(:)%m_NA)) &
                                               /(2.D0*sum(ClustersKind(2*I)%m_Atoms(:)%m_NA))
                        else
                            Concent(INode,I) = Concent(INode,2*I - 1)
                        end if
                    END DO

                    Concent(INode,(CKind -1)/2+2:CKind) = 0.D0

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


            if(Record%GetSimuTimes() .GT. Host_SimuCtrlParam%TermTValue) then
                exit
            end if

            TSTEP = TSTEP*1.5

          END DO
        END Associate

    end subroutine NucleationSimu_SpaceDist_Transient_GrubDumplicate_GPU

    !----------------------------------------------------------------------
    subroutine SolveTridag(IKind,MatrixA,MatrixB,MatrixC,MatrixD,Solver,MatrixSize,w,h)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::IKind
        real(kind=KMCDF),dimension(:),allocatable::MatrixA
        real(kind=KMCDF),dimension(:),allocatable::MatrixB
        real(kind=KMCDF),dimension(:),allocatable::MatrixC
        real(kind=KMCDF),dimension(:),allocatable::MatrixD
        real(kind=KMCDF),dimension(:,:)::Solver
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
        Solver(MatrixSize,IKind) = h(MatrixSize)

        DO I = MatrixSize-1,1,-1
            Solver(I,IKind) = h(I) - MatrixC(I)*Solver(I+1,IKind)/w(I)
        END DO


        return
    end subroutine SolveTridag




end module NUCLEATION_SPACEDIST_GPU
