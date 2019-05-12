module NUCLEATION_SPACEDIST
    use MFLIB_GLOBAL
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


        call AvailableIOUnit(m_StatisticFile)

        fileName = Host_SimuCtrlParam%OutFilePath(1:len_trim(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"Statistic.dat"

        open(Unit=m_StatisticFile,file=fileName(1:len_trim(fileName)))

        write(m_StatisticFile,FMT="(15(A15,1x))") "Step","Time","TStep","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"

        return
    end subroutine InitSimu_SpaceDist


    !***************************************************
    subroutine NucleationSimu_SpaceDist_Balance_GrubDumplicate(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,TheImplantSection)
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

    end subroutine NucleationSimu_SpaceDist_Balance_GrubDumplicate

    !***************************************************
    subroutine NucleationSimu_SpaceDist_Transient_GrubDumplicate(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,TheImplantSection)
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

    end subroutine NucleationSimu_SpaceDist_Transient_GrubDumplicate

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

        Associate(ClustersKind=>Host_SimBoxes%m_ClustersInfo_CPU%ClustersKindArray,Concent=>Host_SimBoxes%m_ClustersInfo_CPU%Concentrate)

          DO INode = 1,NNodes

            NPOWER0Ave(INode) = Dumplicate*sum(Concent(INode,1:CKind))

            Temp = 0.D0
            DO I = 1,CKind
                Temp = Temp + (sum(ClustersKind(I)%m_Atoms(:)%m_NA)**(0.5D0))*Concent(INode,I)
            END DO
            NPOWER1DIV2Ave(INode) = Dumplicate*Temp

            Temp = 0.D0
            DO I = 1,CKind
                Temp = Temp + sum(ClustersKind(I)%m_Atoms(:)%m_NA)*Concent(INode,I)
            END DO
            NPOWER1Ave(INode) = Dumplicate*Temp

            Temp = 0.D0
            DO I = 1,CKind
                Temp = Temp + (sum(ClustersKind(I)%m_Atoms(:)%m_NA)**(1.5D0))*Concent(INode,I)
            END DO
            NPOWER3DIV2Ave(INode) = Dumplicate*Temp
          END DO

        END Associate
        return
    end subroutine

    !---------------------------------------------
    subroutine OutPutCurrent(Host_SimBoxes,Host_SimuCtrlParam,Record)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        type(MigCoalClusterRecord)::Record
        !---Local Vars---
        logical::OutIntegralBoxStatistic
        logical::OutEachBoxStatistic
        !---Body---

        OutIntegralBoxStatistic = .false.
        OutEachBoxStatistic = .true.

        OutIntegralBoxStatistic = Record%WhetherOutSizeDist_IntegralBox(Host_SimuCtrlParam)

        OutEachBoxStatistic = Record%WhetherOutSizeDist_EachBox(Host_SimuCtrlParam)

        if(OutIntegralBoxStatistic .eq. .true. .or. OutEachBoxStatistic .eq. .true.) then
!            call GetBoxesMigCoaleStat_Used_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatisticInfo,Record)
!
!            if(OutIntegralBoxStatistic .eq. .true.) then
!                call PutOut_Instance_Statistic_IntegralBox(Host_Boxes,Host_SimuCtrlParam,TheMigCoaleStatisticInfo,Record,Model=0)
!
!                if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalSteps) then
!                    call Record%SetLastOutSizeDistTime_IntegralBox(dble(Record%GetSimuSteps()))
!                else if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalRealTime) then
!                    call Record%SetLastOutSizeDistTime_IntegralBox(Record%GetSimuTimes())
!                else if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalTimeMagnification) then
!                    call Record%SetLastOutSizeDistTime_IntegralBox(Record%GetSimuTimes())
!                end if
!            end if
!
!            if(OutEachBoxStatistic .eq. .true.) then
!                call PutOut_Instance_Statistic_EachBox(Host_Boxes,Host_SimuCtrlParam,TheMigCoaleStatisticInfo,Record)
!
!                if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalSteps) then
!                    call Record%SetLastOutSizeDistTime_EachBox(dble(Record%GetSimuSteps()))
!                else if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalRealTime) then
!                    call Record%SetLastOutSizeDistTime_EachBox(Record%GetSimuTimes())
!                else if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalTimeMagnification) then
!                    call Record%SetLastOutSizeDistTime_EachBox(Record%GetSimuTimes())
!                end if
!            end if
        end if

        ! check if need to output intermediate configure
        if(Host_SimuCtrlParam%OutPutConfFlag .eq. mp_OutTimeFlag_ByIntervalSteps) then
            if((Record%GetSimuSteps() - Record%GetLastRecordOutConfigTime()) .GE. Host_SimuCtrlParam%OutPutConfValue .OR. &
                Record%GetSimuTimes() .GE. Host_SimuCtrlParam%TermTValue) then

                call Host_SimBoxes%PutoutCfg(Host_SimuCtrlParam,Record,LayerFirst=.true.)

                call Host_SimBoxes%PutoutCfg(Host_SimuCtrlParam,Record,LayerFirst=.false.)

                call Record%SetLastRecordOutConfigTime(dble(Record%GetSimuSteps()))

                call Record%IncreaseOneOutPutIndex()

            end if
        else if(Host_SimuCtrlParam%OutPutConfFlag .eq. mp_OutTimeFlag_ByIntervalRealTime) then
            if((Record%GetSimuTimes() - Record%GetLastRecordOutConfigTime()) .GE. Host_SimuCtrlParam%OutPutConfValue .OR. &
                Record%GetSimuTimes() .GE. Host_SimuCtrlParam%TermTValue) then

                call Host_SimBoxes%PutoutCfg(Host_SimuCtrlParam,Record,LayerFirst=.true.)

                call Host_SimBoxes%PutoutCfg(Host_SimuCtrlParam,Record,LayerFirst=.false.)

                call Record%SetLastRecordOutConfigTime(Record%GetSimuTimes())

                call Record%IncreaseOneOutPutIndex()
            end if

        else if(Host_SimuCtrlParam%OutPutConfFlag .eq. mp_OutTimeFlag_ByIntervalTimeMagnification) then
            if((Record%GetSimuTimes()/Host_SimuCtrlParam%OutPutConfValue) .GE. Record%GetLastRecordOutConfigTime() .OR. &
                Record%GetSimuTimes() .GE. Host_SimuCtrlParam%TermTValue) then

                call Host_SimBoxes%PutoutCfg(Host_SimuCtrlParam,Record,LayerFirst=.true.)

                call Host_SimBoxes%PutoutCfg(Host_SimuCtrlParam,Record,LayerFirst=.false.)

                call Record%SetLastRecordOutConfigTime(Record%GetSimuTimes())

                call Record%IncreaseOneOutPutIndex()
            end if
        end if

        return
    end subroutine OutPutCurrent

    !---------------------------------------------
    subroutine Put_Out_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,Step,TTime,TStep,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        integer,intent(in)::Step
        real(kind=KMCDF),intent(in)::TTime
        real(kind=KMCDF),intent(in)::TStep
        real(kind=KMCDF),dimension(:),allocatable::NPOWER0Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER1Ave
        real(kind=KMCDF),dimension(:),allocatable::NPOWER3DIV2Ave
        real(kind=KMCDF),dimension(:),allocatable::N1
        real(kind=KMCDF),dimension(:),allocatable::N2
        real(kind=KMCDF),dimension(:),allocatable::N3
        real(kind=KMCDF),dimension(:),allocatable::Rave
        !---Local Vars---
        integer::CKind
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

        Associate(ClustersKind=>Host_SimBoxes%m_ClustersInfo_CPU%ClustersKindArray,Concent=>Host_SimBoxes%m_ClustersInfo_CPU%Concentrate)

          DO INode = 1,NNodes

            if(NPOWER0Ave(INode) .GT. 0) then
                N1(INode) = (NPOWER1DIV2Ave((INode))/NPOWER0Ave(INode))**2

                N2(INode) = NPOWER1Ave(INode)/NPOWER0Ave(INode)

                N3(INode) = (NPOWER3DIV2Ave(INode)/NPOWER0Ave(INode))**(dble(2)/dble(3))

                Rave(INode) = DSQRT(N1(INode)/m_RNFACTOR)
            end if
          END DO

          write(6,FMT="(15(A15,1x))") "Step","Time","TStep","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"

          write(6,FMT="(I15,1x,15(E15.8,1x))") Step,TTime,TStep,sum(NPOWER0Ave)/NNodes,sum(NPOWER1DIV2Ave)/NNodes,sum(NPOWER1Ave)/NNodes,sum(NPOWER3DIV2Ave)/NNodes, &
                                             sum(N1)/NNodes,sum(N2)/NNodes,sum(N3)/NNodes,sum(Rave)*C_CM2NM/NNodes

          write(m_StatisticFile,FMT="(I15,1x,15(E15.8,1x))") Step,TTime,TStep,sum(NPOWER0Ave)/NNodes,sum(NPOWER1DIV2Ave)/NNodes,sum(NPOWER1Ave)/NNodes,sum(NPOWER3DIV2Ave)/NNodes, &
                                             sum(N1)/NNodes,sum(N2)/NNodes,sum(N3)/NNodes,sum(Rave)*C_CM2NM/NNodes

          Flush(m_StatisticFile)


         END associate

        return
    end subroutine

end module NUCLEATION_SPACEDIST
