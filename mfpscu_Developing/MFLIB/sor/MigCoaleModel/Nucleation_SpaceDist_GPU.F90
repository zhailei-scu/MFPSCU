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
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_MatrixA
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_MatrixB
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_MatrixC
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_MatrixD
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_w
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_h

    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_ImplantedRate

    real(kind=KMCDF),device,dimension(:),allocatable::dm_Reduced_MaxChangeRate
    real(kind=KMCDF),device,dimension(:),allocatable::dm_Reduced_MaxConcent

    integer::Dumplicate = 1

    integer::LastOutPutConifgIndex = 0

    integer::m_StatisticFile

    !--------
    type(ACluster),device,dimension(:),allocatable::dm_ClustersKindArray
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_Concentrate
    real(kind=KMCDF),device,dimension(:),allocatable::dm_NodeSpace
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_tempNBPVChangeRate

    contains

    subroutine InitSimu_SpaceDist_GPU(Host_Boxes,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer::I
        integer::NNodes
        integer::CKind
        integer::NReduceSize
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

        if(allocated(dm_tempNBPVChangeRate)) then
            deallocate(dm_tempNBPVChangeRate)
        end if
        allocate(dm_tempNBPVChangeRate(NNodes,CKind))
        dm_tempNBPVChangeRate = 0.D0

        if(allocated(dm_MatrixA)) then
            deallocate(dm_MatrixA)
        end if
        allocate(dm_MatrixA(NNodes,CKind))

        if(allocated(dm_MatrixB)) then
            deallocate(dm_MatrixB)
        end if
        allocate(dm_MatrixB(NNodes,CKind))

        if(allocated(dm_MatrixC)) then
            deallocate(dm_MatrixC)
        end if
        allocate(dm_MatrixC(NNodes,CKind))

        if(allocated(dm_MatrixD)) then
            deallocate(dm_MatrixD)
        end if
        allocate(dm_MatrixD(NNodes,CKind))

        if(allocated(dm_w)) then
            deallocate(dm_w)
        end if
        allocate(dm_w(NNodes,CKind))

        if(allocated(dm_h)) then
            deallocate(dm_h)
        end if
        allocate(dm_h(NNodes,CKind))

        if(allocated(dm_ImplantedRate)) then
            deallocate(dm_ImplantedRate)
        end if
        allocate(dm_ImplantedRate(NNodes,CKind))

        NReduceSize = (NNodes*CKind - 1)/p_BLOCKSIZE + 1
        if(allocated(dm_Reduced_MaxChangeRate)) then
            deallocate(dm_Reduced_MaxChangeRate)
        end if
        allocate(dm_Reduced_MaxChangeRate(NReduceSize))

        if(allocated(dm_Reduced_MaxConcent)) then
            deallocate(dm_Reduced_MaxConcent)
        end if
        allocate(dm_Reduced_MaxConcent(NReduceSize))


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

        !call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

        Associate(ClustersKind=>Host_SimBoxes%m_ClustersInfo_CPU%ClustersKindArray,Concent=>Host_SimBoxes%m_ClustersInfo_CPU%Concentrate,NodeSpace=>Host_SimBoxes%NodeSpace)

          ConCentrat0 = sum(Concent)

          ImplantedRate = 0.D0

          call TheImplantSection%Cal_ImplantClustersRate(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,ImplantedRate)

          DO While(.true.)

            call Record%IncreaseOneSimuStep()

            tempNBPVChangeRate = 0.D0

            FSurfEachStep = 0.D0

            FOutEachStep = 0.D0

            CSurfEachStep = 0.D0

            COutEachStep = 0.D0

            call Calc_ReactionRate(Host_SimBoxes,Host_SimuCtrlParam,dm_Concentrate,dm_ClustersKindArray,dm_tempNBPVChangeRate,dm_Reduced_MaxConcent,dm_Reduced_MaxChangeRate)

            TSTEP = Host_SimuCtrlParam%MaxReactChangeRate*maxval(dm_Reduced_MaxConcent)/maxval(dm_Reduced_MaxChangeRate)

            write(*,*) "TSTEP",TSTEP

            call DoReactionAndDiffusion(Host_SimBoxes,Host_SimuCtrlParam,dm_Concentrate,dm_NodeSpace,dm_ClustersKindArray,dm_tempNBPVChangeRate, &
                                        dm_ImplantedRate,dm_MatrixA,dm_MatrixB,dm_MatrixC,dm_MatrixD,dm_w,dm_h,TSTEP)

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

!                if(IKind .eq. 1) then
!
!                    write(*,*) "Accumulated out flux from up surface",FSurfAccum(IKind)
!                    write(*,*) "Accumulated out concentrate from up surface",CSurfAccum(IKind)
!                    write(*,*) "Out flux from up surface in current step",FSurfEachStep(IKind)
!                    write(*,*) "Out concentrate from up surface in current step",CSurfEachStep(IKind)
!                    DO INode = 1,NNodes
!                        write(*,*) "INode",INode,"Concent(INode,1)",Concent(INode,1)
!                    END DO
!                    write(*,*) "Accumulated out flux from lower surface",FOutAccum(IKind)
!                    write(*,*) "Accumulated out concentrate from lower surface",COutAccum(IKind)
!                    write(*,*) "Out flux from lower surface in current step",FOutEachStep(IKind)
!                    write(*,*) "Out concentrate from lower surface in current step",COutEachStep(IKind)
!
!                    write(*,*) "----------------------------------------"
!                    write(*,*) "sum(Concent(:,1))",sum(Concent(:,1))
!                    write(*,*) "sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)",sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)
!                    write(*,*) "(sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1))/ConCentrat0",(ConCentrat0 - (sum(Concent(:,1)) + CSurfAccum(1) + COutAccum(1)))/ConCentrat0
!                    write(*,*) "----------------------------------------"
!                end if

            call Record%AddSimuTimes(TSTEP)

            !call OutPutCurrent(Host_SimBoxes,Host_SimuCtrlParam,Record)

            if(mod(Record%GetSimuSteps(),1) .eq. 0) then

                !call Cal_Statistic_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave)

                !call Put_Out_IMPLANT(Host_SimBoxes,Host_SimuCtrlParam,Record%GetSimuSteps(),Record%GetSimuTimes(),TSTEP,NPOWER0Ave,NPOWER1DIV2Ave,NPOWER1Ave,NPOWER3DIV2Ave,N1,N2,N3,Rave)
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
    subroutine Calc_ReactionRate(Host_SimBoxes,Host_SimuCtrlParam,Dev_Concent,Dev_ClusterKindArray,Dev_tempNBPVChangeRate,Reduced_MaxConcent,Reduced_MaxChangeRate)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_Concent
        type(ACluster),device,dimension(:),allocatable::Dev_ClusterKindArray
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_tempNBPVChangeRate
        real(kind=KMCDF),device,dimension(:),allocatable::Reduced_MaxConcent
        real(kind=KMCDF),device,dimension(:),allocatable::Reduced_MaxChangeRate
        !---Local Vars---
        integer::NNodes
        integer::CKind
        integer::NT
        integer::NB
        integer::BX
        integer::BY
        type(dim3)::blocks
        type(dim3)::threads
        !---Body---
        NNodes = Host_SimBoxes%NNodes
        CKind = Host_SimBoxes%CKind

        NT = NNodes*CKind
        NB = (NT - 1)/p_BLOCKSIZE + 1
        BX = p_BLOCKSIZE
        BY = 1
        blocks = dim3(NB,1,1)
        threads = dim3(BX,BY,1)

        call Kernel_CalReaction_Vanish<<<blocks,threads>>>(NNodes,CKind,Dev_Concent,Dev_ClusterKindArray,Dev_tempNBPVChangeRate)

        call Kernel_CalReaction_Generate<<<blocks,threads>>>(NNodes,CKind,Host_SimuCtrlParam%MaxReactChangeRate,Dev_Concent,Dev_ClusterKindArray,&
                                                             Dev_tempNBPVChangeRate,Reduced_MaxConcent,Reduced_MaxChangeRate)

        return
    end subroutine Calc_ReactionRate


    !***************************************************
    subroutine DoReactionAndDiffusion(Host_SimBoxes,Host_SimuCtrlParam,Dev_Concent,Dev_NodeSpace,Dev_ClusterKindArray,Dev_tempNBPVChangeRate,&
                                    Dev_ImplantedRate,Dev_MatrixA,Dev_MatrixB,Dev_MatrixC,Dev_MatrixD,Dev_w,Dev_h,TSTEP)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_Concent
        real(kind=KMCDF),device,dimension(:),allocatable::Dev_NodeSpace
        type(ACluster),device,dimension(:),allocatable::Dev_ClusterKindArray
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_tempNBPVChangeRate(:,:)
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_ImplantedRate(:,:)
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_MatrixA(:,:)
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_MatrixB(:,:)
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_MatrixC(:,:)
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_MatrixD(:,:)
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_w(:,:)
        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_h(:,:)
        real(kind=KMCDF)::TSTEP
        !---Local Vars---
        integer::NNodes
        integer::CKind
        integer::NB
        integer::BX
        integer::BY
        type(dim3)::blocks
        type(dim3)::threads
        !---Body---
        NNodes = Host_SimBoxes%NNodes
        CKind = Host_SimBoxes%CKind

        NB = (CKind -1)/p_BLOCKSIZE + 1
        BX = p_BLOCKSIZE
        BY = 1
        blocks = dim3(NB,1,1)
        threads = dim3(BX,BY,1)

        call Kernel_DoReaction_And_NodeDiffusion_Balance<<<blocks,threads>>>(NNodes,CKind,TSTEP,Dev_Concent,Dev_NodeSpace,&
                                                         Dev_ClusterKindArray,Dev_tempNBPVChangeRate,Dev_ImplantedRate,&
                                                         Dev_MatrixA,Dev_MatrixB,Dev_MatrixC,Dev_MatrixD,Dev_w,Dev_h)

        return
    end subroutine DoReactionAndDiffusion

    !***************************************************
    attributes(global) subroutine Kernel_CalReaction_Vanish(NNodes,CKind,Dev_Concent,Dev_ClusterKindArray,Dev_tempNBPVChangeRate)
        implicit none
        !---Dummy Vars---
        integer,value::NNodes
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

        IKind = (cid - 1)/NNodes + 1
        INode = cid - (IKind-1)*NNodes

        if(IKind .LE. CKind) then

            Dev_tempNBPVChangeRate(INode,IKind) = 0.D0

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
        end if

        return
    end subroutine Kernel_CalReaction_Vanish

    !***************************************************
    attributes(global) subroutine Kernel_CalReaction_Generate(NNodes,CKind,MaxReactChangeRate,Dev_Concent,Dev_ClusterKindArray, &
                                                              Dev_tempNBPVChangeRate,Reduced_MaxConcent,Reduced_MaxChangeRate)
        implicit none
        !---Dummy Vars---
        integer,value::NNodes
        integer,value::CKind
        real(kind=KMCDF),value::MaxReactChangeRate
        real(kind=KMCDF),device::Dev_Concent(:,:)
        type(ACluster),device::Dev_ClusterKindArray(:)
        real(kind=KMCDF),device::Dev_tempNBPVChangeRate(:,:)
        real(kind=KMCDF),device::Reduced_MaxConcent(:)
        real(kind=KMCDF),device::Reduced_MaxChangeRate(:)
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
        integer::I
        real(kind=KMCDF),shared::Shared_MaxChangeRate(p_BLOCKSIZE)
        real(kind=KMCDF),shared::Shared_MaxConcent(p_BLOCKSIZE)
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x +threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*p_BlockSize + tid

        LKind = (cid - 1)/NNodes + 1
        INode = cid - (LKind-1)*NNodes

        Shared_MaxConcent(tid) = 0.D0

        Shared_MaxChangeRate(tid) = 0.D0

        if(LKind .LE. CKind) then
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

            Shared_MaxConcent(tid) = Dev_Concent(INode,LKind)
            Shared_MaxChangeRate(tid) = DABS(Dev_tempNBPVChangeRate(INode,LKind))
        end if

        call syncthreads()

        I = p_BLOCKSIZE/2

        DO While(I .GT. 0)

            if(tid .LE. I) then
                if(Shared_MaxChangeRate(tid) .LT. Shared_MaxChangeRate(tid+I)) then
                    Shared_MaxChangeRate(tid) = Shared_MaxChangeRate(tid+I)
                end if

                if(Shared_MaxConcent(tid) .LT. Shared_MaxConcent(tid+I)) then
                    Shared_MaxConcent(tid) = Shared_MaxConcent(tid+I)
                end if
            end if

            call syncthreads()

            I = I/2
        END DO

        if(tid .eq. 1) then
            Reduced_MaxConcent(bid) = Shared_MaxConcent(1)
            Reduced_MaxChangeRate(bid) = Shared_MaxChangeRate(1)
        end if

        return
    end subroutine Kernel_CalReaction_Generate

    !***************************************************
    attributes(global) subroutine Kernel_DoReaction_And_NodeDiffusion_Balance(NNodes,CKind,TSTEP,Dev_Concent,Dev_NodeSpace,&
                                                                              Dev_ClusterKindArray,Dev_tempNBPVChangeRate,Dev_ImplantedRate,&
                                                                              Dev_MatrixA,Dev_MatrixB,Dev_MatrixC,Dev_MatrixD,Dev_w,Dev_h)
        implicit none
        !---Dummy Vars---
        integer,value::NNodes
        integer,value::CKind
        real(kind=KMCDF),value::TSTEP
        real(kind=KMCDF),device::Dev_Concent(NNodes,*)
        real(kind=KMCDF),device::Dev_NodeSpace(:)
        type(ACluster),device::Dev_ClusterKindArray(:)
        real(kind=KMCDF),device::Dev_tempNBPVChangeRate(:,:)
        real(kind=KMCDF),device::Dev_ImplantedRate(:,:)
        real(kind=KMCDF),device::Dev_MatrixA(NNodes,*)
        real(kind=KMCDF),device::Dev_MatrixB(NNodes,*)
        real(kind=KMCDF),device::Dev_MatrixC(NNodes,*)
        real(kind=KMCDF),device::Dev_MatrixD(NNodes,*)
        real(kind=KMCDF),device::Dev_w(NNodes,*)
        real(kind=KMCDF),device::Dev_h(NNodes,*)
        !---Local Vars---
        integer::tid
        integer::bid
        integer::cid
        integer::IKind
        integer::INode
        real(kind=KMCDF)::DiffGradient1
        real(kind=KMCDF)::DiffGradient2
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x +threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*p_BlockSize + tid

        IKind = (cid - 1)/NNodes + 1

        if(IKind .LE. CKind) then

            Dev_MatrixA(:,IKind) = 0.D0
            Dev_MatrixB(:,IKind) = 0.D0
            Dev_MatrixC(:,IKind) = 0.D0
            Dev_MatrixD(:,IKind) = 0.D0

            DO INode = 1,NNodes
                if(INode .eq. 1) then  ! upper surface

                    DiffGradient1 = Dev_ClusterKindArray(IKind)%m_DiffCoeff/Dev_NodeSpace(INode)

                    Dev_MatrixA(INode,IKind) = 0.D0
                    if(NNodes .LE. 1) then
                        DiffGradient2 = Dev_ClusterKindArray(IKind)%m_DiffCoeff/Dev_NodeSpace(INode)

                        select case(dm_BDCTYPE(3,1))
                            case(p_Dirichlet_BDC)
                                Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                            case(p_Neumann_BDC)
                                Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP
                        end select

                        Dev_MatrixC(INode,IKind) = 0.D0
                    else
                        DiffGradient2 = (Dev_ClusterKindArray(IKind)%m_DiffCoeff + Dev_ClusterKindArray(IKind)%m_DiffCoeff)/(Dev_NodeSpace(INode) + Dev_NodeSpace(INode+1))

                        select case(dm_BDCTYPE(3,1))
                            case(p_Dirichlet_BDC)
                                Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                            case(p_Neumann_BDC)
                                Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP + DiffGradient2
                        end select

                        Dev_MatrixC(INode,IKind) = -DiffGradient2
                    end if

                    Dev_MatrixD(INode,IKind) = Dev_Concent(INode,IKind)*Dev_NodeSpace(INode)/TSTEP + Dev_tempNBPVChangeRate(INode,IKind)*Dev_NodeSpace(INode) + Dev_ImplantedRate(INode,IKind)*Dev_NodeSpace(INode)
                else if(INode .eq. NNodes) then  ! Low surface

                    DiffGradient2 = Dev_ClusterKindArray(IKind)%m_DiffCoeff/Dev_NodeSpace(INode)

                    if(NNodes .LE. 1) then
                        DiffGradient1 = Dev_ClusterKindArray(IKind)%m_DiffCoeff/Dev_NodeSpace(INode)

                        Dev_MatrixA(INode,IKind) = 0.D0

                        select case(dm_BDCTYPE(3,2))
                            case(p_Dirichlet_BDC)
                                Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                            case(p_Neumann_BDC)
                                Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP
                        end select
                    else
                        DiffGradient1 = (Dev_ClusterKindArray(IKind)%m_DiffCoeff + Dev_ClusterKindArray(IKind)%m_DiffCoeff)/(Dev_NodeSpace(INode-1) + Dev_NodeSpace(INode))

                        Dev_MatrixA(INode,IKind) = -DiffGradient1

                        select case(dm_BDCTYPE(3,2))
                            case(p_Dirichlet_BDC)
                                Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP + DiffGradient1 + DiffGradient2
                            case(p_Neumann_BDC)
                                Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP + DiffGradient1
                        end select

                    end if

                    Dev_MatrixC(INode,IKind) = 0.D0
                    Dev_MatrixD(INode,IKind) = Dev_Concent(INode,IKind)*Dev_NodeSpace(INode)/TSTEP + Dev_tempNBPVChangeRate(INode,IKind)*Dev_NodeSpace(INode) + Dev_ImplantedRate(INode,IKind)*Dev_NodeSpace(INode)
                else
                    DiffGradient1 = (Dev_ClusterKindArray(IKind)%m_DiffCoeff + Dev_ClusterKindArray(IKind)%m_DiffCoeff)/(Dev_NodeSpace(INode-1) + Dev_NodeSpace(INode))
                    DiffGradient2 = (Dev_ClusterKindArray(IKind)%m_DiffCoeff + Dev_ClusterKindArray(IKind)%m_DiffCoeff)/(Dev_NodeSpace(INode) + Dev_NodeSpace(INode+1))
                    Dev_MatrixA(INode,IKind) = -DiffGradient1
                    Dev_MatrixB(INode,IKind) = Dev_NodeSpace(INode)/TSTEP + (DiffGradient1 + DiffGradient2)
                    Dev_MatrixC(INode,IKind) = -DiffGradient2
                    Dev_MatrixD(INode,IKind) = Dev_Concent(INode,IKind)*Dev_NodeSpace(INode)/TSTEP + Dev_tempNBPVChangeRate(INode,IKind)*Dev_NodeSpace(INode) + Dev_ImplantedRate(INode,IKind)*Dev_NodeSpace(INode)
                end if

            END DO

            call Dev_SolveTridag(NNodes,IKind,Dev_MatrixA,Dev_MatrixB,Dev_MatrixC,Dev_MatrixD,Dev_Concent,Dev_w,Dev_h)

        end if

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

    !----------------------------------------------------------------------
    attributes(device) subroutine Dev_SolveTridag(NNodes,IKind,Dev_MatrixA,Dev_MatrixB,Dev_MatrixC,Dev_MatrixD,Dev_Solver,Dev_w,Dev_h)
        implicit none
        !---Dummy Vars---
        integer,value::NNodes
        integer,value::IKind
        real(kind=KMCDF),device::Dev_MatrixA(NNodes,*)
        real(kind=KMCDF),device::Dev_MatrixB(NNodes,*)
        real(kind=KMCDF),device::Dev_MatrixC(NNodes,*)
        real(kind=KMCDF),device::Dev_MatrixD(NNodes,*)
        real(kind=KMCDF),device::Dev_Solver(NNodes,*)
        real(kind=KMCDF),device::Dev_w(NNodes,*)
        real(kind=KMCDF),device::Dev_h(NNodes,*)
        !---Local Vars---
        integer::I
        !---Body---
        Dev_w(:,IKind) = 0.D0
        Dev_h(:,IKind) = 0.D0

        Dev_w(1,IKind) = Dev_MatrixB(1,IKind)
        Dev_h(1,IKind) = Dev_MatrixD(1,IKind)/Dev_w(1,IKind)
        DO I = 2,NNodes
            Dev_w(I,IKind) = Dev_MatrixB(I,IKind) - Dev_MatrixC(I,IKind)*Dev_MatrixA(I,IKind)/Dev_w(I-1,IKind)
            Dev_h(I,IKind) = (Dev_MatrixD(I,IKind) - Dev_MatrixA(I,IKind)*Dev_h(I-1,IKind))/Dev_w(I,IKind)
        END DO
        Dev_Solver(NNodes,IKind) = Dev_h(NNodes,IKind)

        DO I = NNodes-1,1,-1
            Dev_Solver(I,IKind) = Dev_h(I,IKind) - Dev_MatrixC(I,IKind)*Dev_Solver(I+1,IKind)/Dev_w(I,IKind)
        END DO

        return
    end subroutine Dev_SolveTridag




end module NUCLEATION_SPACEDIST_GPU
