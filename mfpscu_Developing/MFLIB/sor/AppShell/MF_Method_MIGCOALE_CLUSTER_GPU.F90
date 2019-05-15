module MF_Method_MIGCOALE_CLUSTER_GPU
    use MF_Method_MIGCOALE_CLUSTER_CPU
    use MIGCOALE_ADDONDATA_DEV
    use MCMF_CONSTANTS_GPU
    use NUCLEATION_SPACEDIST_GPU
    implicit none

    contains

    !*****************************************************************
    subroutine For_One_Test_GPU(Host_SimBoxes,Host_SimuCtrlParam,JobIndex)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        integer,intent(in)::JobIndex
        !---Local Vars---
        type(SimulationCtrlParam),pointer::PSimCtrlParam=>null()
        type(ImplantSection),pointer::PImplantSection=>null()
        integer::ITIME0
        real::TIME0
        integer::ITEST
        integer::TestLoop0,TestLoop1
        integer::TimeSection0
        !---Body---
        if(Host_SimuCtrlParam%INDEPBOX) then
            TestLoop0 = 1
            TestLoop1 = 1
        else
            TestLoop0 = 1
            TestLoop1 = Host_SimuCtrlParam%TOTALBOX/Host_SimuCtrlParam%MultiBox
        end if

        ITIME0 = 0
        TIME0 = ZERO
        TimeSection0 = 1

        PSimCtrlParam=>Host_SimuCtrlParam

        DO ITEST = TestLoop0,TestLoop1

            if(Host_SimuCtrlParam%INDEPBOX) then
                call m_MigCoalClusterRecord%InitMigCoalClusterRecord(PSimCtrlParam%MultiBox,SimuSteps=ITIME0,SimuTimes=TIME0,SimuPatchs=JobIndex,TimeSection=TimeSection0)
            else
                call m_MigCoalClusterRecord%InitMigCoalClusterRecord(PSimCtrlParam%MultiBox,SimuSteps=ITIME0,SimuTimes=TIME0,SimuPatchs=ITEST,TimeSection=TimeSection0)
            end if

            if(ITEST .eq. 1) then

                call m_MigCoaleStatInfoWrap%Init(Host_SimuCtrlParam%MultiBox)

                call m_ImplantList%Init(Host_SimBoxes,PSimCtrlParam)

                call resolveAddOnData(Host_SimBoxes,PSimCtrlParam)
                call InitSimulationBoxesConfig(Host_SimBoxes,PSimCtrlParam,m_InitBoxSimCfgList,m_MigCoaleStatInfoWrap,m_MigCoalClusterRecord)

            end if

            DO While(.true.)

                write(*,*) "Start to evolution for time section ",m_MigCoalClusterRecord%GetTimeSections()

                call m_MigCoalClusterRecord%SetStartImplantTime(m_MigCoalClusterRecord%GetSimuTimes())

                call copyInPhyParamsConstant(PSimCtrlParam)

                call resolveAddOnData(Host_SimBoxes,PSimCtrlParam)

                call CopyAddOnDataToDev()

                PImplantSection=>m_ImplantList%Get_P(PSimCtrlParam%ImplantSectID)

                call For_One_TimeSect_GPU(Host_SimBoxes,PSimCtrlParam,PImplantSection,m_MigCoaleStatInfoWrap,m_MigCoalClusterRecord)

                call m_MigCoalClusterRecord%SetLastRecordImplantNum(m_MigCoalClusterRecord%GetImplantedEntitiesNum())

                PSimCtrlParam=>PSimCtrlParam%next

                if(.not. associated(PSimCtrlParam)) then
                    exit
                end if

                call m_MigCoalClusterRecord%IncreaseOneTimeSection()

            END DO

        END DO

        return

    end subroutine For_One_Test_GPU

    !****************************************************************
    subroutine For_One_TimeSect_GPU(Host_SimBoxes,Host_SimuCtrlParam,TheImplantSection,TheMigCoaleStatInfoWrap,Record)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(ImplantSection)::TheImplantSection
        type(MigCoaleStatInfoWrap)::TheMigCoaleStatInfoWrap
        type(MigCoalClusterRecord)::Record
        !---Local Vars---
        integer::TotalSize
        !---Body---

        call InitSimu_SpaceDist_GPU(Host_SimBoxes,Host_SimuCtrlParam)

        call NucleationSimu_SpaceDist_Balance_GrubDumplicate_GPU(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,TheImplantSection)


!        Associate(Host_ClustesInfo=>Host_Boxes%m_ClustersInfo_CPU,Dev_ClustesInfo=>Dev_Boxes%dm_ClusterInfo_GPU, &
!              TBasicInfo=>Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral)
!
!            DO WHILE(.TRUE.)
!
!                call For_One_Step(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,TheImplantSection,TheMigCoaleStatInfoWrap,Record,TSTEP)
!
!                if(Host_SimuCtrlParam%TUpdateStatisFlag .eq. mp_UpdateStatisFlag_ByIntervalSteps) then
!                    if ((Record%GetSimuSteps() - Record%GetLastUpdateStatisTime()) .GE. Host_SimuCtrlParam%TUpdateStatisValue) then
!                        call GetBoxesMigCoaleStat_Expd_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd,Record)
!
!                        call Record%SetLastUpdateStatisTime(dble(Record%GetSimuSteps()))
!
!                    end if
!
!                else if(Host_SimuCtrlParam%TUpdateStatisFlag .eq. mp_UpdateStatisFlag_ByIntervalRealTime) then
!                    if((Record%GetSimuTimes() - Record%GetLastUpdateStatisTime()) .GE. Host_SimuCtrlParam%TUpdateStatisValue) then
!
!                        call GetBoxesMigCoaleStat_Expd_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd,Record)
!
!                        call Record%SetLastUpdateStatisTime(Record%GetSimuTimes())
!                    end if
!                end if
!
!                call OutPutCurrent(Host_Boxes,Dev_Boxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Used,Record)
!
!
!                if(Record%GetSimuTimes() .GE. Host_SimuCtrlParam%TermTValue) then
!                    exit
!                end if
!
!                !Check if need to duplicate the box
!                if(present(NCUT)) then
!                    if(NAct .LE. NCUT) then
!                        exit
!                    end if
!                end if
!
!            END DO
!
!        END Associate
!
!        call PutOut_Instance_Statistic_IntegralBox(Host_SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Used,Record,Model=1)

        return
    end subroutine For_One_TimeSect_GPU
!
!    !*********************************************************************
!    subroutine For_One_Step_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,TheImplantSection,TheMigCoaleStatInfoWrap,Record,TSTEP)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        type(ImplantSection)::TheImplantSection
!        type(MigCoaleStatInfoWrap)::TheMigCoaleStatInfoWrap
!        type(MigCoalClusterRecord)::Record
!        real(kind=KMCDF)::TSTEP
!        !---Local Vars---
!
!        call UpdateTimeStep_MigCoal(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd,Record,TSTEP)
!
!        if(TheImplantSection%ImplantFlux .GT. 0.D0) then
!            call TheImplantSection%ImplantClusters_FastStrategy(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,TheMigCoaleStatInfoWrap,Record,TSTEP)
!        end if
!
!        call EvoluteOneStep_DepthDist(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,TSTEP)
!
!        call Record%IncreaseOneSimuStep()
!
!        call Record%AddSimuTimes(TSTEP)
!
!        return
!    end subroutine For_One_Step_GPU

end module MF_Method_MIGCOALE_CLUSTER_GPU
