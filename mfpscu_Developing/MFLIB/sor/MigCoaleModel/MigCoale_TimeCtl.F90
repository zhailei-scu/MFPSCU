module MIGCOALE_TIMECTL
    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
    use MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
    use MIGCOALE_TYPEDEF_STATISTICINFO
    implicit none


    real(kind=KMCDF), private, parameter::p_ZeroNCStep = 1.D-10

    contains

    !**************************************************************
!    subroutine UpdateTimeStep_MigCoal(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatisticInfo,Record,TSTEP)
!        ! To automatically determine the time step
!        !  Host_Boxes: the boxes info in host
!        !       INPUT: DIF - the diffusion coefficient
!        !      OUTPUT: TSTEP - the time step
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoaleStatisticInfo_Expd)::TheMigCoaleStatisticInfo
!        CLASS(SimulationRecord)::Record
!        real(kind=KMCDF)::TSTEP
!        !---Local Vars---
!        real(kind=KMCDF)::SEP, DIF
!        integer::NCActFree
!        integer::NCActGB
!        real(kind=KMCDF)::TSTEPFREE,TSTEPGB
!        !---Body---
!
!        TSTEPFREE = 1.D32
!        TSTEPGB = 1.D32
!
!        ASSOCIATE(TBasicInfo=>Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral,TMigStatInfo=>TheMigCoaleStatisticInfo%statistic_IntegralBox)
!
!        select case(Host_SimuCtrlParam%UPDATETSTEPSTRATEGY)
!            case(mp_SelfAdjustlStep_NearestSep)
!                if(Dev_Boxes%dm_ClusterInfo_GPU%GetNLUpdateCount_Dev() .LE. 0) then
!                    call Cal_Neighbor_List_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Record,IfDirectly=.true.,RMAX= &
!                                                max(TMigStatInfo%RMAX(p_ACTIVEFREE_STATU),TMigStatInfo%RMAX(p_ACTIVEINGB_STATU)))
!                end if
!
!                if(TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .LE. 0.D0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .LE. 0.D0) then
!                    TSTEP = p_ZeroNCStep
!                    return
!                end if
!
!                NCActFree = TBasicInfo%NC(p_ACTIVEFREE_STATU)
!                if(NCActFree .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .GT. 0.D0) then
!                    TSTEPFREE = (TBasicInfo%AveNearestSpeFreeClusters**2)/(6.D0*TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU))*Host_SimuCtrlParam%EnlageTStepScale
!                end if
!
!                NCActGB = TBasicInfo%NC(p_ACTIVEINGB_STATU)
!                if(NCActGB .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .GT. 0.D0) then
!                    TSTEPGB = (TBasicInfo%AveNearestSpeGBClusters**2)/(4.D0*TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU))*Host_SimuCtrlParam%EnlageTStepScale
!                end if
!
!                if((NCActFree .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .GT. 0.D0) &
!                   .or. (NCActGB .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .GT. 0.D0)) then
!                    TSTEP = min(TSTEPFREE,TSTEPGB)
!                else
!                    TSTEP = p_ZeroNCStep
!                end if
!
!            case(mp_SelfAdjustlStep_AveSep)
!                if(TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .LE. 0.D0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .LE. 0.D0) then
!                    TSTEP = p_ZeroNCStep
!                    return
!                end if
!
!                NCActFree = TBasicInfo%NC(p_ACTIVEFREE_STATU)
!                if(NCActFree .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .GT. 0.D0) then
!                    SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(NCActFree)
!
!                    SEP = SEP**(0.33333333333333D0)
!                    SEP = SEP - 2.D0*TMigStatInfo%RAVA(p_ACTIVEFREE_STATU)
!
!                    TSTEPFREE = SEP*SEP/(6.D0*TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU))*Host_SimuCtrlParam%EnlageTStepScale
!                end if
!
!                NCActGB = TBasicInfo%NC(p_ACTIVEINGB_STATU)
!                if(NCActGB .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .GT. 0.D0) then
!                    SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(Host_Boxes%m_GrainBoundary%GrainNum)
!
!                    SEP = 6*(SEP**(0.66666666667D0))*Host_Boxes%m_GrainBoundary%GrainNum
!                    SEP = SEP/dble(NCActGB)
!                    SEP = SEP**(0.5D0)
!                    SEP = SEP - 2.D0*TMigStatInfo%RAVA(p_ACTIVEINGB_STATU)
!
!                    TSTEPGB = SEP*SEP/(4.D0*TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU))*Host_SimuCtrlParam%EnlageTStepScale
!                end if
!
!                if((NCActFree .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .GT. 0.D0) &
!                   .or. (NCActGB .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .GT. 0.D0)) then
!                    TSTEP = min(TSTEPFREE,TSTEPGB)
!                else
!                    TSTEP = p_ZeroNCStep
!                end if
!
!            case(mp_FixedTimeStep)
!                TSTEP = Host_SimuCtrlParam%FixedTimeStepValue
!
!            case default
!                write(*,*) "MFPSCUERROR: Unknown strategy to update time step :",Host_SimuCtrlParam%UPDATETSTEPSTRATEGY
!                pause
!                stop
!        end select
!
!        END ASSOCIATE
!
!        return
!    end subroutine UpdateTimeStep_MigCoal
!
!    !*********************************************************
!    function Cal_VerifyTime_Implant(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatisticInfo,Record,ImplantedNumEachBox) result(TheVerifyTime)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoaleStatisticInfo_Virtual)::TheMigCoaleStatisticInfo
!        CLASS(SimulationRecord)::Record
!        integer,intent(in)::ImplantedNumEachBox
!        real(kind=KMCDF),intent(out)::TheVerifyTime
!        !---Local Vars---
!        real(kind=KMCDF)::SEP, DIF
!        integer::NCActFree,NCActGB
!        real(kind=KMCDF)::RAVA
!        real(kind=KMCDF)::TSTEPFREE,TSTEPGB
!        integer::MultiBox
!        !---Body---
!
!        TSTEPFREE = 1.D32
!        TSTEPGB = 1.D32
!
!        MultiBox = Host_SimuCtrlParam%MultiBox
!
!        ASSOCIATE(TBasicInfo=>Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral,TMigStatInfo=>TheMigCoaleStatisticInfo%statistic_IntegralBox)
!
!        select case(Host_SimuCtrlParam%UPDATETSTEPSTRATEGY)
!            case(mp_SelfAdjustlStep_NearestSep)
!                if(Dev_Boxes%dm_ClusterInfo_GPU%GetNLUpdateCount_Dev() .LE. 0) then
!                    call Cal_Neighbor_List_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Record,IfDirectly=.true.,RMAX= &
!                                                max(TMigStatInfo%RMAX(p_ACTIVEFREE_STATU),TMigStatInfo%RMAX(p_ACTIVEINGB_STATU)))
!                end if
!
!                if(TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .LE. 0.D0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .LE. 0.D0) then
!                    TheVerifyTime = p_ZeroNCStep
!                    return
!                end if
!
!                NCActFree = TBasicInfo%NC(p_ACTIVEFREE_STATU) + MultiBox*ImplantedNumEachBox
!                if(NCActFree .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .GT. 0.D0) then
!                    TSTEPFREE = (TBasicInfo%AveNearestSpeFreeClusters**2)/(6.D0*TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU))*Host_SimuCtrlParam%EnlageTStepScale
!                end if
!
!                NCActGB = TBasicInfo%NC(p_ACTIVEINGB_STATU) + MultiBox*ImplantedNumEachBox
!                if(NCActGB .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .GT. 0.D0) then
!                    TSTEPGB = (TBasicInfo%AveNearestSpeGBClusters**2)/(4.D0*TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU))*Host_SimuCtrlParam%EnlageTStepScale
!                end if
!
!                if((NCActFree .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .GT. 0.D0) &
!                   .or. (NCActGB .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .GT. 0.D0)) then
!                    TheVerifyTime = min(TSTEPFREE,TSTEPGB)
!                else
!                    TheVerifyTime = p_ZeroNCStep
!                end if
!
!            case(mp_SelfAdjustlStep_AveSep)
!                if(TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .LE. 0.D0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .LE. 0.D0) then
!                    TheVerifyTime = p_ZeroNCStep
!                    return
!                end if
!
!                NCActFree = TBasicInfo%NC(p_ACTIVEFREE_STATU) + MultiBox*ImplantedNumEachBox
!                if(NCActFree .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .GT. 0.D0) then
!                    SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(NCActFree)
!
!                    SEP = SEP**(0.33333333333333D0)
!                    SEP = SEP - 2.D0*TMigStatInfo%RAVA(p_ACTIVEFREE_STATU)
!
!                    TSTEPFREE = SEP*SEP/(6.D0*TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU))*Host_SimuCtrlParam%EnlageTStepScale
!                end if
!
!                NCActGB = TBasicInfo%NC(p_ACTIVEINGB_STATU) + MultiBox*ImplantedNumEachBox
!                if(NCActGB .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .GT. 0.D0) then
!                    SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(Host_Boxes%m_GrainBoundary%GrainNum)
!
!                    SEP = 6*(SEP**(0.66666666667D0))*Host_Boxes%m_GrainBoundary%GrainNum
!                    SEP = SEP/dble(NCActGB)
!                    SEP = SEP**(0.5D0)
!                    SEP = SEP - 2.D0*TMigStatInfo%RAVA(p_ACTIVEINGB_STATU)
!
!                    TSTEPGB = SEP*SEP/(4.D0*TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU))*Host_SimuCtrlParam%EnlageTStepScale
!                end if
!
!                if((NCActFree .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU) .GT. 0.D0) &
!                   .or. (NCActGB .GT. 0 .AND. TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU) .GT. 0.D0)) then
!                    TheVerifyTime = min(TSTEPFREE,TSTEPGB)
!                else
!                    TheVerifyTime = p_ZeroNCStep
!                end if
!
!            case(mp_FixedTimeStep)
!                TheVerifyTime = Host_SimuCtrlParam%FixedTimeStepValue
!
!            case default
!                write(*,*) "MFPSCUERROR: Unknown strategy to update time step :",Host_SimuCtrlParam%UPDATETSTEPSTRATEGY
!                pause
!                stop
!            end select
!
!        END ASSOCIATE
!
!
!    end function Cal_VerifyTime_Implant


end module MIGCOALE_TIMECTL
