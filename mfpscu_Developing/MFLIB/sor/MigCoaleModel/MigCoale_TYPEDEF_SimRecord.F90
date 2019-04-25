module MIGCOALE_TYPEDEF_SIMRECORD
    use MCMF_TYPEDEF_USUAL
    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
    implicit none

    type,public,extends(SimulationRecord)::MigCoalClusterRecord
        real(kind=KMCDF),private::StartImplantTime = 0.D0
        integer,private::ImplantedEntities = 0
        integer,private::LastRecordImplantNum = 0
        integer,private::NCUT = 0

        integer::RecordNCBeforeSweepOut_Integal(p_NUMBER_OF_STATU)  = 0
        integer,dimension(:,:),allocatable::RecordNCBeforeSweepOut_SingleBox

        real(kind=KMCDF),private::LastUpdateAveSepTime = 0.D0

        integer,private::rescaleCount = 0

        integer::HSizeStatistic_TotalBox = 0
        integer::HSizeStatistic_EachBox = 0
        real(kind=KMCDF),private::LastOutSizeDistTime_IntegralBox = 0.D0
        real(kind=KMCDF),private::LastOutSizeDistTime_EachBox = 0.D0

        contains
        procedure,NON_OVERRIDABLE,public,pass::InitMigCoalClusterRecord
        procedure,non_overridable,public,pass::SetStartImplantTime=>Set_StartImplantTime
        procedure,non_overridable,public,pass::GetStartImplantTime=>Get_StartImplantTime
        procedure,non_overridable,public,pass::GetLastUpdateAveSepTime=>Get_LastUpdateAveSepTime
        procedure,non_overridable,public,pass::SetLastUpdateAveSepTime=>Set_LastUpdateAveSepTime
        procedure,non_overridable,public,pass::IncreaseOneRescaleCount=>Increase_OneRescaleCount
        procedure,non_overridable,public,pass::GetRescaleCount=>Get_RescaleCount
        procedure,NON_OVERRIDABLE,public,pass::AddImplantedEntitiesNum=>Add_ImplantedEntitiesNum
        procedure,NON_OVERRIDABLE,public,pass::GetImplantedEntitiesNum=>Get_ImplantedEntitiesNum
        procedure,NON_OVERRIDABLE,public,pass::SetImplantedEntitiesNum=>Set_ImplantedEntitiesNum
        procedure,NON_OVERRIDABLE,public,pass::SetLastRecordImplantNum=>Set_LastRecordImplantNum
        procedure,NON_OVERRIDABLE,public,pass::GetLastRecordImplantNum=>Get_LastRecordImplantNum
        procedure,NON_OVERRIDABLE,public,pass::SetLastOutSizeDistTime_IntegralBox
        procedure,NON_OVERRIDABLE,public,pass::GetLastOutSizeDistTime_IntegralBox
        procedure,NON_OVERRIDABLE,public,pass::SetLastOutSizeDistTime_EachBox
        procedure,NON_OVERRIDABLE,public,pass::GetLastOutSizeDistTime_EachBox
        procedure,non_overridable,public,pass::WhetherOutSizeDist_IntegralBox
        procedure,non_overridable,public,pass::WhetherOutSizeDist_EachBox
        procedure,non_overridable,pass,private::TheDefProc=>TheDefProc_MigCoalClusterRecord

    end type MigCoalClusterRecord

    private::InitMigCoalClusterRecord
    private::Set_StartImplantTime
    private::Get_StartImplantTime
    private::Set_LastUpdateAveSepTime
    private::Get_LastUpdateAveSepTime
    private::Increase_OneRescaleCount
    private::Get_RescaleCount
    private::Add_ImplantedEntitiesNum
    private::Get_ImplantedEntitiesNum
    private::Set_ImplantedEntitiesNum
    private::Set_LastRecordImplantNum
    private::Get_LastRecordImplantNum
    private::SetLastOutSizeDistTime_IntegralBox
    private::GetLastOutSizeDistTime_IntegralBox
    private::SetLastOutSizeDistTime_EachBox
    private::GetLastOutSizeDistTime_EachBox
    private::WhetherOutSizeDist_IntegralBox
    private::WhetherOutSizeDist_EachBox
    private::TheDefProc_MigCoalClusterRecord

    contains

    !***********type MigCoalClusterRecord *****************
    subroutine InitMigCoalClusterRecord(this,MultiBox,SimuSteps,SimuTimes,SimuPatchs,TimeSection)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        integer,intent(in)::MultiBox
        integer,optional::SimuSteps
        real,optional::SimuTimes
        integer,optional::SimuPatchs
        integer,optional::TimeSection
        !---Local Vars---
        integer::Steps
        real(kind=KMCDF)::Times
        integer::Patchs
        integer::TheTimeSection
        !---Body-- -::
        Steps = 0
        Times = 0.D0
        Patchs = 1
        TheTimeSection = 1

        if(present(SimuSteps)) then
            Steps = SimuSteps
        end if

        if(present(SimuTimes)) then
            Times = SimuTimes
        end if

        if(present(SimuPatchs)) then
            Patchs = SimuPatchs
        end if

        if(present(TimeSection)) then
            TheTimeSection = TimeSection
        end if

        call AllocateArray_Host(this%RecordNCBeforeSweepOut_SingleBox,MultiBox,p_NUMBER_OF_STATU,"RecordNCBeforeSweepOut_SingleBox")
        this%RecordNCBeforeSweepOut_SingleBox = 0
        this%RecordNCBeforeSweepOut_Integal = 0

        this%LastUpdateAveSepTime = 0.D0

        this%rescaleCount = 0

        this%LastOutSizeDistTime_IntegralBox = 0.D0
        this%LastOutSizeDistTime_EachBox = 0.D0

        call this%InitSimulationRecord(SimuSteps=Steps,SimuTimes=Times,SimuPatchs=Patchs,TimeSections=TheTimeSection)

        this%StartImplantTime = 0

        this%ImplantedEntities = 0

        this%LastRecordImplantNum = 0

        this%NCUT = 0

        return
    end subroutine InitMigCoalClusterRecord

    subroutine Set_StartImplantTime(this,TheTime)
        implicit none
        !---Dummy Vars---
        CLass(MigCoalClusterRecord)::this
        real(kind=KMCDF),intent(in)::TheTime
        !---Body---
        this%StartImplantTime = TheTime
        return
    end subroutine Set_StartImplantTime

    function Get_StartImplantTime(this) result(TheTime)
        implicit none
        !---Dummy Vars---
        CLass(MigCoalClusterRecord)::this
        real(kind=KMCDF),intent(out)::TheTime
        !---Body---
        TheTime = this%StartImplantTime
        return
    end function Get_StartImplantTime

    subroutine Add_ImplantedEntitiesNum(this,AddNum)
        implicit none
        !---Dummy Vars---
        CLass(MigCoalClusterRecord)::this
        integer, intent(in)::AddNum
        !---Body---
        if(AddNum .LT. 0) then
            write(*,*) "MCPSCUERROR: The new implanted clusters number is not possible less than 0 :",AddNum
            pause
            stop
        end if

        this%ImplantedEntities = this%ImplantedEntities + AddNum

        return
    end subroutine Add_ImplantedEntitiesNum

    integer function Get_ImplantedEntitiesNum(this)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this

        Get_ImplantedEntitiesNum = this%ImplantedEntities
        return
    end function Get_ImplantedEntitiesNum

    subroutine Set_ImplantedEntitiesNum(this,TheNum)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        integer,intent(in)::TheNum
        !---Body---
        this%ImplantedEntities = TheNum

        return
    end subroutine

    integer function Get_LastRecordImplantNum(this)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this

        Get_LastRecordImplantNum = this%LastRecordImplantNum
        return
    end function Get_LastRecordImplantNum

    subroutine Set_LastRecordImplantNum(this,TheNum)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        integer,intent(in)::TheNum
        !---Body---
        this%LastRecordImplantNum = TheNum

        return
    end subroutine Set_LastRecordImplantNum

    integer function Get_NCUT(this)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this

        Get_NCUT = this%NCUT
        return
    end function Get_NCUT

    subroutine Set_NCUT(this,TheNCUT)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        integer,intent(in)::TheNCUT
        !---Body---
        this%NCUT = TheNCUT

        return
    end subroutine Set_NCUT


    subroutine SetLastOutSizeDistTime_IntegralBox(this,TheTime)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        real(kind=KMCDF),intent(in)::TheTime
        !---Body---
        this%LastOutSizeDistTime_IntegralBox = TheTime

        return
    end subroutine SetLastOutSizeDistTime_IntegralBox

    function GetLastOutSizeDistTime_IntegralBox(this) result(TheTime)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        real(kind=KMCDF),intent(out)::TheTime
        !---Body---
        TheTime = this%LastOutSizeDistTime_IntegralBox

        return
    end function GetLastOutSizeDistTime_IntegralBox

    subroutine SetLastOutSizeDistTime_EachBox(this,TheTime)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        real(kind=KMCDF),intent(in)::TheTime
        !---Body---
        this%LastOutSizeDistTime_EachBox = TheTime

        return
    end subroutine SetLastOutSizeDistTime_EachBox

    function GetLastOutSizeDistTime_EachBox(this) result(TheTime)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        real(kind=KMCDF),intent(out)::TheTime
        !---Body---
        TheTime = this%LastOutSizeDistTime_EachBox

        return
    end function GetLastOutSizeDistTime_EachBox

    !******************************************
    function WhetherOutSizeDist_IntegralBox(this,Host_SimuCtrlParam) result(TheResult)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        logical,intent(inout)::TheResult
        !---Body---
        TheResult = .false.

        if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalSteps) then
            if((this%GetSimuSteps() - this%GetLastOutSizeDistTime_IntegralBox()) .GE. Host_SimuCtrlParam%OutPutSCValue_IntegralBox) then
                TheResult = .true.
            end if

        else if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalRealTime) then
            if((this%GetSimuTimes() - this%GetLastOutSizeDistTime_IntegralBox()) .GE. Host_SimuCtrlParam%OutPutSCValue_IntegralBox) then
                TheResult = .true.
            end if

        else if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalTimeMagnification) then
            if((this%GetSimuTimes()/Host_SimuCtrlParam%OutPutSCValue_IntegralBox) .GE. this%GetLastOutSizeDistTime_IntegralBox()) then
                TheResult = .true.
            end if
        end if

        return
    end function WhetherOutSizeDist_IntegralBox


    !******************************************
    function WhetherOutSizeDist_EachBox(this,Host_SimuCtrlParam) result(TheResult)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoalClusterRecord)::this
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        logical,intent(inout)::TheResult
        !---Body---
        TheResult = .false.

        if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalSteps) then
            if((this%GetSimuSteps() - this%GetLastOutSizeDistTime_EachBox()) .GE. Host_SimuCtrlParam%OutPutSCValue_EachBox) then
                TheResult = .true.
            end if

        else if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalRealTime) then
            if((this%GetSimuTimes() - this%GetLastOutSizeDistTime_EachBox()) .GE. Host_SimuCtrlParam%OutPutSCValue_EachBox) then
                TheResult = .true.
            end if

        else if(Host_SimuCtrlParam%OutPutSCFlag .eq. mp_OutTimeFlag_ByIntervalTimeMagnification) then
            if((this%GetSimuTimes()/Host_SimuCtrlParam%OutPutSCValue_EachBox) .GE. this%GetLastOutSizeDistTime_EachBox()) then
                TheResult = .true.
            end if
        end if

        return
    end function WhetherOutSizeDist_EachBox



    !*******************************************
    subroutine TheDefProc_MigCoalClusterRecord(this)
        implicit none
        CLASS(MigCoalClusterRecord)::this
    end subroutine


    !**************************************************************
    subroutine Increase_OneRescaleCount(this)
        implicit none
        Class(MigCoalClusterRecord)::this

        this%rescaleCount = this%rescaleCount + 1
        return
    end subroutine Increase_OneRescaleCount

    !**************************************************************
    function Get_RescaleCount(this) result(rescaleCount)
        implicit none
        Class(MigCoalClusterRecord)::this
        integer::rescaleCount

        rescaleCount = this%rescaleCount
        return
    end function Get_RescaleCount

    !*******************************************************
    subroutine Set_LastUpdateAveSepTime(this,TheTime)
        implicit none
        !---Dummy Vars---
        CLass(MigCoalClusterRecord)::this
        real(kind=KMCDF),intent(in)::TheTime
        !---Body---
        this%LastUpdateAveSepTime = TheTime
        return
    end subroutine Set_LastUpdateAveSepTime

    function Get_LastUpdateAveSepTime(this) result(TheTime)
        implicit none
        !---Dummy Vars---
        CLass(MigCoalClusterRecord)::this
        real(kind=KMCDF),intent(out)::TheTime
        !---Body---
        TheTime = this%LastUpdateAveSepTime
        return
    end function Get_LastUpdateAveSepTime

end module
