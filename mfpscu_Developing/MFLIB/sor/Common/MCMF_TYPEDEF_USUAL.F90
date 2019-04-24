module MCMF_TYPEDEF_USUAL

    USE MCMF_CONSTANTS
    USE MCMF_UTILITIES

    implicit none

    !------------------
    type,public::RunningRecord

        real::LastRecordOutprofileTime = 0.D0

        logical,private::stopRunningFlag = .false.
        character(LEN = 12)::Start_Clock(3), End_Clock(3)
        integer::Start_DateTime(8),End_DateTime(8)

        contains
        procedure,NON_OVERRIDABLE,public,pass::InitRunningRecord
        procedure,NON_OVERRIDABLE,public,pass::StopRunning=>Stop_Running
        procedure,NON_OVERRIDABLE,public,pass::IsStoppedRunning=>Is_StoppedRunning

    end type RunningRecord


    type,abstract,public::SimulationRecord
        type(RunningRecord)::Running_Record

        integer,private::SimulaitonSteps = 0
        integer,private::SimulationPatch = 1
        real(kind=KMCDF),private::SimulationTimes = 0.D0
        integer,private::TimeSections = 1

        real(kind=KMCDF),private::LastUpdateStatisTime = 0.D0

        real(kind=KMCDF),private::LastUpdateNLTime = 0.D0
        integer,private::LastUpdateNL_NC0 = 1

        real(kind=KMCDF),private::LastRecordOutConfigTime = 0.D0
        integer,private::OutPutIndex = 0

        contains
        procedure,NON_OVERRIDABLE,public,pass::InitSimulationRecord
        procedure,NON_OVERRIDABLE,public,pass::SetSimuSteps=>Set_SimuSteps
        procedure,NON_OVERRIDABLE,public,pass::GetSimuSteps=>Get_SimuSteps
        procedure,NON_OVERRIDABLE,public,pass::IncreaseOneSimuStep=>Increase_OneSimuStep

        procedure,NON_OVERRIDABLE,public,pass::SetSimuTimes=>Set_SimuTimes
        procedure,NON_OVERRIDABLE,public,pass::GetSimuTimes=>Get_SimuTimes
        procedure,NON_OVERRIDABLE,public,pass::AddSimuTimes=>Add_SimuTimes

        procedure,NON_OVERRIDABLE,public,pass::SetSimuPatch=>Set_SimuPatch
        procedure,NON_OVERRIDABLE,public,pass::GetSimuPatch=>Get_SimuPatch

        procedure,NON_OVERRIDABLE,public,pass::SetTimeSections=>Set_TimeSections
        procedure,NON_OVERRIDABLE,public,pass::GetTimeSections=>Get_TimeSections
        procedure,NON_OVERRIDABLE,public,pass::IncreaseOneTimeSection=>Increase_OneTimeSection

        procedure,NON_OVERRIDABLE,public,pass::GetLastUpdateStatisTime=>Get_LastUpdateStatisTime
        procedure,NON_OVERRIDABLE,public,pass::SetLastUpdateStatisTime=>Set_LastUpdateStatisTime

        procedure,NON_OVERRIDABLE,public,pass::GetLastUpdateNLTime=>Get_LastUpdateNLTime
        procedure,NON_OVERRIDABLE,public,pass::SetLastUpdateNLTime=>Set_LastUpdateNLTime

        procedure,NON_OVERRIDABLE,public,pass::GetLastUpdateNLNC0=>Get_LastUpdateNLNC0
        procedure,NON_OVERRIDABLE,public,pass::SetLastUpdateNLNC0=>Set_LastUpdateNLNC0

        procedure,NON_OVERRIDABLE,public,pass::GetLastRecordOutConfigTime=>Get_LastRecordOutConfigTime
        procedure,NON_OVERRIDABLE,public,pass::SetLastRecordOutConfigTime=>Set_LastRecordOutConfigTime

        procedure,NON_OVERRIDABLE,public,pass::GetOutPutIndex=>Get_OutPutIndex
        procedure,NON_OVERRIDABLE,public,pass::SetOutPutIndex=>Set_OutPutIndex
        procedure,NON_OVERRIDABLE,public,pass::IncreaseOneOutPutIndex=>Increase_OneOutPutIndex

        !---abstract method---
        procedure(DefProc),pass,deferred,private::TheDefProc

    end type SimulationRecord

    private::InitRunningRecord
    private::Stop_Running
    private::Is_StoppedRunning
    private::InitSimulationRecord
    private::Set_SimuSteps
    private::Get_SimuSteps
    private::Increase_OneSimuStep
    private::Set_SimuTimes
    private::Get_SimuTimes
    private::Add_SimuTimes
    private::Set_SimuPatch
    private::Get_SimuPatch
    private::Set_TimeSections
    private::Get_TimeSections
    private::Increase_OneTimeSection
    private::Get_LastUpdateStatisTime
    private::Set_LastUpdateStatisTime
    private::Get_LastUpdateNLTime
    private::Set_LastUpdateNLTime
    private::Get_LastUpdateNLNC0
    private::Set_LastUpdateNLNC0
    private::Get_LastRecordOutConfigTime
    private::Set_LastRecordOutConfigTime
    private::Get_OutPutIndex
    private::Set_OutPutIndex
    private::Increase_OneOutPutIndex

    contains

    !************type RunningRecord*******************
    subroutine InitRunningRecord(this)
        implicit none
        CLASS(RunningRecord)::this

        this%LastRecordOutprofileTime = 0.D0

        this%stopRunningFlag= .false.
        this%Start_Clock = ''
        this%End_Clock = ''
        this%Start_DateTime = 0
        this%End_DateTime = 0

    end subroutine


    subroutine Stop_Running(this)
        implicit none
        CLASS(RunningRecord)::this

        this%stopRunningFlag = .true.
        return
    end subroutine

    logical function Is_StoppedRunning(this)
        implicit none
        CLASS(RunningRecord)::this

        Is_StoppedRunning = this%stopRunningFlag
        return
    end function Is_StoppedRunning

    !****abstract type SimulationRecord*****************
    subroutine InitSimulationRecord(this,SimuSteps,SimuTimes,SimuPatchs,TimeSections)
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        integer,optional::SimuSteps
        real(kind=KMCDF),optional::SimuTimes
        integer,optional::SimuPatchs
        integer,optional::TimeSections
        !---Local Vars---
        integer::Steps
        real(kind=KMCDF)::Times
        integer::Patchs
        integer::TheTimeSection
        !---Body---
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

        if(present(TimeSections)) then
            TheTimeSection = TimeSections
        end if

        call this%SetSimuSteps(Steps)

        call this%SetSimuTimes(Times)

        call this%SetSimuPatch(Patchs)

        call this%SetTimeSections(TheTimeSection)

        call this%Running_Record%InitRunningRecord()

        this%LastUpdateStatisTime = 0.D0

        this%LastUpdateNLTime = 0.D0
        this%LastUpdateNL_NC0 = 1

        this%LastRecordOutConfigTime = 0.D0
        this%OutPutIndex = 0

        return
    end subroutine InitSimulationRecord

    !***************************************************
    subroutine Set_SimuSteps(this,Steps)
            implicit none
            CLASS(SimulationRecord)::this
            integer::Steps
            this%SimulaitonSteps = Steps
    end subroutine

    !***************************************************
    integer function Get_SimuSteps(this)
            implicit none
            CLASS(SimulationRecord)::this
            Get_SimuSteps = this%SimulaitonSteps
    end function


    !***************************************************
    subroutine Increase_OneSimuStep(this)
        implicit none
        CLASS(SimulationRecord)::this

        this%SimulaitonSteps = this%SimulaitonSteps + 1

        return
    end subroutine Increase_OneSimuStep

    !***************************************************
    subroutine Set_SimuTimes(this,Times)
            implicit none
            CLASS(SimulationRecord)::this
            real(kind=KMCDF)::Times
            this%SimulationTimes = Times
    end subroutine

    !***************************************************
    real(kind=KMCDF) function Get_SimuTimes(this)
            implicit none
            CLASS(SimulationRecord)::this
            Get_SimuTimes = this%SimulationTimes
    end function

    !*****************************************************
    subroutine Add_SimuTimes(this,increaseTime)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        real(kind=KMCDF),intent(in)::increaseTime
        !---Body---
        this%SimulationTimes = this%SimulationTimes + increaseTime

        return
    end subroutine Add_SimuTimes

   !********************************************************
   subroutine Set_SimuPatch(this,SimPath)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        integer,intent(in)::SimPath
        !---Body---

        this%SimulationPatch = SimPath

        return
   end subroutine Set_SimuPatch

   !********************************************************
   function Get_SimuPatch(this) result(SimPath)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        integer,intent(out)::SimPath
        !---Body---

        SimPath = this%SimulationPatch

        return
   end function Get_SimuPatch

    !***************************************************
    subroutine Set_TimeSections(this,TimeSection)
        implicit none
        CLASS(SimulationRecord)::this
        integer::TimeSection
        this%TimeSections = TimeSection
    end subroutine Set_TimeSections

    !***************************************************
    integer function Get_TimeSections(this)
        implicit none
        CLASS(SimulationRecord)::this
        Get_TimeSections = this%TimeSections
    end function Get_TimeSections


    !***************************************************
    subroutine Increase_OneTimeSection(this)
        implicit none
        CLASS(SimulationRecord)::this

        this%TimeSections = this%TimeSections + 1

        return
    end subroutine Increase_OneTimeSection

    !****************************************************
    real(kind=KMCDF) function Get_LastUpdateStatisTime(this)
        implicit none
        CLASS(SimulationRecord)::this

        Get_LastUpdateStatisTime = this%LastUpdateStatisTime
        return
    end function Get_LastUpdateStatisTime

    !***************************************************
    subroutine Set_LastUpdateStatisTime(this,TIME)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        real(kind=KMCDF),intent(in)::TIME
        !---Body---
        this%LastUpdateStatisTime = TIME

        return
    end subroutine Set_LastUpdateStatisTime

    !****************************************************
    function Get_LastUpdateNLTime(this) result(TheTime)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        real(kind=KMCDF),intent(out)::TheTime
        !---Body---
        TheTime = this%LastUpdateNLTime
        return
    end function Get_LastUpdateNLTime

    !****************************************************
    subroutine Set_LastUpdateNLTime(this,TheTime)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        real(kind=KMCDF),intent(in)::TheTime
        !---Body---
        this%LastUpdateNLTime = TheTime
        return
    end subroutine Set_LastUpdateNLTime

    !****************************************************
    function Get_LastUpdateNLNC0(this) result(NC0)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        integer,intent(out)::NC0
        !---Body---
        NC0 = this%LastUpdateNL_NC0
        return
    end function Get_LastUpdateNLNC0

    !****************************************************
    subroutine Set_LastUpdateNLNC0(this,NC0)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        integer,intent(in)::NC0
        !---Body---
        this%LastUpdateNL_NC0 = NC0

        return
    end subroutine Set_LastUpdateNLNC0

    !****************************************************
    real function Get_LastRecordOutConfigTime(this)
        implicit none
        CLASS(SimulationRecord)::this

        Get_LastRecordOutConfigTime = this%LastRecordOutConfigTime
        return
    end function Get_LastRecordOutConfigTime

    !****************************************************
    subroutine Set_LastRecordOutConfigTime(this,TIME)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        real(kind=KMCDF),intent(in)::TIME
        !---Body---
        this%LastRecordOutConfigTime = TIME

        return
    end subroutine Set_LastRecordOutConfigTime

    !****************************************************
    integer function Get_OutPutIndex(this)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this

        Get_OutPutIndex = this%OutPutIndex

        return
    end function Get_OutPutIndex

    !****************************************************
    subroutine Set_OutPutIndex(this,OutIndex)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this
        integer,intent(in)::OutIndex

        this%OutPutIndex = OutIndex

        return
    end subroutine Set_OutPutIndex

    !****************************************************
    subroutine Increase_OneOutPutIndex(this)
        implicit none
        !---Dummy Vars---
        CLASS(SimulationRecord)::this

        this%OutPutIndex = this%OutPutIndex + 1
        return
    end subroutine
    !*****The declare for abstract method***************
    !*****Note: this is not same with the "Fortran 95/2003 For Scientists and Engineers, Third Edition (chinese version)(P668)"
    !*****Because the abstract useage way in this book is not depended on PGFORTRAN compiler. The  following is based on our test
    !*****and verified to suit for pgfortran.
    subroutine DefProc(this)
            implicit none
            CLASS(SimulationRecord)::this
    end subroutine

end module MCMF_TYPEDEF_USUAL
