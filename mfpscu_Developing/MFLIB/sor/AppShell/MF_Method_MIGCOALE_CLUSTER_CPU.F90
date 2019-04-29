module MF_Method_MIGCOALE_CLUSTER_CPU
    use MCMF_CONSTANTS
    use MIGCOALE_TYPEDEF_SIMRECORD
    use MIGCOALE_TYPEDEF_STATISTICINFO
    use MIGCOALE_ADDONDATA_HOST
    use MIGCOALE_IMPLANTATION
    use NUCLEATION_SPACEDIST
    implicit none

    integer, parameter, private::p_ClusterIniConfig_Simple = 0
    integer, parameter, private::p_ClusterIniConfig_SpecialDistFromFile = 1
    integer, parameter, private::p_ClusterIniConfig_SpecialDistFromExteFunc = 2

    character(len=25),  private:: m_INIFSTARTFLAG = "&INITINPUTF"

    !--------------------
    type,public::InitBoxSimCfg

        integer::InitType = -1

        integer::InitCKind = 200

        integer::CKind_Range_From = 1
        integer::CKind_Range_To = 200
        integer::CKind_Range_Interval = 1

        integer::InitNNodes = 1

        character*256::InitCfgFileName = ""

        character*10::Elemets(p_ATOMS_GROUPS_NUMBER) = ""

        real(kind=KMCDF)::NAINI = 0.D0

        real(kind=KMCDF)::NASDINI = 0.D0

        real(kind=KMCDF)::NACUT(2) = 0.D0

        real(kind=KMCDF)::CompositWeight(p_ATOMS_GROUPS_NUMBER) = 0.D0

        integer::InitDepthDistType = -1

        real(kind=KMCDF)::InitConcentrate = 0.D0

        real(kind=KMCDF),dimension(:),allocatable::LayerThick

        real(kind=KMCDF),dimension(:),allocatable::PNCLayers

        real(kind=KMCDF)::DepthINI = 0.D0

        real(kind=KMCDF)::DepthSDINI = 0.D0

        contains
        procedure,non_overridable,public,pass::CopyInitBoxSimCfgFromOther
        procedure,non_overridable,public,pass::Clean_InitBoxSimCfg
        Generic::Assignment(=)=>CopyInitBoxSimCfgFromOther
        Final::CleanInitBoxSimCfg
    end type

    !--------------------
    type,public::InitBoxSimCfgList
        type(InitBoxSimCfg)::TheValue
        integer::ListCount = 0
        type(InitBoxSimCfgList),pointer::next=>null()

        contains
        procedure,non_overridable,public,pass::AppendOne_InintSimBoxCfg
        procedure,non_overridable,public,pass::GetList_Count=>GetInitBoxSimCfgList_Count
        procedure,non_overridable,public,pass::Clean_InitBoxSimCfgList
        Final::CleanInitBoxSimCfgList
    end type

    type(InitBoxSimCfgList),target::m_InitBoxSimCfgList
    type(ImplantList)::m_ImplantList
    type(MigCoalClusterRecord)::m_MigCoalClusterRecord
    type(MigCoaleStatInfoWrap),private::m_MigCoaleStatInfoWrap
    !type(MigCoale_GVarsDev),private::m_MigCoale_GVarsDev

    private::CopyInitBoxSimCfgFromOther
    private::Clean_InitBoxSimCfg
    private::CleanInitBoxSimCfg
    private::AppendOne_InintSimBoxCfg
    private::GetInitBoxSimCfgList_Count
    private::Clean_InitBoxSimCfgList
    private::CleanInitBoxSimCfgList

    contains

    !*****************************************************************
    subroutine For_One_Test(Host_SimBoxes,Host_SimuCtrlParam,JobIndex)
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

                call resolveAddOnData(Host_SimBoxes,PSimCtrlParam)

                PImplantSection=>m_ImplantList%Get_P(PSimCtrlParam%ImplantSectID)

                call For_One_TimeSect(Host_SimBoxes,PSimCtrlParam,PImplantSection,m_MigCoaleStatInfoWrap,m_MigCoalClusterRecord)

                call m_MigCoalClusterRecord%SetLastRecordImplantNum(m_MigCoalClusterRecord%GetImplantedEntitiesNum())

                PSimCtrlParam=>PSimCtrlParam%next

                if(.not. associated(PSimCtrlParam)) then
                    exit
                end if

                call m_MigCoalClusterRecord%IncreaseOneTimeSection()

            END DO

        END DO

        return

    end subroutine For_One_Test

    !****************************************************************
    subroutine For_One_TimeSect(Host_SimBoxes,Host_SimuCtrlParam,TheImplantSection,TheMigCoaleStatInfoWrap,Record)
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

        call InitSimu_SpaceDist(Host_SimBoxes,Host_SimuCtrlParam)

        call NucleationSimu_SpaceDist(Host_SimBoxes,Host_SimuCtrlParam)


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
    end subroutine For_One_TimeSect
!
!    !*********************************************************************
!    subroutine For_One_Step(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,TheImplantSection,TheMigCoaleStatInfoWrap,Record,TSTEP)
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
!    end subroutine For_One_Step

    !*************************************************************
    subroutine CopyInitBoxSimCfgFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(InitBoxSimCfg),intent(out)::this
        type(InitBoxSimCfg),intent(in)::other
        !---Local Vars---
        integer::I
        !---Body---
        this%InitType = other%InitType

        this%InitCKind = other%InitCKind

        this%CKind_Range_From = other%CKind_Range_From
        this%CKind_Range_To = other%CKind_Range_To
        this%CKind_Range_Interval = other%CKind_Range_Interval

        this%InitNNodes = other%InitNNodes

        this%InitCfgFileName = other%InitCfgFileName

        DO I = 1,p_ATOMS_GROUPS_NUMBER
            this%Elemets(I) = other%Elemets(I)
            this%CompositWeight(I) = other%CompositWeight(I)
        END DO

        this%NAINI = other%NAINI

        this%NASDINI = other%NASDINI

        this%NACUT = other%NACUT

        this%InitDepthDistType = other%InitDepthDistType

        this%InitConcentrate = other%InitConcentrate

        this%DepthINI = other%DepthINI

        this%DepthSDINI = other%DepthSDINI

        call DeAllocateArray_Host(this%LayerThick,"LayerThick")
        call AllocateArray_Host(this%LayerThick,size(other%LayerThick),"LayerThick")
        this%LayerThick = other%LayerThick

        call DeAllocateArray_Host(this%PNCLayers,"PNCLayers")
        call AllocateArray_Host(this%PNCLayers,size(other%PNCLayers),"PNCLayers")
        this%PNCLayers = other%PNCLayers

        return
    end subroutine CopyInitBoxSimCfgFromOther

    !*************************************************************
    subroutine Clean_InitBoxSimCfg(this)
        implicit none
        !---Dummy Vars---
        CLASS(InitBoxSimCfg),intent(out)::this
        !---Body---
        this%InitType = -1

        this%InitCKind = 200

        this%CKind_Range_From = 1
        this%CKind_Range_To = 200
        this%CKind_Range_Interval = 1

        this%InitNNodes = 1

        this%InitCfgFileName = ""

        this%Elemets = ""

        this%NAINI = 0.D0

        this%NASDINI = 0.D0

        this%NACUT = 0.D0

        this%CompositWeight = 0.D0

        this%InitDepthDistType = -1

        this%InitConcentrate = 0.D0

        call DeAllocateArray_Host(this%LayerThick,"LayerThick")

        call DeAllocateArray_Host(this%PNCLayers,"PNCLayers")

        this%DepthINI = 0.D0

        this%DepthSDINI = 0.D0

        return
    end subroutine

    !*************************************************
    subroutine CleanInitBoxSimCfg(this)
        implicit none
        !---Dummy Vars---
        type(InitBoxSimCfg)::this
        !---Body---
        call this%Clean_InitBoxSimCfg()

        return
    end subroutine CleanInitBoxSimCfg

    !*************************************************************
    subroutine AppendOne_InintSimBoxCfg(this,newOne)
        implicit none
        !---Dummy Vars---
        CLASS(InitBoxSimCfgList),target::this
        type(InitBoxSimCfg)::newOne
        !---Local Vars---
        type(InitBoxSimCfgList),pointer::cursor=>null(),cursorP=>null()
        !---Body---
        if(this%GetList_Count() .LE. 0) then
            this%ListCount = 1
            this%TheValue = newOne
        else
            cursor=>this%next
            cursorP=>this

            DO while(associated(cursor))
                cursor=>cursor%next
                cursorP=>cursorP%next

            END DO

            this%ListCount = this%ListCount + 1

            allocate(cursor)
            NUllify(cursor%next)
            ! The assignment(=) had been overrided
            cursor%TheValue = newOne
            cursorP%next=>cursor
        end if

        Nullify(cursorP)
        cursorP=>null()
        Nullify(cursor)
        cursor=>null()
        return
    end subroutine AppendOne_InintSimBoxCfg

    !**************************************
    integer function GetInitBoxSimCfgList_Count(this)
        implicit none
        !---Dummy Vars---
        CLASS(InitBoxSimCfgList)::this
        !---Body---
        GetInitBoxSimCfgList_Count = this%ListCount

        return
    end function

    !**************************************
    subroutine Clean_InitBoxSimCfgList(this)
        implicit none
        !---Dummy Vars---
        CLASS(InitBoxSimCfgList),target::this
        !---Local Vars---
        type(InitBoxSimCfgList),pointer::cursor=>null()
        type(InitBoxSimCfgList),pointer::next=>null()
        !---Body---
        cursor=>this%next

        call this%TheValue%Clean_InitBoxSimCfg()

        DO While(associated(cursor))
            next=>cursor%next
            call cursor%TheValue%Clean_InitBoxSimCfg()
            deallocate(cursor)
            Nullify(cursor)
            cursor=>next
        END DO


        this%next=>null()

        this%ListCount = 0

        Nullify(cursor)
        Nullify(next)
        cursor=>null()
        next=>null()

        return
    end subroutine Clean_InitBoxSimCfgList

    !************************************
    subroutine CleanInitBoxSimCfgList(this)
        implicit none
        !---Dummy Vars---
        type(InitBoxSimCfgList)::this
        !---Body---

        call this%Clean_InitBoxSimCfgList()

        return
    end subroutine CleanInitBoxSimCfgList

    !*****************************************************************
    subroutine InitSimulationBoxesConfig(SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,TheMigCoaleStatInfoWrap,Record)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfgList)::InitBoxCfgList
        type(MigCoaleStatInfoWrap)::TheMigCoaleStatInfoWrap
        type(MigCoalClusterRecord)::Record
        !---Local Vars---
        logical::existed
        integer::hFile
        character*256::STR
        character*32::KEYWORD
        character*32::STRTMP(10)
        integer::LINE
        integer::N
        character*256::path
        character*18,dimension(:),allocatable::CRMin
        character*18,dimension(:),allocatable::CRMax
        character*18,dimension(:),allocatable::CCNum
        character*18,dimension(:),allocatable::CAcumNum
        integer::I
        integer::length,trueLength
        !---Body---
        existed = .false.

        LINE = 0

        call InitBoxCfgList%Clean_InitBoxSimCfgList()

        INQUIRE(File=SimBoxes%IniConfig(1:LENTRIM(SimBoxes%IniConfig)),exist=existed)

        if(.not. existed) then
            write(*,*) "MFPSCUERROR: The box initial file do not existed!"
            write(*,*) SimBoxes%IniConfig
            pause
            stop
        end if

        hFile = openExistedFile(SimBoxes%IniConfig)

        call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)
        if(KEYWORD(1:LENTRIM(KEYWORD)) .ne. m_INIFSTARTFLAG) then
            write(*,*) "MFPSCUERROR: The Start Flag of Init box Parameters is Illegal: ",KEYWORD(1:LENTRIM(KEYWORD))
            pause
            stop
        end if

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE, "!", *100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDINITINPUTF")
                    exit
                case("&SIMPLEDISTSUBCTL")
                    InitBoxCfgList%TheValue%InitType = p_ClusterIniConfig_Simple
                    call ReadInitSimulationBoxesConfig_Simple(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,LINE)
                case("&FILEDISTSUBCTL")
                    InitBoxCfgList%TheValue%InitType = p_ClusterIniConfig_SpecialDistFromFile
                    call ReadInitSimulationBoxesConfig_SpecialDistFromFile(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,LINE)
                case("&FUNCDISTSUBCTL")
                    InitBoxCfgList%TheValue%InitType = p_ClusterIniConfig_SpecialDistFromExteFunc
                    call ReadInitSimulationBoxesConfig_SpecialDistFromExteFunc(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,LINE)
                case default
                    write(*,*) "MFPSCUERROR: Unknown strategy to initial configuration."
                    write(*,*) "At LINE: ",LINE
                    write(*,*) KEYWORD(1:LENTRIM(KEYWORD))
                    pause
                    stop
            end select
        END DO

        call CheckInitSimulationBoxesConfig(SimBoxes,Host_SimuCtrlParam,InitBoxCfgList)

        call DOInitSimulationBoxesConfig(SimBoxes,Host_SimuCtrlParam,Record,InitBoxCfgList)

!        call GetBoxesMigCoaleStat_Used_CPU(SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Used)
!
!        call GetBoxesMigCoaleStat_Expd_CPU(SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd)
!
!        !---The final one is used for the total boxes
!        allocate(CRMin(p_NUMBER_OF_STATU))
!        allocate(CRMax(p_NUMBER_OF_STATU))
!        allocate(CCNum(p_NUMBER_OF_STATU))
!        allocate(CAcumNum(p_NUMBER_OF_STATU))
!
!        DO I = 1,p_NUMBER_OF_STATU
!            CRMin(I) = "RMIN_"//trim(p_CStatu(I))
!            length = len(CRMin(I))
!            trueLength = LENTRIM(CRMin(I))
!            CRMin(I)(length-trueLength+1:length) = CRMin(I)(1:trueLength)
!            CRMin(I)(1:length-trueLength) = ""
!
!            CRMax(I) = "RMAX_"//trim(p_CStatu(I))
!            length = len(CRMax(I))
!            trueLength = LENTRIM(CRMax(I))
!            CRMax(I)(length-trueLength+1:length) = CRMax(I)(1:trueLength)
!            CRMax(I)(1:length-trueLength) = ""
!
!            CCNum(I) = "NUM_"//trim(p_CStatu(I))
!            length = len(CCNum(I))
!            trueLength = LENTRIM(CCNum(I))
!            CCNum(I)(length-trueLength+1:length) = CCNum(I)(1:trueLength)
!            CCNum(I)(1:length-trueLength) = ""
!
!            CAcumNum(I) = "ACUM_"//trim(p_CStatu(I))
!            length = len(CAcumNum(I))
!            trueLength = LENTRIM(CAcumNum(I))
!            CAcumNum(I)(length-trueLength+1:length) = CAcumNum(I)(1:trueLength)
!            CAcumNum(I)(1:length-trueLength) = ""
!
!        END DO
!
!        path = Host_SimuCtrlParam%OutFilePath(1:LENTRIM(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"RTStatistic_EachBox_.out"
!
!        Record%HSizeStatistic_EachBox = CreateNewFile(path)
!
!        write(Record%HSizeStatistic_EachBox, fmt="(130(A20,1x))") "Step",                              &
!                                                                  "IBox",                              &
!                                                                  "Time(s)",                           &
!                                                                  "NACTClusters",                      &
!                                                                  "TotalCluster",                      &
!                                                                  "TotalAtoms",                        &
!                                                                  CCNum(1:p_NUMBER_OF_STATU),          &
!                                                                  CRMin(1:p_NUMBER_OF_STATU),          &
!                                                                  CRMax(1:p_NUMBER_OF_STATU),          &
!                                                                  "RAVE",                              &
!                                                                  "NAVA",                              &
!                                                                  "Concentrate",                       &
!                                                                  "TotalImplant",                      &
!                                                                  CAcumNum(1:p_NUMBER_OF_STATU)
!
!        FLUSH(Record%HSizeStatistic_EachBox)
!
!        path = Host_SimuCtrlParam%OutFilePath(1:LENTRIM(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"RTStatistic_TotalBox_.out"
!
!        Record%HSizeStatistic_TotalBox = CreateNewFile(path)
!
!        write(Record%HSizeStatistic_TotalBox, fmt="(130(A20,1x))")  "Step",                              &
!                                                                    "Time(s)",                           &
!                                                                    "NACTClusters",                      &
!                                                                    "TotalCluster",                      &
!                                                                    "TotalAtoms",                        &
!                                                                    CCNum(1:p_NUMBER_OF_STATU),          &
!                                                                    CRMin(1:p_NUMBER_OF_STATU),          &
!                                                                    CRMax(1:p_NUMBER_OF_STATU),          &
!                                                                    "RAVE",                              &
!                                                                    "NAVA",                              &
!                                                                    "Concentrate",                       &
!                                                                    "TotalImplant",                      &
!                                                                    CAcumNum(1:p_NUMBER_OF_STATU)
!
!
!        FLUSH(Record%HSizeStatistic_TotalBox)
!
!        deallocate(CRMin)
!        deallocate(CRMax)
!        deallocate(CCNum)
!        deallocate(CAcumNum)
!
!        call SimBoxes%PutoutCfg(Host_SimuCtrlParam,Record)
!
!        call PutOut_Instance_Statistic_IntegralBox(SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Used,Record,Model=0)
!        call PutOut_Instance_Statistic_EachBox(SimBoxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Used,Record)

        return

        100 write(*,*) "MFPSCUERROR : Load init config file"//SimBoxes%IniConfig(1:LENTRIM(SimBoxes%IniConfig))//"failed !"
            write(*,*) "At line :",LINE
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine InitSimulationBoxesConfig

    !*****************************************************************
    subroutine CheckInitSimulationBoxesConfig(SimBoxes,Host_SimuCtrlParam,InitBoxCfgList)
        !---Dummy Vars---
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfgList),target::InitBoxCfgList
        !---Local Vars---
        type(InitBoxSimCfgList),pointer::PInitBoxCfgList=>null()
        type(InitBoxSimCfgList),pointer::NPInitBoxCfgList=>null()
        !---Local Vars---
        integer::I
        integer::AIndex
        integer::BIndex
        integer::GroupIndex
        character*10::Symbol1
        character*10::Symbol2
        !---Body---
        PInitBoxCfgList=>InitBoxCfgList

        GroupIndex = 0

        Do While(associated(PInitBoxCfgList))
            GroupIndex = GroupIndex + 1
            AIndex = 0
            BIndex = 0

            DO I = 1,p_ATOMS_GROUPS_NUMBER
                if(PInitBoxCfgList%TheValue%CompositWeight(I) .GT. 0.D0) then
                    if(AIndex .ne. 0) then
                        call SimBoxes%Atoms_list%GetSymbolByIndex(AIndex,Symbol1)
                        call SimBoxes%Atoms_list%GetSymbolByIndex(I,Symbol2)
                        write(*,*) "MCPSCUERROR: The current application only support one element ."
                        write(*,*) "By you seems specialed element",Symbol1
                        write(*,*) "And ",Symbol2, " simultaneously."
                        write(*,*) "At group : ",GroupIndex
                        pause
                        stop
                    end if
                    AIndex = I
                end if
            END DO

            NPInitBoxCfgList=>PInitBoxCfgList%next
            if(associated(NPInitBoxCfgList)) then
                DO I = 1,p_ATOMS_GROUPS_NUMBER
                    if(NPInitBoxCfgList%TheValue%CompositWeight(I) .GT. 0.D0) then
                        if(BIndex .ne. 0) then
                            call SimBoxes%Atoms_list%GetSymbolByIndex(BIndex,Symbol1)
                            call SimBoxes%Atoms_list%GetSymbolByIndex(I,Symbol2)

                            write(*,*) "MCPSCUERROR: The current application only support one element ."
                            write(*,*) "By you seems specialed element", Symbol1
                            write(*,*) "And ",Symbol2, " simultaneously."
                            write(*,*) "At group : ",GroupIndex + 1
                            pause
                            stop
                        end if
                        BIndex = I
                    end if
                END DO

                if(AIndex .ne. BIndex) then
                    call SimBoxes%Atoms_list%GetSymbolByIndex(AIndex,Symbol1)
                    call SimBoxes%Atoms_list%GetSymbolByIndex(BIndex,Symbol2)

                    write(*,*) "MCPSCUERROR: The current application only support one element ."
                    write(*,*) "By you seems specialed element", Symbol1 ," in group ",GroupIndex
                    write(*,*) "And ",Symbol2," in group ",GroupIndex+1, " simultaneously."
                    pause
                    stop
                end if

            end if

            PInitBoxCfgList=>PInitBoxCfgList%next

        End Do

        return
    end subroutine

    !*****************************************************************
    subroutine ReadInitSimulationBoxesConfig_Simple(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,LINE)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfgList),target::InitBoxCfgList
        integer::LINE
        !---Local Vars---
        character*256::STR
        character*32::KEYWORD
        character*32::STRTMP(10)
        integer::N
        integer::IGroup
        type(InitBoxSimCfgList),pointer::PInitBoxCfgList=>null()
        type(InitBoxSimCfg)::tempInitBoxSimCfg
        !---Body---
        IGroup = 0

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE, "!", *100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSUBCTL")
                    exit
                case("&SPACEDISTSUBCTL")
                    call ReadInitBoxSimCfg_Simple_SpaceDist(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,LINE,*100)
                    PInitBoxCfgList=>InitBoxCfgList%next

                    DO While(associated(PInitBoxCfgList))

                        PInitBoxCfgList%TheValue%InitNNodes = InitBoxCfgList%TheValue%InitNNodes

                        call AllocateArray_Host(PInitBoxCfgList%TheValue%LayerThick,InitBoxCfgList%TheValue%InitNNodes,"LayerThick")

                        PInitBoxCfgList%TheValue%LayerThick = InitBoxCfgList%TheValue%LayerThick

                        PInitBoxCfgList=>PInitBoxCfgList%next
                    END DO

                case("&GROUPSUBCTL")
                    if(IGroup .eq. 0) then
                        call ReadInitBoxSimCfg_Simple_OneGroup(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList%TheValue,LINE,*100)
                    else
                        call tempInitBoxSimCfg%Clean_InitBoxSimCfg()
                        tempInitBoxSimCfg = InitBoxCfgList%TheValue
                        call ReadInitBoxSimCfg_Simple_OneGroup(hFile,SimBoxes,Host_SimuCtrlParam,tempInitBoxSimCfg,LINE,*100)
                        call InitBoxCfgList%AppendOne_InintSimBoxCfg(tempInitBoxSimCfg)
                    end if

                case default
                    write(*,*) "MFPSCUERROR: The Illegal flag: ",KEYWORD
                    pause
                    stop
            end select

        END DO

        return

        100 write(*,*) "MFPSCUERROR : Load init config file failed !"
            write(*,*) "At line :",LINE
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine ReadInitSimulationBoxesConfig_Simple


    !****************************************************************
    subroutine ReadInitBoxSimCfg_Simple_SpaceDist(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,LINE,*)
        !---Dummy Vars---
        integer,intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfgList),target::InitBoxCfgList
        integer::LINE
        !---Local Vars---
        character*256::STR
        character*32::KEYWORD
        character*32::STRNUMB(10)
        integer::N
        integer::INode
        !---Body---

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE, "!", *100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)


            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSBCTL")
                    exit
                case("&NNODES")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCUERROR : Too few parameters for number of nodes in mean field"
                        write(*,*) "You should special: '&NNODES The spatical nodes number = ' "
                        pause
                        stop
                    end if
                    InitBoxCfgList%TheValue%InitNNodes = ISTR(STRNUMB(1))

                case("&NODESPACE")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCUERROR : Too few parameters for node spaces in mean field"
                        write(*,*) "You should special: 'NODESPACE The space for each node = ' "
                        pause
                        stop
                    end if

                    if(InitBoxCfgList%TheValue%InitNNodes .LE. 0) then
                        write(*,*) "MFPSCUERROR: The node number must greater than 0"
                        write(*,*) InitBoxCfgList%TheValue%InitNNodes
                        pause
                        stop
                    end if

                    call AllocateArray_Host(InitBoxCfgList%TheValue%LayerThick,InitBoxCfgList%TheValue%InitNNodes,"LayerThick")

                    if(N .eq. 1) then
                        InitBoxCfgList%TheValue%LayerThick = DRSTR(STRNUMB(1))*SimBoxes%LatticeLength
                    else if(N .GT. 1) then
                        if(N .LT. InitBoxCfgList%TheValue%InitNNodes) then
                            write(*,*) "MFPSCUERROR: If you special the node thickness one by one, you should special all of them"
                            write(*,*) "The nodes number set is : ",InitBoxCfgList%TheValue%InitNNodes
                            write(*,*) "However, you only give : ",N, " node space"
                            write(*,*) "Or, you can give the nodes a same thickness"
                            pause
                            stop
                        end if

                        DO INode = 1,InitBoxCfgList%TheValue%InitNNodes
                            InitBoxCfgList%TheValue%LayerThick(INode) = DRSTR(STRNUMB(INode))*SimBoxes%LatticeLength
                        END DO
                    end if

                    if(sum(InitBoxCfgList%TheValue%LayerThick) .ne. SimBoxes%BOXSIZE(3)) then
                        write(*,*) "MFPSCUERROR: The total nodes thicknesses given is not equal with the box depth."
                        write(*,*) "The total nodes thicknesses is : ",sum(InitBoxCfgList%TheValue%LayerThick)
                        write(*,*) "The box depth is : ",SimBoxes%BOXSIZE(3)
                        pause
                        stop
                    end if
                case default
                    write(*,*) "MFPSCUERROR: Unknown key word :",KEYWORD
                    write(*,*) STR
                    write(*,*) "At Line: ",LINE
            end select
        END DO

        return
        100 return 1
    end subroutine ReadInitBoxSimCfg_Simple_SpaceDist


    !*****************************************************************
    subroutine ReadInitBoxSimCfg_Simple_OneGroup(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfg,LINE,*)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfg)::InitBoxCfg
        integer::LINE
        !---Local Vars---
        type(InitBoxSimCfg)::tempInitBoxSimCfg
        character*256::STR
        character*32::KEYWORD
        character*32::STRTMP(10)
        integer::N
        !---Body---

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE, "!", *100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSUBCTL")
                    exit
                case("&CLUSTERRANGESUBCTL")
                    call ReadClusterRange_Simple(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfg,LINE,*100)
                case("&SIZESUBCTL")
                    call ReadClusterSizeDist_Simple(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfg,LINE,*100)
                case("&CONCENTDISTSUBCTL")
                    call ReadClusterConcentrateDist_Simple(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfg,LINE,*100)
                case default
                    write(*,*) "MFPSCUERROR: Unknown Keyword :",KEYWORD
                    write(*,*) "At LINE",LINE
                    pause
                    stop
            end select
        END DO

        return
        100 return 1
    end subroutine ReadInitBoxSimCfg_Simple_OneGroup

    !****************************************************************
    subroutine ReadInitSimulationBoxesConfig_SpecialDistFromFile(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,LINE)
        !---Dummy Vars---
        integer,intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfgList)::InitBoxCfgList
        integer::LINE
        !---Local Vars---
        character*256::STR
        character*32::KEYWORD
        character*32::STRTMP(10)
        integer::N
        !---Body---

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE, "!", *100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSUBCTL")
                    exit
                case("&INITFILE")
                    call EXTRACT_SUBSTR(STR,1,N,STRTMP)

                    if(N .LT. 1) then
                        write(*,*) "MFPSCUERROR: Too few parameters for the clusters initialize file path"
                        write(*,*) "You should special: &INITFILE The initialize file path = "
                        pause
                        stop
                    end if

                    InitBoxCfgList%TheValue%InitCfgFileName = INQUIREFILE(STRTMP(1),Host_SimuCtrlParam%InputFilePath)

                case default
                    write(*,*) "MFPSCUERROR: The Illegal flag: ",KEYWORD
                    pause
                    stop
            end select
        END DO

        return

        100 write(*,*) "MFPSCUERROR : Load init config file failed !"
            write(*,*) "At line :",LINE
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine ReadInitSimulationBoxesConfig_SpecialDistFromFile

    !*****************************************************************
    subroutine ReadInitSimulationBoxesConfig_SpecialDistFromExteFunc(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfgList,LINE)
        !---Dummy Vars---
        integer,intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfgList)::InitBoxCfgList
        integer::LINE

        ! @todo (zhail#1#):


        return
    end subroutine ReadInitSimulationBoxesConfig_SpecialDistFromExteFunc

    !*****************************************************************
    subroutine ReadClusterRange_Simple(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfg,LINE,*)
        !---Dummy Vars---
        integer,intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfg)::InitBoxCfg
        integer::LINE
        !---Local Vars---
        character*256::STR
        character*30::KEYWORD
        character*32::STRTMP(10)
        integer::N
        !---Body---

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSBUCTL")
                    exit
                case("&CLUSTERSGROUP")
                    call EXTRACT_NUMB(STR,1,N,STRTMP)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCUERROR: Too few parameters for clusters group number ."
                        write(*,*) "You should special: '&CLUSTERSGROUP The max clusters group = ' "
                        pause
                        stop
                    end if
                    InitBoxCfg%InitCKind = ISTR(STRTMP(1))

                    if(InitBoxCfg%InitCKind .LE. 0) then
                        write(*,*) "MFPSCUERROR: The clusters kind number must greater than 0"
                        write(*,*) InitBoxCfg%InitCKind
                        pause
                        stop
                    end if

                case("&RANGE")
                    call EXTRACT_NUMB(STR,3,N,STRTMP)
                    if(N .LT. 2) then
                        write(*,*) "MFPSCUERROR: You must special the clusters range beginning and ending."
                        write(*,*) "You should special like that : ' &RANGE The clusters atoms start from = , end to = , the interval = '"
                        write(*,*) "At LINE :",LINE
                        pause
                        stop
                    end if

                    if(N .GE. 2) then
                        InitBoxCfg%CKind_Range_From = ISTR(STRTMP(1))
                        InitBoxCfg%CKind_Range_To = ISTR(STRTMP(2))
                    end if

                    if(N .GE. 3) then
                        InitBoxCfg%CKind_Range_Interval = ISTR(STRTMP(3))
                    end if

                    if((InitBoxCfg%CKind_Range_To - InitBoxCfg%CKind_Range_From) .LT. InitBoxCfg%CKind_Range_Interval) then
                        write(*,*) "MFPSCUERROR: The clusters range define is not true"
                        write(*,*) "At line: ",LINE
                        write(*,*) InitBoxCfg%CKind_Range_From, InitBoxCfg%CKind_Range_To, InitBoxCfg%CKind_Range_Interval
                        pause
                        stop
                    end if

                    if(((InitBoxCfg%CKind_Range_To - InitBoxCfg%CKind_Range_From)/InitBoxCfg%CKind_Range_Interval +2 ) .GT. InitBoxCfg%InitCKind) then
                        write(*,*) "MFPSCUERROR: The total clusters kind you defined is : ",InitBoxCfg%InitCKind
                        write(*,*) "However, the ranges you define is: ",InitBoxCfg%CKind_Range_From, InitBoxCfg%CKind_Range_To, InitBoxCfg%CKind_Range_Interval
                        write(*,*) "Which had exceed the clusters kind."
                        pause
                    end if



                case default
                    write(*,*) "MFPSCUERROR: Unknown keyword: ",KEYWORD
                    pause
                    stop
            end select

        END DO

        return
        100 return 1
    end subroutine

    !***************************************************************
    subroutine ReadClusterSizeDist_Simple(hFile,SimBoxes,Host_SimuCtrlParam,InitBoxCfg,LINE,*)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfg)::InitBoxCfg
        integer::LINE
        !---Local Vars---
        character*512::STR
        character*32::KEYWORD
        character*32::STRTMP(10)
        integer::N
        integer::NElements
        integer::I
        integer::TheIndex
        !---Body---
        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            SELECT CASE(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSUBCTL")
                    exit
                CASE("&NATOMDIST")
                    call EXTRACT_NUMB(STR,4,N,STRTMP)
                    if(N .LT. 4) then
                        write(*,*) "MFPSCUERROR: Too few parameters for the cluster size distribution"
                        write(*,*) "You should special: &NATOMDIST The atoms number in each cluster distribution as Gauss that central = , distribution half width = , left cut = ,right cut ="
                        pause
                        stop
                    end if
                    InitBoxCfg%NAINI = DRSTR(STRTMP(1))
                    InitBoxCfg%NASDINI  = DRSTR(STRTMP(2))
                    InitBoxCfg%NACUT(1) = DRSTR(STRTMP(3))
                    InitBoxCfg%NACUT(2) = DRSTR(STRTMP(4))

                    if(InitBoxCfg%NACUT(1) .GE. InitBoxCfg%NACUT(2)) then
                        write(*,*) "MFPSCUERROR: The right cut cannot less than left cut."
                        write(*,*) "LCut",InitBoxCfg%NACUT(1)
                        write(*,*) "RCut",InitBoxCfg%NACUT(2)
                        pause
                        stop
                    end if
                case("&ELEMENTCOMPOSIT")
                    call EXTRACT_SUBSTR(STR,p_ATOMS_GROUPS_NUMBER,NElements,STRTMP)
                    if(NElements .LE. 0) then
                        write(*,*) "MFPSCUERROR: None of atoms kind (Elements) are specialized "
                        write(*,*) "You should special like that : &ELEMENTCOMPOSIT The included element = 'A', 'B' ."
                        pause
                        stop
                    else if(NElements .GT. p_ATOMS_GROUPS_NUMBER) then
                        write(*,*) "MFPSCUERROR: the specialized elements kinds is : ",N
                        write(*,*) "which is great than the max permitted elements kinds :",p_ATOMS_GROUPS_NUMBER
                        pause
                        stop
                    else if(NElements .GT. 1) then
                        write(*,*) "MFPSCU Info: Currently, the program only support one element."
                        write(*,*) "So the first element would be chosen"
                    else
                        DO I = 1,NElements
                            InitBoxCfg%Elemets(I) = adjustl(trim(STRTMP(I)))
                            call UPCASE(InitBoxCfg%Elemets(I))
                        END DO
                    end if
                    ! Current, we only support one element
                    InitBoxCfg%CompositWeight = 0
                    TheIndex = SimBoxes%Atoms_list%FindIndexBySymbol(InitBoxCfg%Elemets(1))
                    InitBoxCfg%CompositWeight(TheIndex) = 1.D0
                    InitBoxCfg%CompositWeight = InitBoxCfg%CompositWeight/sum(InitBoxCfg%CompositWeight)

!                    call EXTRACT_NUMB(STR,p_ATOMS_GROUPS_NUMBER,N,STRTMP)
!                    if(N .ne. NElements) then
!                        write(*,*) "MFPSCUERROR: The elements weights number is not equal with the elements kinds which given."
!                        write(*,*) "The elements kinds number is :",NElements
!                        write(*,*) "But the weights number is :",N
!                        pause
!                        stop
!                    else
!                        InitBoxCfg%CompositWeight = 0.D0
!
!                        DO I = 1,N
!                            TheIndex = SimBoxes%Atoms_list%FindIndexBySymbol(InitBoxCfg%Elemets(I))
!                            InitBoxCfg%CompositWeight(TheIndex) = DRSTR(STRTMP(I))
!                        END DO
!
!                        if(sum(InitBoxCfg%CompositWeight) .LE. 0.D0) then
!                            write(*,*) "MFPSCUERROR: The sum of elements weights must great than 0 ."
!                            write(*,*) STR
!                            write(*,*) "At Line :",LINE
!                            pause
!                            stop
!                        end if
!
!                        InitBoxCfg%CompositWeight = InitBoxCfg%CompositWeight/sum(InitBoxCfg%CompositWeight)
!                    end if
                CASE default
                    write(*,*) "MFPSCUERROR: Illegal Symbol: ", KEYWORD
                    pause
                    stop
            END SELECT

        END DO

        return

        100 return 1
    end subroutine ReadClusterSizeDist_Simple

    !***************************************************************
    subroutine ReadClusterConcentrateDist_Simple(hFile,Host_SimBoxes,Host_SimuCtrlParam,InitBoxCfg,LINE,*)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::hFile
        type(SimulationBoxes)::Host_SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfg)::InitBoxCfg
        integer::LINE
        !---Local Vars---
        character*256::STR
        character*32::KEYWORD
        character*32::STRTMP(10)
        integer::N
        integer::LayerStart
        integer::LayerEnd
        integer::LayerNum
        integer::I
        real(kind=KMCDF)::TotalLayerThick
        real(kind=KMCDF)::TotalPNC
        !---Body---

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            SELECT CASE(KEYWORD(1:LENTRIM(KEYWORD)))
                CASE("&ENDSUBCTL")
                    exit
                CASE("&CONCENTRATE")
                    call EXTRACT_NUMB(STR,1,N,STRTMP)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCUERROR: Too few parameters for concentration setting"
                        write(*,*) "You should special: '&CONCENTRATE THE INITIAL CONCENTRATE OF ENTRIES =   ! (1/(cm^3))'"
                        pause
                        stop
                    end if
                    InitBoxCfg%InitConcentrate = DRSTR(STRTMP(1))

                CASE("&DEPTH_LAYER")
                    InitBoxCfg%InitDepthDistType = p_DEPT_DIS_Layer

                    call EXTRACT_NUMB(STR,2,N,STRTMP)
                    if(N .LT. 2) then
                        write(*,*) "MFPSCUERROR: Too few parameters for the bubble depth distribution layer type"
                        write(*,*) "You should special: &DEPTH_LAYER The Layer from =  , end =   , THE CONCENTRATION DISTRIBUTION = "
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if

                    LayerStart = ISTR(STRTMP(1))
                    LayerEnd = ISTR(STRTMP(2))

                    LayerNum =  LayerEnd - LayerStart + 1

                    if(LayerStart .LT. 0) then
                        write(*,*) "MFPSCUERROR: The layer start layer should greater than 1"
                        write(*,*) LayerStart
                        write(*,*) "At line :",LINE
                        write(*,*) STR
                        pause
                        stop
                    end if

                    if(LayerEnd .GT. InitBoxCfg%InitNNodes) then
                        write(*,*) "MFPSCUERROR: The initial layer had exceed the max layer number"
                        write(*,*) LayerEnd
                        write(*,*) "At line :",LINE
                        write(*,*) STR
                        pause
                        stop
                    end if

                    if(LayerNum .LT. 1) then
                        write(*,*) "MFPSCUERROR: The layer number should greater than 1"
                        write(*,*) "At line :",LINE
                        write(*,*) STR
                        pause
                        stop
                    end if

                    call AllocateArray_Host(InitBoxCfg%PNCLayers,InitBoxCfg%InitNNodes,"PNCLayers")
                    InitBoxCfg%PNCLayers = 0.D0

                    if(N .eq. 3) then
                        InitBoxCfg%PNCLayers = DRSTR(STRTMP(3))
                    else if(N .GT. 3) then
                        if((N-2) .LT. LayerNum) then
                            write(*,*) "MFPSCUERROR: If you special the initial concentrates for layers(nodes) one by one,"
                            write(*,*) "You should give each of them."
                            pause
                            stop
                        end if

                        DO I = LayerStart,LayerEnd
                            InitBoxCfg%PNCLayers(I) = DRSTR(STRTMP(I - LayerStart + 1 + 2))
                        END DO
                    end if

                    TotalPNC = sum(InitBoxCfg%PNCLayers)
                    if(TotalPNC .LE. 0) then
                        write(*,*) "MFPSCUERROR: The total percent cannot less than 0"
                        pause
                        stop
                    end if
                    InitBoxCfg%PNCLayers = InitBoxCfg%PNCLayers/TotalPNC

                CASE("&DEPTH_GAUSS")

                    InitBoxCfg%InitDepthDistType = p_DEPT_DIS_GAS

                    call EXTRACT_NUMB(STR,2,N,STRTMP)
                    if(N .LT. 2) then
                        write(*,*) "MFPSCUERROR: Too few parameters for the bubble depth distribution gauss type"
                        write(*,*) "You shoud special: &DEPTH_GAUSS THE GAUSS DISTRIBUTION CENTRAL = , THE HALF WIDTH = "
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if
                    InitBoxCfg%DepthINI = DRSTR(STRTMP(1))*C_NM2CM
                    InitBoxCfg%DepthSDINI = DRSTR(STRTMP(2))*C_NM2CM
                CASE default
                    write(*,*) "MFPSCUERROR: Illegal Symbol: ", KEYWORD
                    pause
                    stop
            END SELECT

        END DO

        return
        100 return 1
    end subroutine ReadClusterConcentrateDist_Simple

    !*****************************************************************
    subroutine DOInitSimulationBoxesConfig(SimBoxes,Host_SimuCtrlParam,Record,InitBoxCfgList)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(MigCoalClusterRecord)::Record
        type(InitBoxSimCfgList),target::InitBoxCfgList
        !---Local Vars---
        integer::MultiBox
        type(InitBoxSimCfgList),pointer::cursor=>null()
        integer::IBox
        integer::TNC0
        !---Body---

        MultiBox = Host_SimuCtrlParam%MultiBox

        cursor=>InitBoxCfgList

        call SimBoxes%m_ClustersInfo_CPU%Clean()

        DO While(associated(cursor))

            select case(cursor%TheValue%InitType)
                case(p_ClusterIniConfig_Simple)
                    call DoInitSimulationBoxesConfig_Simple(SimBoxes,Host_SimuCtrlParam,cursor%TheValue)
                case(p_ClusterIniConfig_SpecialDistFromFile)
                    call DoInitSimulationBoxesConfig_SpecialDistFromFile(SimBoxes,Host_SimuCtrlParam,Record,cursor%TheValue)
                case(p_ClusterIniConfig_SpecialDistFromExteFunc)
                    call DoInitSimulationBoxesConfig_SpecialDistFromExteFunc(SimBoxes,Host_SimuCtrlParam,cursor%TheValue)
                case default
                    write(*,*) "MFPSCUERROR: Unknown strategy for the box initialization:",cursor%TheValue%InitType
                    pause
                    stop
            end select

            cursor=>cursor%next
        END DO

        Nullify(cursor)

        return
    end subroutine

    !****************************************************************
    subroutine DoInitSimulationBoxesConfig_Simple(SimBoxes,Host_SimuCtrlParam,InitBoxCfg)
        !---Dummy Vars---
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfg)::InitBoxCfg
        !---Body--
        select case(InitBoxCfg%InitDepthDistType)
            case(p_DEPT_DIS_Layer)
                call Init_Depth_Dis_LAY(SimBoxes,Host_SimuCtrlParam,InitBoxCfg)
            case(p_DEPT_DIS_BOX)
                call Init_Depth_Dis_SubBox(SimBoxes,Host_SimuCtrlParam,InitBoxCfg)
            case(p_DEPT_DIS_GAS)
                call Init_Depth_Dis_Gauss(SimBoxes,Host_SimuCtrlParam,InitBoxCfg)
            case default
                write(*,*) "MFPSCUERROR : Unknown way to initial the simulation box conifuration :",InitBoxCfg%InitDepthDistType
                pause
                stop
        end select

        return
    end subroutine DoInitSimulationBoxesConfig_Simple

    !****************************************************************
    subroutine DoInitSimulationBoxesConfig_SpecialDistFromFile(SimBoxes,Host_SimuCtrlParam,Record,InitBoxCfg)
        !---Dummy Vars---
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(MigCoalClusterRecord)::Record
        type(InitBoxSimCfg)::InitBoxCfg
        !---Body---

        call SimBoxes%PutinCfg(Host_SimuCtrlParam,Record,InitBoxCfg%InitCfgFileName,m_RNFACTOR,m_FREESURDIFPRE,m_GBSURDIFPRE)

        return
    end subroutine DoInitSimulationBoxesConfig_SpecialDistFromFile

    !*****************************************************************
    subroutine DoInitSimulationBoxesConfig_SpecialDistFromExteFunc(SimBoxes,Host_SimuCtrlParam,InitBoxCfg)
        !---Dummy Vars---
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfg)::InitBoxCfg
        ! @todo (zhail#1#):


        return
    end subroutine DoInitSimulationBoxesConfig_SpecialDistFromExteFunc

    !**************************************************************
    subroutine Init_Depth_Dis_LAY(Host_Boxes,Host_SimuCtrlParam,InitBoxCfg)
      !*** Purpose: To initialize the system (clusters distributed as the form of layer)
      ! Host_Boxes: the boxes information in host
      use RAND32_MODULE
      implicit none
      !---Dummy Vars---
      type(SimulationBoxes)::Host_Boxes
      type(SimulationCtrlParam)::Host_SimuCtrlParam
      type(InitBoxSimCfg)::InitBoxCfg
      !-----local variables---
      integer::MultiBox
      real(kind=KMCDF)::POS(3)
      real(kind=KMCDF)::Z0
      real(kind=KMCDF)::BOXBOUNDARY(3,2)
      real(kind=KMCDF)::BOXSIZE(3)
      integer::IBox, II, IC, LAY, PNC
      integer::J
      integer::SNC0,SNC
      integer::LayerNum
      integer::NAtoms
      type(DiffusorValue)::TheDiffusorValue
      !---Body---
!      MultiBox = Host_SimuCtrlParam%MultiBox
!
!      BOXBOUNDARY = Host_Boxes%BOXBOUNDARY
!      BOXSIZE = Host_Boxes%BOXSIZE
!
!      LayerNum = size(InitBoxCfg%LayerThick)
!
!      call Host_Boxes%ExpandClustersInfor_CPU(Host_SimuCtrlParam,InitBoxCfg%NClusters)
!
!      DO IBox = 1,MultiBox
!
!        SNC0 = InitBoxCfg%NClusters(IBox)
!        SNC  = 0
!        Z0  = 0.D0
!
!        IC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - SNC0
!
!        DO LAY=1, LayerNum
!
!            if(LAY .EQ. LayerNum) THEN
!                PNC = SNC0 - SNC
!            else
!                PNC = SNC0*InitBoxCfg%PNCLayers(LAY)
!            end if
!
!            SNC = SNC + PNC
!            Z0 = BOXBOUNDARY(3,1) + sum(InitBoxCfg%LayerThick(1:LAY-1))
!
!            DO II = 1, PNC
!
!                IC = IC + 1
!                !Initialize the position of clusters
!                POS(1) = DRAND32()*BOXSIZE(1)+BOXBOUNDARY(1,1)
!                POS(2) = DRAND32()*BOXSIZE(2)+BOXBOUNDARY(2,1)
!                POS(3) = DRAND32()*InitBoxCfg%LayerThick(LAY) + Z0
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_POS = POS
!
!                !Give the cluster an type(layer) ID for the convenience of visualization
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Layer = LAY
!
!                !*** Initialize the size of the clusters
!                NAtoms = RGAUSS0_WithCut(InitBoxCfg%NAINI, InitBoxCfg%NASDINI,InitBoxCfg%NACUT(1),InitBoxCfg%NACUT(2))
!
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA = NAtoms*InitBoxCfg%CompositWeight
!
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) = Host_Boxes%m_GrainBoundary%GrainBelongsTo(POS,Host_Boxes%HBOXSIZE,Host_Boxes%BOXSIZE,Host_SimuCtrlParam)
!
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu = p_ACTIVEFREE_STATU
!
!                TheDiffusorValue = Host_Boxes%m_DiffusorTypesMap%Get(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC))
!
!                !-- In Current application, the simple init distribution is only considered in free matrix, if you want to init the clusters in GB---
!                !---you should init the distribution by external file---
!                select case(TheDiffusorValue%ECRValueType_Free)
!                    case(p_ECR_ByValue)
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = TheDiffusorValue%ECR_Free
!                    case(p_ECR_ByBCluster)
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = DSQRT(sum(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA)/m_RNFACTOR)
!                end select
!
!                select case(TheDiffusorValue%DiffusorValueType_Free)
!                    case(p_DiffuseCoefficient_ByValue)
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
!                    case(p_DiffuseCoefficient_ByArrhenius)
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
!                    case(p_DiffuseCoefficient_ByBCluster)
!                        ! Here we adopt a model that D=D0*(1/R)**Gama
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = m_FREESURDIFPRE*(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD**(-p_GAMMA))
!                end select
!
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(p_ACTIVEFREE_STATU) + 1
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEFREE_STATU) + 1
!
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC0(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC0(p_ACTIVEFREE_STATU) + 1
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC0(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC0(p_ACTIVEFREE_STATU) + 1
!
!                Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) + 1
!                Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) = Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) + 1
!
!            END DO
!        END DO
!
!      END DO

      return
    end subroutine Init_Depth_Dis_LAY

    !**************************************************************
    subroutine Init_Depth_Dis_SubBox(Host_Boxes,Host_SimuCtrlParam,InitBoxCfg)
      !*** Purpose: To initialize the system (clusters distributed as the form of layer)
      ! Host_Boxes: the boxes information in host
      use RAND32_MODULE
      implicit none
      !---Dummy Vars---
      type(SimulationBoxes)::Host_Boxes
      type(SimulationCtrlParam)::Host_SimuCtrlParam
      type(InitBoxSimCfg)::InitBoxCfg
      !-----local variables---
      integer::MultiBox
      real(kind=KMCDF)::POS(3)
      real(kind=KMCDF)::SUBBOXSIZE(3)
      integer::IBox, II, IC
      integer::SNC0
      integer::I
      integer::NAtoms
      type(DiffusorValue)::TheDiffusorValue
      !---Body---
!      MultiBox = Host_SimuCtrlParam%MultiBox
!
!      DO I = 1,3
!        SUBBOXSIZE(I) = InitBoxCfg%SUBBOXBOUNDARY(I,2) - InitBoxCfg%SUBBOXBOUNDARY(I,1)
!      END DO
!
!      call Host_Boxes%ExpandClustersInfor_CPU(Host_SimuCtrlParam,InitBoxCfg%NClusters)
!
!      DO IBox = 1,MultiBox
!
!        SNC0 = InitBoxCfg%NClusters(IBox)
!
!        IC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - SNC0
!
!        DO II = 1, SNC0
!
!            IC = IC + 1
!            !Initialize the position of clusters
!            DO I = 1,3
!                POS(I) = DRAND32()*SUBBOXSIZE(I) + InitBoxCfg%SUBBOXBOUNDARY(I,1)
!            END DO
!
!            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_POS = POS
!            !Give the cluster an type(layer) ID for the convenience of visualization
!            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Layer = 1
!
!            NAtoms = RGAUSS0_WithCut(InitBoxCfg%NAINI, InitBoxCfg%NASDINI,InitBoxCfg%NACUT(1),InitBoxCfg%NACUT(2))
!
!            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA = NAtoms*InitBoxCfg%CompositWeight
!
!            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) = Host_Boxes%m_GrainBoundary%GrainBelongsTo(POS,Host_Boxes%HBOXSIZE,Host_Boxes%BOXSIZE,Host_SimuCtrlParam)
!
!            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu = p_ACTIVEFREE_STATU
!
!            TheDiffusorValue = Host_Boxes%m_DiffusorTypesMap%Get(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC))
!
!            !-- In Current application, the simple init distribution is only considered in free matrix, if you want to init the clusters in GB---
!            !---you should init the distribution by external file---
!            select case(TheDiffusorValue%ECRValueType_Free)
!                case(p_ECR_ByValue)
!                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = TheDiffusorValue%ECR_Free
!                case(p_ECR_ByBCluster)
!                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = DSQRT(sum(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA)/m_RNFACTOR)
!            end select
!
!            select case(TheDiffusorValue%DiffusorValueType_Free)
!                case(p_DiffuseCoefficient_ByValue)
!                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
!                case(p_DiffuseCoefficient_ByArrhenius)
!                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
!                case(p_DiffuseCoefficient_ByBCluster)
!                    ! Here we adopt a model that D=D0*(1/R)**Gama
!                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = m_FREESURDIFPRE*(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD**(-p_GAMMA))
!            end select
!
!            Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(p_ACTIVEFREE_STATU) + 1
!            Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEFREE_STATU) + 1
!
!            Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC0(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC0(p_ACTIVEFREE_STATU) + 1
!            Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC0(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC0(p_ACTIVEFREE_STATU) + 1
!
!            Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) + 1
!            Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) = Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) + 1
!        END DO
!
!      END DO

      return
    end subroutine Init_Depth_Dis_SubBox

    !**************************************************************
    subroutine Init_Depth_Dis_Gauss(Host_Boxes,Host_SimuCtrlParam,InitBoxCfg)
        !*** Purpose: To initialize the system (clusters distributed as the form of gauss in depth)
        ! Host_Boxes: the boxes information in host
        use RAND32_MODULE
        implicit none
        !-------Dummy Vars------
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(InitBoxSimCfg)::InitBoxCfg
        !---local variables---
        integer::MultiBox
        real(kind=KMCDF)::POS(3)
        real(kind=KMCDF)::SEP
        real(kind=KMCDF)::BOXBOUNDARY(3,2)
        real(kind=KMCDF)::BOXSIZE(3)
        integer::IBox,II,IC
        integer::SNC0
        integer::NAtoms
        type(DiffusorValue)::TheDiffusorValue
        !---Body---
!        MultiBox = Host_SimuCtrlParam%MultiBox
!        BOXBOUNDARY = Host_Boxes%BOXBOUNDARY
!        BOXSIZE = Host_Boxes%BOXSIZE
!
!        call Host_Boxes%ExpandClustersInfor_CPU(Host_SimuCtrlParam,InitBoxCfg%NClusters)
!
!        DO IBox = 1,MultiBox
!
!            SNC0 = InitBoxCfg%NClusters(IBox)
!
!            IC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - SNC0
!
!            DO II = 1, SNC0
!
!                IC = IC + 1
!                !Initialize the position of clusters
!                POS(1) = DRAND32()*BOXSIZE(1) + BOXBOUNDARY(1,1)
!                POS(2) = DRAND32()*BOXSIZE(2) + BOXBOUNDARY(2,1)
!                POS(3) = RGAUSS0_WithCut(InitBoxCfg%DepthINI, InitBoxCfg%DepthSDINI,BOXBOUNDARY(3,1),BOXBOUNDARY(3,2))
!
!                if(POS(3) .LT. BOXBOUNDARY(3,1)) then
!                    POS(3) = BOXBOUNDARY(3,1)
!                end if
!
!                if(POS(3) .GT. BOXBOUNDARY(3,2)) then
!                    POS(3) = BOXBOUNDARY(3,2)
!                end if
!
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_POS = POS
!
!                !Give the cluster an type(layper) ID for the convenience of visualization
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Layer = 1
!
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) = Host_Boxes%m_GrainBoundary%GrainBelongsTo(POS,Host_Boxes%HBOXSIZE,Host_Boxes%BOXSIZE,Host_SimuCtrlParam)
!
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu = p_ACTIVEFREE_STATU
!
!                !*** Initialize the size of the clusters
!                NAtoms = RGAUSS0_WithCut(InitBoxCfg%NAINI, InitBoxCfg%NASDINI,InitBoxCfg%NACUT(1),InitBoxCfg%NACUT(2))
!
!                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA = NAtoms*InitBoxCfg%CompositWeight
!
!                TheDiffusorValue = Host_Boxes%m_DiffusorTypesMap%Get(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC))
!
!                !-- In Current application, the simple init distribution is only considered in free matrix, if you want to init the clusters in GB---
!                !---you should init the distribution by external file---
!                select case(TheDiffusorValue%ECRValueType_Free)
!                    case(p_ECR_ByValue)
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = TheDiffusorValue%ECR_Free
!                    case(p_ECR_ByBCluster)
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = DSQRT(sum(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA)/m_RNFACTOR)
!                end select
!
!                select case(TheDiffusorValue%DiffusorValueType_Free)
!                    case(p_DiffuseCoefficient_ByValue)
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
!                    case(p_DiffuseCoefficient_ByArrhenius)
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
!                    case(p_DiffuseCoefficient_ByBCluster)
!                        ! Here we adopt a model that D=D0*(1/R)**Gama
!                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = m_FREESURDIFPRE*(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD**(-p_GAMMA))
!                end select
!
!
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(p_ACTIVEFREE_STATU) + 1
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEFREE_STATU) + 1
!
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC0(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC0(p_ACTIVEFREE_STATU) + 1
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC0(p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC0(p_ACTIVEFREE_STATU) + 1
!
!                Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) + 1
!                Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) = Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) + 1
!            END DO
!
!        END DO

        return
   end subroutine Init_Depth_Dis_Gauss

!    !*****************************************************************
!    subroutine OutPutCurrent(Host_Boxes,Dev_Boxes,Host_SimuCtrlParam,TheMigCoaleStatisticInfo,Record)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(MigCoaleStatisticInfo_Used)::TheMigCoaleStatisticInfo
!        type(MigCoalClusterRecord)::Record
!        !---Local Vars---
!        logical::OutIntegralBoxStatistic
!        logical::OutEachBoxStatistic
!        integer::NC0
!        !---Body---
!
!        OutIntegralBoxStatistic = .false.
!        OutEachBoxStatistic = .true.
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(Host_SimuCtrlParam%MultiBox,2) .GT. 0) then
!            NC0 = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(Host_SimuCtrlParam%MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            NC0 = 0
!        end if
!
!        OutIntegralBoxStatistic = Record%WhetherOutSizeDist_IntegralBox(Host_SimuCtrlParam)
!
!        OutEachBoxStatistic = Record%WhetherOutSizeDist_EachBox(Host_SimuCtrlParam)
!
!
!        if(OutIntegralBoxStatistic .eq. .true. .or. OutEachBoxStatistic .eq. .true.) then
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
!
!        end if
!
!        ! check if need to output intermediate configure
!        if(Host_SimuCtrlParam%OutPutConfFlag .eq. mp_OutTimeFlag_ByIntervalSteps) then
!            if((Record%GetSimuSteps() - Record%GetLastRecordOutConfigTime()) .GE. Host_SimuCtrlParam%OutPutConfValue .OR. &
!                Record%GetSimuTimes() .GE. Host_SimuCtrlParam%TermTValue) then
!
!                call Dev_Boxes%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,NC0,IfCpyNL=.false.)
!
!                call Host_Boxes%PutoutCfg(Host_SimuCtrlParam,Record)
!
!                call Record%SetLastRecordOutConfigTime(dble(Record%GetSimuSteps()))
!
!            end if
!        else if(Host_SimuCtrlParam%OutPutConfFlag .eq. mp_OutTimeFlag_ByIntervalRealTime) then
!            if((Record%GetSimuTimes() - Record%GetLastRecordOutConfigTime()) .GE. Host_SimuCtrlParam%OutPutConfValue .OR. &
!                Record%GetSimuTimes() .GE. Host_SimuCtrlParam%TermTValue) then
!
!                call Dev_Boxes%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,NC0,IfCpyNL=.false.)
!
!                call Host_Boxes%PutoutCfg(Host_SimuCtrlParam,Record)
!
!                call Record%SetLastRecordOutConfigTime(Record%GetSimuTimes())
!            end if
!
!        else if(Host_SimuCtrlParam%OutPutConfFlag .eq. mp_OutTimeFlag_ByIntervalTimeMagnification) then
!            if((Record%GetSimuTimes()/Host_SimuCtrlParam%OutPutConfValue) .GE. Record%GetLastRecordOutConfigTime() .OR. &
!                Record%GetSimuTimes() .GE. Host_SimuCtrlParam%TermTValue) then
!
!                call Dev_Boxes%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,NC0,IfCpyNL=.false.)
!
!                call Host_Boxes%PutoutCfg(Host_SimuCtrlParam,Record)
!
!                call Record%SetLastRecordOutConfigTime(Record%GetSimuTimes())
!            end if
!        end if
!
!    end subroutine
!
!    !*****************************************************************
!    subroutine PutOut_Instance_Statistic_IntegralBox(Host_Boxes,Host_SimuCtrlParam,TheMigCoaleStatisticInfo,Record,Model)
!        implicit none
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(MigCoaleStatisticInfo_Used)::TheMigCoaleStatisticInfo
!        type(MigCoalClusterRecord)::Record
!        integer,intent(in)::Model                           ! 1 would output ot file
!        !---Local Vars---
!        integer::MultiBox
!        real(kind=KMCDF)::RMIN
!        real(kind=KMCDF)::Concentrate
!        integer::NCAct
!        real(kind=KMCDF)::RAVA
!        real(kind=KMCDF)::NAVA
!        !---Body---
!        MultiBox = Host_SimuCtrlParam%MultiBox
!
!        ASSOCIATE(TBasicInfo=>Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral,TMigStatInfo=>TheMigCoaleStatisticInfo%statistic_IntegralBox)
!
!            if(Model .eq. 0) then
!
!                NCAct = TBasicInfo%NC(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU)
!
!                RAVA = (TBasicInfo%NC(p_ACTIVEFREE_STATU)*TMigStatInfo%RAVA(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU)*TMigStatInfo%RAVA(p_ACTIVEINGB_STATU))/NCAct
!
!                if(NCAct .GT. 0) then
!                    NAVA = dble(TBasicInfo%NA(p_ACTIVEFREE_STATU) +  TBasicInfo%NA(p_ACTIVEINGB_STATU))/dble(NCAct)
!                else
!                    NAVA = 0.D0
!                end if
!
!                Concentrate = NCAct/(MultiBox*Host_Boxes%BOXVOLUM)
!
!                write(Record%HSizeStatistic_TotalBox, fmt="(I20,1x,1PE20.4,1x,9(I20,1x),15(1PE20.4,1x),7(I20,1x))") Record%GetSimuSteps(),                 &
!                                                                                                           Record%GetSimuTimes(),                          &
!                                                                                                           NCAct,                                          &
!                                                                                                           sum(TBasicInfo%NC),                             &
!                                                                                                           sum(TBasicInfo%NA),                             &
!                                                                                                           TBasicInfo%NC(1:p_NUMBER_OF_STATU),             &
!                                                                                                           TMigStatInfo%RMIN(1:p_NUMBER_OF_STATU)*C_CM2NM, &
!                                                                                                           TMigStatInfo%RMAX(1:p_NUMBER_OF_STATU)*C_CM2NM, &
!                                                                                                           RAVA*C_CM2NM,                                   &
!                                                                                                           NAVA,                                           &
!                                                                                                           Concentrate,                                    &
!                                                                                                           Record%GetImplantedEntitiesNum(),                &
!                                                                                                           TBasicInfo%NC(p_ACTIVEFREE_STATU),               &
!                                                                                                           TBasicInfo%NC(p_ACTIVEINGB_STATU),               &
!                                                                                                           TBasicInfo%NC(p_OUT_DESTROY_STATU:p_ABSORBED_STATU) + &
!                                                                                                           Record%RecordNCBeforeSweepOut_Integal(p_OUT_DESTROY_STATU:p_ABSORBED_STATU)
!                call flush(Record%HSizeStatistic_TotalBox)
!
!                if((TBasicInfo%NC(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU) + &
!                   sum(TBasicInfo%NC(p_OUT_DESTROY_STATU:p_ABSORBED_STATU) + Record%RecordNCBeforeSweepOut_Integal(p_OUT_DESTROY_STATU:p_ABSORBED_STATU)) - &
!                   Record%GetImplantedEntitiesNum() - sum(TBasicInfo%NC0) - TBasicInfo%NCDumpAdded) .ne. 0) then
!
!                   write(*,*) "MFPSCUERROR: The clusters number is not conservation."
!                   write(*,*) "The accumulated clusters for all kinds =",TBasicInfo%NC(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU) + &
!                               sum(TBasicInfo%NC(p_OUT_DESTROY_STATU:p_ABSORBED_STATU) + Record%RecordNCBeforeSweepOut_Integal(p_OUT_DESTROY_STATU:p_ABSORBED_STATU))
!                   write(*,*) "The total implanted cluster number = ",Record%GetImplantedEntitiesNum()
!                   write(*,*) "The initial cluster number plus rescale added number = ",TBasicInfo%NC0(p_ACTIVEFREE_STATU) + TBasicInfo%NC0(p_ACTIVEINGB_STATU) + TBasicInfo%NCDumpAdded
!                   pause
!                end if
!
!            end if
!
!            NCAct = TBasicInfo%NC(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU)
!
!            RAVA = (TBasicInfo%NC(p_ACTIVEFREE_STATU)*TMigStatInfo%RAVA(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU)*TMigStatInfo%RAVA(p_ACTIVEINGB_STATU))/NCAct
!
!            if(NCAct .GT. 0) then
!                NAVA = dble(TBasicInfo%NA(p_ACTIVEFREE_STATU) +  TBasicInfo%NA(p_ACTIVEINGB_STATU))/dble(NCAct)
!            else
!                NAVA = 0.D0
!            end if
!
!            Concentrate = NCAct/(MultiBox*Host_Boxes%BOXVOLUM)
!
!            write(6, fmt= "(130(A15,1x))")   "Step",            &
!                                             "Time",            &
!                                             "NC(ACTIVEFREE)",  &
!                                             "NC(ACTIVEINGB)",  &
!                                             "TNCAct",          &
!                                             "TNC",             &
!                                             "RMin",            &
!                                             "RMax",            &
!                                             "Rava",            &
!                                              "NAVA",           &
!                                             "MaxDiff",         &
!                                             "Concentrate"
!            RMIN = 1.D32
!            if(TBasicInfo%NC(p_ACTIVEFREE_STATU) .GT. 0) then
!                RMIN = min(RMIN,TMigStatInfo%RMIN(p_ACTIVEFREE_STATU))
!            else if(TBasicInfo%NC(p_ACTIVEINGB_STATU) .GT. 0) then
!                RMIN = min(RMIN,TMigStatInfo%RMIN(p_ACTIVEINGB_STATU))
!            end if
!
!            write(6, fmt= "(I15,1x,1PE15.4,1x,3(I15,1x),I15,1x,130(1PE15.4,1x))")   Record%GetSimuSteps(),                                                                   &
!                                                                                    Record%GetSimuTimes(),                                                                   &
!                                                                                    TBasicInfo%NC(p_ACTIVEFREE_STATU),                                                       &
!                                                                                    TBasicInfo%NC(p_ACTIVEINGB_STATU),                                                       &
!                                                                                    NCAct,                                                                                   &
!                                                                                    sum(TBasicInfo%NC),                                                                      &
!                                                                                    RMIN*C_CM2NM,                                                                            &
!                                                                                    max(TMigStatInfo%RMAX(p_ACTIVEFREE_STATU),TMigStatInfo%RMAX(p_ACTIVEINGB_STATU))*C_CM2NM,&
!                                                                                    RAVA*C_CM2NM,                                                                            &
!                                                                                    NAVA,                                                                                    &
!                                                                                    max(TMigStatInfo%DiffusorValueMax(p_ACTIVEFREE_STATU),TMigStatInfo%DiffusorValueMax(p_ACTIVEINGB_STATU)), &
!                                                                                    Concentrate
!
!
!        END ASSOCIATE
!
!        return
!    end subroutine PutOut_Instance_Statistic_IntegralBox
!
!    !*****************************************************************
!    subroutine PutOut_Instance_Statistic_EachBox(Host_Boxes,Host_SimuCtrlParam,TheMigCoaleStatisticInfo,Record)
!        implicit none
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(MigCoaleStatisticInfo_Used)::TheMigCoaleStatisticInfo
!        type(MigCoalClusterRecord)::Record
!        !---Local Vars---
!        integer::IBox
!        integer::MultiBox
!        real(kind=KMCDF)::RMIN
!        real(kind=KMCDF)::Concentrate
!        integer::NCAct
!        real(kind=KMCDF)::RAVA
!        real(kind=KMCDF)::NAVA
!        !---Body---
!        MultiBox = Host_SimuCtrlParam%MultiBox
!
!        DO IBox = 1,MultiBox
!
!            ASSOCIATE(SBasicInfo=>Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox),SMigStatInfo=>TheMigCoaleStatisticInfo%statistic_SingleBoxes(IBox))
!
!                NCAct = SBasicInfo%NC(p_ACTIVEFREE_STATU) + SBasicInfo%NC(p_ACTIVEINGB_STATU)
!
!                RAVA = (SBasicInfo%NC(p_ACTIVEFREE_STATU)*SMigStatInfo%RAVA(p_ACTIVEFREE_STATU) + SBasicInfo%NC(p_ACTIVEINGB_STATU)*SMigStatInfo%RAVA(p_ACTIVEINGB_STATU))/NCAct
!
!                if(NCAct .GT. 0) then
!                    NAVA = dble(SBasicInfo%NA(p_ACTIVEFREE_STATU) +  SBasicInfo%NA(p_ACTIVEINGB_STATU))/dble(NCAct)
!                else
!                    NAVA = 0.D0
!                end if
!
!                Concentrate = NCAct/Host_Boxes%BOXVOLUM
!
!                write(Record%HSizeStatistic_EachBox,fmt="(2(I20,1x),1PE20.4,1x,9(I20,1x),15(1PE20.4,1x),7(I20,1x))") Record%GetSimuSteps(),                   &
!                                                                                                             IBox,                                            &
!                                                                                                             Record%GetSimuTimes(),                           &
!                                                                                                             NCAct,                                           &
!                                                                                                             sum(SBasicInfo%NC),                              &
!                                                                                                             sum(SBasicInfo%NA(1:p_NUMBER_OF_STATU)),         &
!                                                                                                             SBasicInfo%NC(1:p_NUMBER_OF_STATU),              &
!                                                                                                             SMigStatInfo%RMIN(1:p_NUMBER_OF_STATU)*C_CM2NM,  &
!                                                                                                             SMigStatInfo%RMAX(1:p_NUMBER_OF_STATU)*C_CM2NM,  &
!                                                                                                             RAVA*C_CM2NM,                                    &
!                                                                                                             NAVA,                                            &
!                                                                                                             Concentrate,                                     &
!                                                                                                             Record%GetImplantedEntitiesNum(),                &
!                                                                                                             SBasicInfo%NC(p_ACTIVEFREE_STATU),               &
!                                                                                                             SBasicInfo%NC(p_ACTIVEINGB_STATU),               &
!                                                                                                             SBasicInfo%NC(p_OUT_DESTROY_STATU:p_ABSORBED_STATU) + &
!                                                                                                             Record%RecordNCBeforeSweepOut_SingleBox(IBox,p_OUT_DESTROY_STATU:p_ABSORBED_STATU)
!
!            END ASSOCIATE
!
!            call flush(Record%HSizeStatistic_EachBox)
!        END DO
!
!        return
!    end subroutine PutOut_Instance_Statistic_EachBox


end module MF_Method_MIGCOALE_CLUSTER_CPU
