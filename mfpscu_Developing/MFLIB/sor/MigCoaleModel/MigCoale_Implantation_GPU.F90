module MIGCOALE_IMPLANTATION_GPU
    use cudafor
    use MCLIB_CONSTANTS_GPU
    use MCLIB_TYPEDEF_USUAL
    use MCLIB_TYPEDEF_ACLUSTER
    use MCLIB_GLOBAL_GPU
    use MCLIB_TYPEDEF_SIMULATIONBOXARRAY_GPU
    use MCLIB_CAL_NEIGHBOR_LIST_GPU
    use MIGCOALE_TIMECTL
    use MIGCOALE_TYPEDEF_STATISTICINFO
    use MIGCOALE_STATISTIC_GPU
    use MIGCOALE_STATISTIC_CPU
    use MIGCOALE_TYPEDEF_SIMRECORD
    use MIGCOALE_ADDONDATA_HOST
    use MIGCOALE_GLOBALVARS_DEV
    implicit none

    character(len=11),private,parameter::m_IMPFINPUTF = "&IMPFINPUTF"

    ! Note: This is DISTOKMC18 different with OKMC_OUTCFG_FORMAT18, the main different is cause that for implant,
    !       what we need for the implant source is a discrete distribution(size and space), not a special configuration.
    !       So, to use the OKMC configuration for implant, we need to convert the OKMC configuration to a discrete distribution.
    !       However, for the mean-filed result, it is just the distribution that we need, so we need not to do any other things.
    character(len=11),parameter::OKMC_DIST_FORMAT18 = "&DISTOKMC18"
    character(len=9),private,parameter::SRIM_DIST = "&DISTSRIM"
    character(len=10),private,parameter::PANDA_DIST = "&DISTPANDA"

    integer,private,parameter::p_Implant_Hunger = 0
    integer,private,parameter::p_Implant_MemSaving = 1

    integer, parameter, private::p_ImplantConfig_Simple = 0
    integer, parameter, private::p_ImplantConfig_SpecialDistFromFile = 1
    integer, parameter, private::p_ImplantConfig_SpecialDistFromExteFunc = 2

    type,public::ImplantInfo_DevPart

        logical::InitFlag = .false.    ! the flag to record if the data structure had been initialization

        real(kind=KMCDF),device,allocatable,dimension(:)::Dev_CompositWeight

        real(kind=KMCDF),device,allocatable,dimension(:,:)::Dev_SUBBOXBOUNDARY

        real(kind=KMCDF),device,dimension(:),allocatable::Dev_LayerThick

        type(ACluster),device,dimension(:,:),allocatable::Dev_ClustersSample

        real(kind=KMCDF),device,dimension(:,:),allocatable::Dev_ClustersSampleRate

        contains
        procedure,non_overridable,public,pass::CopyImplantInfo_DevPartFromOther
        procedure,non_overridable,public,pass::Clean=>Clean_ImplantInfo_DevPart
        Generic::ASSIGNMENT(=)=>CopyImplantInfo_DevPartFromOther
        Final::CleanImplantInfo_DevPart
    end type ImplantInfo_DevPart

    type,public::ImplantSection
        integer::MemoryOccupyFactor = 10
        integer::ExpandFactor = 10
        real(kind=KMCDF)::ImplantFlux = 0

        integer::ImplantConfigType = -1

        character*256::ImplantCfgFileName = ""

        character*20::ImplantCfgFileType = ""

        character*10::Elemets(p_ATOMS_GROUPS_NUMBER) = ""

        real(kind=KMCDF)::NAINI = 0.D0

        real(kind=KMCDF)::NASDINI = 0.D0

        real(kind=KMCDF)::NACUT(2) = 0.D0

        real(kind=KMCDF)::CompositWeight(p_ATOMS_GROUPS_NUMBER) = 0.D0

        integer::ImplantDepthDistType = -1

        real(kind=KMCDF)::SUBBOXBOUNDARY(3,2) = 0.D0

        real(kind=KMCDF)::DepthINI = 0.D0

        real(kind=KMCDF)::DepthSDINI = 0.D0

        real(kind=KMCDF),dimension(:),allocatable::LayerThick

        real(kind=KMCDF),dimension(:,:),allocatable::ClustersSampleRate

        type(ACluster),dimension(:,:),allocatable::ClustersSample

        type(ImplantInfo_DevPart)::dm_ImplantInfo_DevPart

        contains
        procedure,non_overridable,public,pass::Load_ImplantSection
        procedure,non_overridable,public,pass::InitImplantInfo_DevPart
        procedure,non_overridable,public,pass::ReadImplantSection
        procedure,non_overridable,public,pass::ReadImplantSection_Simple
        procedure,non_overridable,public,pass::ReadImplantClusterSizeDist_Simple
        procedure,non_overridable,public,pass::ReadImplantClusterDepthDist_Simple
        procedure,non_overridable,public,pass::ReadImplantSection_SpecialDistFromFile
        procedure,non_overridable,public,pass::Putin_PANDA_OUTCFG_Distribution
        procedure,non_overridable,public,pass::Putin_SRIM2003_OUTCFG_Distribution
        procedure,non_overridable,public,pass::Putin_OKMC_FORMAT18_Distribution
        procedure,non_overridable,public,pass::ReadImplantSection_SpecialDistFromExteFunc
        procedure,non_overridable,public,pass::ImplantClusters_FastStrategy
        procedure,non_overridable,public,pass::AdjustTimeStep_Implant
        procedure,non_overridable,private,pass::Cal_ImplantANDExpandSize
        procedure,non_overridable,private,pass::DoImplantTillVirtualBoundary_CPUTOGPU
        procedure,non_overridable,private,pass::DoImplantTillVirtualBoundary_CPU
        procedure,non_overridable,private,pass::FillVirtualBoundary_CPU_Simple
        procedure,non_overridable,private,pass::FillVirtualBoundary_CPU_Depth_LAY
        procedure,non_overridable,private,pass::FillVirtualBoundary_CPU_Depth_SubBox
        procedure,non_overridable,private,pass::FillVirtualBoundary_CPU_Depth_Gauss
        procedure,non_overridable,private,pass::FillVirtualBoundary_CPU_FromFile
        procedure,non_overridable,private,pass::FillVirtualBoundary_CPU_FromExteFunc
        procedure,non_overridable,private,pass::DoImplantTillVirtualBoundary_GPUTOCPU
        procedure,non_overridable,private,pass::DoImplantTillVirtualBoundary_GPU
        procedure,non_overridable,private,pass::FillVirtualBoundary_GPU_Simple
        procedure,non_overridable,private,pass::FillVirtualBoundary_GPU_Depth_LAY
        procedure,non_overridable,private,pass::FillVirtualBoundary_GPU_Depth_SubBox
        procedure,non_overridable,private,pass::FillVirtualBoundary_GPU_Depth_Gauss
        procedure,non_overridable,private,pass::FillVirtualBoundary_GPU_FromFile
        procedure,non_overridable,private,pass::FillVirtualBoundary_GPU_FromExteFunc
        procedure,non_overridable,public,pass::CopyImplantSectionFromOther
        procedure,non_overridable,public,pass::Clean=>Clean_ImplantSection
        Generic::ASSIGNMENT(=)=>CopyImplantSectionFromOther
        Final::CleanImplantSection
    end type ImplantSection

    type,public::ImplantList
        type(ImplantSection)::TheImplantSection

        type(ImplantList),pointer::next=>null()

        integer::ListCount = 0

        contains
        procedure,non_overridable,public,pass::Init=>Init_ImplantList
        procedure,non_overridable,private,pass::Load_ImplantList
        procedure,non_overridable,private,pass::CheckImplantList
        procedure,non_overridable,public,pass::AppendOneSection=>AppendOne_ImplantSection
        procedure,non_overridable,public,pass::Get_P=>GetImplantSection_P
        procedure,non_overridable,public,pass::Clean=>Clean_ImplantList
        Final::CleanImplantList
    end type

    private::CopyImplantInfo_DevPartFromOther
    private::Clean_ImplantInfo_DevPart
    private::CleanImplantInfo_DevPart
    private::Load_ImplantSection
    private::InitImplantInfo_DevPart
    private::ReadImplantSection
    private::ReadImplantSection_Simple
    private::ReadImplantClusterSizeDist_Simple
    private::ReadImplantClusterDepthDist_Simple
    private::ReadImplantSection_SpecialDistFromFile
    private::Putin_PANDA_OUTCFG_Distribution
    private::Putin_SRIM2003_OUTCFG_Distribution
    private::Putin_OKMC_FORMAT18_Distribution
    private::ReadImplantSection_SpecialDistFromExteFunc
    private::ImplantClusters_FastStrategy
    private::AdjustTimeStep_Implant
    private::Cal_ImplantANDExpandSize
    private::DoImplantTillVirtualBoundary_CPUTOGPU
    private::DoImplantTillVirtualBoundary_CPU
    private::FillVirtualBoundary_CPU_Simple
    private::FillVirtualBoundary_CPU_Depth_LAY
    private::FillVirtualBoundary_CPU_Depth_SubBox
    private::FillVirtualBoundary_CPU_Depth_Gauss
    private::FillVirtualBoundary_CPU_FromFile
    private::FillVirtualBoundary_CPU_FromExteFunc
    private::DoImplantTillVirtualBoundary_GPUTOCPU
    private::DoImplantTillVirtualBoundary_GPU
    private::FillVirtualBoundary_GPU_Simple
    private::FillVirtualBoundary_GPU_Depth_LAY
    private::FillVirtualBoundary_GPU_Depth_SubBox
    private::FillVirtualBoundary_GPU_Depth_Gauss
    private::FillVirtualBoundary_GPU_FromFile
    private::FillVirtualBoundary_GPU_FromExteFunc
    private::CopyImplantSectionFromOther
    private::Clean_ImplantSection
    private::CleanImplantSection
    private::Init_ImplantList
    private::Load_ImplantList
    private::CheckImplantList
    private::AppendOne_ImplantSection
    private::GetImplantSection_P
    private::Clean_ImplantList
    private::CleanImplantList

    contains

    !********************For type ImplantInfo_DevPart**********************
    subroutine CopyImplantInfo_DevPartFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantInfo_DevPart),intent(out)::this
        TYPE(ImplantInfo_DevPart),intent(in)::other
        !---Local Vars--
        integer::err
        !---Body---

        call DeAllocateArray_GPU(this%Dev_CompositWeight,"Dev_CompositWeight")
        call AllocateArray_GPU(this%Dev_CompositWeight,p_ATOMS_GROUPS_NUMBER,"Dev_CompositWeight")
        err = cudaMemcpy(this%Dev_CompositWeight,other%Dev_CompositWeight,size(other%Dev_CompositWeight),cudaMemcpyDeviceToDevice)

        call DeAllocateArray_GPU(this%Dev_SUBBOXBOUNDARY,"Dev_SUBBOXBOUNDARY")
        call AllocateArray_GPU(this%Dev_SUBBOXBOUNDARY,3,2,"Dev_SUBBOXBOUNDARY")
        err = cudaMemcpy(this%Dev_SUBBOXBOUNDARY,other%Dev_SUBBOXBOUNDARY,size(other%Dev_SUBBOXBOUNDARY),cudaMemcpyDeviceToDevice)

        call DeAllocateArray_GPU(this%Dev_LayerThick,"Dev_LayerThick")
        call AllocateArray_GPU(this%Dev_LayerThick,size(other%Dev_LayerThick),"Dev_LayerThick")
        err = cudaMemcpy(this%Dev_LayerThick,other%Dev_LayerThick,size(other%Dev_LayerThick),cudaMemcpyDeviceToDevice)

        call DeAllocateArray_GPU(this%Dev_ClustersSample,"Dev_ClustersSample")
        call AllocateArray_GPU(this%Dev_ClustersSample,size(other%Dev_ClustersSample,dim=1),size(other%Dev_ClustersSample,dim=2),"Dev_ClustersSample")
        call copyClustersDevToDevSync2D(other%Dev_ClustersSample,this%Dev_ClustersSample,size(other%Dev_ClustersSample))

        call DeAllocateArray_GPU(this%Dev_ClustersSampleRate,"Dev_ClustersSampleRate")
        call AllocateArray_GPU(this%Dev_ClustersSampleRate,size(other%Dev_ClustersSampleRate,dim=1),size(other%Dev_ClustersSampleRate,dim=2),"Dev_ClustersSampleRate")
        err = cudaMemcpy(this%Dev_ClustersSampleRate,other%Dev_ClustersSampleRate,size(other%Dev_ClustersSampleRate),cudaMemcpyDeviceToDevice)

        return
    end subroutine

    !**********************************************************************
    subroutine Clean_ImplantInfo_DevPart(this)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantInfo_DevPart)::this
        !---Body---

        call DeAllocateArray_GPU(this%Dev_CompositWeight,"Dev_CompositWeight")

        call DeAllocateArray_GPU(this%Dev_SUBBOXBOUNDARY,"Dev_SUBBOXBOUNDARY")

        call DeAllocateArray_GPU(this%Dev_LayerThick,"Dev_LayerThick")

        call DeAllocateArray_GPU(this%Dev_ClustersSample,"Dev_ClustersSample")

        call DeAllocateArray_GPU(this%Dev_ClustersSampleRate,"Dev_ClustersSampleRate")

    end subroutine

    !**********************************************************************
    subroutine CleanImplantInfo_DevPart(this)
        implicit none
        !---Dummy Vars---
        TYPE(ImplantInfo_DevPart)::this
        !---Body---

        call this%Clean()
        return
    end subroutine

    !***************For type ImplantList************************************
    subroutine Init_ImplantList(this,Host_Boxes,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantList)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Body---

        call this%Clean()

        call this%Load_ImplantList(Host_Boxes,Host_SimuCtrlParam)

        call this%CheckImplantList(Host_SimuCtrlParam)

        return
    end subroutine

    !***********************************************************************
    subroutine CheckImplantList(this,Host_SimuCtrlParam)
         implicit none
        !---Dummy Vars---
        CLASS(ImplantList),target::this
        type(SimulationCtrlParam),target::Host_SimuCtrlParam
        !---Local Vars---
        type(SimulationCtrlParam),pointer::PSimuCtrlParamCursor=>Null()
        type(ImplantSection),pointer::PImplantSection=>null()
        integer::ICount
        !---Body---
        PSimuCtrlParamCursor=>Host_SimuCtrlParam

        ICount = 0
        DO While(associated(PSimuCtrlParamCursor))
            ICount = ICount + 1

            if(PSimuCtrlParamCursor%ImplantSectID .GE. 1) then
                PImplantSection=>this%Get_P(PSimuCtrlParamCursor%ImplantSectID)

                if(.not. associated(PImplantSection)) then
                    write(*,*) "MCPSCUERROR: The implantation section is not special :",PSimuCtrlParamCursor%ImplantSectID
                    write(*,*) "For the simulation section :",ICount
                    pause
                    stop
                end if

                if(PImplantSection%ImplantFlux .GT. 0.D0 .AND. PSimuCtrlParamCursor%NEIGHBORUPDATESTRATEGY .eq. mp_NEIGHBORUPDATEBYNCREMIND) then
                    write(*,*) "MCPSCUERROR: You cannot use the neighbor-list update strategy by clusters number remind percent when the implantation"
                    write(*,*) "flux exist."
                    write(*,*) "For the simulation section :",ICount
                    pause
                    stop
                end if

            end if

            PSimuCtrlParamCursor=>PSimuCtrlParamCursor%next
        END DO

        return
    end subroutine

    !***********************************************************************
    subroutine Load_ImplantList(this,Host_Boxes,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantList)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Local Vars---
        type(ImplantSection)::tempImplantSection
        character*256::truePath
        character*256::STR
        character*32::KEYWORD
        integer::hFile
        integer::LINE
        !---Body---
        LINE = 0

        truePath = INQUIREFILE(Host_Boxes%ImpFile)

        hFile = OpenExistedFile(truePath)

        call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")

        STR = adjustl(STR)

        call GETKEYWORD("&",STR,KEYWORD)

        call UPCASE(KEYWORD)

        if(.not. IsStrEqual(KEYWORD,m_IMPFINPUTF)) then
            write(*,*) "MCPSCUERROR: Unknown file header: ",KEYWORD
            write(*,*) "In file: ",truePath
            pause
            stop
        end if

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
            call RemoveComments(STR,"!")

            STR = adjustl(STR)

            call GETKEYWORD("&",STR,KEYWORD)

            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDIMPFINPUTF")
                    exit
                case("&GROUPSUBCTL")
                    call tempImplantSection%Clean()
                    call tempImplantSection%Load_ImplantSection(hFile,Host_Boxes,Host_SimuCtrlParam,LINE)

                    call this%AppendOneSection(tempImplantSection)
                case default
                    write(*,*) "MCPSCUERROR: Unknown flag: ",KEYWORD
                    write(*,*) "At Line: ",LINE
                    pause
                    stop
            end select

        END DO

        return

        100 write(*,*) "MCPSCUERROR: Fail to read the file: ",truePath
            write(*,*) "At Line: ",LINE
            pause
            stop
    end subroutine Load_ImplantList

    !***********************************************************************
    subroutine AppendOne_ImplantSection(this,TheImplantSection)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantList),target::this
        type(ImplantSection)::TheImplantSection
        !---Local Vars---
        type(ImplantList),pointer::cursor=>null()
        type(ImplantList),pointer::next=>null()
        !---Body---
        cursor=>this

        if(.not. associated(cursor)) then
            write(*,*) "MCPSCUERROR: You should allocate the ImplantList first!"
            pause
            stop
        end if

        if(this%ListCount .eq. 0) then
            ! The assignment(=) had been override
            this%TheImplantSection = TheImplantSection
        else
            cursor=>this
            next=>cursor%next

            Do While(associated(next))
                cursor=>next
                next=>cursor%next
            End Do

            allocate(next)
            ! The assignment(=) had been override
            next%TheImplantSection = TheImplantSection
            Nullify(next%next)
            cursor%next=>next
        end if

        this%ListCount = this%ListCount + 1

        return
    end subroutine

    !***********************************************************************
    function GetImplantSection_P(this,TheIndex) result(TheResult)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantList),target::this
        integer,intent(in)::TheIndex
        type(ImplantSection),intent(out),pointer::TheResult
        !---Local Vars---
        type(ImplantList),pointer::cursor=>null()
        integer::CountTemp
        !---Body---

        TheResult=>null()

        cursor=>this

        CountTemp = 0

        DO While(associated(cursor))


            CountTemp = CountTemp + 1

            if(CountTemp .eq. TheIndex) then
                TheResult=>cursor%TheImplantSection
                exit
            end if

            cursor=>cursor%next
        END DO

        Nullify(cursor)

        if(.not. associated(TheResult)) then
            write(*,*) "MCPSCUERROR: Cannot find the Implantation section by the id: ",TheIndex
            pause
            stop
        end if

        return
    end function GetImplantSection_P


    !***********************************************************************
    subroutine Clean_ImplantList(this)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantList),target::this
        !---Local Vars---
        type(ImplantList),pointer::cursor=>null()
        type(ImplantList),pointer::next=>null()
        !---Body---
        cursor=>this

        if(.not. associated(cursor)) then
            return
        end if

        cursor=>this%next

        call this%TheImplantSection%Clean()

        Do while(associated(cursor))
            next=>cursor%next
            call cursor%TheImplantSection%Clean()
            deallocate(cursor)
            Nullify(cursor)
            cursor=>next
        End Do

        this%next=>null()

        this%ListCount = 0
        Nullify(cursor)
        cursor=>null()
        Nullify(next)
        next=>null()

        return
    end subroutine Clean_ImplantList

    !***********************************************************************
    subroutine CleanImplantList(this)
        implicit none
        !---Dummy Vars---
        type(ImplantList)::this
        !---Body---

        call this%Clean()
        return
    end subroutine

    !*****************For Type ImplantSection****************************
    subroutine Load_ImplantSection(this,hFile,SimBoxes,Host_SimuCtrlParam,LINE)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        integer, intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer::LINE
        !---Local Vars---
        character*256::STR
        character*32::KEYWORD
        character*20::STRTMP(10)
        integer::N
        real(kind=KMCDF)::ReflectRatio
        !---Body---
        Do While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSUBCTL")
                    exit
                case("&TYPE")
                    call EXTRACT_NUMB(STR,1,N,STRTMP)
                    if(N .LT. 1) then
                        write(*,*) "MCPSCUERROR: Too few parameters for implantation distribution type."
                        write(*,*) "At Line :", LINE
                        write(*,*) "You should special by the way : &TYPE The implantation cluster distribution type =  "
                        pause
                        stop
                    end if
                    this%ImplantConfigType = ISTR(STRTMP(1))
                    exit
                case default
                    write(*,*) "MCPSCUERROR: You must special the implantation distribution type first!"
                    write(*,*) "By the way: &TYPE The implantation cluster distribution type = "
                    write(*,*) "However, the words you input is: ",STR
                    pause
                    stop
            end select
        End Do

        Do While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSUBCTL")
                    exit
                case("&FLUX")
                    call EXTRACT_NUMB(STR,2,N,STRTMP)
                    if(N .LT. 2) then
                        write(*,*) "MCPSCUERROR: Too few parameters for the implantation flux ."
                        write(*,*) "At Line :", LINE
                        write(*,*) "You should special by the way: &FLUX THE implantation flux = , the reflect ratio = "
                        pause
                        stop
                    end if
                    this%ImplantFlux = DRSTR(STRTMP(1))
                    ReflectRatio = DRSTR(STRTMP(2))

                    if(ReflectRatio .LT. 0 .or. ReflectRatio .GT. 1.D0) then
                        write(*,*) "MCPSCUERROR: The reflect ratio should between 0 and 1"
                        write(*,*) "However, the current reflect ratio is: ",ReflectRatio
                        pause
                        stop
                    end if
                    this%ImplantFlux = this%ImplantFlux*(1.D0-ReflectRatio)

                case("&FEXPAND")
                    call EXTRACT_NUMB(STR,1,N,STRTMP)
                    if(N .LT. 1) then
                        write(*,*) "MCPSCUERROR: Too few parameters for the implantation expand factor ."
                        write(*,*) "At Line :", LINE
                        write(*,*) "You should special by the way: &FEXPAND The expand size factor = "
                        pause
                        stop
                    end if
                    this%ExpandFactor = ISTR(STRTMP(1))

                case("&FMEMOCCUP")
                    call EXTRACT_NUMB(STR,1,N,STRTMP)
                    if(N .LT. 1) then
                        write(*,*) "MCPSCUERROR: Too few parameters for the memory occupy factor."
                        write(*,*) "At Line :", LINE
                        write(*,*) "You should special by the way: &FMEMOCCUP TThe memory occupied factor ="
                        pause
                        stop
                    end if
                    this%MemoryOccupyFactor = ISTR(STRTMP(1))

                    if(this%MemoryOccupyFactor .LE. 1) then
                        write(*,*) "MCPSCUERROR: The MemoryOccupyFactor cannot less than 1"
                        pause
                        stop
                    end if

                case("&SIZESUBCTL","&DEPTHSUBCTL","&EXTFSUBCTL")
                    call this%ReadImplantSection(hFile,KEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)

                case default
                    write(*,*) "MCPSCUERROR: Unknown Flag: ",KEYWORD
                    write(*,*) "At LINE: ",LINE
                    pause
                    stop
            end select
        End Do

        return

        100 write(*,*) "MCPSCUERROR : Load implantation configuration file failed !"
            write(*,*) "At line :",LINE
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine

    !****************************************************************
    subroutine InitImplantInfo_DevPart(this)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        !---Body---
        if(this%dm_ImplantInfo_DevPart%InitFlag .eq. .false.) then

            this%dm_ImplantInfo_DevPart%InitFlag = .true.

            call AllocateArray_GPU(this%dm_ImplantInfo_DevPart%Dev_CompositWeight,size(this%CompositWeight),"Dev_CompositWeight")
            this%dm_ImplantInfo_DevPart%Dev_CompositWeight = this%CompositWeight

            call AllocateArray_GPU(this%dm_ImplantInfo_DevPart%Dev_SUBBOXBOUNDARY,size(this%SUBBOXBOUNDARY,DIM=1),size(this%SUBBOXBOUNDARY,DIM=2),"Dev_SUBBOXBOUNDARY")
            this%dm_ImplantInfo_DevPart%Dev_SUBBOXBOUNDARY = this%SUBBOXBOUNDARY

            call AllocateArray_GPU(this%dm_ImplantInfo_DevPart%Dev_LayerThick,size(this%LayerThick),"Dev_LayerThick")
            this%dm_ImplantInfo_DevPart%Dev_LayerThick = this%LayerThick

            call AllocateArray_GPU(this%dm_ImplantInfo_DevPart%Dev_ClustersSample,size(this%ClustersSample,dim=1),size(this%ClustersSample,dim=2),"Dev_ClustersSample")
            call copyInClustersSync2D(this%ClustersSample,this%dm_ImplantInfo_DevPart%Dev_ClustersSample,size(this%ClustersSample))

            call AllocateArray_GPU(this%dm_ImplantInfo_DevPart%Dev_ClustersSampleRate,size(this%ClustersSampleRate,dim=1),size(this%ClustersSampleRate,dim=2),"Dev_ClustersSampleRate")
            this%dm_ImplantInfo_DevPart%Dev_ClustersSampleRate = this%ClustersSampleRate

        end if

        return
    end subroutine InitImplantInfo_DevPart

    !***************************************************************
    subroutine ReadImplantSection(this,hFile,KEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        integer,intent(in)::hFile
        character*(*)::KEYWORD
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer::LINE
        !--Body---
        select case(this%ImplantConfigType)
            case(p_ImplantConfig_Simple)
                call this%ReadImplantSection_Simple(hFile,KEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)
            case(p_ImplantConfig_SpecialDistFromFile)
                call this%ReadImplantSection_SpecialDistFromFile(hFile,KEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)
            case(p_ImplantConfig_SpecialDistFromExteFunc)
                call this%ReadImplantSection_SpecialDistFromExteFunc(hFile,KEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)
            case default
                write(*,*) "MCPSCUERROR: Unknown strategy for the implantation configuration:",this%ImplantConfigType
                pause
                stop
        end select

        return
    end subroutine ReadImplantSection

    !*****************************************************************
    subroutine ReadImplantSection_Simple(this,hFile,KEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        integer,intent(in)::hFile
        character*(*)::KEYWORD
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer::LINE
        !---Body---

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                return
            case("&SIZESUBCTL")
                    call this%ReadImplantClusterSizeDist_Simple(hFile,SimBoxes,LINE)
            case("&DEPTHSUBCTL")
                    call this%ReadImplantClusterDepthDist_Simple(hFile,SimBoxes,LINE)
            case default
                write(*,*) "MCPSCUERROR: The Illegal flag: ",KEYWORD
                pause
                stop
        end select

        return

        100 write(*,*) "MCPSCUERROR : Load implantation configuration file failed !"
            write(*,*) "At line :",LINE
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine ReadImplantSection_Simple

    !****************************************************************
    subroutine ReadImplantSection_SpecialDistFromFile(this,hFile,PreKEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        integer,intent(in)::hFile
        character*(*)::PreKEYWORD
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer::LINE
        !---Local Vars---
        character*256::STR
        character*32::KEYWORD
        character*256::STRTEMP(1)
        integer::N
        type(MigCoalClusterRecord)::tempRecord
        real(kind=KMCDF)::TotalSampleRate
        integer::LayerNum
        !---Body---

        if(.not. IsStrEqual(PreKEYWORD,"&EXTFSUBCTL")) then
            write(*,*) "MCPSCUERROR: You must special the &EXTFSUBCTL when the implant strategy is chosen by outer file ."
            write(*,*) "However, you had special the key word :",KEYWORD
            write(*,*) "At line: ",LINE
            pause
            stop
        end if

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
            call RemoveComments(STR,"!")
            STR = adjustl(STR)
            call GETKEYWORD("&",STR,KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:LENTRIM(KEYWORD)))
                case("&ENDSUBCTL")
                    exit
                case("&DISTFILETYPE")
                    call EXTRACT_SUBSTR(STR,1,N,STRTEMP)
                    if(N .LT. 1) then
                        write(*,*) "MCPSCUERROR: You must special the implantation configuration file type if you had chosen the file model."
                        write(*,*) "By the way: &DISTFILETYPE The distribution file type = ! 'DISTOKMC18','BOXMF18','BOXSPMF18','DISTSRIM','DISTPANDA' "
                        pause
                        stop
                    end if

                    if(LENTRIM(STRTEMP(1)) .LE. 0) then
                        write(*,*) "MCPSCUERROR: The implant configuration file type is null."
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if

                    call UPCASE(STRTEMP(1))

                    this%ImplantCfgFileType = adjustl(trim(KEYWORD_HEAD))//adjustl(trim(STRTEMP(1)))

                    if(IsStrEqual(trim(this%ImplantCfgFileType),SRIM_DIST) .or.  IsStrEqual(trim(this%ImplantCfgFileType),PANDA_DIST)) then
                        call EXTRACT_SUBSTR(STR,2,N,STRTEMP)

                        if(N .LT. 2) then
                            write(*,*) "MCPSCUERROR: when the specialized distribution type is 'DISTSRIM' or 'DISTPANDA'"
                            write(*,*) "MCPSCUERROR: you must special the implantation ion type."
                            pause
                            stop
                        end if

                        ! For both of PANDA and SRIM, only one kind of ion are injected to matrix
                        call UPCASE(STRTEMP(2))
                        this%Elemets(1) = adjustl(trim((STRTEMP(2))))
                    end if

                    if(IsStrEqual(trim(this%ImplantCfgFileType),SRIM_DIST) .or. IsStrEqual( trim(this%ImplantCfgFileType),OKMC_DIST_FORMAT18)) then
                        call EXTRACT_NUMB(STR,1,N,STRTEMP)

                        if(N .LT. 1) then
                            write(*,*) "MCPSCUERROR: when the specialized distribution type is 'DISTSRIM' or 'DISTOKMC18'"
                            write(*,*) "MCPSCUERROR: you must special the layer number that you want to divide."
                            pause
                            stop
                        end if

                        LayerNum = ISTR(STRTEMP(1))

                        if(LayerNum .LE. 0) then
                            write(*,*) "MCPSCUERROR: the total layer number cannot be less than 0 when it is set for SRIM or OKMC18 distribution"
                            pause
                            stop
                        end if
                    end if


                case("&DISTFILE")
                    call EXTRACT_SUBSTR(STR,1,N,STRTEMP)
                    if(N .LT. 1) then
                        write(*,*) "MCPSCUERROR: You must special the implantation configuration file if you had chosen the file model."
                        write(*,*) "By the way: &DISTFILE The distribution file path = "
                        pause
                        stop
                    end if

                    if(LENTRIM(STRTEMP(1)) .LE. 0) then
                        write(*,*) "MCPSCUERROR: The implant configuration file name is null."
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if

                    this%ImplantCfgFileName = INQUIREFILE(STRTEMP(1),Host_SimuCtrlParam%InputFilePath)

                case default
                    write(*,*) "MCPSCUERROR: Illegal flag: ",KEYWORD
                    pause
                    stop
            end select
        END DO

        select case(adjustl(trim(this%ImplantCfgFileType)))
            case(OKMC_DIST_FORMAT18)
                call this%Putin_OKMC_FORMAT18_Distribution(LayerNum)

            case(MF_OUTCFG_FORMAT18)
                call SimBoxes%Putin_MF_OUTCFG_FORMAT18_Distribution(Host_SimuCtrlParam,this%ImplantCfgFileName,this%LayerThick,this%ClustersSampleRate,this%ClustersSample,tempRecord,m_RNFACTOR,m_FREESURDIFPRE)

            case(SPMF_OUTCFG_FORMAT18)
                call SimBoxes%Putin_SPMF_OUTCFG_FORMAT18_Distribution(Host_SimuCtrlParam,this%ImplantCfgFileName,this%LayerThick,this%ClustersSampleRate,this%ClustersSample,tempRecord,m_RNFACTOR,m_FREESURDIFPRE,m_GBSURDIFPRE)

            case(SRIM_DIST)
                call this%Putin_SRIM2003_OUTCFG_Distribution(SimBoxes,Host_SimuCtrlParam,LayerNum)

            case(PANDA_DIST)
                call this%Putin_PANDA_OUTCFG_Distribution(SimBoxes,Host_SimuCtrlParam)

            case default
                write(*,*) "MCPSCUERROR: Unknown Implant Configuration file type : ",this%ImplantCfgFileType
                write(*,*) "In current version, only the ", &
                            OKMC_DIST_FORMAT18," ",         &
                            MF_OUTCFG_FORMAT18," ",         &
                            SPMF_OUTCFG_FORMAT18," ",       &
                            SRIM_DIST," ",                  &
                            PANDA_DIST," ",                 &
                            "are supported. "
                pause
                stop
        end select

        ! Note, the out put SampleRate may be the concentrate, we need to convert it to rate now.
        TotalSampleRate = sum(this%ClustersSampleRate)
        if(TotalSampleRate .LE. 0) then
            write(*,*) "MCPSCUERROR: The total concentrate cannot less equal with 0"
            pause
            stop
        end if
        this%ClustersSampleRate = this%ClustersSampleRate/TotalSampleRate

        return
        100 write(*,*) "MCPSCUERROR : Load implantation configuration file failed !"
            write(*,*) "At line :",LINE
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine ReadImplantSection_SpecialDistFromFile

    !*****************************************************************
    subroutine Putin_PANDA_OUTCFG_Distribution(this,SimBoxes,Host_SimuCtrlParam)
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes),intent(in)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Local Vars---
        integer::hFile
        integer::LINE
        integer::N
        character*256::STR
        character*32::KEYWORD
        character*32::STRTMP(20)
        integer::IElement
        integer::LayerNum
        integer::ILayer
        integer::Layer
        real(kind=KMCDF)::SumOfThick
        type(ACluster)::ImplantIon
        type(DiffusorValue)::TheDiffusorValue
        !---Body---

        LINE = 0

        hFile = OpenExistedFile(this%ImplantCfgFileName)

        LayerNum = 0

        LINE = 0

        DO While(.not. GETINPUTSTRLINE_New(hFile,STR,LINE,"!"))
            call RemoveComments(STR,"!")

            if(LENTRIM(adjustl(STR)) .LE. 0) then
                cycle
            end if

            LINE = LINE + 1

            LayerNum = LayerNum + 1

        END DO

        if(LayerNum .LE. 0) then
            write(*,*) "MCPSCUERROR: The layers number in panda distribution is less than 1, that is impossible."
            pause
            stop
        end if

        call AllocateArray_Host(this%LayerThick,LayerNum,"LayersThick")
        call AllocateArray_Host(this%ClustersSampleRate,LayerNum,1,"ClustersSampleRate")
        call AllocateArray_Host(this%ClustersSample,LayerNum,1,"ClustersSample")

        !---For SRIM2003 and PANDA, only one kind of ions would be implanted
        DO IElement = 1,p_ATOMS_GROUPS_NUMBER
            ImplantIon%m_Atoms(IElement)%m_ID = IElement
            ImplantIon%m_Atoms(IElement)%m_NA = 0
        END DO
        IElement = SimBoxes%Atoms_list%FindIndexBySymbol(this%Elemets(1))
        ImplantIon%m_Atoms(IElement)%m_ID = IElement
        ImplantIon%m_Atoms(IElement)%m_NA = 1

        TheDiffusorValue = SimBoxes%m_DiffusorTypesMap%Get(ImplantIon)

        !---In PANDA, the matrix is amorphous---
        select case(TheDiffusorValue%ECRValueType_Free)
            case(p_ECR_ByValue)
                ImplantIon%m_RAD = TheDiffusorValue%ECR_Free
            case(p_ECR_ByBCluster)
                ! Convert the number of atoms to radius
                ! Ref. Modelling Simul. Mater. Sci. Eng.16(2008)055003
                ImplantIon%m_RAD = DSQRT(sum(ImplantIon%m_Atoms(:)%m_NA)/m_RNFACTOR)
        end select

        select case(TheDiffusorValue%DiffusorValueType_Free)
            case(p_DiffuseCoefficient_ByValue)
                ImplantIon%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
            case(p_DiffuseCoefficient_ByArrhenius)
                ImplantIon%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
            case(p_DiffuseCoefficient_ByBCluster)
                ! Here we adopt a model that D=D0*(1/R)**Gama
                ImplantIon%m_DiffCoeff = m_FREESURDIFPRE*(ImplantIon%m_RAD**(-p_GAMMA))
        end select

        ImplantIon%m_Statu = p_ACTIVEFREE_STATU

        DO Layer = 1,LayerNum
            this%ClustersSample(Layer,1) = ImplantIon
        END DO


        ReWind(hFile)

        ILayer = 1

        SumOfThick = 0.D0

        LINE = 0

        DO While(.not. GETINPUTSTRLINE_New(hFile,STR,LINE,"!"))
            call RemoveComments(STR,"!")

            if(LENTRIM(adjustl(STR)) .LE. 0) then
                cycle
            end if

            LINE = LINE + 1

            call EXTRACT_NUMB(STR,2,N,STRTMP)

            if(N .LT. 2) then
                write(*,*) "MCPSCUERROR: The panda distribution file data cannot be recognized in line: ",LINE
                write(*,*) STR
                write(*,*) "At file: ",this%ImplantCfgFileName
                pause
                stop
            end if

            this%LayerThick(ILayer) = 2*(DRSTR(STRTMP(1))*C_UM2CM - SumOfThick)
            SumOfThick = SumOfThick + this%LayerThick(ILayer)

            if(SumOfThick .GT. SimBoxes%BOXSIZE(3)) then
                write(*,*) "MCPSCUERROR: The PANDA depth distribution is greater than simulation box depth."
                write(*,*) "The PANDA depth is: ",SumOfThick
                write(*,*) "The simulation box depth is : ",SimBoxes%BOXSIZE(3)
                pause
                stop
            end if

            this%ClustersSampleRate(ILayer,1) = DRSTR(STRTMP(2))

            this%ClustersSample(ILayer,1)%m_Statu = p_ACTIVEFREE_STATU

            this%ClustersSample(ILayer,1)%m_Layer = ILayer

            ILayer =  ILayer + 1

        END DO

        close(hFile)

        return
    end subroutine Putin_PANDA_OUTCFG_Distribution

    !******************************************************************
    subroutine Putin_SRIM2003_OUTCFG_Distribution(this,SimBoxes,Host_SimuCtrlParam,LayerNum)
                !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes),intent(in)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer,intent(in)::LayerNum
        !---Local Vars---
        integer::hFile
        integer::LINE
        integer::N
        character*256::STR
        character*32::KEYWORD
        character*32::STRTMP(20)
        integer::IElement
        integer::StoppedNum
        integer::IIon
        integer::ILayer
        real(kind=KMCDF)::SumOfThick
        type(ACluster)::ImplantIon
        type(DiffusorValue)::TheDiffusorValue
        real(kind=KMCDF),dimension(:,:),allocatable::StoppedPosition
        real(kind=KMCDF)::Thickness
        !---Body---
        LINE = 0

        hFile = OpenExistedFile(this%ImplantCfgFileName)

        StoppedNum = 0

        DO While(.not. GETINPUTSTRLINE_New(hFile,STR,LINE,"!"))
            call RemoveComments(STR,"!")

            if(LENTRIM(adjustl(STR)) .LE. 0) then
                cycle
            end if

            LINE = LINE + 1

            if((iachar(STR(1:1)) .GE. iachar('0') .AND. iachar(STR(1:1)) .LE. iachar('9'))) then
                StoppedNum = StoppedNum + 1
            end if

        END DO

        if(StoppedNum .LE. 0) then
            write(*,*) "MCPSCUERROR: The bins number in SRIM2003 distribution is less than 1, that is impossible."
            pause
            stop
        end if

        call AllocateArray_Host(StoppedPosition,StoppedNum,3,"StoppedPosition")

        ReWind(hFile)

        LINE = 0

        IIon = 1

        DO While(.not. GETINPUTSTRLINE_New(hFile,STR,LINE,"!"))
            call RemoveComments(STR,"!")

            if(LENTRIM(adjustl(STR)) .LE. 0) then
                cycle
            end if

            LINE = LINE + 1

            if((iachar(STR(1:1)) .GE. iachar('0') .AND. iachar(STR(1:1)) .LE. iachar('9'))) then

                call EXTRACT_NUMB(STR,4,N,STRTMP)

                if(N .LT. 4) then
                    write(*,*) "MCPSCUERROR: The SRIM2003 distribution file data cannot be recognized in line: ",LINE
                    write(*,*) STR
                    write(*,*) "At file: ",this%ImplantCfgFileName
                    pause
                    stop
                end if

                StoppedPosition(IIon,1) = DRSTR(STRTMP(2))*C_AM2CM ! depth   X
                StoppedPosition(IIon,2) = DRSTR(STRTMP(3))*C_AM2CM ! lateral Y
                StoppedPosition(IIon,3) = DRSTR(STRTMP(4))*C_AM2CM ! lateral Z

                if(StoppedPosition(IIon,2) .LT. SimBoxes%BOXBOUNDARY(1,1) .or. StoppedPosition(IIon,2) .GT. SimBoxes%BOXBOUNDARY(1,2)) then
                    write(*,*) "MCPSCUERROR: The SRIM2003 distribution is out of the simulation box in lateral X."
                    write(*,*) STR
                    write(*,*) "Current position in lateral X is (cm) : ",StoppedPosition(IIon,2)
                    write(*,*) "However, the box boundary range from ",SimBoxes%BOXBOUNDARY(1,1)," To ",SimBoxes%BOXBOUNDARY(1,2)
                    pause
                    stop
                end if

                if(StoppedPosition(IIon,3) .LT. SimBoxes%BOXBOUNDARY(2,1) .or. StoppedPosition(IIon,3) .GT. SimBoxes%BOXBOUNDARY(2,2)) then
                    write(*,*) "MCPSCUERROR: The SRIM2003 distribution is out of the simulation box in lateral Y."
                    write(*,*) STR
                    write(*,*) "Current position in lateral Y is (cm) : ",StoppedPosition(IIon,3)
                    write(*,*) "However, the box boundary range from ",SimBoxes%BOXBOUNDARY(2,1)," To ",SimBoxes%BOXBOUNDARY(2,2)
                    pause
                    stop
                end if

                if(StoppedPosition(IIon,1) .LT. SimBoxes%BOXBOUNDARY(3,1) .or. StoppedPosition(IIon,1) .GT. SimBoxes%BOXBOUNDARY(3,2)) then
                    write(*,*) "MCPSCUERROR: The SRIM2003 distribution is out of the simulation box in depth."
                    write(*,*) STR
                    write(*,*) "Current depth is (cm) : ",StoppedPosition(IIon,1)
                    write(*,*) "However, the box boundary range from ",SimBoxes%BOXBOUNDARY(3,1)," To ",SimBoxes%BOXBOUNDARY(3,2)
                    pause
                    stop
                end if

                IIon = IIon + 1
            end if

        END DO

        call AllocateArray_Host(this%LayerThick,LayerNum,"LayersThick")
        call AllocateArray_Host(this%ClustersSampleRate,LayerNum,1,"ClustersSampleRate")
        call AllocateArray_Host(this%ClustersSample,LayerNum,1,"ClustersSample")

        Thickness = (maxval(StoppedPosition(:,1)) - minval(StoppedPosition(:,1)))/LayerNum

        this%LayerThick = Thickness

        this%ClustersSampleRate = 0

        DO IIon = 1,StoppedNum
            ILayer = max(floor(StoppedPosition(IIon,1)/Thickness),1)
            this%ClustersSampleRate(ILayer,1) = this%ClustersSampleRate(ILayer,1) + 1
        END DO

        !---For SRIM2003 and PANDA, only one kind of ions would be implanted
        DO IElement = 1,p_ATOMS_GROUPS_NUMBER
            ImplantIon%m_Atoms(IElement)%m_ID = IElement
            ImplantIon%m_Atoms(IElement)%m_NA = 0
        END DO
        IElement = SimBoxes%Atoms_list%FindIndexBySymbol(this%Elemets(1))
        ImplantIon%m_Atoms(IElement)%m_ID = IElement
        ImplantIon%m_Atoms(IElement)%m_NA = 1

        TheDiffusorValue = SimBoxes%m_DiffusorTypesMap%Get(ImplantIon)

        !---In SRIM2003, the matrix is amorphous---
        select case(TheDiffusorValue%ECRValueType_Free)
            case(p_ECR_ByValue)
                ImplantIon%m_RAD = TheDiffusorValue%ECR_Free
            case(p_ECR_ByBCluster)
                ! Convert the number of atoms to radius
                ! Ref. Modelling Simul. Mater. Sci. Eng.16(2008)055003
                ImplantIon%m_RAD = DSQRT(sum(ImplantIon%m_Atoms(:)%m_NA)/m_RNFACTOR)
        end select

        select case(TheDiffusorValue%DiffusorValueType_Free)
            case(p_DiffuseCoefficient_ByValue)
                ImplantIon%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
            case(p_DiffuseCoefficient_ByArrhenius)
                ImplantIon%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
            case(p_DiffuseCoefficient_ByBCluster)
                ! Here we adopt a model that D=D0*(1/R)**Gama
                ImplantIon%m_DiffCoeff = m_FREESURDIFPRE*(ImplantIon%m_RAD**(-p_GAMMA))
        end select

        ImplantIon%m_Statu = p_ACTIVEFREE_STATU

        this%ClustersSample(:,:) = ImplantIon

        call DeAllocateArray_Host(StoppedPosition,"StoppedPosition")

        close(hFile)

        return
    end subroutine Putin_SRIM2003_OUTCFG_Distribution

    !******************************************************************
    subroutine Putin_OKMC_FORMAT18_Distribution(this,LayerNum)
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        integer,intent(in)::LayerNum

    end subroutine Putin_OKMC_FORMAT18_Distribution

    !*****************************************************************
    subroutine ReadImplantSection_SpecialDistFromExteFunc(this,hFile,KEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        integer,intent(in)::hFile
        character*(*)::KEYWORD
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer::LINE

        ! @todo (zhail#1#):


        return
    end subroutine ReadImplantSection_SpecialDistFromExteFunc

    !***************************************************************
    subroutine PrePareImplantSection(this,SimBoxes,Host_SimuCtrlParam)
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Local Vars---
        integer::hFile
        character*256::STR
        character*32::KEYWORD
        integer::LINE
        type(MigCoalClusterRecord)::tempRecord
        real(kind=KMCDF)::TotalSampleRate
        !---Body---



        return
        100 write(*,*) "MCPSCUERROR: Fail to load the implant distribution at file: ",this%ImplantCfgFileName
            write(*,*) "At line: ",LINE
            write(*,*) STR
            pause
            stop
    end subroutine PrePareImplantSection

    !***************************************************************
    subroutine ReadImplantClusterSizeDist_Simple(this,hFile,SimBoxes,LINE)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        integer,intent(in)::hFile
        type(SimulationBoxes)::SimBoxes
        integer::LINE
        !---Local Vars---
        character*512::STR
        character*32::KEYWORD
        character*32::STRTMP(10)
        integer::N
        integer::NElements
        integer::I
        integer::TheIndex
        real(kind=KMCDF)::TotalSampleRate
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
                case("&NATOMDIST")
                    call EXTRACT_NUMB(STR,4,N,STRTMP)
                    if(N .LT. 4) then
                        write(*,*) "MCPSCUERROR: Too few parameters for the cluster size distribution"
                        write(*,*) "You should special: &NATOMDIST The atoms number in each cluster distribution as Gauss that central = , distribution half width = , left cut = ,right cut = "
                        pause
                        stop
                    end if
                    this%NAINI = DRSTR(STRTMP(1))
                    this%NASDINI  = DRSTR(STRTMP(2))
                    this%NACUT(1) = DRSTR(STRTMP(3))
                    this%NACUT(2) = DRSTR(STRTMP(4))
                    if(this%NACUT(1) .GE. this%NACUT(2)) then
                        write(*,*) "MCPSCUERROR: The right cut cannot less than left cut."
                        write(*,*) "LCut",this%NACUT(1)
                        write(*,*) "RCut",this%NACUT(2)
                        pause
                        stop
                    end if

                case("&ELEMENTCOMPOSIT")
                    call EXTRACT_SUBSTR(STR,p_ATOMS_GROUPS_NUMBER,NElements,STRTMP)
                    if(NElements .LE. 0) then
                        write(*,*) "MCPSCUERROR: None of atoms kind (Elements) are specialized "
                        write(*,*) "You should special like that : &ELEMENTCOMPOSIT The included element = 'A', 'B' ."
                        pause
                        stop
                    else if(NElements .GT. p_ATOMS_GROUPS_NUMBER) then
                        write(*,*) "MCPSCUERROR: the specialized elements kinds is : ",N
                        write(*,*) "which is great than the max permitted elements kinds :",p_ATOMS_GROUPS_NUMBER
                        pause
                        stop
                    else
                        DO I = 1,NElements
                            this%Elemets(I) = adjustl(trim(STRTMP(I)))
                            call UPCASE(this%Elemets(I))
                        END DO
                    end if

                    call EXTRACT_NUMB(STR,p_ATOMS_GROUPS_NUMBER,N,STRTMP)
                    if(N .ne. NElements) then
                        write(*,*) "MCPSCUERROR: The elements weights number is not equal with the elements kinds which given."
                        write(*,*) "The elements kinds number is :",NElements
                        write(*,*) "But the weights number is :",N
                        pause
                        stop
                    else
                        this%CompositWeight = 0.D0

                        DO I = 1,N
                            TheIndex = SimBoxes%Atoms_list%FindIndexBySymbol(this%Elemets(I))
                            this%CompositWeight(TheIndex) = DRSTR(STRTMP(I))
                        END DO

                        if(sum(this%CompositWeight) .LE. 0.D0) then
                            write(*,*) "MCPSCUERROR: The sum of elements weights must great than 0 ."
                            write(*,*) STR
                            write(*,*) "At Line :",LINE
                            pause
                            stop
                        end if

                        this%CompositWeight = this%CompositWeight/sum(this%CompositWeight)
                    end if
                CASE default
                    write(*,*) "MCPSCUERROR: Illegal Symbol: ", KEYWORD
                    pause
                    stop
            END SELECT

        END DO

        return

        100 write(*,*) "MCPSCUERROR : Load implantation configuration file failed for cluster size!"
            write(*,*) "At line :",LINE
            write(*,*) STR
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine ReadImplantClusterSizeDist_Simple

    !***************************************************************
    subroutine ReadImplantClusterDepthDist_Simple(this,hFile,Host_SimBoxes,LINE)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        integer,intent(in)::hFile
        type(SimulationBoxes)::Host_SimBoxes
        integer::LINE
        !---Local Vars---
        character*256::STR
        character*32::KEYWORD
        character*32::STRTMP(10)
        integer::N
        integer::I
        integer::LayerNum
        real(kind=KMCDF)::SumOfLayer
        real(kind=KMCDF)::TotalSampleRate
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
                CASE("&DEPTH_LAYER")
                    this%ImplantDepthDistType = p_DEPT_DIS_Layer

                    call EXTRACT_NUMB(STR,1,N,STRTMP)
                    if(N .LT. 1) then
                        write(*,*) "MCPSCUERROR: Too few parameters for the clusters depth distribution layer type"
                        write(*,*) "You should special: &DEPTH_LAYER THE NUMBER OF DEPTH DISTRIBUTION LAYER = , THE ENTRIES DISTRIBUTION ="
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if

                    LayerNum = ISTR(STRTMP(1))
                    if(LayerNum .LT. 1) then
                        write(*,*) "MCPSCUERROR: The layer number should greater than 1"
                        write(*,*) "At line :",LINE
                        write(*,*) STR
                        pause
                        stop
                    end if

                    call EXTRACT_NUMB(STR,LayerNum*2+1,N,STRTMP)

                    if((N-1) .NE. LayerNum*2) then
                        write(*,*) "MCPSCUERROR: the specialeld layer is not equal with your setting"
                        write(*,*) STR
                        write(*,*) "At line :", LINE
                        pause
                        stop
                    end if

                    call AllocateArray_Host(this%LayerThick,LayerNum,"LayerThick")
                    this%LayerThick = 0.D0
                    call AllocateArray_Host(this%ClustersSampleRate,LayerNum,1,"ClustersSampleRate")
                    this%ClustersSampleRate = 0.D0

                    DO I = 1,LayerNum
                        this%LayerThick(I) = DRSTR(STRTMP(I+1))
                    END DO

                    if(sum(this%LayerThick) .LE. 0) then
                        this%LayerThick(1) = 1.D0
                    end if

                    SumOfLayer = sum(this%LayerThick)

                    DO I = 1,LayerNum
                        this%LayerThick(I) = Host_SimBoxes%BOXSIZE(3)*this%LayerThick(I)/SumOfLayer
                    END DO

                    DO I = 1,LayerNum
                        this%ClustersSampleRate(I,1) = DRSTR(STRTMP(I + LayerNum + 1))
                    END DO

                    ! Note, the out put SampleRate may be the concentrate, we need to convert it to rate now.
                    TotalSampleRate = sum(this%ClustersSampleRate)
                    if(TotalSampleRate .LE. 0) then
                        write(*,*) "MCPSCUERROR: The total concentrate cannot less equal with 0"
                        pause
                        stop
                    end if
                    this%ClustersSampleRate = this%ClustersSampleRate/TotalSampleRate


                CASE("&DEPTH_SUBBOX")

                    this%ImplantDepthDistType = p_DEPT_DIS_BOX

                    call EXTRACT_NUMB(STR,3,N,STRTMP)
                    if(N .LT. 3) then
                        write(*,*) "MCPSCUERROR: Too few parameters for the clusters depth distribution subbox type"
                        write(*,*) "You shoud special: &DEPTH_SUBBOX THE SUBOX SHAPE IS THAT: X =, Y =, Z ="
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if
                    DO I=1, 3
                        this%SUBBOXBOUNDARY(I,1) = Host_SimBoxes%BOXBOUNDARY(I,1) - DRSTR(STRTMP(I))*C_NM2CM/2
                        this%SUBBOXBOUNDARY(I,2) = Host_SimBoxes%BOXBOUNDARY(I,2) + DRSTR(STRTMP(I))*C_NM2CM/2
                    END DO

                CASE("&DEPTH_GAUSS")

                    this%ImplantDepthDistType = p_DEPT_DIS_GAS

                    call EXTRACT_NUMB(STR,2,N,STRTMP)
                    if(N .LT. 2) then
                        write(*,*) "MCPSCUERROR: Too few parameters for the clusters depth distribution gauss type"
                        write(*,*) "You shoud special: &DEPTH_GAUSS THE GAUSS DISTRIBUTION CENTRAL = , THE HALF WIDTH = "
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if
                    this%DepthINI = DRSTR(STRTMP(1))*C_NM2CM
                    this%DepthSDINI = DRSTR(STRTMP(2))*C_NM2CM
                CASE default
                    write(*,*) "MCPSCUERROR: Illegal Symbol: ", KEYWORD
                    pause
                    stop
            END SELECT

        END DO

        return

        100 write(*,*) "MCPSCUERROR : Load Implantation configuration file failed for clusters depth distribution!"
            write(*,*) "At line :",LINE
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine ReadImplantClusterDepthDist_Simple

    !********************************************************************
    subroutine CopyImplantSectionFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection),intent(out)::this
        type(ImplantSection),intent(in)::other
        !---Local Vars---
        integer::I
        !---Body---

        this%MemoryOccupyFactor = other%MemoryOccupyFactor

        this%ExpandFactor = other%ExpandFactor

        this%ImplantFlux = other%ImplantFlux

        this%ImplantConfigType =other%ImplantConfigType

        this%ImplantCfgFileName = other%ImplantCfgFileName

        this%ImplantCfgFileType = other%ImplantCfgFileType

        DO I = 1,p_ATOMS_GROUPS_NUMBER
            this%Elemets(I) = other%Elemets(I)
        END DO

        this%NAINI = other%NAINI

        this%NASDINI = other%NASDINI

        this%NACUT = other%NACUT

        this%CompositWeight = other%CompositWeight

        this%ImplantDepthDistType = other%ImplantDepthDistType

        this%SUBBOXBOUNDARY = other%SUBBOXBOUNDARY

        this%DepthINI = other%DepthINI

        this%DepthSDINI = other%DepthSDINI

        call DeAllocateArray_Host(this%LayerThick,"LayerThick")
        call AllocateArray_Host(this%LayerThick,size(other%LayerThick),"LayerThick")
        this%LayerThick = other%LayerThick

        !---The assignment(=) had been override
        call DeAllocateArray_Host(this%ClustersSample,"ClustersSample")
        call AllocateArray_Host(this%ClustersSample,size(other%ClustersSample,dim=1),size(other%ClustersSample,dim=2),"ClustersSample")
        this%ClustersSample = other%ClustersSample

        call DeAllocateArray_Host(this%ClustersSampleRate,"ClustersSampleRate")
        call AllocateArray_Host(this%ClustersSampleRate,size(other%ClustersSampleRate,dim=1),size(other%ClustersSampleRate,dim=2),"ClustersSampleRate")
        this%ClustersSampleRate = other%ClustersSampleRate

        !---The assignment(=) had been override
        this%dm_ImplantInfo_DevPart = other%dm_ImplantInfo_DevPart

        return
    end subroutine CopyImplantSectionFromOther

    !*********************************************************************
    subroutine Clean_ImplantSection(this)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        !---Body---
        this%MemoryOccupyFactor = 100
        this%ExpandFactor = 10
        this%ImplantFlux = 0

        this%ImplantConfigType = -1

        this%ImplantCfgFileName = ""

        this%ImplantCfgFileType = ""

        this%Elemets = ""

        this%NAINI = 0.D0

        this%NASDINI = 0.D0

        this%NACUT = 0.D0

        this%CompositWeight = 0.D0

        this%ImplantDepthDistType = -1

        this%SUBBOXBOUNDARY = 0.D0

        this%DepthINI = 0.D0

        this%DepthSDINI = 0.D0

        call DeAllocateArray_Host(this%LayerThick,"LayerThick")

        call DeAllocateArray_Host(this%ClustersSample,"ClustersSample")

        call DeAllocateArray_Host(this%ClustersSampleRate,"ClustersSampleRate")

        call this%dm_ImplantInfo_DevPart%Clean()
        return
    end subroutine Clean_ImplantSection

    !*********************************************************************
    subroutine CleanImplantSection(this)
        implicit none
        !---Dummy Vars---
        type(ImplantSection)::this
        !---Body---
        call this%Clean()

        return
    end subroutine CleanImplantSection


    !*********************************************************************
    subroutine ImplantClusters_FastStrategy(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,TheMigCoaleStatInfoWrap,Record,TSTEP)
        use RAND32_MODULE
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        type(MigCoaleStatInfoWrap)::TheMigCoaleStatInfoWrap
        type(MigCoalClusterRecord)::Record
        real(kind=KMCDF)::TSTEP
        !---Local Vars---
        integer::err
        integer::MultiBox
        integer::IBox
        integer::ImplantNumEachBox_Ceiling
        logical::NeedAddVirtualRange
        logical::NeedAddExpdRange
        integer::NewAllocateNCEachBox
        integer::NewTotalSize
        real(kind=KMCDF)::tempTSTEP
        integer::tempImplantNumEachBox
        integer::NSIZE
        integer::ImplantNumEachBox
        integer::TotalImplantNum
        !integer::NCUsed
        !integer::NCAct
        !logical::SweepOut
        !---Body---

        TotalImplantNum = 0

        NeedAddVirtualRange = .false.

        NeedAddExpdRange = .false.

        MultiBox = Host_SimuCtrlParam%MultiBox

        call this%AdjustTimeStep_Implant(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap,Record,TSTEP,ImplantNumEachBox_Ceiling)

        DO IBox = 1,MultiBox
            if((Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) + ImplantNumEachBox_Ceiling) .GT. Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2)) then
                NeedAddVirtualRange = .true.
                exit
            end if
        END DO

        if(NeedAddVirtualRange .eq. .true.) then

            call Dev_Boxes%GetBoxesBasicStatistic_AllStatu_GPU(Host_Boxes,Host_SimuCtrlParam)
            call Record%RecordNC_ForSweepOut(MultiBox,Host_Boxes%m_BoxesBasicStatistic)

            call Dev_Boxes%SweepUnActiveMemory_GPUToCPU(Host_Boxes,Host_SimuCtrlParam)

            call Dev_Boxes%GetBoxesBasicStatistic_AllStatu_GPU(Host_Boxes,Host_SimuCtrlParam)

            if(this%Cal_ImplantANDExpandSize(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Record,TSTEP,ImplantNumEachBox_Ceiling,NewAllocateNCEachBox) .eq. .false.) then
                write(*,*) "MCPSCUInfo: There are no enough memory to do the future implant job, so the time step is deduced."
            end if

            write(*,*) "NewAllocateNCEachBox",NewAllocateNCEachBox

            if(NewAllocateNCEachBox .GT. 0) then ! need to allocate new memory for new added cluster

                write(*,*) ".....Add virtual range...."

                call Dev_Boxes%ExpandClustersInfo_GPUToCPU_EqualNum(Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)

                if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
                    NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
                else
                    NewTotalSize = 0
                end if

                call Dev_MigCoaleGVars%dm_MigCoale_RandDev%ReSizeWalkRandNum(NewTotalSize)
                call Dev_MigCoaleGVars%dm_MigCoale_RandDev%ReSizeImplantRandNum(MultiBox*NewAllocateNCEachBox)

                call this%DoImplantTillVirtualBoundary_CPUTOGPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,NewAllocateNCEachBox)

                call GetBoxesMigCoaleStat_Virtual_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Virtual,Record)

                DO IBox = 1,MultiBox
                    write(*,*) "The virtual range for box ",IBox, " is ",Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2)
                END DO

            end if

            DO IBox = 1,MultiBox
                Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) = min(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) + ImplantNumEachBox_Ceiling*this%ExpandFactor,Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2))
            END DO
            Dev_Boxes%dm_SEExpdIndexBox = Host_Boxes%m_BoxesInfo%SEExpdIndexBox

            call GetBoxesMigCoaleStat_Expd_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd,Record)

            if(Host_SimuCtrlParam%TUpdateStatisFlag .eq. mp_UpdateStatisFlag_ByIntervalSteps) then
                call Record%SetLastUpdateStatisTime(Record%GetSimuSteps() + 1.D0)
            else if(Host_SimuCtrlParam%TUpdateStatisFlag .eq. mp_UpdateStatisFlag_ByIntervalRealTime) then
                call Record%SetLastUpdateStatisTime(Record%GetSimuTimes() + TSTEP)
            end if

            call Cal_Neighbor_List_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Record,IfDirectly=.true.,RMAX= &
                                      max(TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd%statistic_IntegralBox%RMAX(p_ACTIVEFREE_STATU), &
                                          TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd%statistic_IntegralBox%RMAX(p_ACTIVEINGB_STATU)))

        else

            DO IBox = 1,MultiBox
                if((Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) + ImplantNumEachBox_Ceiling) .GT. Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2)) then
                    NeedAddExpdRange = .true.
                    exit
                end if
            END DO

            if(NeedAddExpdRange .eq. .true.) then

                write(*,*) ".....Add expand range...."

                DO IBox = 1,MultiBox
                    Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) = min(Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2) + ImplantNumEachBox_Ceiling*this%ExpandFactor,Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2))

                    !if(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) .GT. 0) then
                    !    NCUsed = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,1) + 1
                    !else
                    !    NCUsed = 0
                    !end if

                    !NCAct = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(p_ACTIVEFREE_STATU) + Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(p_ACTIVEINGB_STATU)

                    !if(NCUsed .GT. this%SweepOutFactor*NCAct) then
                    !    SweepOut = .true.
                    !end if

                    write(*,*) "The expanded range for box ",IBox, " is ",Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2)
                END DO
                Dev_Boxes%dm_SEExpdIndexBox = Host_Boxes%m_BoxesInfo%SEExpdIndexBox

                call GetBoxesMigCoaleStat_Expd_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd,Record)

                !if(SweepOut) then
                    !write(*,*) "...Sweep out inactive cluster memory...."

                    !call Dev_Boxes%SweepUnActiveMemory_GPUToCPU(Host_Boxes,Host_SimuCtrlParam)

                !end if

                if(Host_SimuCtrlParam%TUpdateStatisFlag .eq. mp_UpdateStatisFlag_ByIntervalSteps) then
                    call Record%SetLastUpdateStatisTime(Record%GetSimuSteps() + 1.D0)
                else if(Host_SimuCtrlParam%TUpdateStatisFlag .eq. mp_UpdateStatisFlag_ByIntervalRealTime) then
                    call Record%SetLastUpdateStatisTime(Record%GetSimuTimes() + TSTEP)
                end if

                call Cal_Neighbor_List_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Record,IfDirectly=.true.,RMAX= &
                                            max(TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd%statistic_IntegralBox%RMAX(p_ACTIVEFREE_STATU), &
                                                TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Expd%statistic_IntegralBox%RMAX(p_ACTIVEINGB_STATU)))

            end if

        end if


        DO IBox = 1,MultiBox

            ImplantNumEachBox = floor(Host_Boxes%BOXSIZE(1)*Host_Boxes%BOXSIZE(2)*this%ImplantFlux*TSTEP)

            if(DRAND32() .LE. (Host_Boxes%BOXSIZE(1)*Host_Boxes%BOXSIZE(2)*this%ImplantFlux*TSTEP - ImplantNumEachBox)) then
                ImplantNumEachBox = ImplantNumEachBox + 1
            end if

            Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) =  Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) + ImplantNumEachBox

            if(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) .GT. Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2)) then
                write(*,*) "MCPSCUERROR: The expand size for box ,",IBox," is not enough!"
                write(*,*) "Implant number is: ",ImplantNumEachBox
                write(*,*) "End of the used clusters index is, ",Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2)
                write(*,*) "End of the expd clusters index is, ",Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2)
                pause
                stop
            end if

            if(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) .GT. Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2)) then
                write(*,*) "MCPSCUERROR: The virtual size for box ,",IBox," is not enough!"
                write(*,*) "Implant number is: ",ImplantNumEachBox
                write(*,*) "End of the used clusters index is, ",Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2)
                write(*,*) "End of the virtual clusters index is, ",Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2)
                pause
                stop
            end if

            TotalImplantNum = TotalImplantNum + ImplantNumEachBox
        END DO

        Dev_Boxes%dm_SEUsedIndexBox = Host_Boxes%m_BoxesInfo%SEUsedIndexBox

        call Record%AddImplantedEntitiesNum(TotalImplantNum)

        return
    end subroutine ImplantClusters_FastStrategy

    !*********************************************
    subroutine AdjustTimeStep_Implant(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap,Record,TSTEP,ImplantNumEachBox_Ceiling)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoaleStatInfoWrap)::TheMigCoaleStatInfoWrap
        type(MigCoalClusterRecord)::Record
        real(kind=KMCDF)::TSTEP
        integer::ImplantNumEachBox_Ceiling
        !---Local Vars---
        integer::ImplantNumEachBox
        real(kind=KMCDF)::VerifyTime
        !---Body---

        DO While(.true.)

            ImplantNumEachBox_Ceiling = ceiling(Host_Boxes%BOXSIZE(1)*Host_Boxes%BOXSIZE(2)*this%ImplantFlux*TSTEP)

            VerifyTime = Cal_VerifyTime_Implant(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TheMigCoaleStatInfoWrap%m_MigCoaleStatisticInfo_Virtual,Record,ImplantNumEachBox_Ceiling)

            if(VerifyTime .LT. TSTEP) then
                TSTEP = TSTEP*0.95
                cycle
            else
                exit
            end if
        END DO

        return
    end subroutine AdjustTimeStep_Implant

    !*********************************************
    function Cal_ImplantANDExpandSize(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Record,TSTEP,ImplantNumEachBox_Ceiling,NewAllocateNCEachBox) result(TheStatu)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoalClusterRecord)::Record
        real(kind=KMCDF), intent(inout)::TSTEP
        integer, intent(inout)::ImplantNumEachBox_Ceiling
        integer, intent(inout)::NewAllocateNCEachBox
        logical, intent(inout)::TheStatu
        !---Local Vars---
        integer::err
        real(kind=KMCDF)::newSimulatedTime
        integer::MultiBox
        integer::IBox
        integer(kind=cuda_count_kind)::FreeMemSize,TotalMemSize
        integer::FreeMemSize2NCEachBox
        integer::TotalMemSize2NCEachBox
        integer::AllocatedFreeNCEachBox
        integer::newFreeNCEachBox
        integer::PreAllocatedNCEachBox
        integer::PreFreeNCEachBox
        integer::RestoreImplantNumEachBox
        real::ImplantPersistTime
        integer,dimension(:),allocatable::NCFree
        !---Body---

        RestoreImplantNumEachBox = ImplantNumEachBox_Ceiling

        FreeMemSize2NCEachBox = 0
        TotalMemSize2NCEachBox = 0
        NewAllocateNCEachBox = 0

        MultiBox = Host_SimuCtrlParam%MultiBox

        ImplantPersistTime = Record%GetSimuTimes() - Record%GetStartImplantTime()

        call AllocateArray_Host(NCFree,MultiBox,"NCFree")

        DO IBox=1,MultiBox
            PreAllocatedNCEachBox = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,1) + 1
            if(PreAllocatedNCEachBox .LE. 0) then
                PreAllocatedNCEachBox = 0
            end if

            PreFreeNCEachBox = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,1) + 1
            if(PreFreeNCEachBox .LE. 0) then
                PreFreeNCEachBox = 0
            end if

            NCFree(IBox) = PreAllocatedNCEachBox - PreFreeNCEachBox

            if(NCFree(IBox) .LE. 0) then
                NCFree(IBox) = 0
            end if
        END DO
        AllocatedFreeNCEachBox = minval(NCFree)

        err = cudaMemGetInfo(FreeMemSize,TotalMemSize)

        !The 8*3 is for the random walking random number
        if(Host_SimuCtrlParam%FreeDiffusion .eq. .false.) then
            FreeMemSize2NCEachBox  = (FreeMemSize/(Dev_Boxes%dm_ClusterInfo_GPU%GetMemConsumOneClusterInfo(Host_SimuCtrlParam%MAXNEIGHBORNUM) + 8*3))/MultiBox
            TotalMemSize2NCEachBox = (TotalMemSize/(Dev_Boxes%dm_ClusterInfo_GPU%GetMemConsumOneClusterInfo(Host_SimuCtrlParam%MAXNEIGHBORNUM)+ 8*3))/MultiBox
        else
            FreeMemSize2NCEachBox  = (FreeMemSize/(Dev_Boxes%dm_ClusterInfo_GPU%GetMemConsumOneClusterInfo(0) + 8*3))/MultiBox
            TotalMemSize2NCEachBox = (TotalMemSize/(Dev_Boxes%dm_ClusterInfo_GPU%GetMemConsumOneClusterInfo(0)+ 8*3))/MultiBox
        end if

        DO while(.true.)

            NewAllocateNCEachBox = 0

            newFreeNCEachBox = AllocatedFreeNCEachBox

            if(newFreeNCEachBox .GE. ImplantNumEachBox_Ceiling*this%ExpandFactor) then
                exit
            end if

            NewAllocateNCEachBox = min(FreeMemSize2NCEachBox,TotalMemSize2NCEachBox/this%MemoryOccupyFactor)

            newFreeNCEachBox = AllocatedFreeNCEachBox + NewAllocateNCEachBox

            if(newFreeNCEachBox .LT. ImplantNumEachBox_Ceiling*this%ExpandFactor .or. (FreeMemSize2NCEachBox - NewAllocateNCEachBox) .LE. 0) then
                TSTEP = TSTEP/2.D0
                ImplantNumEachBox_Ceiling = ceiling(Host_Boxes%BOXSIZE(1)*Host_Boxes%BOXSIZE(2)*this%ImplantFlux*TSTEP)
            else
                exit
            end if

            if(ImplantNumEachBox_Ceiling .LE. 0) then
                exit
            end if
        END DO

        if(RestoreImplantNumEachBox .GT. 0 .AND. ImplantNumEachBox_Ceiling .LE. 0) then
            TheStatu = .false.
        else
            TheStatu = .true.
        end if

        call DeAllocateArray_Host(NCFree,"NCFree")

        return
    end function Cal_ImplantANDExpandSize

    !*************************************************************
    subroutine DoImplantTillVirtualBoundary_CPU(this,Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer,intent(in)::NewAllocateNCEachBox
        !---Body---
        select case(this%ImplantConfigType)
            case(p_ImplantConfig_Simple)
                call this%FillVirtualBoundary_CPU_Simple(Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
            case(p_ImplantConfig_SpecialDistFromFile)
                call this%FillVirtualBoundary_CPU_FromFile(Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
            case(p_ImplantConfig_SpecialDistFromExteFunc)
                call this%FillVirtualBoundary_CPU_FromExteFunc(Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
            case default
                write(*,*) "MCPSCUERROR: Unknown strategy for the implantation configuration:",this%ImplantConfigType
                pause
                stop
        end select
        return
    end subroutine DoImplantTillVirtualBoundary_CPU

    !*************************************************************
    subroutine FillVirtualBoundary_CPU_Simple(this,SimBoxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer,intent(in)::NewAllocateNCEachBox
        !---Body--
        select case(this%ImplantDepthDistType)
            case(p_DEPT_DIS_Layer)
                call this%FillVirtualBoundary_CPU_Depth_LAY(SimBoxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
            case(p_DEPT_DIS_BOX)
                call this%FillVirtualBoundary_CPU_Depth_SubBox(SimBoxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
            case(p_DEPT_DIS_GAS)
                call this%FillVirtualBoundary_CPU_Depth_Gauss(SimBoxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
            case default
                write(*,*) "MCPSCUERROR : Unknown way to Unknown strategy for the simple implantation configuration: ",this%ImplantDepthDistType
                pause
                stop
        end select

        return
    end subroutine FillVirtualBoundary_CPU_Simple

    !*************************************************************
    subroutine FillVirtualBoundary_CPU_FromFile(this,Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
        use RAND32_MODULE
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer,intent(in)::NewAllocateNCEachBox
        !---Local Vars---
        integer::MultiBox
        integer::IBox
        integer::IC
        integer::ICFROM
        integer::ICTO
        type(DiffusorValue)::TheDiffusorValue
        integer::NC
        integer::NCUsed
        real(kind=KMCDF)::POS(3)
        integer::MaxGroups
        real(kind=KMCDF)::tempRand
        real(kind=KMCDF)::GroupRateTemp
        integer::ILayer
        integer::IGroup
        logical::exitFlag
        integer::LayerNum
        !---Body---
        MultiBox = Host_SimuCtrlParam%MultiBox

        LayerNum = size(this%LayerThick)

        DO IBox = 1,MultiBox

            if(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) .LE. 0) then
                NCUsed = 0
            else
                NCUsed = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,1) + 1
            end if

            if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) .LE. 0) then
                NC = 0
            else
                NC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,1) + 1
            end if

            if((NC - NCUsed) .LT. NewAllocateNCEachBox) then
                write(*,*) "MCPSCUERROR: The allocated  memory space are to implant the clusters"
                write(*,*) "For box :",IBox
                write(*,*) "The free of allocated allocated  memory space is: ",NC - NCUsed
                write(*,*) "The waiting to be implanted clusters number is:",NewAllocateNCEachBox
                pause
                stop
            end if

            ICFROM = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1
            ICTO = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2)

            if(NewAllocateNCEachBox .LE. 0) then
                cycle
            end if

            if(ICFROM .LT. Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2)) then
                write(*,*) "MCPSCUERROR: The allocated  memory space are too small to implant the clusters"
                write(*,*) "For box :",IBox
                write(*,*) "It would occupy other free clusters for id: ",ICFROM
                pause
                stop
            end if

            if(size(this%ClustersSample) .LE. 0) then
                write(*,*) "MCPSCUERROR: The number of sample clusters is less than 1: "
                write(*,*) "There would be no clusters to be implanted from the sample distribution."
                pause
                stop
            end if

            MaxGroups = size(this%ClustersSampleRate,dim=2)

            DO IC = ICFROM,ICTO

                tempRand = DRAND32()

                GroupRateTemp = 0.D0

                exitFlag = .false.
                DO ILayer = 1,LayerNum

                    if(exitFlag .eq. .true.) then
                        exit
                    end if

                     DO IGroup = 1,MaxGroups
                        GroupRateTemp = GroupRateTemp + this%ClustersSampleRate(ILayer,IGroup)
                        if(GroupRateTemp .GE. tempRand) then
                            !---The assignment(=) had been override
                            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms = this%ClustersSample(ILayer,IGroup)%m_Atoms

                            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Layer = ILayer

                            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu = this%ClustersSample(ILayer,IGroup)%m_Statu

                            POS(1) = DRAND32()*Host_Boxes%BOXSIZE(1) + Host_Boxes%BOXBOUNDARY(1,1)
                            POS(2) = DRAND32()*Host_Boxes%BOXSIZE(2) + Host_Boxes%BOXBOUNDARY(2,1)
                            POS(3) = DRAND32()*this%LayerThick(ILayer) + sum(this%LayerThick(1:ILayer-1)) + Host_Boxes%BOXBOUNDARY(3,1)
                            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_POS = POS

                            if(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD .GT. 0.D0) then
                                write(*,*) "MCPSCUERROR: the implant position had been occupied in memory",IC
                                pause
                            end if

                            TheDiffusorValue = Host_Boxes%m_DiffusorTypesMap%Get(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC))

                            if(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu .eq. p_ACTIVEFREE_STATU) then
                                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) = Host_Boxes%m_GrainBoundary%GrainBelongsTo(POS,Host_Boxes%HBOXSIZE,Host_Boxes%BOXSIZE,Host_SimuCtrlParam)

                                select case(TheDiffusorValue%ECRValueType_Free)
                                    case(p_ECR_ByValue)
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = TheDiffusorValue%ECR_Free
                                    case(p_ECR_ByBCluster)
                                        ! Convert the number of atoms to radius
                                        ! Ref. Modelling Simul. Mater. Sci. Eng.16(2008)055003
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = DSQRT(sum(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA)/m_RNFACTOR)
                                end select

                                select case(TheDiffusorValue%DiffusorValueType_Free)
                                    case(p_DiffuseCoefficient_ByValue)
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                                    case(p_DiffuseCoefficient_ByArrhenius)
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
                                    case(p_DiffuseCoefficient_ByBCluster)
                                        ! Here we adopt a model that D=D0*(1/R)**Gama
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = m_FREESURDIFPRE*(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD**(-p_GAMMA))
                                end select

                            else if(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu .eq. p_ACTIVEINGB_STATU) then

                                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) = this%ClustersSample(ILayer,IGroup)%m_GrainID(1)

                                if(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) .GT. Host_Boxes%m_GrainBoundary%GrainNum) then
                                    write(*,*) "MCPSCUERROR: The grain number is greater than the seeds number in system."
                                    write(*,*) Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1)
                                    pause
                                    stop
                                end if

                                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(2) = this%ClustersSample(ILayer,IGroup)%m_GrainID(2)

                                if(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(2) .GT. Host_Boxes%m_GrainBoundary%GrainNum) then
                                    write(*,*) "MCPSCUERROR: The grain number is greater than the seeds number in system."
                                    write(*,*) Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(2)
                                    pause
                                    stop
                                end if

                                select case(TheDiffusorValue%ECRValueType_InGB)
                                    case(p_ECR_ByValue)
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = TheDiffusorValue%ECR_InGB
                                    case(p_ECR_ByBCluster)
                                        ! Convert the number of atoms to radius
                                        ! Ref. Modelling Simul. Mater. Sci. Eng.16(2008)055003
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = DSQRT(sum(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA)/m_RNFACTOR)
                                end select

                                select case(TheDiffusorValue%DiffusorValueType_InGB)
                                    case(p_DiffuseCoefficient_ByValue)
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_InGB_Value
                                    case(p_DiffuseCoefficient_ByArrhenius)
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%PreFactor_InGB*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_InGB/Host_SimuCtrlParam%TKB)
                                    case(p_DiffuseCoefficient_ByBCluster)
                                        ! Here we adopt a model that D=D0*(1/R)**Gama
                                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = m_GBSURDIFPRE*(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD**(-p_GAMMA))
                                end select
                            end if

                            exitFlag = .true.

                            exit

                        end if

                    END DO
                END DO

            END DO

        END DO

        return
    end subroutine FillVirtualBoundary_CPU_FromFile

    !*************************************************************
    subroutine FillVirtualBoundary_CPU_FromExteFunc(this,Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer,intent(in)::NewAllocateNCEachBox
        !---Body---

    end subroutine FillVirtualBoundary_CPU_FromExteFunc


    !*************************************************************
    subroutine DoImplantTillVirtualBoundary_CPUTOGPU(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        integer,intent(in)::NewAllocateNCEachBox
        !---Local Vars---
        integer::MultiBox
        integer::NSIZE
        !---Body---
        MultiBox = Host_SimuCtrlParam%MultiBox

        call this%DoImplantTillVirtualBoundary_CPU(Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)

        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
            NSIZE = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
        else
            NSIZE = 0
        end if

        if(NSIZE .ne. size(Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters)) then

            call Dev_Boxes%dm_ClusterInfo_GPU%ReleaseClustersInfo_GPU()

            if(Host_SimuCtrlParam%FreeDiffusion .eq. .false.) then
                call Dev_Boxes%dm_ClusterInfo_GPU%AllocateClustersInfo_GPU(NSIZE,Host_SimuCtrlParam%MAXNEIGHBORNUM)
            else
                call Dev_Boxes%dm_ClusterInfo_GPU%AllocateClustersInfo_GPU(NSIZE,0)
            end if
        end if

        call Dev_Boxes%dm_ClusterInfo_GPU%CopyInFromHost(Host_Boxes%m_ClustersInfo_CPU,NSIZE,IfCpyNL=.false.)

        return
    end subroutine DoImplantTillVirtualBoundary_CPUTOGPU

    !**************************************************************
    subroutine FillVirtualBoundary_CPU_Depth_LAY(this,Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
      !*** Purpose: To initialize the system (clusters distributed as the form of layer)
      ! Host_Boxes: the boxes information in host
      use RAND32_MODULE
      implicit none
      !---Dummy Vars---
      CLASS(ImplantSection)::this
      type(SimulationBoxes)::Host_Boxes
      type(SimulationCtrlParam)::Host_SimuCtrlParam
      integer,intent(in)::NewAllocateNCEachBox
      !-----local variables---
      integer::MultiBox
      real(kind=KMCDF)::POS(3)
      integer::IBox
      integer::NC
      integer::NCUsed
      integer::IC
      integer::ICFROM
      integer::ICTO
      integer::MaxGroups
      real(kind=KMCDF)::tempRand
      real(kind=KMCDF)::GroupRateTemp
      integer::ILayer
      integer::IGroup
      integer::IElement
      integer::NAtoms
      type(DiffusorValue)::TheDiffusorValue
      logical::exitFlag
      integer::LayerNum
      !---Body---
      MultiBox = Host_SimuCtrlParam%MultiBox

      LayerNum = size(this%LayerThick)

      DO IBox = 1,MultiBox

        if(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) .LE. 0) then
            NCUsed = 0
        else
            NCUsed = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,1) + 1
        end if

        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) .LE. 0) then
            NC = 0
        else
            NC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,1) + 1
        end if

        if((NC - NCUsed) .LT. NewAllocateNCEachBox) then
            write(*,*) "MCPSCUERROR: The allocated  memory space are to implant the clusters"
            write(*,*) "For box :",IBox
            write(*,*) "The free of allocated allocated  memory space is: ",NC - NCUsed
            write(*,*) "The waiting to be implanted clusters number is:",NewAllocateNCEachBox
            pause
            stop
        end if

        ICFROM = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1
        ICTO = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2)

        if(ICFROM .LT. Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2)) then
            write(*,*) "MCPSCUERROR: The allocated  memory space are too small to implant the clusters"
            write(*,*) "For box :",IBox
            write(*,*) "It would occupy other free clusters for id: ",IC
            pause
            stop
        end if

        MaxGroups = size(this%ClustersSampleRate,dim=2)

        DO IC=ICFROM,ICTO

            tempRand = DRAND32()

            GroupRateTemp = 0.D0

            exitFlag = .false.
            DO ILayer = 1,LayerNum
                if(exitFlag .eq. .true.) then
                    exit
                end if

                DO IGroup = 1,MaxGroups
                    GroupRateTemp = GroupRateTemp + this%ClustersSampleRate(ILayer,IGroup)
                    if(GroupRateTemp .GE. tempRand) then

                        !*** Initialize the size of the clusters
                        NAtoms = RGAUSS0_WithCut(this%NAINI, this%NASDINI,this%NACUT(1),this%NACUT(2))

                        DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(IElement)%m_NA = floor(NAtoms*this%CompositWeight(IElement)+0.5D0)
                            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(IElement)%m_ID = IElement
                        END DO

                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu = p_ACTIVEFREE_STATU

                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Layer = ILayer

                        POS(1) = DRAND32()*Host_Boxes%BOXSIZE(1) + Host_Boxes%BOXBOUNDARY(1,1)
                        POS(2) = DRAND32()*Host_Boxes%BOXSIZE(2) + Host_Boxes%BOXBOUNDARY(2,1)
                        POS(3) = DRAND32()*this%LayerThick(ILayer) + sum(this%LayerThick(1:ILayer-1)) + Host_Boxes%BOXBOUNDARY(3,1)
                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_POS = POS

                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) = Host_Boxes%m_GrainBoundary%GrainBelongsTo(POS,Host_Boxes%HBOXSIZE,Host_Boxes%BOXSIZE,Host_SimuCtrlParam)

                        TheDiffusorValue = Host_Boxes%m_DiffusorTypesMap%Get(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC))

                        !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
                        !---you should init the distribution by external file---
                        select case(TheDiffusorValue%ECRValueType_Free)
                            case(p_ECR_ByValue)
                                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = TheDiffusorValue%ECR_Free
                            case(p_ECR_ByBCluster)
                                ! Convert the number of atoms to radius
                                ! Ref. Modelling Simul. Mater. Sci. Eng.16(2008)055003
                                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = DSQRT(sum(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA)/m_RNFACTOR)
                        end select

                        select case(TheDiffusorValue%DiffusorValueType_Free)
                            case(p_DiffuseCoefficient_ByValue)
                                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                            case(p_DiffuseCoefficient_ByArrhenius)
                                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
                            case(p_DiffuseCoefficient_ByBCluster)
                                    ! Here we adopt a model that D=D0*(1/R)**Gama
                                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = m_FREESURDIFPRE*(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD**(-p_GAMMA))
                        end select

                        exitFlag = .true.
                        exit

                    end if
                END DO
            END DO

        END DO

      END DO

      return
    end subroutine FillVirtualBoundary_CPU_Depth_LAY

    !**************************************************************
    subroutine FillVirtualBoundary_CPU_Depth_SubBox(this,Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
      !*** Purpose: To initialize the system (clusters distributed as the form of layer)
      ! Host_Boxes: the boxes information in host
      use RAND32_MODULE
      implicit none
      !---Dummy Vars---
      CLASS(ImplantSection)::this
      type(SimulationBoxes)::Host_Boxes
      type(SimulationCtrlParam)::Host_SimuCtrlParam
      integer,intent(in)::NewAllocateNCEachBox
      !-----local variables---
      integer::MultiBox
      real(kind=KMCDF)::POS(3)
      real(kind=KMCDF)::SUBBOXSIZE(3)
      integer::IBox, II, IC
      integer::I
      integer::NC
      integer::NCUsed
      integer::NAtoms
      integer::IElement
      type(DiffusorValue)::TheDiffusorValue
      !---Body---
      MultiBox = Host_SimuCtrlParam%MultiBox

      DO I = 1,3
        SUBBOXSIZE(I) = this%SUBBOXBOUNDARY(I,2) - this%SUBBOXBOUNDARY(I,1)
      END DO

      DO IBox = 1,MultiBox

        if(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) .LE. 0) then
            NCUsed = 0
        else
            NCUsed = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,1) + 1
        end if

        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) .LE. 0) then
            NC = 0
        else
            NC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,1) + 1
        end if

        if((NC - NCUsed) .LT. NewAllocateNCEachBox) then
            write(*,*) "MCPSCUERROR: The allocated  memory space are to implant the clusters"
            write(*,*) "For box :",IBox
            write(*,*) "The free of allocated allocated  memory space is: ",NC - NCUsed
            write(*,*) "The waiting to be implanted clusters number is:",NewAllocateNCEachBox
            pause
            stop
        end if

        IC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox

        if(IC .LT. Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2)) then
            write(*,*) "MCPSCUERROR: The allocated  memory space are too small to implant the clusters"
            write(*,*) "For box :",IBox
            write(*,*) "It would occupy other free clusters for id: ",IC
            pause
            stop
        end if

        DO II = 1, NewAllocateNCEachBox

            IC = IC + 1
            !Initialize the position of clusters
            DO I = 1,3
                POS(I) = DRAND32()*SUBBOXSIZE(I) + this%SUBBOXBOUNDARY(I,1)
            END DO

            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_POS = POS
            !Give the cluster an type(layer) ID for the convenience of visualization
            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Layer = 1

            NAtoms = RGAUSS0_WithCut(this%NAINI, this%NASDINI,this%NACUT(1),this%NACUT(2))

            DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(IElement)%m_NA = floor(NAtoms*this%CompositWeight(IElement)+0.5D0)
                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(IElement)%m_ID = IElement
            END DO

            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) = Host_Boxes%m_GrainBoundary%GrainBelongsTo(POS,Host_Boxes%HBOXSIZE,Host_Boxes%BOXSIZE,Host_SimuCtrlParam)

            Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu = p_ACTIVEFREE_STATU

            TheDiffusorValue = Host_Boxes%m_DiffusorTypesMap%Get(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC))

            !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
            !---you should init the distribution by external file---
            select case(TheDiffusorValue%ECRValueType_Free)
                case(p_ECR_ByValue)
                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = TheDiffusorValue%ECR_Free
                case(p_ECR_ByBCluster)
                    ! Convert the number of atoms to radius
                    ! Ref. Modelling Simul. Mater. Sci. Eng.16(2008)055003
                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = DSQRT(sum(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA)/m_RNFACTOR)
             end select

            select case(TheDiffusorValue%DiffusorValueType_Free)
                case(p_DiffuseCoefficient_ByValue)
                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                case(p_DiffuseCoefficient_ByArrhenius)
                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
                case(p_DiffuseCoefficient_ByBCluster)
                    ! Here we adopt a model that D=D0*(1/R)**Gama
                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = m_FREESURDIFPRE*(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD**(-p_GAMMA))
            end select

        END DO

      END DO

      return
    end subroutine FillVirtualBoundary_CPU_Depth_SubBox

    !**************************************************************
    subroutine FillVirtualBoundary_CPU_Depth_Gauss(this,Host_Boxes,Host_SimuCtrlParam,NewAllocateNCEachBox)
        !*** Purpose: To initialize the system (clusters distributed as the form of gauss in depth)
        ! Host_Boxes: the boxes information in host
        use RAND32_MODULE
        implicit none
        !-------Dummy Vars------
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer,intent(in)::NewAllocateNCEachBox
        !---local variables---
        integer::MultiBox
        real(kind=KMCDF)::POS(3)
        real(kind=KMCDF)::SEP
        real(kind=KMCDF)::BOXBOUNDARY(3,2)
        real(kind=KMCDF)::BOXSIZE(3)
        integer::IBox,II,IC
        integer::NAtoms
        integer::NC
        integer::NCUsed
        integer::IElement
        type(DiffusorValue)::TheDiffusorValue
        !---Body---
        MultiBox = Host_SimuCtrlParam%MultiBox
        BOXBOUNDARY = Host_Boxes%BOXBOUNDARY
        BOXSIZE = Host_Boxes%BOXSIZE

        DO IBox = 1,MultiBox

            if(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) .LE. 0) then
                NCUsed = 0
            else
                NCUsed = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,1) + 1
            end if

            if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) .LE. 0) then
                NC = 0
            else
                NC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,1) + 1
            end if

            if((NC - NCUsed) .LT. NewAllocateNCEachBox) then
                write(*,*) "MCPSCUERROR: The allocated  memory space are to implant the clusters"
                write(*,*) "For box :",IBox
                write(*,*) "The free of allocated allocated  memory space is: ",NC - NCUsed
                write(*,*) "The waiting to be implanted clusters number is:",NewAllocateNCEachBox
                pause
                stop
            end if

            IC = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox

            if(IC .LT. Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2)) then
                write(*,*) "MCPSCUERROR: The allocated  memory space are too small to implant the clusters"
                write(*,*) "For box :",IBox
                write(*,*) "It would occupy other free clusters for id: ",IC
                pause
                stop
            end if

            DO II = 1, NewAllocateNCEachBox

                IC = IC + 1
                !Initialize the position of clusters
                POS(1) = DRAND32()*BOXSIZE(1) + BOXBOUNDARY(1,1)
                POS(2) = DRAND32()*BOXSIZE(2) + BOXBOUNDARY(2,1)
                POS(3) = RGAUSS0_WithCut(this%DepthINI, this%DepthSDINI,BOXBOUNDARY(3,1),BOXBOUNDARY(3,2))

                if(POS(3) .LT. BOXBOUNDARY(3,1)) then
                    POS(3) = BOXBOUNDARY(3,1)
                end if

                if(POS(3) .GT. BOXBOUNDARY(3,2)) then
                    POS(3) = BOXBOUNDARY(3,2)
                end if

                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_POS = POS

                !Give the cluster an type(layper) ID for the convenience of visualization
                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Layer = 1

                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_GrainID(1) = Host_Boxes%m_GrainBoundary%GrainBelongsTo(POS,Host_Boxes%HBOXSIZE,Host_Boxes%BOXSIZE,Host_SimuCtrlParam)

                Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu = p_ACTIVEFREE_STATU

                !*** Initialize the size of the clusters
                NAtoms = RGAUSS0_WithCut(this%NAINI, this%NASDINI,this%NACUT(1),this%NACUT(2))

                DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(IElement)%m_NA = floor(NAtoms*this%CompositWeight(IElement)+0.5D0)
                    Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(IElement)%m_ID = IElement
                END DO

                TheDiffusorValue = Host_Boxes%m_DiffusorTypesMap%Get(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC))

                !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
                !---you should init the distribution by external file---
                select case(TheDiffusorValue%ECRValueType_Free)
                    case(p_ECR_ByValue)
                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = TheDiffusorValue%ECR_Free
                    case(p_ECR_ByBCluster)
                        ! Convert the number of atoms to radius
                        ! Ref. Modelling Simul. Mater. Sci. Eng.16(2008)055003
                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD = DSQRT(sum(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Atoms(:)%m_NA)/m_RNFACTOR)
                end select

                select case(TheDiffusorValue%DiffusorValueType_Free)
                    case(p_DiffuseCoefficient_ByValue)
                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                    case(p_DiffuseCoefficient_ByArrhenius)
                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/Host_SimuCtrlParam%TKB)
                    case(p_DiffuseCoefficient_ByBCluster)
                        ! Here we adopt a model that D=D0*(1/R)**Gama
                        Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff = m_FREESURDIFPRE*(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD**(-p_GAMMA))
                end select

            END DO

        END DO

        return
    end subroutine FillVirtualBoundary_CPU_Depth_Gauss


    !*************************************************************
    subroutine DoImplantTillVirtualBoundary_GPUTOCPU(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        integer,intent(in)::NewAllocateNCEachBox
        !---Local Vars---
        integer::MultiBox
        integer::NSIZE
        !---Body---

        MultiBox = Host_SimuCtrlParam%MultiBox

        call this%DoImplantTillVirtualBoundary_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)

        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
            NSIZE = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,1) + 1
        else
            NSIZE = 0
        end if

        call Host_Boxes%Clean()

        call Dev_Boxes%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,NSIZE,IfCpyNL=.false.)

        return
    end subroutine DoImplantTillVirtualBoundary_GPUTOCPU


    !*************************************************************
    subroutine DoImplantTillVirtualBoundary_GPU(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        integer,intent(in)::NewAllocateNCEachBox
        !---Body---

        call this%InitImplantInfo_DevPart()

        select case(this%ImplantConfigType)
            case(p_ImplantConfig_Simple)
                call this%FillVirtualBoundary_GPU_Simple(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
            case(p_ImplantConfig_SpecialDistFromFile)
                call this%FillVirtualBoundary_GPU_FromFile(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
            case(p_ImplantConfig_SpecialDistFromExteFunc)
                call this%FillVirtualBoundary_GPU_FromExteFunc(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
            case default
                write(*,*) "MCPSCUERROR: Unknown strategy for the implantation configuration:",this%ImplantConfigType
                pause
                stop
        end select

        return
    end subroutine DoImplantTillVirtualBoundary_GPU


    !*************************************************************
    subroutine FillVirtualBoundary_GPU_Simple(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        integer,intent(in)::NewAllocateNCEachBox
        !---Body---
        select case(this%ImplantDepthDistType)
            case(p_DEPT_DIS_Layer)
                call this%FillVirtualBoundary_GPU_Depth_LAY(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
            case(p_DEPT_DIS_BOX)
                call this%FillVirtualBoundary_GPU_Depth_SubBox(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
            case(p_DEPT_DIS_GAS)
                call this%FillVirtualBoundary_GPU_Depth_Gauss(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
            case default
                write(*,*) "MCPSCUERROR : Unknown way to Unknown strategy for the simple implantation configuration: ",this%ImplantDepthDistType
                pause
                stop
        end select

        return
    end subroutine FillVirtualBoundary_GPU_Simple

    !*************************************************************
    subroutine FillVirtualBoundary_GPU_Depth_LAY(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        integer,intent(in)::NewAllocateNCEachBox
        !---Local Vars---
        integer::MultiBox
        integer::TotalAllocateNC
        type(dim3)::blocks
        type(dim3)::threads
        integer::NB
        integer::NBX,NBY
        Integer::BX
        integer::BY
        integer::err
        !---Body---

        ASSOCIATE(ImplantRand=>Dev_MigCoaleGVars%dm_MigCoale_RandDev)

            if(NewAllocateNCEachBox .GT. 0) then

                MultiBox = Host_SimuCtrlParam%MultiBox

                TotalAllocateNC = MultiBox*NewAllocateNCEachBox

                NB = (TotalAllocateNC - 1)/p_BLOCKSIZE + 1
                NBX  = min(NB,p_BLOCKDIMX)
                NBY = (NB - 1)/NBX + 1

                !*** to determine the block size
                BX = p_BLOCKSIZE
                BY = 1
                !*** to determine the dimension of blocks

                blocks  = dim3(NBX, NBY, 1)
                threads = dim3(BX,  BY,  1)

                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Layer,ImplantRand%dm_SpaceDist_Implant(1:TotalAllocateNC),TotalAllocateNC)
                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_X,ImplantRand%dm_SpaceDist_Implant(TotalAllocateNC+1:2*TotalAllocateNC),TotalAllocateNC)
                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Y,ImplantRand%dm_SpaceDist_Implant(2*TotalAllocateNC+1:3*TotalAllocateNC),TotalAllocateNC)
                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Z,ImplantRand%dm_SpaceDist_Implant(3*TotalAllocateNC+1:4*TotalAllocateNC),TotalAllocateNC)

                err = curandGenerateNormal(ImplantRand%m_ranGen_ClustersSizeDist,ImplantRand%dm_SizeDist_Implant,TotalAllocateNC,this%NAINI,this%NASDINI)

                call Kernel_ImplantClusters_Depth_Layer<<<blocks,threads>>>(TotalAllocateNC,                                          &
                                                                            NewAllocateNCEachBox,                                     &
                                                                            Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,                 &
                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_TypesEntities,          &
                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_SingleAtomsDivideArrays,&
                                                                            Host_Boxes%m_GrainBoundary%GrainNum,                      &
                                                                            Dev_Boxes%dm_GrainBoundary%dm_GrainSeeds,                 &
                                                                            ImplantRand%dm_SpaceDist_Implant,                         &
                                                                            ImplantRand%dm_SizeDist_Implant,                          &
                                                                            this%NACUT(1),                                            &
                                                                            this%NACUT(2),                                            &
                                                                            Dev_Boxes%dm_SEVirtualIndexBox,                           &
                                                                            this%dm_ImplantInfo_DevPart%Dev_LayerThick,               &
                                                                            this%dm_ImplantInfo_DevPart%Dev_ClustersSampleRate,       &
                                                                            this%dm_ImplantInfo_DevPart%Dev_CompositWeight)
            end if

        END ASSOCIATE

        return
    end subroutine FillVirtualBoundary_GPU_Depth_LAY

    !**********************************************
    attributes(global) subroutine Kernel_ImplantClusters_Depth_Layer(TotalAllocateNC,            &
                                                                    NewAllocateNCEachBox,        &
                                                                    Dev_Clusters,                &
                                                                    Dev_TypesEntities,           &
                                                                    Dev_SingleAtomsDivideArrays, &
                                                                    Nseeds,                      &
                                                                    Dev_GrainSeeds,              &
                                                                    Dev_RandArray_SpaceDist,     &
                                                                    Dev_RandArray_SizeDist,      &
                                                                    LNACUT,                      &
                                                                    RNACUT,                      &
                                                                    Dev_SEVirtualIndexBox,       &
                                                                    Dev_LayerThick,              &
                                                                    Dev_ClustersSampleRate,      &
                                                                    Dev_CompositWeight)
        implicit none
        !---Dummy Vars---
        integer, value::TotalAllocateNC
        integer, value::NewAllocateNCEachBox
        type(ACluster), device::Dev_Clusters(:)
        type(DiffusorTypeEntity),device::Dev_TypesEntities(:)
        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! If the two dimension array would be delivered to attributes(device), the first dimension must be known
        integer,value::Nseeds
        type(GrainSeed),device::Dev_GrainSeeds(:)
        real(kind=KMCDF),device::Dev_RandArray_SpaceDist(:)
        real(kind=KMCDF),device::Dev_RandArray_SizeDist(:)
        real(kind=KMCDF),value::LNACUT
        real(kind=KMCDF),value::RNACUT
        integer, device::Dev_SEVirtualIndexBox(:,:)
        integer, device::Dev_SEAddedClustersBoxes(:,:)
        real(kind=KMCDF),device::Dev_LayerThick(:)
        real(kind=KMCDF),device::Dev_ClustersSampleRate(:,:)
        real(kind=KMCDF),device::Dev_CompositWeight(:)
        !---Local Vars---
        integer::tid
        integer::bid
        integer::cid
        integer::IBox
        integer::cid0
        integer::ICTRUE
        real(kind=KMCDF)::POS(3)
        integer::NLayer
        integer::MaxGroups
        integer::ILayer
        integer::IGroup
        logical::exitFlag
        real(kind=KMCDF)::tempRand
        real(kind=KMCDF)::GroupRateTemp
        integer::IElement
        type(DiffusorValue)::TheDiffusorValue
        real(kind=KMCDF)::randSize
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*blockdim%x*blockdim%y + tid

        IBox = (cid - 1)/NewAllocateNCEachBox + 1
        cid0 = (IBox - 1)*NewAllocateNCEachBox + 1

        if(cid .LE. TotalAllocateNC) then
            NLayer = size(Dev_ClustersSampleRate,dim=1)
            MaxGroups = size(Dev_ClustersSampleRate,dim=2)

            tempRand = Dev_RandArray_SpaceDist(cid)

            GroupRateTemp = 0.D0
            exitFlag = .false.
            DO ILayer = 1,NLayer

                if(exitFlag .eq. .true.) then
                    exit
                end if

                DO IGroup = 1,MaxGroups
                    GroupRateTemp = GroupRateTemp + Dev_ClustersSampleRate(ILayer,IGroup)
                    if(GroupRateTemp .GE. tempRand) then

                        ICTRUE = Dev_SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1 + (cid - cid0)

                        !Initialize the position of clusters
                        POS(1) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC)*dm_BOXSIZE(1)+dm_BOXBOUNDARY(1,1)
                        POS(2) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*2)*dm_BOXSIZE(2)+dm_BOXBOUNDARY(2,1)
                        POS(3) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*3)*Dev_LayerThick(ILayer) + sum(Dev_LayerThick(1:ILayer-1),dim=1) + dm_BOXBOUNDARY(3,1)
                        Dev_Clusters(ICTRUE)%m_POS = POS

                        !Give the cluster an type(layer) ID for the convenience of visualization
                        Dev_Clusters(ICTRUE)%m_Layer = ILayer

                        !*** Initialize the size of the clusters
                        randSize = Dev_RandArray_SizeDist(cid)
                        if(randSize .GT. RNACUT) then
                            randSize = 2*RNACUT - randSize - (RNACUT-LNACUT)*floor((randSize-RNACUT)/(RNACUT-LNACUT))
                        else if(randSize .LT. LNACUT) then
                            randSize = 2*LNACUT - randSize - (RNACUT-LNACUT)*floor((LNACUT - randSize)/(RNACUT-LNACUT))
                        end if
                        DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                            Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_NA = floor(randSize*Dev_CompositWeight(IElement)+0.5D0)
                            Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_ID = IElement
                        END DO

                        Dev_Clusters(ICTRUE)%m_GrainID(1) = GrainBelongsTo_Dev(Nseeds,Dev_GrainSeeds,POS)

                        Dev_Clusters(ICTRUE)%m_Statu = p_ACTIVEFREE_STATU

                        call Dev_GetValueFromDiffusorsMap(Dev_Clusters(ICTRUE),Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheDiffusorValue)

                        !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
                        !---you should init the distribution by external file---
                        select case(TheDiffusorValue%ECRValueType_Free)
                            case(p_ECR_ByValue)
                                Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_Free
                            case(p_ECR_ByBCluster)
                                Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
                        end select

                        select case(TheDiffusorValue%DiffusorValueType_Free)
                            case(p_DiffuseCoefficient_ByValue)
                                Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                            case(p_DiffuseCoefficient_ByArrhenius)
                                Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/dm_TKB)
                            case(p_DiffuseCoefficient_ByBCluster)
                                ! Here we adopt a model that D=D0*(1/R)**Gama
                                Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_FREESURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
                        end select

                        exitFlag = .true.
                        exit


                    end if
                END DO
            END DO

        end if

        return
    end subroutine Kernel_ImplantClusters_Depth_Layer

    !*************************************************************
    subroutine FillVirtualBoundary_GPU_Depth_SubBox(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        integer,intent(in)::NewAllocateNCEachBox
        !---Local Vars---
        integer::MultiBox
        integer::TotalAllocateNC
        type(dim3)::blocks
        type(dim3)::threads
        integer::NB
        integer::NBX,NBY
        Integer::BX
        integer::BY
        integer::err
        !---Body---

        ASSOCIATE(ImplantRand=>Dev_MigCoaleGVars%dm_MigCoale_RandDev)

            if(NewAllocateNCEachBox .GT. 0) then

                MultiBox = Host_SimuCtrlParam%MultiBox

                TotalAllocateNC = MultiBox*NewAllocateNCEachBox

                NB = (TotalAllocateNC - 1)/p_BLOCKSIZE + 1
                NBX  = min(NB,p_BLOCKDIMX)
                NBY = (NB - 1)/NBX + 1

                !*** to determine the block size
                BX = p_BLOCKSIZE
                BY = 1
                !*** to determine the dimension of blocks

                blocks  = dim3(NBX, NBY, 1)
                threads = dim3(BX,  BY,  1)

                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_X,ImplantRand%dm_SpaceDist_Implant(1:TotalAllocateNC),TotalAllocateNC)
                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Y,ImplantRand%dm_SpaceDist_Implant(TotalAllocateNC+1:2*TotalAllocateNC),TotalAllocateNC)
                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Z,ImplantRand%dm_SpaceDist_Implant(2*TotalAllocateNC+1:3*TotalAllocateNC),TotalAllocateNC)

                err = curandGenerateNormal(ImplantRand%m_ranGen_ClustersSizeDist,ImplantRand%dm_SizeDist_Implant,TotalAllocateNC,this%NAINI,this%NASDINI)

                call Kernel_ImplantClusters_Depth_SubBox<<<blocks,threads>>>(TotalAllocateNC,                                         &
                                                                            NewAllocateNCEachBox,                                     &
                                                                            Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,                 &
                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_TypesEntities,          &
                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_SingleAtomsDivideArrays,&
                                                                            Host_Boxes%m_GrainBoundary%GrainNum,                      &
                                                                            Dev_Boxes%dm_GrainBoundary%dm_GrainSeeds,                 &
                                                                            ImplantRand%dm_SpaceDist_Implant,                         &
                                                                            ImplantRand%dm_SizeDist_Implant,                          &
                                                                            this%NACUT(1),                                            &
                                                                            this%NACUT(2),                                            &
                                                                            Dev_Boxes%dm_SEVirtualIndexBox,                           &
                                                                            this%dm_ImplantInfo_DevPart%Dev_SUBBOXBOUNDARY,           &
                                                                            this%dm_ImplantInfo_DevPart%Dev_CompositWeight)
            end if

        END ASSOCIATE


        return
    end subroutine FillVirtualBoundary_GPU_Depth_SubBox

    !**********************************************
    attributes(global) subroutine Kernel_ImplantClusters_Depth_SubBox(TotalAllocateNC,           &
                                                                    NewAllocateNCEachBox,        &
                                                                    Dev_Clusters,                &
                                                                    Dev_TypesEntities,           &
                                                                    Dev_SingleAtomsDivideArrays, &
                                                                    Nseeds,                      &
                                                                    Dev_GrainSeeds,              &
                                                                    Dev_RandArray_SpaceDist,     &
                                                                    Dev_RandArray_SizeDist,      &
                                                                    LNACUT,                      &
                                                                    RNACUT,                      &
                                                                    Dev_SEVirtualIndexBox,       &
                                                                    Dev_SUBBOXBOUNDARY,          &
                                                                    Dev_CompositWeight)
        implicit none
        !---Dummy Vars---
        integer, value::TotalAllocateNC
        integer, value::NewAllocateNCEachBox
        type(ACluster), device::Dev_Clusters(:)
        type(DiffusorTypeEntity),device::Dev_TypesEntities(:)
        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! If the two dimension array would be delivered to attributes(device), the first dimension must be known
        integer,value::Nseeds
        type(GrainSeed),device::Dev_GrainSeeds(:)
        real(kind=KMCDF),device::Dev_RandArray_SpaceDist(:)
        real(kind=KMCDF),device::Dev_RandArray_SizeDist(:)
        real(kind=KMCDF),value::LNACUT
        real(kind=KMCDF),value::RNACUT
        integer, device::Dev_SEVirtualIndexBox(:,:)
        real(kind=KMCDF),device::Dev_SUBBOXBOUNDARY(:,:)
        real(kind=KMCDF),device::Dev_CompositWeight(:)
        !---Local Vars---
        integer::tid
        integer::bid
        integer::cid
        integer::IBox
        integer::cid0
        integer::ICTRUE
        integer::I
        real(kind=KMCDF)::POS(3)
        integer::IElement
        type(DiffusorValue)::TheDiffusorValue
        real(kind=KMCDF)::randSize
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*blockdim%x*blockdim%y + tid

        IBox = (cid - 1)/NewAllocateNCEachBox + 1
        cid0 = (IBox - 1)*NewAllocateNCEachBox + 1

        if(cid .LE. TotalAllocateNC) then
            ICTRUE = Dev_SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1 + (cid - cid0)

            DO I = 1,3
                POS(I) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*(I - 1))*(Dev_SUBBOXBOUNDARY(I,2) - Dev_SUBBOXBOUNDARY(I,1)) + Dev_SUBBOXBOUNDARY(I,1)
            END DO
            !Initialize the position of clusters
            Dev_Clusters(ICTRUE)%m_POS = POS

            !Give the cluster an type(layer) ID for the convenience of visualization
            Dev_Clusters(ICTRUE)%m_Layer = 1

            !*** Initialize the size of the clusters
            randSize = Dev_RandArray_SizeDist(cid)
            if(randSize .GT. RNACUT) then
                randSize = 2*RNACUT - randSize - (RNACUT-LNACUT)*floor((randSize-RNACUT)/(RNACUT-LNACUT))
            else if(randSize .LT. LNACUT) then
                randSize = 2*LNACUT - randSize - (RNACUT-LNACUT)*floor((LNACUT - randSize)/(RNACUT-LNACUT))
            end if
            DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_NA = floor(randSize*Dev_CompositWeight(IElement)+0.5D0)
                Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_ID = IElement
            END DO

            Dev_Clusters(ICTRUE)%m_GrainID(1) = GrainBelongsTo_Dev(Nseeds,Dev_GrainSeeds,POS)

            Dev_Clusters(ICTRUE)%m_Statu = p_ACTIVEFREE_STATU

            call Dev_GetValueFromDiffusorsMap(Dev_Clusters(ICTRUE),Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheDiffusorValue)

            !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
            !---you should init the distribution by external file---
            select case(TheDiffusorValue%ECRValueType_Free)
                case(p_ECR_ByValue)
                    Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_Free
                case(p_ECR_ByBCluster)
                    Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
            end select

            select case(TheDiffusorValue%DiffusorValueType_Free)
                case(p_DiffuseCoefficient_ByValue)
                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                case(p_DiffuseCoefficient_ByArrhenius)
                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/dm_TKB)
                case(p_DiffuseCoefficient_ByBCluster)
                    ! Here we adopt a model that D=D0*(1/R)**Gama
                    Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_FREESURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
            end select

        end if

        return
    end subroutine Kernel_ImplantClusters_Depth_SubBox

    !*************************************************************
    subroutine FillVirtualBoundary_GPU_Depth_Gauss(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        integer,intent(in)::NewAllocateNCEachBox
        !---Local Vars---
        integer::MultiBox
        integer::TotalAllocateNC
        type(dim3)::blocks
        type(dim3)::threads
        integer::NB
        integer::NBX,NBY
        Integer::BX
        integer::BY
        integer::err
        !---Body---
        ASSOCIATE(ImplantRand=>Dev_MigCoaleGVars%dm_MigCoale_RandDev)

            if(NewAllocateNCEachBox .GT. 0) then

                MultiBox = Host_SimuCtrlParam%MultiBox

                TotalAllocateNC = MultiBox*NewAllocateNCEachBox

                NB = (TotalAllocateNC - 1)/p_BLOCKSIZE + 1
                NBX  = min(NB,p_BLOCKDIMX)
                NBY = (NB - 1)/NBX + 1

                !*** to determine the block size
                BX = p_BLOCKSIZE
                BY = 1
                !*** to determine the dimension of blocks

                blocks  = dim3(NBX, NBY, 1)
                threads = dim3(BX,  BY,  1)

                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_X,ImplantRand%dm_SpaceDist_Implant(1:TotalAllocateNC),TotalAllocateNC)
                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Y,ImplantRand%dm_SpaceDist_Implant(TotalAllocateNC+1:2*TotalAllocateNC),TotalAllocateNC)
                err = curandGenerateNormal(ImplantRand%m_ranGen_ClustersSpaceDist_Z,ImplantRand%dm_SpaceDist_Implant(2*TotalAllocateNC+1:3*TotalAllocateNC),TotalAllocateNC,this%DepthINI,this%DepthSDINI)

                err = curandGenerateNormal(ImplantRand%m_ranGen_ClustersSizeDist,ImplantRand%dm_SizeDist_Implant,TotalAllocateNC,this%NAINI,this%NASDINI)

                call Kernel_ImplantClusters_Depth_Gauss<<<blocks,threads>>>(TotalAllocateNC,                                          &
                                                                            NewAllocateNCEachBox,                                     &
                                                                            Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,                 &
                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_TypesEntities,          &
                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_SingleAtomsDivideArrays,&
                                                                            Host_Boxes%m_GrainBoundary%GrainNum,                      &
                                                                            Dev_Boxes%dm_GrainBoundary%dm_GrainSeeds,                 &
                                                                            ImplantRand%dm_SpaceDist_Implant,                         &
                                                                            ImplantRand%dm_SizeDist_Implant,                          &
                                                                            this%NACUT(1),                                            &
                                                                            this%NACUT(2),                                            &
                                                                            Dev_Boxes%dm_SEVirtualIndexBox,                           &
                                                                            this%dm_ImplantInfo_DevPart%Dev_CompositWeight)
            end if

        END ASSOCIATE

        return
    end subroutine FillVirtualBoundary_GPU_Depth_Gauss

    !**********************************************
    attributes(global) subroutine Kernel_ImplantClusters_Depth_Gauss(TotalAllocateNC,            &
                                                                    NewAllocateNCEachBox,        &
                                                                    Dev_Clusters,                &
                                                                    Dev_TypesEntities,           &
                                                                    Dev_SingleAtomsDivideArrays, &
                                                                    Nseeds,                      &
                                                                    Dev_GrainSeeds,              &
                                                                    Dev_RandArray_SpaceDist,     &
                                                                    Dev_RandArray_SizeDist,      &
                                                                    LNACUT,                      &
                                                                    RNACUT,                      &
                                                                    Dev_SEVirtualIndexBox,       &
                                                                    Dev_CompositWeight)
        implicit none
        !---Dummy Vars---
        integer, value::TotalAllocateNC
        integer, value::NewAllocateNCEachBox
        type(ACluster), device::Dev_Clusters(:)
        type(DiffusorTypeEntity),device::Dev_TypesEntities(:)
        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! If the two dimension array would be delivered to attributes(device), the first dimension must be known
        integer,value::Nseeds
        type(GrainSeed),device::Dev_GrainSeeds(:)
        real(kind=KMCDF),device::Dev_RandArray_SpaceDist(:)
        real(kind=KMCDF),device::Dev_RandArray_SizeDist(:)
        real(kind=KMCDF),value::LNACUT
        real(kind=KMCDF),value::RNACUT
        integer, device::Dev_SEVirtualIndexBox(:,:)
        real(kind=KMCDF),device::Dev_CompositWeight(:)
        !---Local Vars---
        integer::tid
        integer::bid
        integer::cid
        integer::IBox
        integer::cid0
        integer::ICTRUE
        integer::I
        real(kind=KMCDF)::POS(3)
        integer::IElement
        type(DiffusorValue)::TheDiffusorValue
        real(kind=KMCDF)::randSize
        real(kind=KMCDF)::randDepth
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*blockdim%x*blockdim%y + tid

        IBox = (cid - 1)/NewAllocateNCEachBox + 1
        cid0 = (IBox - 1)*NewAllocateNCEachBox + 1

        if(cid .LE. TotalAllocateNC) then
            ICTRUE = Dev_SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1 + (cid - cid0)

            POS(1) = Dev_RandArray_SpaceDist(cid)*dm_BOXSIZE(1)+dm_BOXBOUNDARY(1,1)
            POS(2) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC)*dm_BOXSIZE(2)+dm_BOXBOUNDARY(2,1)

            randDepth = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*2)
            if(randDepth .GT. dm_BOXBOUNDARY(3,2)) then
                randDepth = 2*dm_BOXBOUNDARY(3,2) - randDepth - (dm_BOXBOUNDARY(3,2)-dm_BOXBOUNDARY(3,1))*floor((randDepth-dm_BOXBOUNDARY(3,2))/(dm_BOXBOUNDARY(3,2)-dm_BOXBOUNDARY(3,1)))
            else if(randDepth .LT. dm_BOXBOUNDARY(3,1)) then
                randDepth = 2*dm_BOXBOUNDARY(3,1) - randDepth - (dm_BOXBOUNDARY(3,2)-dm_BOXBOUNDARY(3,1))*floor((dm_BOXBOUNDARY(3,1)-randDepth)/(dm_BOXBOUNDARY(3,2)-dm_BOXBOUNDARY(3,1)))
            end if

            POS(3) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*2)*dm_BOXSIZE(3) + dm_BOXBOUNDARY(3,1)
            !Initialize the position of clusters
            Dev_Clusters(ICTRUE)%m_POS = POS

            !Give the cluster an type(layer) ID for the convenience of visualization
            Dev_Clusters(ICTRUE)%m_Layer = 1

            !*** Initialize the size of the clusters
            randSize = Dev_RandArray_SizeDist(cid)
            if(randSize .GT. RNACUT) then
                randSize = 2*RNACUT - randSize - (RNACUT-LNACUT)*floor((randSize-RNACUT)/(RNACUT-LNACUT))
            else if(randSize .LT. LNACUT) then
                randSize = 2*LNACUT - randSize - (RNACUT-LNACUT)*floor((LNACUT - randSize)/(RNACUT-LNACUT))
            end if
            DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_NA = floor(randSize*Dev_CompositWeight(IElement)+0.5D0)
                Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_ID = IElement
            END DO

            Dev_Clusters(ICTRUE)%m_GrainID(1) = GrainBelongsTo_Dev(Nseeds,Dev_GrainSeeds,POS)

            Dev_Clusters(ICTRUE)%m_Statu = p_ACTIVEFREE_STATU

            call Dev_GetValueFromDiffusorsMap(Dev_Clusters(ICTRUE),Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheDiffusorValue)

            !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
            !---you should init the distribution by external file---
            select case(TheDiffusorValue%ECRValueType_Free)
                case(p_ECR_ByValue)
                    Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_Free
                case(p_ECR_ByBCluster)
                    Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
            end select

            select case(TheDiffusorValue%DiffusorValueType_Free)
                case(p_DiffuseCoefficient_ByValue)
                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                case(p_DiffuseCoefficient_ByArrhenius)
                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/dm_TKB)
                case(p_DiffuseCoefficient_ByBCluster)
                    ! Here we adopt a model that D=D0*(1/R)**Gama
                    Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_FREESURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
            end select

        end if

        return
    end subroutine Kernel_ImplantClusters_Depth_Gauss

    !*************************************************************
    subroutine FillVirtualBoundary_GPU_FromFile(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        integer,intent(in)::NewAllocateNCEachBox
        !---Local Vars---
        integer::MultiBox
        integer::TotalAllocateNC
        type(dim3)::blocks
        type(dim3)::threads
        integer::NB
        integer::NBX,NBY
        Integer::BX
        integer::BY
        integer::err
        !---Body---
        ASSOCIATE(ImplantRand=>Dev_MigCoaleGVars%dm_MigCoale_RandDev)

            if(NewAllocateNCEachBox .GT. 0) then

                MultiBox = Host_SimuCtrlParam%MultiBox

                TotalAllocateNC = MultiBox*NewAllocateNCEachBox

                NB = (TotalAllocateNC - 1)/p_BLOCKSIZE + 1
                NBX  = min(NB,p_BLOCKDIMX)
                NBY = (NB - 1)/NBX + 1

                !*** to determine the block size
                BX = p_BLOCKSIZE
                BY = 1
                !*** to determine the dimension of blocks

                blocks  = dim3(NBX, NBY, 1)
                threads = dim3(BX,  BY,  1)

                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Layer,ImplantRand%dm_SpaceDist_Implant(1:TotalAllocateNC),TotalAllocateNC)
                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_X,ImplantRand%dm_SpaceDist_Implant(TotalAllocateNC+1:2*TotalAllocateNC),TotalAllocateNC)

                call Kernel_ImplantClusters_FromFile<<<blocks,threads>>>(TotalAllocateNC,                                         &
                                                                        NewAllocateNCEachBox,                                     &
                                                                        Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,                 &
                                                                        Dev_Boxes%dm_DiffusorTypesMap%Dev_TypesEntities,          &
                                                                        Dev_Boxes%dm_DiffusorTypesMap%Dev_SingleAtomsDivideArrays,&
                                                                        Host_Boxes%m_GrainBoundary%GrainNum,                      &
                                                                        Dev_Boxes%dm_GrainBoundary%dm_GrainSeeds,                 &
                                                                        ImplantRand%dm_SpaceDist_Implant,                         &
                                                                        Dev_Boxes%dm_SEVirtualIndexBox,                           &
                                                                        this%dm_ImplantInfo_DevPart%Dev_LayerThick,               &
                                                                        this%dm_ImplantInfo_DevPart%Dev_ClustersSample,           &
                                                                        this%dm_ImplantInfo_DevPart%Dev_ClustersSampleRate)
            end if

        END ASSOCIATE

        return
    end subroutine FillVirtualBoundary_GPU_FromFile


    !**********************************************
    attributes(global) subroutine Kernel_ImplantClusters_FromFile(TotalAllocateNC,             &
                                                                  NewAllocateNCEachBox,        &
                                                                  Dev_Clusters,                &
                                                                  Dev_TypesEntities,                &
                                                                  Dev_SingleAtomsDivideArrays, &
                                                                  Nseeds,                      &
                                                                  Dev_GrainSeeds,              &
                                                                  Dev_RandArray_SpaceDist,     &
                                                                  Dev_SEVirtualIndexBox,       &
                                                                  Dev_LayerThick,              &
                                                                  Dev_ClustersSample,          &
                                                                  Dev_ClustersSampleRate)
        implicit none
        !---Dummy Vars---
        integer, value::TotalAllocateNC
        integer, value::NewAllocateNCEachBox
        type(ACluster), device::Dev_Clusters(:)
        type(DiffusorTypeEntity),device::Dev_TypesEntities(:)
        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! If the two dimension array would be delivered to attributes(device), the first dimension must be known
        integer,value::Nseeds
        type(GrainSeed),device::Dev_GrainSeeds(:)
        real(kind=KMCDF),device::Dev_RandArray_SpaceDist(:)
        integer, device::Dev_SEVirtualIndexBox(:,:)
        real(kind=KMCDF),device::Dev_LayerThick(:)
        type(ACluster),device::Dev_ClustersSample(:,:)
        real(kind=KMCDF),device::Dev_ClustersSampleRate(:,:)
        !---Local Vars---
        integer::tid
        integer::bid
        integer::cid
        integer::IBox
        integer::cid0
        integer::ICTRUE
        real(kind=KMCDF)::POS(3)
        integer::NLayer
        integer::MaxGroups
        integer::ILayer
        integer::IGroup
        logical::exitFlag
        real(kind=KMCDF)::tempRand
        real(kind=KMCDF)::GroupRateTemp
        type(DiffusorValue)::TheDiffusorValue
        !---Body---
        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
        bid = (blockidx%y - 1)*griddim%x + blockidx%x
        cid = (bid - 1)*blockdim%x*blockdim%y + tid

        IBox = (cid - 1)/NewAllocateNCEachBox + 1
        cid0 = (IBox - 1)*NewAllocateNCEachBox + 1

        if(cid .LE. TotalAllocateNC) then
            NLayer = size(Dev_ClustersSampleRate,dim=1)
            MaxGroups = size(Dev_ClustersSampleRate,dim=2)

            tempRand = Dev_RandArray_SpaceDist(cid)

            GroupRateTemp = 0.D0
            exitFlag = .false.
            DO ILayer = 1,NLayer

                if(exitFlag .eq. .true.) then
                    exit
                end if

                DO IGroup = 1,MaxGroups
                    GroupRateTemp = GroupRateTemp + Dev_ClustersSampleRate(ILayer,IGroup)
                    if(GroupRateTemp .GE. tempRand) then

                        ICTRUE = Dev_SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1 + (cid - cid0)

                        Dev_Clusters(ICTRUE)%m_Atoms = Dev_ClustersSample(ILayer,IGroup)%m_Atoms

                        !Initialize the position of clusters
                        POS(1) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC)*dm_BOXSIZE(1)+dm_BOXBOUNDARY(1,1)
                        POS(2) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*2)*dm_BOXSIZE(2)+dm_BOXBOUNDARY(2,1)
                        POS(3) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*3)*Dev_LayerThick(ILayer) + sum(Dev_LayerThick(1:ILayer-1),dim=1) + dm_BOXBOUNDARY(3,1)
                        Dev_Clusters(ICTRUE)%m_POS = POS

                        !Give the cluster an type(layer) ID for the convenience of visualization
                        Dev_Clusters(ICTRUE)%m_Layer = ILayer

                        Dev_Clusters(ICTRUE)%m_Statu = Dev_ClustersSample(ILayer,IGroup)%m_Statu

                        call Dev_GetValueFromDiffusorsMap(Dev_Clusters(ICTRUE),Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheDiffusorValue)

                        if(Dev_Clusters(ICTRUE)%m_Statu .eq. p_ACTIVEFREE_STATU) then

                            Dev_Clusters(ICTRUE)%m_GrainID(1) = GrainBelongsTo_Dev(Nseeds,Dev_GrainSeeds,POS)

                            select case(TheDiffusorValue%ECRValueType_Free)
                                case(p_ECR_ByValue)
                                    Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_Free
                                case(p_ECR_ByBCluster)
                                    Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
                            end select

                            select case(TheDiffusorValue%DiffusorValueType_Free)
                                case(p_DiffuseCoefficient_ByValue)
                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
                                case(p_DiffuseCoefficient_ByArrhenius)
                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/dm_TKB)
                                case(p_DiffuseCoefficient_ByBCluster)
                                    ! Here we adopt a model that D=D0*(1/R)**Gama
                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_FREESURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
                            end select

                        else if(Dev_Clusters(ICTRUE)%m_Statu .eq. p_ACTIVEINGB_STATU) then

                            Dev_Clusters(ICTRUE)%m_GrainID = Dev_ClustersSample(ILayer,IGroup)%m_GrainID

                            select case(TheDiffusorValue%ECRValueType_InGB)
                                case(p_ECR_ByValue)
                                    Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_InGB
                                case(p_ECR_ByBCluster)
                                    Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
                            end select

                            select case(TheDiffusorValue%DiffusorValueType_InGB)
                                case(p_DiffuseCoefficient_ByValue)
                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_InGB_Value
                                case(p_DiffuseCoefficient_ByArrhenius)
                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_InGB*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_InGB/dm_TKB)
                                case(p_DiffuseCoefficient_ByBCluster)
                                    ! Here we adopt a model that D=D0*(1/R)**Gama
                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_GBSURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
                            end select
                        end if

                        exitFlag = .true.
                        exit

                    end if
                END DO
            END DO

        end if

        return
    end subroutine Kernel_ImplantClusters_FromFile

    !*************************************************************
    subroutine FillVirtualBoundary_GPU_FromExteFunc(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(SimulationBoxes_GPU)::Dev_Boxes
        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
        integer,intent(in)::NewAllocateNCEachBox
        !---Body---

        return
    end subroutine FillVirtualBoundary_GPU_FromExteFunc


!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!    !*********************************************************************
!    subroutine ImplantClusters(Host_Boxes,Host_SimuCtrlParam,Record,Dev_Boxes,TSTEP,ImplantClusterType,Flux)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(ImplantRecord)::Record
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        real(kind=KMCDF)::TSTEP
!        type(ACluster)::ImplantClusterType
!        real(kind=KMCDF)::Flux
!        !---Local Vars---
!        integer::err
!        integer::MultiBox
!        integer::IBox
!        integer::AddedTotalClustersNum
!        integer::NB,NBX,NBY,BX,BY
!        type(dim3)::blocks,threads
!        integer::NewTotalSize
!        integer::OldActSize
!        integer::ITYPE
!        !---Body---
!
!        !Associate(Dev_ClustersInfo=>Dev_Boxes%dm_ClusterInfo_GPU)
!
!        MultiBox = Host_SimuCtrlParam%MultiBox
!
!
!        dm_Layer  = ImplantClusterType%m_Layer
!        dm_RAD = ImplantClusterType%m_RAD
!        dm_TYPE   = ImplantClusterType%m_TYPE
!        dm_Atoms = ImplantClusterType%m_Atoms
!
!
!        OldActSize = Host_Boxes%m_BoxesInfo%SEActIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEActIndexBox(1,1) + 1
!
!        call Cal_ImplantClusters_Num(Host_Boxes,Host_SimuCtrlParam,Record,Dev_Boxes,TSTEP,Flux)
!
!        AddedTotalClustersNum = Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes(1,1) + 1
!
!        if(AddedTotalClustersNum .GE. 1 .AND. Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes(1,1) .GT. 0) then
!
!            err = curandGenerateUniformDouble(m_ranGen_ClustersSpaceDist,Dev_Boxes%dm_ClusterInfo_GPU%dm_SpaceDist_Implant,AddedTotalClustersNum*3)
!
!            BX = p_BLOCKSIZE
!            BY = 1
!            NB = (AddedTotalClustersNum-1)/p_BLOCKSIZE + 1
!            NBX = min(NB, p_BLOCKDIMX)
!            NBY = (NB - 1)/NBX + 1
!            NB = NBX*NBY
!
!            threads = dim3(BX,BY,1)
!            blocks = dim3(NBX,NBY,1)
!
!            NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!
!            call Kernel_ImplantClusters_New<<<blocks,threads>>>(MultiBox,                                &
!                                                                NewTotalSize,                            &
!                                                                Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,            &
!                                                                Dev_Boxes%dm_ClusterInfo_GPU%dm_ActiveIndex,         &
!                                                                Dev_Boxes%dm_ClusterInfo_GPU%dm_SpaceDist_Implant,   &
!                                                                AddedTotalClustersNum,                   &
!                                                                Dev_Boxes%dm_SEAddedClustersBoxes,       &
!                                                                Dev_Boxes%dm_SEActIndexBox,              &
!                                                                Dev_Boxes%dm_SEVirtualIndexBox,           &
!                                                                Dev_Boxes%dm_ClusterInfo_GPU%dm_ActiveStatus)
!
!
!            DO IBox = 1,MultiBox
!
!                OldActSize = Host_Boxes%m_BoxesInfo%SEActIndexBox(IBox,2) - Host_Boxes%m_BoxesInfo%SEActIndexBox(IBox,1) + 1
!
!                if(IBox .eq. 1) then
!                    Host_Boxes%m_BoxesInfo%SEActIndexBox(IBox,1) = 1
!                else
!                    Host_Boxes%m_BoxesInfo%SEActIndexBox(IBox,1) = Host_Boxes%m_BoxesInfo%SEActIndexBox(IBox-1,2) + 1
!                end if
!
!                Host_Boxes%m_BoxesInfo%SEActIndexBox(IBox,2) = Host_Boxes%m_BoxesInfo%SEActIndexBox(IBox,1) + OldActSize + AddedTotalClustersNum - 1
!
!                ITYPE = ImplantClusterType%m_TYPE
!
!                Host_Boxes%m_StatisticInfo_Integral%NC(ITYPE,p_ACTIVEFREE_STATU) = Host_Boxes%m_StatisticInfo_Integral%NC(ITYPE,p_ACTIVEFREE_STATU) &
!                                                                                 + AddedTotalClustersNum
!
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(ITYPE,p_ACTIVEFREE_STATU) = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(ITYPE,p_ACTIVEFREE_STATU) &
!                                                                                     + AddedTotalClustersNum
!            END DO
!
!            Dev_Boxes%dm_SEActIndexBox = Host_Boxes%m_BoxesInfo%SEActIndexBox
!
!        end if
!
!
!        !END Associate
!
!        return
!    end subroutine ImplantClusters
!
!    !*********************************************
!    subroutine Cal_ImplantClusters_Num(Host_Boxes,Host_SimuCtrlParam,Record,Dev_Boxes,TSTEP,Flux)
!        use cudafor
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(ImplantRecord)::Record
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        real(kind=KMCDF), intent(inout)::TSTEP
!        real(kind=KMCDF)::Flux
!        !---Local Vars---
!        integer::err
!        real(kind=KMCDF)::Host_SurfaceArea
!        real(kind=KMCDF)::newSimulatedTime
!        integer::MultiBox
!        integer::IBox
!        integer(kind=cuda_count_kind)::FreeMemSize,TotalMemSize
!        integer::FreeMemSize2NC
!        integer::TotalMemSize2NC
!        integer::newAllocatedNCSize
!        integer::tempInjectNumber
!        integer::tempMaxNCActive,tempNC
!        integer::AllocatedFreeNCSize,newFreeNCSize
!        real(kind=KMCDF)::verifyTime
!        integer::PreAllocatedSystemSize
!        real(kind=KMCDF)::MaxDifCoeffes
!        integer,dimension(:),allocatable::tempArray
!        !---Body---
!
!        FreeMemSize2NC = 0
!        TotalMemSize2NC = 0
!        newAllocatedNCSize = 0
!        tempInjectNumber = 0
!
!        MultiBox = Host_SimuCtrlParam%MultiBox
!
!        Host_SurfaceArea = Host_Boxes%BOXSIZE(1)*Host_Boxes%BOXSIZE(2)
!
!        !---Each box own same memory size
!        PreAllocatedSystemSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!
!        allocate(tempArray(MultiBox))
!
!        FORALL(IBox=1:MultiBox)
!            tempArray(IBox) = sum(Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(:,p_ACTIVEFREE_STATU) + Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(:,p_ACTIVEINGB_STATU))
!        END FORALL
!
!        tempMaxNCActive = maxval(tempArray)
!
!        AllocatedFreeNCSize = PreAllocatedSystemSize - tempMaxNCActive
!
!        tempNC = tempMaxNCActive
!
!        call Get_MaxDiffusion(Host_Boxes,Host_SimuCtrlParam,DIFCOES(1),MaxDifCoeffes)
!
!        !call Cal_TimeStep_ByDensityChanging(Host_Boxes,m_Implant_DensityDrop,m_LC_DiffCoeffesTable(p_LCluster),p_LC_CaptureRadium,TSTEP)
!
!        !write(*,*) "Cal_TimeStep_ByDensityChanging",TSTEP
!
!        !call Cal_SegmentsTimeStep(Host_Boxes,0,MaxDifCoeffes,TSTEP)
!
!        !write(*,*) "Cal_SegmentsTimeStep",TSTEP
!
!        call Cal_TimeStep_New(Host_Boxes,Host_SimuCtrlParam,0,MaxDifCoeffes,TSTEP)
!
!        !write(*,*) "Cal_TimeStep",TSTEP
!
!        DO while(.true.)
!
!            newSimulatedTime = Record%GetSimuTimes() + TSTEP
!
!            newAllocatedNCSize = 0
!
!            newFreeNCSize = AllocatedFreeNCSize
!
!            tempInjectNumber = max(floor(newSimulatedTime*Flux*Host_SurfaceArea) - Record%Get(),0)
!
!            tempNC = tempMaxNCActive + tempInjectNumber
!
!            verifyTime = 1.D0/(50*MaxDifCoeffes*((tempNC/Host_Boxes%BOXVOLUM)**0.666666666666666667))
!
!            if(verifyTime .GT. TSTEP) then
!                if(newFreeNCSize .GE. tempInjectNumber) then
!                    exit
!                else
!
!                    call mE_MEMORYMAN_EVOLUTION%GetDeviceMemInfo(FreeMemSize,TotalMemSize)
!
!                    FreeMemSize2NC  = FreeMemSize/mE_MEMORYMAN_EVOLUTION%GetOneClusterMemConsume_GPU()
!                    TotalMemSize2NC = TotalMemSize/mE_MEMORYMAN_EVOLUTION%GetOneClusterMemConsume_GPU()
!
!                    DO while(.true.)
!
!                            if((FreeMemSize2NC - newAllocatedNCSize) .LE. 0) then
!                                exit
!                            end if
!
!                            newAllocatedNCSize = newAllocatedNCSize + min(FreeMemSize2NC,TotalMemSize2NC/10)
!
!                            newFreeNCSize = AllocatedFreeNCSize + newAllocatedNCSize
!
!                            if(newFreeNCSize .GE. tempInjectNumber) then
!                                exit
!                            end if
!
!                    END DO
!
!                    if(newFreeNCSize .GE. tempInjectNumber) then
!                        exit
!                    else
!                        TSTEP = TSTEP/2.D0
!                    end if
!
!                end if
!
!            else
!                TSTEP = TSTEP/2.D0
!            end if
!
!        END DO
!
!        call Record%AddImplantedEntitiesNum(tempInjectNumber)
!
!        if(newAllocatedNCSize .GT. 0) then  ! need to allocate new memory for new added cluster
!
!            call Dev_Boxes%ExpandClustersInfo_GPUToCPU(Host_Boxes,Host_SimuCtrlParam,newAllocatedNCSize)
!
!        end if
!
!
!        DO IBox = 1,MultiBox
!
!            if(IBox .eq. 1) then
!                Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes(IBox,1) = 1
!            else
!                Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes(IBox,1) = Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes(IBox-1,2) + 1
!            end if
!            Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes(IBox,2) = Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes(IBox,1) + tempInjectNumber - 1
!
!        END DO
!
!        Dev_Boxes%dm_SEAddedClustersBoxes = Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes
!
!        deallocate(tempArray)
!
!        return
!    end subroutine Cal_ImplantClusters_Num
!
!
!    !**********************************************
!    subroutine Cal_TimeStep_New(Host_Boxes,Host_SimuCtrlParam,AddedNC,MaxDifCoeffes,TSTEP)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes),intent(in)::Host_Boxes
!        type(SimulationCtrlParam),intent(in)::Host_SimuCtrlParam
!        integer, intent(in)::AddedNC
!        real(kind=KMCDF),intent(in)::MaxDifCoeffes
!        real(kind=KMCDF)::TSTEP
!        !---Local Vars---
!        real(kind=KMCDF)::SEP
!        integer::TotalActiveNC
!        integer::TotalActiveNCBCluster
!        integer::NewTotalActiveNC
!        real(kind=KMCDF)::RAVA
!        !---Body---
!
!        ASSOCIATE(TBasicInfo=>Host_Boxes%m_StatisticInfo_Integral)
!            if(Host_SimuCtrlParam%UPDATETSTEPSTRATEGY .eq. mp_SelfAdjustlStep) then
!                TotalActiveNC = sum(TBasicInfo%NC(:,p_ACTIVEFREE_STATU) + TBasicInfo%NC(:,p_ACTIVEINGB_STATU))
!
!                TotalActiveNCBCluster = TBasicInfo%NC(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU)
!
!                NewTotalActiveNC =  TotalActiveNC + AddedNC
!
!                if(NewTotalActiveNC .LE. 0) then
!                    NewTotalActiveNC = 1
!
!                    SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(NewTotalActiveNC)
!                    SEP = SEP**(0.33333333333333D0)
!                    TSTEP = SEP*SEP/(6.D0*MaxDifCoeffes)
!                else
!                    if(TotalActiveNCBCluster .eq. 0) then
!
!                        SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(NewTotalActiveNC)
!                        SEP = SEP**(0.33333333333333D0)
!                        TSTEP = SEP*SEP/(6.D0*MaxDifCoeffes)
!                    else
!                        RAVA = (TBasicInfo%RAVA(p_ACTIVEFREE_STATU)*TBasicInfo%NC(p_ACTIVEFREE_STATU) + TBasicInfo%RAVA(p_ACTIVEINGB_STATU)*TBasicInfo%NC(p_ACTIVEINGB_STATU))/TotalActiveNCBCluster
!
!                        SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(NewTotalActiveNC)
!                        SEP = SEP**(0.33333333333333D0)
!                        SEP = SEP - 2.D0*RAVA
!
!                        TSTEP = SEP*SEP/(6.D0*MaxDifCoeffes)*Host_SimuCtrlParam%TStepValue
!                    end if
!
!                end if
!            else
!                TSTEP = Host_SimuCtrlParam%TStepValue
!            end if
!
!        END ASSOCIATE
!
!        return
!    end subroutine Cal_TimeStep_New
!
!    !**********************************************
!    subroutine Cal_TimeStep(Host_Boxes,Host_SimuCtrlParam,AddedNC,MaxDifCoeffes,TSTEP)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        integer::AddedNC
!        real(kind=KMCDF)::MaxDifCoeffes
!        real(kind=KMCDF)::TSTEP
!        !---Local Vars---
!        real(kind=KMCDF)::SEP
!        integer::NCAct
!        integer::NCActBCLuster
!        real(kind=KMCDF)::RAVA
!        !---Body---
!        ASSOCIATE(TBasicInfo=>Host_Boxes%m_StatisticInfo_Integral)
!            if(Host_SimuCtrlParam%UPDATETSTEPSTRATEGY .eq. mp_SelfAdjustlStep) then
!                NCAct = sum(TBasicInfo%NC(:,p_ACTIVEFREE_STATU) + TBasicInfo%NC(:,p_ACTIVEINGB_STATU))
!
!                NCActBCLuster = TBasicInfo%NC(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU)
!
!                RAVA = (TBasicInfo%NC(p_ACTIVEFREE_STATU)*TBasicInfo%RAVA(p_ACTIVEFREE_STATU) + TBasicInfo%NC(p_ACTIVEINGB_STATU)*TBasicInfo%RAVA(p_ACTIVEINGB_STATU))/dble(NCActBCLuster)
!
!                SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(NCAct + AddedNC)
!                SEP = SEP**(0.33333333333333D0)
!                SEP = SEP - 2.D0*RAVA
!
!                TSTEP = SEP*SEP/(6.D0*MaxDifCoeffes)*Host_SimuCtrlParam%TStepValue
!            else
!                TSTEP = Host_SimuCtrlParam%TStepValue
!            end if
!
!        END ASSOCIATE
!
!        return
!    end subroutine Cal_TimeStep
!
!    !**********************************************
!    subroutine Cal_TimeStep_ByDensityChanging(Host_Boxes,Host_SimuCtrlParam,DensityDrop,SinglePartDifCoeffes,SinglePartECRS,TSTEP)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes),intent(in)::Host_Boxes
!        type(SimulationCtrlParam),intent(in)::Host_SimuCtrlParam
!        real(kind=KMCDF),intent(in)::DensityDrop
!        real(kind=KMCDF),intent(in)::SinglePartDifCoeffes
!        real(kind=KMCDF),intent(in)::SinglePartECRS
!        real(kind=KMCDF),intent(out)::TSTEP
!        !---Local Vars---
!        real(kind=KMCDF)::OriginalDensity
!        integer::NAct
!        !---Body---
!        NAct = sum(Host_Boxes%m_StatisticInfo_Integral%NC(:,p_ACTIVEFREE_STATU) + Host_Boxes%m_StatisticInfo_Integral%NC(:,p_ACTIVEINGB_STATU))
!        OriginalDensity = NAct/dble(Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM)
!        !***Ref: [1] and deduce the following relationship ***
!        TSTEP = ((1-DensityDrop)**(-0.5D0) - 1)/(4*PI*SinglePartDifCoeffes*SinglePartECRS*OriginalDensity)
!
!        return
!    end subroutine Cal_TimeStep_ByDensityChanging
!
!    !**********************************************
!    subroutine Cal_SegmentsTimeStep(Host_Boxes,Host_SimuCtrlParam,AddedNC,MaxDifCoeffes,SegmentsTSTEP)
!        implicit none
!        !---Dummy Vars---
!        type(SimulationBoxes),intent(in)::Host_Boxes
!        type(SimulationCtrlParam),intent(in)::Host_SimuCtrlParam
!        integer, intent(in)::AddedNC
!        real(kind=KMCDF),intent(in)::MaxDifCoeffes
!        real(kind=KMCDF)::SegmentsTSTEP
!        !---Local Vars---
!        real(kind=KMCDF)::SEP
!        integer::NAct
!        !---Body---
!        NAct = sum(Host_Boxes%m_StatisticInfo_Integral%NC(:,p_ACTIVEFREE_STATU) + Host_Boxes%m_StatisticInfo_Integral%NC(:,p_ACTIVEINGB_STATU))
!        SEP = Host_SimuCtrlParam%MultiBox*Host_Boxes%BOXVOLUM/dble(NAct + AddedNC)
!        SEP = SEP**(0.33333333333333D0)
!
!        SegmentsTSTEP = SEP*SEP/(6.D0*MaxDifCoeffes)
!
!        return
!    end subroutine Cal_SegmentsTimeStep
!
!

!    !**********************************************
!    attributes(global) subroutine Kernel_ImplantClusters_New(MultiBox,TNC,Dev_Clusters,Dev_ActiveIndex,Dev_SpaceDist, &
!                                                             AddedTotalClustersNum,Dev_SEAddedClustersBoxes,Dev_SEActIndexBox,Dev_SEIndexBox,Dev_ActiveStatu)
!        implicit none
!        !---Dummy Vars---
!        integer, value::MultiBox
!        integer, value::TNC
!        type(ACluster), device, dimension(TNC)::Dev_Clusters
!        integer, device, dimension(TNC)::Dev_ActiveIndex
!        real(kind=KMCDF), dimension(TNC*3)::Dev_SpaceDist
!        integer, value::AddedTotalClustersNum
!        integer, device, dimension(MultiBox,2)::Dev_SEAddedClustersBoxes
!        integer, device, dimension(MultiBox,2)::Dev_SEActIndexBox
!        integer, device, dimension(MultiBox,2)::Dev_SEIndexBox
!        integer,device,dimension(TNC)::Dev_ActiveStatu
!        !---Local Vars---
!        integer::IT,IB,IC
!        integer::I,IBox
!        integer::SBoxClusterNumber
!        integer::tempIndex,mappedIndex,relativeToFirstdex
!        !---Body---
!        IT = (threadidx%y - 1)*blockdim%x + threadidx%x
!        IB = (blockidx%y - 1)*griddim%x + blockidx%x
!        IC = (IB - 1)*blockdim%x*blockdim%y + IT
!
!        if(IC .LE. AddedTotalClustersNum) then
!
!            ! binary way to find the box
!            call BinarySearchBoxIndex(MultiBox,IC,Dev_SEAddedClustersBoxes,IBox)
!
!            relativeToFirstdex = IC - Dev_SEAddedClustersBoxes(IBox,1) + 1
!            tempIndex = Dev_SEActIndexBox(IBox,2) + relativeToFirstdex
!            mappedIndex = Dev_ActiveIndex(tempIndex)
!
!            !---Construct each clusters--
!            Dev_ActiveStatu(mappedIndex) = p_ACTIVEFREE_STATU
!
!            Dev_Clusters(mappedIndex)%m_Pos(1) = Dev_SpaceDist(IC)*dm_BOXSIZE(1) + dm_BOXBOUNDARY(1,1)
!            Dev_Clusters(mappedIndex)%m_Pos(2) = Dev_SpaceDist(IC + AddedTotalClustersNum)*dm_BOXSIZE(2) + dm_BOXBOUNDARY(2,1)
!            Dev_Clusters(mappedIndex)%m_Pos(3) = Dev_SpaceDist(IC + 2*AddedTotalClustersNum)*dm_BOXSIZE(3) + dm_BOXBOUNDARY(3,1)
!            Dev_Clusters(mappedIndex)%m_Layer = dm_Layer
!            Dev_Clusters(mappedIndex)%m_Atoms = dm_Atoms
!            Dev_Clusters(mappedIndex)%m_RAD = dm_RAD
!            Dev_Clusters(mappedIndex)%m_Statu = p_ACTIVEFREE_STATU
!            Dev_Clusters(mappedIndex)%m_TYPE = dm_TYPE
!
!        end if
!
!        return
!    end subroutine Kernel_ImplantClusters_New
!
!
!  !************************************************************************
!    attributes(device) subroutine Cluster_TypeDefine_GPU(NAtom,Atoms,ClusterType,Radium)
!        implicit none
!        !---Dummy Vars---
!        integer,intent(in)::NAtom
!        type(Single_AtomsSet),dimension(p_ATOMS_GROUPS_NUMBER)::Atoms
!        integer, intent(inout)::ClusterType
!        real,intent(inout)::Radium
!        !---Body---
!
!        if(NAtom .LT. p_BCNA) then
!            ! @todo (zhail#1#): The LCluster types need to be completed.
!            ClusterType = p_LCluster
!
!            Radium = p_LC_CaptureRadium
!
!        else
!            ClusterType = p_BCluster
!
!            ! Convert number of atoms to bubble radius
!            ! Ref. Modelling Simul. Mater. Sci. Eng.16(2008)055003
!            Radium = SQRT(NAtom/dm_RNFACTOR)
!        end if
!
!        return
!    end subroutine Cluster_TypeDefine_GPU

end module MIGCOALE_IMPLANTATION_GPU
