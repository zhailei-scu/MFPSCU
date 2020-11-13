module MIGCOALE_IMPLANTATION
    use MCMF_TYPEDEF_USUAL
    use MCMF_TYPEDEF_ACLUSTER
    use MFLIB_GLOBAL
    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
    use MIGCOALE_TIMECTL
    use MIGCOALE_TYPEDEF_STATISTICINFO
    use MIGCOALE_STATISTIC_CPU
    use MIGCOALE_TYPEDEF_SIMRECORD
    use MIGCOALE_ADDONDATA_HOST
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


    type,public::ImplantSection

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

        real(kind=KMCDF)::DepthINI = 0.D0

        real(kind=KMCDF)::DepthSDINI = 0.D0

        real(kind=KMCDF),dimension(:),allocatable::LayerThick

        real(kind=KMCDF),dimension(:,:),allocatable::ClustersSampleRate

        type(ACluster),dimension(:,:),allocatable::ClustersSample

        contains
        procedure,non_overridable,public,pass::Load_ImplantSection
        procedure,non_overridable,public,pass::ReadImplantSection
        procedure,non_overridable,public,pass::ReadImplantSection_Simple
        procedure,non_overridable,public,pass::ReadImplantClusterSizeDist_Simple
        procedure,non_overridable,public,pass::ReadImplantClusterDepthDist_Simple
        procedure,non_overridable,public,pass::ReadImplantSection_SpecialDistFromFile
        procedure,non_overridable,public,pass::Putin_PANDA_OUTCFG_Distribution
        procedure,non_overridable,public,pass::Putin_SRIM2003_OUTCFG_Distribution
        procedure,non_overridable,public,pass::Putin_OKMC_FORMAT18_Distribution
        procedure,non_overridable,public,pass::ReadImplantSection_SpecialDistFromExteFunc
        procedure,non_overridable,public,pass::Cal_ImplantClustersRate
        procedure,non_overridable,private,pass::Cal_ImplantClustersRate_Simple
        procedure,non_overridable,private,pass::Cal_ImplantClustersRate_Depth_LAY
        procedure,non_overridable,private,pass::Cal_ImplantClustersRate_Depth_Gauss
        procedure,non_overridable,private,pass::Cal_ImplantClustersRate_FromFile
        procedure,non_overridable,private,pass::Cal_ImplantClustersRate_FromExteFunc
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
        procedure,non_overridable,private,pass::Check_ImplantList
        procedure,non_overridable,public,pass::AppendOneSection=>AppendOne_ImplantSection
        procedure,non_overridable,public,pass::Get_P=>GetImplantSection_P
        procedure,non_overridable,public,pass::Clean=>Clean_ImplantList
        Final::CleanImplantList
    end type

    private::Load_ImplantSection
    private::ReadImplantSection
    private::ReadImplantSection_Simple
    private::ReadImplantClusterSizeDist_Simple
    private::ReadImplantClusterDepthDist_Simple
    private::ReadImplantSection_SpecialDistFromFile
    private::Putin_PANDA_OUTCFG_Distribution
    private::Putin_SRIM2003_OUTCFG_Distribution
    private::Putin_OKMC_FORMAT18_Distribution
    private::ReadImplantSection_SpecialDistFromExteFunc
    private::Cal_ImplantClustersRate
    private::Cal_ImplantClustersRate_Simple
    private::Cal_ImplantClustersRate_Depth_LAY
    private::Cal_ImplantClustersRate_Depth_Gauss
    private::Cal_ImplantClustersRate_FromFile
    private::Cal_ImplantClustersRate_FromExteFunc
    private::CopyImplantSectionFromOther
    private::Clean_ImplantSection
    private::CleanImplantSection
    private::Init_ImplantList
    private::Load_ImplantList
    private::Check_ImplantList
    private::AppendOne_ImplantSection
    private::GetImplantSection_P
    private::Clean_ImplantList
    private::CleanImplantList

    contains


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
            write(*,*) "MFPSCUERROR: Unknown file header: ",KEYWORD
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
                    write(*,*) "MFPSCUERROR: Unknown flag: ",KEYWORD
                    write(*,*) "At Line: ",LINE
                    pause
                    stop
            end select

        END DO

        call this%Check_ImplantList(Host_Boxes,Host_SimuCtrlParam)

        return

        100 write(*,*) "MFPSCUERROR: Fail to read the file: ",truePath
            write(*,*) "At Line: ",LINE
            pause
            stop
    end subroutine Load_ImplantList

    !**********************************************************************
    subroutine Check_ImplantList(this,SimBoxes,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantList),target::this
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Local Vars---
        type(ImplantList),pointer::cursor=>null()
        integer::ISection
        !---Body---
        cursor=>this

        ISection = 0

        DO While(associated(cursor))
            ISection = ISection + 1

            if(allocated(cursor%TheImplantSection%LayerThick)) then
                if(sum(cursor%TheImplantSection%LayerThick) .GT. SimBoxes%BOXSIZE(3)) then
                    write(*,*) "MFPSCUERROR: The implant depth is greater than the simulation box."
                    write(*,*) "In section: ",ISection
                    write(*,*) sum(cursor%TheImplantSection%LayerThick),SimBoxes%BOXSIZE(3)
                    pause
                    stop
                end if
            end if

            cursor=>cursor%next
        END DO
    end subroutine

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
            write(*,*) "MFPSCUERROR: You should allocate the ImplantList first!"
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
            write(*,*) "MFPSCUERROR: Cannot find the Implantation section by the id: ",TheIndex
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
                        write(*,*) "MFPSCUERROR: Too few parameters for implantation distribution type."
                        write(*,*) "At Line :", LINE
                        write(*,*) "You should special by the way : &TYPE The implantation cluster distribution type =  "
                        pause
                        stop
                    end if
                    this%ImplantConfigType = ISTR(STRTMP(1))
                    exit
                case default
                    write(*,*) "MFPSCUERROR: You must special the implantation distribution type first!"
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
                        write(*,*) "MFPSCUERROR: Too few parameters for the implantation flux ."
                        write(*,*) "At Line :", LINE
                        write(*,*) "You should special by the way: &FLUX THE implantation flux = , the reflect ratio = "
                        pause
                        stop
                    end if
                    this%ImplantFlux = DRSTR(STRTMP(1))
                    ReflectRatio = DRSTR(STRTMP(2))

                    if(ReflectRatio .LT. 0 .or. ReflectRatio .GT. 1.D0) then
                        write(*,*) "MFPSCUERROR: The reflect ratio should between 0 and 1"
                        write(*,*) "However, the current reflect ratio is: ",ReflectRatio
                        pause
                        stop
                    end if
                    this%ImplantFlux = this%ImplantFlux*(1.D0-ReflectRatio)

                case("&SIZESUBCTL","&DEPTHSUBCTL","&EXTFSUBCTL")
                    call this%ReadImplantSection(hFile,KEYWORD,SimBoxes,Host_SimuCtrlParam,LINE)

                case default
                    write(*,*) "MFPSCUERROR: Unknown Flag: ",KEYWORD
                    write(*,*) "At LINE: ",LINE
                    pause
                    stop
            end select
        End Do

        return

        100 write(*,*) "MFPSCUERROR : Load implantation configuration file failed !"
            write(*,*) "At line :",LINE
            write(*,*) "The program would stop."
            pause
            stop
    end subroutine

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
                write(*,*) "MFPSCUERROR: Unknown strategy for the implantation configuration:",this%ImplantConfigType
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
                write(*,*) "MFPSCUERROR: The Illegal flag: ",KEYWORD
                pause
                stop
        end select

        return

        100 write(*,*) "MFPSCUERROR : Load implantation configuration file failed !"
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
            write(*,*) "MFPSCUERROR: You must special the &EXTFSUBCTL when the implant strategy is chosen by outer file ."
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
                        write(*,*) "MFPSCUERROR: You must special the implantation configuration file type if you had chosen the file model."
                        write(*,*) "By the way: &DISTFILETYPE The distribution file type = ! 'DISTOKMC18','BOXMF18','BOXSPMF18','DISTSRIM','DISTPANDA' "
                        pause
                        stop
                    end if

                    if(LENTRIM(STRTEMP(1)) .LE. 0) then
                        write(*,*) "MFPSCUERROR: The implant configuration file type is null."
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if

                    call UPCASE(STRTEMP(1))

                    this%ImplantCfgFileType = adjustl(trim(KEYWORD_HEAD))//adjustl(trim(STRTEMP(1)))

                    if(IsStrEqual(trim(this%ImplantCfgFileType),SRIM_DIST) .or.  IsStrEqual(trim(this%ImplantCfgFileType),PANDA_DIST)) then
                        call EXTRACT_SUBSTR(STR,2,N,STRTEMP)

                        if(N .LT. 2) then
                            write(*,*) "MFPSCUERROR: when the specialized distribution type is 'DISTSRIM' or 'DISTPANDA'"
                            write(*,*) "MFPSCUERROR: you must special the implantation ion type."
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
                            write(*,*) "MFPSCUERROR: when the specialized distribution type is 'DISTSRIM' or 'DISTOKMC18'"
                            write(*,*) "MFPSCUERROR: you must special the layer number that you want to divide."
                            pause
                            stop
                        end if

                        LayerNum = ISTR(STRTEMP(1))

                        if(LayerNum .LE. 0) then
                            write(*,*) "MFPSCUERROR: the total layer number cannot be less than 0 when it is set for SRIM or OKMC18 distribution"
                            pause
                            stop
                        end if
                    end if


                case("&DISTFILE")
                    call EXTRACT_SUBSTR(STR,1,N,STRTEMP)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCUERROR: You must special the implantation configuration file if you had chosen the file model."
                        write(*,*) "By the way: &DISTFILE The distribution file path = "
                        pause
                        stop
                    end if

                    if(LENTRIM(STRTEMP(1)) .LE. 0) then
                        write(*,*) "MFPSCUERROR: The implant configuration file name is null."
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if

                    this%ImplantCfgFileName = INQUIREFILE(STRTEMP(1),Host_SimuCtrlParam%InputFilePath)

                case default
                    write(*,*) "MFPSCUERROR: Illegal flag: ",KEYWORD
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
                write(*,*) "MFPSCUERROR: Unknown Implant Configuration file type : ",this%ImplantCfgFileType
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
            write(*,*) "MFPSCUERROR: The total concentrate cannot less equal with 0"
            pause
            stop
        end if
        this%ClustersSampleRate = this%ClustersSampleRate/TotalSampleRate

        return
        100 write(*,*) "MFPSCUERROR : Load implantation configuration file failed !"
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
            write(*,*) "MFPSCUERROR: The layers number in panda distribution is less than 1, that is impossible."
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
                write(*,*) "MFPSCUERROR: The panda distribution file data cannot be recognized in line: ",LINE
                write(*,*) STR
                write(*,*) "At file: ",this%ImplantCfgFileName
                pause
                stop
            end if

            this%LayerThick(ILayer) = 2*(DRSTR(STRTMP(1))*C_UM2CM - SumOfThick)
            SumOfThick = SumOfThick + this%LayerThick(ILayer)

            if(SumOfThick .GT. SimBoxes%BOXSIZE(3)) then
                write(*,*) "MFPSCUERROR: The PANDA depth distribution is greater than simulation box depth."
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
            write(*,*) "MFPSCUERROR: The bins number in SRIM2003 distribution is less than 1, that is impossible."
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
                    write(*,*) "MFPSCUERROR: The SRIM2003 distribution file data cannot be recognized in line: ",LINE
                    write(*,*) STR
                    write(*,*) "At file: ",this%ImplantCfgFileName
                    pause
                    stop
                end if

                StoppedPosition(IIon,1) = DRSTR(STRTMP(2))*C_AM2CM ! depth   X
                StoppedPosition(IIon,2) = DRSTR(STRTMP(3))*C_AM2CM ! lateral Y
                StoppedPosition(IIon,3) = DRSTR(STRTMP(4))*C_AM2CM ! lateral Z

                if(StoppedPosition(IIon,2) .LT. SimBoxes%BOXBOUNDARY(1,1) .or. StoppedPosition(IIon,2) .GT. SimBoxes%BOXBOUNDARY(1,2)) then
                    write(*,*) "MFPSCUERROR: The SRIM2003 distribution is out of the simulation box in lateral X."
                    write(*,*) STR
                    write(*,*) "Current position in lateral X is (cm) : ",StoppedPosition(IIon,2)
                    write(*,*) "However, the box boundary range from ",SimBoxes%BOXBOUNDARY(1,1)," To ",SimBoxes%BOXBOUNDARY(1,2)
                    pause
                    stop
                end if

                if(StoppedPosition(IIon,3) .LT. SimBoxes%BOXBOUNDARY(2,1) .or. StoppedPosition(IIon,3) .GT. SimBoxes%BOXBOUNDARY(2,2)) then
                    write(*,*) "MFPSCUERROR: The SRIM2003 distribution is out of the simulation box in lateral Y."
                    write(*,*) STR
                    write(*,*) "Current position in lateral Y is (cm) : ",StoppedPosition(IIon,3)
                    write(*,*) "However, the box boundary range from ",SimBoxes%BOXBOUNDARY(2,1)," To ",SimBoxes%BOXBOUNDARY(2,2)
                    pause
                    stop
                end if

                if(StoppedPosition(IIon,1) .LT. SimBoxes%BOXBOUNDARY(3,1) .or. StoppedPosition(IIon,1) .GT. SimBoxes%BOXBOUNDARY(3,2)) then
                    write(*,*) "MFPSCUERROR: The SRIM2003 distribution is out of the simulation box in depth."
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
                        write(*,*) "MFPSCUERROR: Too few parameters for the cluster size distribution"
                        write(*,*) "You should special: &NATOMDIST The atoms number in each cluster distribution as Gauss that central = , distribution half width = , left cut = ,right cut = "
                        pause
                        stop
                    end if
                    this%NAINI = DRSTR(STRTMP(1))
                    this%NASDINI  = DRSTR(STRTMP(2))
                    this%NACUT(1) = DRSTR(STRTMP(3))
                    this%NACUT(2) = DRSTR(STRTMP(4))
                    if(this%NACUT(1) .GE. this%NACUT(2)) then
                        write(*,*) "MFPSCUERROR: The right cut cannot less than left cut."
                        write(*,*) "LCut",this%NACUT(1)
                        write(*,*) "RCut",this%NACUT(2)
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
                    else
                        if(NElements .GT. 1) then
                            write(*,*) "MFPSCU Info: Currently, the program only support one element."
                            write(*,*) "So the first element would be chosen"
                        end if

                        DO I = 1,NElements
                            this%Elemets(I) = adjustl(trim(STRTMP(I)))
                            call UPCASE(this%Elemets(I))
                        END DO
                    end if
                    ! Current, we only support one element
                    this%CompositWeight = 0
                    TheIndex = SimBoxes%Atoms_list%FindIndexBySymbol(this%Elemets(1))
                    this%CompositWeight(TheIndex) = 1.D0
                    this%CompositWeight = this%CompositWeight/sum(this%CompositWeight)
                CASE default
                    write(*,*) "MFPSCUERROR: Illegal Symbol: ", KEYWORD
                    pause
                    stop
            END SELECT

        END DO

        return

        100 write(*,*) "MFPSCUERROR : Load implantation configuration file failed for cluster size!"
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
                        write(*,*) "MFPSCUERROR: Too few parameters for the clusters depth distribution layer type"
                        write(*,*) "You should special: &DEPTH_LAYER THE NUMBER OF DEPTH DISTRIBUTION LAYER = , THE ENTRIES DISTRIBUTION ="
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if

                    LayerNum = ISTR(STRTMP(1))
                    if(LayerNum .LT. 1) then
                        write(*,*) "MFPSCUERROR: The layer number should greater than 1"
                        write(*,*) "At line :",LINE
                        write(*,*) STR
                        pause
                        stop
                    end if

                    call EXTRACT_NUMB(STR,LayerNum*2+1,N,STRTMP)

                    if((N-1) .NE. LayerNum*2) then
                        write(*,*) "MFPSCUERROR: the specialeld layer is not equal with your setting"
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
                        write(*,*) "MFPSCUERROR: The total concentrate cannot less equal with 0"
                        pause
                        stop
                    end if
                    this%ClustersSampleRate = this%ClustersSampleRate/TotalSampleRate

                CASE("&DEPTH_GAUSS")

                    this%ImplantDepthDistType = p_DEPT_DIS_GAS

                    call EXTRACT_NUMB(STR,2,N,STRTMP)
                    if(N .LT. 2) then
                        write(*,*) "MFPSCUERROR: Too few parameters for the clusters depth distribution gauss type"
                        write(*,*) "You shoud special: &DEPTH_GAUSS THE GAUSS DISTRIBUTION CENTRAL = , THE HALF WIDTH = "
                        write(*,*) "At line: ",LINE
                        pause
                        stop
                    end if
                    this%DepthINI = DRSTR(STRTMP(1))*C_NM2CM
                    this%DepthSDINI = DRSTR(STRTMP(2))*C_NM2CM
                CASE default
                    write(*,*) "MFPSCUERROR: Illegal Symbol: ", KEYWORD
                    pause
                    stop
            END SELECT

        END DO

        return

        100 write(*,*) "MFPSCUERROR : Load Implantation configuration file failed for clusters depth distribution!"
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

        this%DepthINI = other%DepthINI

        this%DepthSDINI = other%DepthSDINI

        call DeAllocateArray_Host(this%LayerThick,"LayerThick")
        if(allocated(other%LayerThick)) then
            call AllocateArray_Host(this%LayerThick,size(other%LayerThick),"LayerThick")
            this%LayerThick = other%LayerThick
        end if

        !---The assignment(=) had been override
        call DeAllocateArray_Host(this%ClustersSample,"ClustersSample")
        if(allocated(other%ClustersSample)) then
            call AllocateArray_Host(this%ClustersSample,size(other%ClustersSample,dim=1),size(other%ClustersSample,dim=2),"ClustersSample")
            this%ClustersSample = other%ClustersSample
        end if

        call DeAllocateArray_Host(this%ClustersSampleRate,"ClustersSampleRate")
        if(allocated(other%ClustersSampleRate)) then
            call AllocateArray_Host(this%ClustersSampleRate,size(other%ClustersSampleRate,dim=1),size(other%ClustersSampleRate,dim=2),"ClustersSampleRate")
            this%ClustersSampleRate = other%ClustersSampleRate
        end if

        return
    end subroutine CopyImplantSectionFromOther

    !*********************************************************************
    subroutine Clean_ImplantSection(this)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        !---Body---

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

        this%DepthINI = 0.D0

        this%DepthSDINI = 0.D0

        call DeAllocateArray_Host(this%LayerThick,"LayerThick")

        call DeAllocateArray_Host(this%ClustersSample,"ClustersSample")

        call DeAllocateArray_Host(this%ClustersSampleRate,"ClustersSampleRate")

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
    subroutine Cal_ImplantClustersRate(this,Host_Boxes,Host_SimuCtrlParam,TheMigCoaleStatInfoWrap,Record,ImplantRate)
        use RAND32_MODULE
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        type(MigCoaleStatInfoWrap)::TheMigCoaleStatInfoWrap
        type(MigCoalClusterRecord)::Record
        real(kind=KMCDF),dimension(:,:),allocatable::ImplantRate
        !---Local Vars---

        select case(this%ImplantConfigType)
            case(p_ImplantConfig_Simple)
                call this%Cal_ImplantClustersRate_Simple(Host_Boxes,Host_SimuCtrlParam,ImplantRate)
            case(p_ImplantConfig_SpecialDistFromFile)
                call this%Cal_ImplantClustersRate_FromFile(Host_Boxes,Host_SimuCtrlParam,ImplantRate)
            case(p_ImplantConfig_SpecialDistFromExteFunc)
                call this%Cal_ImplantClustersRate_FromExteFunc(Host_Boxes,Host_SimuCtrlParam,ImplantRate)
            case default
                write(*,*) "MFPSCUERROR: Unknown strategy for the implantation configuration:",this%ImplantConfigType
                pause
                stop
        end select


        return
    end subroutine Cal_ImplantClustersRate

    !*************************************************************
    subroutine Cal_ImplantClustersRate_Simple(this,SimBoxes,Host_SimuCtrlParam,ImplantRate)
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::SimBoxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        real(kind=KMCDF),dimension(:,:),allocatable::ImplantRate
        !---Body--
        select case(this%ImplantDepthDistType)
            case(p_DEPT_DIS_Layer)
                call this%Cal_ImplantClustersRate_Depth_LAY(SimBoxes,Host_SimuCtrlParam,ImplantRate)
            case(p_DEPT_DIS_GAS)
                call this%Cal_ImplantClustersRate_Depth_Gauss(SimBoxes,Host_SimuCtrlParam,ImplantRate)
            case default
                write(*,*) "MFPSCUERROR : Unknown way to Unknown strategy for the simple implantation configuration: ",this%ImplantDepthDistType
                pause
                stop
        end select

        return
    end subroutine Cal_ImplantClustersRate_Simple

    !*************************************************************
    subroutine Cal_ImplantClustersRate_FromFile(this,Host_Boxes,Host_SimuCtrlParam,ImplantRate)
        use RAND32_MODULE
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        real(kind=KMCDF),dimension(:,:),allocatable::ImplantRate
        !-----local variables---
        integer::MaxGroups
        integer::INode
        integer::ILayer
        integer::JLayer
        integer::LayerStart
        integer::IGroup
        integer::IElement
        integer::LayerNum
        integer::ClusterIndex
        real(kind=KMCDF)::AccumThick
        real(kind=KMCDF)::AccumThickTemp
        real(kind=KMCDF)::AccumThick_Implant
        integer::OccupyLayer_Start
        integer::OccupyLayer_End
        real(kind=KMCDF)::TotalInterval
        real(kind=KMCDF)::Percent
        !---Body---
        AccumThick = 0.D0

        AccumThick_Implant = 0.D0

        ImplantRate = 0.D0

        if(this%ImplantFlux .LE. 0.D0) then
            return
        end if

        LayerStart = 1

        LayerNum = size(this%LayerThick)

        MaxGroups = size(this%ClustersSampleRate,dim=2)

        DO ILayer = 1,LayerNum

            AccumThick_Implant = AccumThick_Implant + this%LayerThick(ILayer)

            DO JLayer = LayerStart,Host_Boxes%NNodes
                AccumThick = AccumThick + Host_Boxes%NodeSpace(JLayer)

                OccupyLayer_End = JLayer

                LayerStart = JLayer + 1

                if(AccumThick .GE. AccumThick_Implant) then
                    exit
                end if
            END DO

            AccumThickTemp = AccumThick

            DO JLayer = OccupyLayer_End,1,-1

                OccupyLayer_Start = JLayer

                AccumThickTemp = AccumThickTemp - Host_Boxes%NodeSpace(JLayer)

                if(AccumThickTemp .LE. AccumThick_Implant) then
                    exit
                end if
            END DO

            AccumThickTemp = AccumThick

            DO JLayer = OccupyLayer_End,OccupyLayer_Start,-1

                if(JLayer .eq. OccupyLayer_Start) then
                    Percent = (AccumThickTemp - (AccumThick_Implant - this%LayerThick(ILayer)) )/this%LayerThick(ILayer)

                    Percent = min(Percent,1.D0)
                else if(JLayer .eq. OccupyLayer_End) then
                    AccumThickTemp = AccumThickTemp - Host_Boxes%NodeSpace(JLayer)
                    Percent = (AccumThick_Implant - AccumThickTemp)/this%LayerThick(ILayer)

                    Percent = min(Percent,1.D0)
                else
                    AccumThickTemp = AccumThickTemp - Host_Boxes%NodeSpace(JLayer)
                    Percent = Host_Boxes%NodeSpace(JLayer)/this%LayerThick(ILayer)
                end if

                DO IGroup = 1,MaxGroups

                    ClusterIndex = Host_Boxes%m_ClustersInfo_CPU%GetKindIndexByCluster(this%ClustersSample(ILayer,IGroup))

                    ImplantRate(JLayer,ClusterIndex) = ImplantRate(JLayer,ClusterIndex) + Percent*this%ImplantFlux*this%ClustersSampleRate(ILayer,IGroup)
                END DO

            END DO

        END DO

        DO INode = 1,Host_Boxes%NNodes
           ImplantRate(INode,:) = ImplantRate(INode,:)/Host_Boxes%NodeSpace(INode)
        END DO

        return
    end subroutine Cal_ImplantClustersRate_FromFile

    !*************************************************************
    subroutine Cal_ImplantClustersRate_FromExteFunc(this,Host_Boxes,Host_SimuCtrlParam,ImplantRate)
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        real(kind=KMCDF),dimension(:,:),allocatable::ImplantRate
        !---Body---

    end subroutine Cal_ImplantClustersRate_FromExteFunc


    !**************************************************************
    subroutine Cal_ImplantClustersRate_Depth_LAY(this,Host_Boxes,Host_SimuCtrlParam,ImplantRate)
      !*** Purpose: To initialize the system (clusters distributed as the form of layer)
      ! Host_Boxes: the boxes information in host
      use RAND32_MODULE
      implicit none
      !---Dummy Vars---
      CLASS(ImplantSection)::this
      type(SimulationBoxes)::Host_Boxes
      type(SimulationCtrlParam)::Host_SimuCtrlParam
      real(kind=KMCDF),dimension(:,:),allocatable::ImplantRate
      !-----local variables---
      integer::NAtoms
      integer::MaxGroups
      integer::INode
      integer::ILayer
      integer::JLayer
      integer::LayerStart
      integer::IGroup
      integer::IElement
      integer::LayerNum
      integer::ClusterIndex
      real(kind=KMCDF)::AccumThick
      real(kind=KMCDF)::AccumThickTemp
      real(kind=KMCDF)::AccumThick_Implant
      integer::OccupyLayer_Start
      integer::OccupyLayer_End
      type(ACluster)::TheCluster
      real(kind=KMCDF)::TotalInterval
      real(kind=KMCDF)::Percent
      !---Body---

      AccumThick = 0.D0

      AccumThick_Implant = 0.D0

      ImplantRate = 0.D0

      if(this%ImplantFlux .LE. 0.D0) then
        return
      end if

      LayerStart = 1

      LayerNum = size(this%LayerThick)

      MaxGroups = size(this%ClustersSampleRate,dim=2)

      DO ILayer = 1,LayerNum

         AccumThick_Implant = AccumThick_Implant + this%LayerThick(ILayer)

         write(*,*) "AccumThick_Implant",AccumThick_Implant

         DO JLayer = LayerStart,Host_Boxes%NNodes
            AccumThick = AccumThick + Host_Boxes%NodeSpace(JLayer)

            OccupyLayer_End = JLayer

            LayerStart = JLayer + 1

            if(AccumThick .GE. AccumThick_Implant) then
                exit
            end if
         END DO

         AccumThickTemp = AccumThick

         DO JLayer = OccupyLayer_End,1,-1

             OccupyLayer_Start = JLayer

             AccumThickTemp = AccumThickTemp - Host_Boxes%NodeSpace(JLayer)

             if(AccumThickTemp .LE. (AccumThick_Implant - this%LayerThick(ILayer))) then
                exit
             end if
         END DO

         AccumThickTemp = AccumThick

         DO JLayer = OccupyLayer_End,OccupyLayer_Start,-1
            if(JLayer .eq. OccupyLayer_End) then
                AccumThickTemp = AccumThickTemp - Host_Boxes%NodeSpace(JLayer)
                Percent = (AccumThick_Implant - AccumThickTemp)/Host_Boxes%NodeSpace(JLayer)
            else if(JLayer .eq. OccupyLayer_Start) then
                Percent = (AccumThickTemp - (AccumThick_Implant - this%LayerThick(ILayer)) )/Host_Boxes%NodeSpace(JLayer)
            else
                AccumThickTemp = AccumThickTemp - Host_Boxes%NodeSpace(JLayer)
                Percent = 1.D0
            end if

            DO IGroup = 1,MaxGroups
                !*** Initialize the size of the clusters
                NAtoms = RGAUSS0_WithCut(this%NAINI, this%NASDINI,this%NACUT(1),this%NACUT(2))

                DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                    TheCluster%m_Atoms(IElement)%m_NA = floor(NAtoms*this%CompositWeight(IElement)+0.5D0)
                    TheCluster%m_Atoms(IElement)%m_ID = IElement
                END DO

                ClusterIndex = Host_Boxes%m_ClustersInfo_CPU%GetKindIndexByCluster(TheCluster)

                ImplantRate(JLayer,ClusterIndex) = ImplantRate(JLayer,ClusterIndex) + Percent*this%ImplantFlux*this%ClustersSampleRate(ILayer,IGroup)
            END DO

         END DO

      END DO

      DO INode = 1,Host_Boxes%NNodes
         ImplantRate(INode,:) = ImplantRate(INode,:)/Host_Boxes%NodeSpace(INode)
      END DO

      return
    end subroutine Cal_ImplantClustersRate_Depth_LAY


    !**************************************************************
    subroutine Cal_ImplantClustersRate_Depth_Gauss(this,Host_Boxes,Host_SimuCtrlParam,ImplantRate)
        !*** Purpose: To initialize the system (clusters distributed as the form of gauss in depth)
        ! Host_Boxes: the boxes information in host
        use RAND32_MODULE
        implicit none
        !---Dummy Vars---
        CLASS(ImplantSection)::this
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        real(kind=KMCDF),dimension(:,:),allocatable::ImplantRate
        !-----local variables---
        integer::NAtoms
        integer::MaxGroups
        integer::IGroup
        integer::IElement
        integer::ClusterIndex
        type(ACluster)::TheCluster
        real(kind=KMCDF)::Depth
        integer::SampleTimeEachLayer
        integer::ISample
        integer::INode
        real(kind=KMCDF)::AccumThick
        integer::TrueLayer
        !---Body---

        SampleTimeEachLayer = 100

        ImplantRate = 0.D0

        if(this%ImplantFlux .LE. 0.D0) then
            return
        end if

        MaxGroups = 1

        DO ISample = 1,SampleTimeEachLayer*Host_Boxes%NNodes

            AccumThick = 0.D0

            Depth = RGAUSS0_WithCut(this%DepthINI, this%DepthSDINI,Host_Boxes%BOXBOUNDARY(3,1),Host_Boxes%BOXBOUNDARY(3,2))

            DO INode = 1,Host_Boxes%NNodes
                AccumThick = AccumThick + Host_Boxes%NodeSpace(INode)

                if(AccumThick .GE. Depth) then
                    TrueLayer = INode
                    exit
                end if
            END DO

            DO IGroup = 1,MaxGroups
                !*** Initialize the size of the clusters
                NAtoms = RGAUSS0_WithCut(this%NAINI, this%NASDINI,this%NACUT(1),this%NACUT(2))

                DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                    TheCluster%m_Atoms(IElement)%m_NA = floor(NAtoms*this%CompositWeight(IElement)+0.5D0)
                    TheCluster%m_Atoms(IElement)%m_ID = IElement
                END DO

                ClusterIndex = Host_Boxes%m_ClustersInfo_CPU%GetKindIndexByCluster(TheCluster)

                ImplantRate(TrueLayer,ClusterIndex) = ImplantRate(TrueLayer,ClusterIndex) + this%ImplantFlux/(SampleTimeEachLayer*MaxGroups)
            END DO

        END DO

        DO INode = 1,Host_Boxes%NNodes
           ImplantRate(INode,:) = ImplantRate(INode,:)/Host_Boxes%NodeSpace(INode)
        END DO

        return
    end subroutine Cal_ImplantClustersRate_Depth_Gauss
!
!
!    !*************************************************************
!    subroutine DoImplantTillVirtualBoundary_GPUTOCPU(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ImplantSection)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        integer,intent(in)::NewAllocateNCEachBox
!        !---Local Vars---
!        integer::MultiBox
!        integer::NSIZE
!        !---Body---
!
!        MultiBox = Host_SimuCtrlParam%MultiBox
!
!        call this%DoImplantTillVirtualBoundary_GPU(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            NSIZE = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,1) + 1
!        else
!            NSIZE = 0
!        end if
!
!        call Host_Boxes%Clean()
!
!        call Dev_Boxes%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,NSIZE,IfCpyNL=.false.)
!
!        return
!    end subroutine DoImplantTillVirtualBoundary_GPUTOCPU
!
!
!    !*************************************************************
!    subroutine DoImplantTillVirtualBoundary_GPU(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ImplantSection)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        integer,intent(in)::NewAllocateNCEachBox
!        !---Body---
!
!        call this%InitImplantInfo_DevPart()
!
!        select case(this%ImplantConfigType)
!            case(p_ImplantConfig_Simple)
!                call this%FillVirtualBoundary_GPU_Simple(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!            case(p_ImplantConfig_SpecialDistFromFile)
!                call this%FillVirtualBoundary_GPU_FromFile(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!            case(p_ImplantConfig_SpecialDistFromExteFunc)
!                call this%FillVirtualBoundary_GPU_FromExteFunc(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!            case default
!                write(*,*) "MFPSCUERROR: Unknown strategy for the implantation configuration:",this%ImplantConfigType
!                pause
!                stop
!        end select
!
!        return
!    end subroutine DoImplantTillVirtualBoundary_GPU
!
!
!    !*************************************************************
!    subroutine FillVirtualBoundary_GPU_Simple(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ImplantSection)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        integer,intent(in)::NewAllocateNCEachBox
!        !---Body---
!        select case(this%ImplantDepthDistType)
!            case(p_DEPT_DIS_Layer)
!                call this%FillVirtualBoundary_GPU_Depth_LAY(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!            case(p_DEPT_DIS_BOX)
!                call this%FillVirtualBoundary_GPU_Depth_SubBox(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!            case(p_DEPT_DIS_GAS)
!                call this%FillVirtualBoundary_GPU_Depth_Gauss(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!            case default
!                write(*,*) "MFPSCUERROR : Unknown way to Unknown strategy for the simple implantation configuration: ",this%ImplantDepthDistType
!                pause
!                stop
!        end select
!
!        return
!    end subroutine FillVirtualBoundary_GPU_Simple
!
!    !*************************************************************
!    subroutine FillVirtualBoundary_GPU_Depth_LAY(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ImplantSection)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        integer,intent(in)::NewAllocateNCEachBox
!        !---Local Vars---
!        integer::MultiBox
!        integer::TotalAllocateNC
!        type(dim3)::blocks
!        type(dim3)::threads
!        integer::NB
!        integer::NBX,NBY
!        Integer::BX
!        integer::BY
!        integer::err
!        !---Body---
!
!        ASSOCIATE(ImplantRand=>Dev_MigCoaleGVars%dm_MigCoale_RandDev)
!
!            if(NewAllocateNCEachBox .GT. 0) then
!
!                MultiBox = Host_SimuCtrlParam%MultiBox
!
!                TotalAllocateNC = MultiBox*NewAllocateNCEachBox
!
!                NB = (TotalAllocateNC - 1)/p_BLOCKSIZE + 1
!                NBX  = min(NB,p_BLOCKDIMX)
!                NBY = (NB - 1)/NBX + 1
!
!                !*** to determine the block size
!                BX = p_BLOCKSIZE
!                BY = 1
!                !*** to determine the dimension of blocks
!
!                blocks  = dim3(NBX, NBY, 1)
!                threads = dim3(BX,  BY,  1)
!
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Layer,ImplantRand%dm_SpaceDist_Implant(1:TotalAllocateNC),TotalAllocateNC)
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_X,ImplantRand%dm_SpaceDist_Implant(TotalAllocateNC+1:2*TotalAllocateNC),TotalAllocateNC)
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Y,ImplantRand%dm_SpaceDist_Implant(2*TotalAllocateNC+1:3*TotalAllocateNC),TotalAllocateNC)
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Z,ImplantRand%dm_SpaceDist_Implant(3*TotalAllocateNC+1:4*TotalAllocateNC),TotalAllocateNC)
!
!                err = curandGenerateNormal(ImplantRand%m_ranGen_ClustersSizeDist,ImplantRand%dm_SizeDist_Implant,TotalAllocateNC,this%NAINI,this%NASDINI)
!
!                call Kernel_ImplantClusters_Depth_Layer<<<blocks,threads>>>(TotalAllocateNC,                                          &
!                                                                            NewAllocateNCEachBox,                                     &
!                                                                            Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,                 &
!                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_TypesEntities,          &
!                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_SingleAtomsDivideArrays,&
!                                                                            Host_Boxes%m_GrainBoundary%GrainNum,                      &
!                                                                            Dev_Boxes%dm_GrainBoundary%dm_GrainSeeds,                 &
!                                                                            ImplantRand%dm_SpaceDist_Implant,                         &
!                                                                            ImplantRand%dm_SizeDist_Implant,                          &
!                                                                            this%NACUT(1),                                            &
!                                                                            this%NACUT(2),                                            &
!                                                                            Dev_Boxes%dm_SEVirtualIndexBox,                           &
!                                                                            this%dm_ImplantInfo_DevPart%Dev_LayerThick,               &
!                                                                            this%dm_ImplantInfo_DevPart%Dev_ClustersSampleRate,       &
!                                                                            this%dm_ImplantInfo_DevPart%Dev_CompositWeight)
!            end if
!
!        END ASSOCIATE
!
!        return
!    end subroutine FillVirtualBoundary_GPU_Depth_LAY
!
!    !**********************************************
!    attributes(global) subroutine Kernel_ImplantClusters_Depth_Layer(TotalAllocateNC,            &
!                                                                    NewAllocateNCEachBox,        &
!                                                                    Dev_Clusters,                &
!                                                                    Dev_TypesEntities,           &
!                                                                    Dev_SingleAtomsDivideArrays, &
!                                                                    Nseeds,                      &
!                                                                    Dev_GrainSeeds,              &
!                                                                    Dev_RandArray_SpaceDist,     &
!                                                                    Dev_RandArray_SizeDist,      &
!                                                                    LNACUT,                      &
!                                                                    RNACUT,                      &
!                                                                    Dev_SEVirtualIndexBox,       &
!                                                                    Dev_LayerThick,              &
!                                                                    Dev_ClustersSampleRate,      &
!                                                                    Dev_CompositWeight)
!        implicit none
!        !---Dummy Vars---
!        integer, value::TotalAllocateNC
!        integer, value::NewAllocateNCEachBox
!        type(ACluster), device::Dev_Clusters(:)
!        type(DiffusorTypeEntity),device::Dev_TypesEntities(:)
!        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! If the two dimension array would be delivered to attributes(device), the first dimension must be known
!        integer,value::Nseeds
!        type(GrainSeed),device::Dev_GrainSeeds(:)
!        real(kind=KMCDF),device::Dev_RandArray_SpaceDist(:)
!        real(kind=KMCDF),device::Dev_RandArray_SizeDist(:)
!        real(kind=KMCDF),value::LNACUT
!        real(kind=KMCDF),value::RNACUT
!        integer, device::Dev_SEVirtualIndexBox(:,:)
!        integer, device::Dev_SEAddedClustersBoxes(:,:)
!        real(kind=KMCDF),device::Dev_LayerThick(:)
!        real(kind=KMCDF),device::Dev_ClustersSampleRate(:,:)
!        real(kind=KMCDF),device::Dev_CompositWeight(:)
!        !---Local Vars---
!        integer::tid
!        integer::bid
!        integer::cid
!        integer::IBox
!        integer::cid0
!        integer::ICTRUE
!        real(kind=KMCDF)::POS(3)
!        integer::NLayer
!        integer::MaxGroups
!        integer::ILayer
!        integer::IGroup
!        logical::exitFlag
!        real(kind=KMCDF)::tempRand
!        real(kind=KMCDF)::GroupRateTemp
!        integer::IElement
!        type(DiffusorValue)::TheDiffusorValue
!        real(kind=KMCDF)::randSize
!        !---Body---
!        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
!        bid = (blockidx%y - 1)*griddim%x + blockidx%x
!        cid = (bid - 1)*blockdim%x*blockdim%y + tid
!
!        IBox = (cid - 1)/NewAllocateNCEachBox + 1
!        cid0 = (IBox - 1)*NewAllocateNCEachBox + 1
!
!        if(cid .LE. TotalAllocateNC) then
!            NLayer = size(Dev_ClustersSampleRate,dim=1)
!            MaxGroups = size(Dev_ClustersSampleRate,dim=2)
!
!            tempRand = Dev_RandArray_SpaceDist(cid)
!
!            GroupRateTemp = 0.D0
!            exitFlag = .false.
!            DO ILayer = 1,NLayer
!
!                if(exitFlag .eq. .true.) then
!                    exit
!                end if
!
!                DO IGroup = 1,MaxGroups
!                    GroupRateTemp = GroupRateTemp + Dev_ClustersSampleRate(ILayer,IGroup)
!                    if(GroupRateTemp .GE. tempRand) then
!
!                        ICTRUE = Dev_SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1 + (cid - cid0)
!
!                        !Initialize the position of clusters
!                        POS(1) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC)*dm_BOXSIZE(1)+dm_BOXBOUNDARY(1,1)
!                        POS(2) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*2)*dm_BOXSIZE(2)+dm_BOXBOUNDARY(2,1)
!                        POS(3) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*3)*Dev_LayerThick(ILayer) + sum(Dev_LayerThick(1:ILayer-1),dim=1) + dm_BOXBOUNDARY(3,1)
!                        Dev_Clusters(ICTRUE)%m_POS = POS
!
!                        !Give the cluster an type(layer) ID for the convenience of visualization
!                        Dev_Clusters(ICTRUE)%m_Layer = ILayer
!
!                        !*** Initialize the size of the clusters
!                        randSize = Dev_RandArray_SizeDist(cid)
!                        if(randSize .GT. RNACUT) then
!                            randSize = 2*RNACUT - randSize - (RNACUT-LNACUT)*floor((randSize-RNACUT)/(RNACUT-LNACUT))
!                        else if(randSize .LT. LNACUT) then
!                            randSize = 2*LNACUT - randSize - (RNACUT-LNACUT)*floor((LNACUT - randSize)/(RNACUT-LNACUT))
!                        end if
!                        DO IElement = 1,p_ATOMS_GROUPS_NUMBER
!                            Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_NA = floor(randSize*Dev_CompositWeight(IElement)+0.5D0)
!                            Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_ID = IElement
!                        END DO
!
!                        Dev_Clusters(ICTRUE)%m_GrainID(1) = GrainBelongsTo_Dev(Nseeds,Dev_GrainSeeds,POS)
!
!                        Dev_Clusters(ICTRUE)%m_Statu = p_ACTIVEFREE_STATU
!
!                        call Dev_GetValueFromDiffusorsMap(Dev_Clusters(ICTRUE),Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheDiffusorValue)
!
!                        !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
!                        !---you should init the distribution by external file---
!                        select case(TheDiffusorValue%ECRValueType_Free)
!                            case(p_ECR_ByValue)
!                                Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_Free
!                            case(p_ECR_ByBCluster)
!                                Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
!                        end select
!
!                        select case(TheDiffusorValue%DiffusorValueType_Free)
!                            case(p_DiffuseCoefficient_ByValue)
!                                Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
!                            case(p_DiffuseCoefficient_ByArrhenius)
!                                Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/dm_TKB)
!                            case(p_DiffuseCoefficient_ByBCluster)
!                                ! Here we adopt a model that D=D0*(1/R)**Gama
!                                Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_FREESURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
!                        end select
!
!                        exitFlag = .true.
!                        exit
!
!
!                    end if
!                END DO
!            END DO
!
!        end if
!
!        return
!    end subroutine Kernel_ImplantClusters_Depth_Layer
!
!    !*************************************************************
!    subroutine FillVirtualBoundary_GPU_Depth_SubBox(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ImplantSection)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        integer,intent(in)::NewAllocateNCEachBox
!        !---Local Vars---
!        integer::MultiBox
!        integer::TotalAllocateNC
!        type(dim3)::blocks
!        type(dim3)::threads
!        integer::NB
!        integer::NBX,NBY
!        Integer::BX
!        integer::BY
!        integer::err
!        !---Body---
!
!        ASSOCIATE(ImplantRand=>Dev_MigCoaleGVars%dm_MigCoale_RandDev)
!
!            if(NewAllocateNCEachBox .GT. 0) then
!
!                MultiBox = Host_SimuCtrlParam%MultiBox
!
!                TotalAllocateNC = MultiBox*NewAllocateNCEachBox
!
!                NB = (TotalAllocateNC - 1)/p_BLOCKSIZE + 1
!                NBX  = min(NB,p_BLOCKDIMX)
!                NBY = (NB - 1)/NBX + 1
!
!                !*** to determine the block size
!                BX = p_BLOCKSIZE
!                BY = 1
!                !*** to determine the dimension of blocks
!
!                blocks  = dim3(NBX, NBY, 1)
!                threads = dim3(BX,  BY,  1)
!
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_X,ImplantRand%dm_SpaceDist_Implant(1:TotalAllocateNC),TotalAllocateNC)
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Y,ImplantRand%dm_SpaceDist_Implant(TotalAllocateNC+1:2*TotalAllocateNC),TotalAllocateNC)
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Z,ImplantRand%dm_SpaceDist_Implant(2*TotalAllocateNC+1:3*TotalAllocateNC),TotalAllocateNC)
!
!                err = curandGenerateNormal(ImplantRand%m_ranGen_ClustersSizeDist,ImplantRand%dm_SizeDist_Implant,TotalAllocateNC,this%NAINI,this%NASDINI)
!
!                call Kernel_ImplantClusters_Depth_SubBox<<<blocks,threads>>>(TotalAllocateNC,                                         &
!                                                                            NewAllocateNCEachBox,                                     &
!                                                                            Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,                 &
!                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_TypesEntities,          &
!                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_SingleAtomsDivideArrays,&
!                                                                            Host_Boxes%m_GrainBoundary%GrainNum,                      &
!                                                                            Dev_Boxes%dm_GrainBoundary%dm_GrainSeeds,                 &
!                                                                            ImplantRand%dm_SpaceDist_Implant,                         &
!                                                                            ImplantRand%dm_SizeDist_Implant,                          &
!                                                                            this%NACUT(1),                                            &
!                                                                            this%NACUT(2),                                            &
!                                                                            Dev_Boxes%dm_SEVirtualIndexBox,                           &
!                                                                            this%dm_ImplantInfo_DevPart%Dev_SUBBOXBOUNDARY,           &
!                                                                            this%dm_ImplantInfo_DevPart%Dev_CompositWeight)
!            end if
!
!        END ASSOCIATE
!
!
!        return
!    end subroutine FillVirtualBoundary_GPU_Depth_SubBox
!
!    !**********************************************
!    attributes(global) subroutine Kernel_ImplantClusters_Depth_SubBox(TotalAllocateNC,           &
!                                                                    NewAllocateNCEachBox,        &
!                                                                    Dev_Clusters,                &
!                                                                    Dev_TypesEntities,           &
!                                                                    Dev_SingleAtomsDivideArrays, &
!                                                                    Nseeds,                      &
!                                                                    Dev_GrainSeeds,              &
!                                                                    Dev_RandArray_SpaceDist,     &
!                                                                    Dev_RandArray_SizeDist,      &
!                                                                    LNACUT,                      &
!                                                                    RNACUT,                      &
!                                                                    Dev_SEVirtualIndexBox,       &
!                                                                    Dev_SUBBOXBOUNDARY,          &
!                                                                    Dev_CompositWeight)
!        implicit none
!        !---Dummy Vars---
!        integer, value::TotalAllocateNC
!        integer, value::NewAllocateNCEachBox
!        type(ACluster), device::Dev_Clusters(:)
!        type(DiffusorTypeEntity),device::Dev_TypesEntities(:)
!        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! If the two dimension array would be delivered to attributes(device), the first dimension must be known
!        integer,value::Nseeds
!        type(GrainSeed),device::Dev_GrainSeeds(:)
!        real(kind=KMCDF),device::Dev_RandArray_SpaceDist(:)
!        real(kind=KMCDF),device::Dev_RandArray_SizeDist(:)
!        real(kind=KMCDF),value::LNACUT
!        real(kind=KMCDF),value::RNACUT
!        integer, device::Dev_SEVirtualIndexBox(:,:)
!        real(kind=KMCDF),device::Dev_SUBBOXBOUNDARY(:,:)
!        real(kind=KMCDF),device::Dev_CompositWeight(:)
!        !---Local Vars---
!        integer::tid
!        integer::bid
!        integer::cid
!        integer::IBox
!        integer::cid0
!        integer::ICTRUE
!        integer::I
!        real(kind=KMCDF)::POS(3)
!        integer::IElement
!        type(DiffusorValue)::TheDiffusorValue
!        real(kind=KMCDF)::randSize
!        !---Body---
!        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
!        bid = (blockidx%y - 1)*griddim%x + blockidx%x
!        cid = (bid - 1)*blockdim%x*blockdim%y + tid
!
!        IBox = (cid - 1)/NewAllocateNCEachBox + 1
!        cid0 = (IBox - 1)*NewAllocateNCEachBox + 1
!
!        if(cid .LE. TotalAllocateNC) then
!            ICTRUE = Dev_SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1 + (cid - cid0)
!
!            DO I = 1,3
!                POS(I) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*(I - 1))*(Dev_SUBBOXBOUNDARY(I,2) - Dev_SUBBOXBOUNDARY(I,1)) + Dev_SUBBOXBOUNDARY(I,1)
!            END DO
!            !Initialize the position of clusters
!            Dev_Clusters(ICTRUE)%m_POS = POS
!
!            !Give the cluster an type(layer) ID for the convenience of visualization
!            Dev_Clusters(ICTRUE)%m_Layer = 1
!
!            !*** Initialize the size of the clusters
!            randSize = Dev_RandArray_SizeDist(cid)
!            if(randSize .GT. RNACUT) then
!                randSize = 2*RNACUT - randSize - (RNACUT-LNACUT)*floor((randSize-RNACUT)/(RNACUT-LNACUT))
!            else if(randSize .LT. LNACUT) then
!                randSize = 2*LNACUT - randSize - (RNACUT-LNACUT)*floor((LNACUT - randSize)/(RNACUT-LNACUT))
!            end if
!            DO IElement = 1,p_ATOMS_GROUPS_NUMBER
!                Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_NA = floor(randSize*Dev_CompositWeight(IElement)+0.5D0)
!                Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_ID = IElement
!            END DO
!
!            Dev_Clusters(ICTRUE)%m_GrainID(1) = GrainBelongsTo_Dev(Nseeds,Dev_GrainSeeds,POS)
!
!            Dev_Clusters(ICTRUE)%m_Statu = p_ACTIVEFREE_STATU
!
!            call Dev_GetValueFromDiffusorsMap(Dev_Clusters(ICTRUE),Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheDiffusorValue)
!
!            !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
!            !---you should init the distribution by external file---
!            select case(TheDiffusorValue%ECRValueType_Free)
!                case(p_ECR_ByValue)
!                    Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_Free
!                case(p_ECR_ByBCluster)
!                    Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
!            end select
!
!            select case(TheDiffusorValue%DiffusorValueType_Free)
!                case(p_DiffuseCoefficient_ByValue)
!                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
!                case(p_DiffuseCoefficient_ByArrhenius)
!                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/dm_TKB)
!                case(p_DiffuseCoefficient_ByBCluster)
!                    ! Here we adopt a model that D=D0*(1/R)**Gama
!                    Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_FREESURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
!            end select
!
!        end if
!
!        return
!    end subroutine Kernel_ImplantClusters_Depth_SubBox
!
!    !*************************************************************
!    subroutine FillVirtualBoundary_GPU_Depth_Gauss(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ImplantSection)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        integer,intent(in)::NewAllocateNCEachBox
!        !---Local Vars---
!        integer::MultiBox
!        integer::TotalAllocateNC
!        type(dim3)::blocks
!        type(dim3)::threads
!        integer::NB
!        integer::NBX,NBY
!        Integer::BX
!        integer::BY
!        integer::err
!        !---Body---
!        ASSOCIATE(ImplantRand=>Dev_MigCoaleGVars%dm_MigCoale_RandDev)
!
!            if(NewAllocateNCEachBox .GT. 0) then
!
!                MultiBox = Host_SimuCtrlParam%MultiBox
!
!                TotalAllocateNC = MultiBox*NewAllocateNCEachBox
!
!                NB = (TotalAllocateNC - 1)/p_BLOCKSIZE + 1
!                NBX  = min(NB,p_BLOCKDIMX)
!                NBY = (NB - 1)/NBX + 1
!
!                !*** to determine the block size
!                BX = p_BLOCKSIZE
!                BY = 1
!                !*** to determine the dimension of blocks
!
!                blocks  = dim3(NBX, NBY, 1)
!                threads = dim3(BX,  BY,  1)
!
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_X,ImplantRand%dm_SpaceDist_Implant(1:TotalAllocateNC),TotalAllocateNC)
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Y,ImplantRand%dm_SpaceDist_Implant(TotalAllocateNC+1:2*TotalAllocateNC),TotalAllocateNC)
!                err = curandGenerateNormal(ImplantRand%m_ranGen_ClustersSpaceDist_Z,ImplantRand%dm_SpaceDist_Implant(2*TotalAllocateNC+1:3*TotalAllocateNC),TotalAllocateNC,this%DepthINI,this%DepthSDINI)
!
!                err = curandGenerateNormal(ImplantRand%m_ranGen_ClustersSizeDist,ImplantRand%dm_SizeDist_Implant,TotalAllocateNC,this%NAINI,this%NASDINI)
!
!                call Kernel_ImplantClusters_Depth_Gauss<<<blocks,threads>>>(TotalAllocateNC,                                          &
!                                                                            NewAllocateNCEachBox,                                     &
!                                                                            Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,                 &
!                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_TypesEntities,          &
!                                                                            Dev_Boxes%dm_DiffusorTypesMap%Dev_SingleAtomsDivideArrays,&
!                                                                            Host_Boxes%m_GrainBoundary%GrainNum,                      &
!                                                                            Dev_Boxes%dm_GrainBoundary%dm_GrainSeeds,                 &
!                                                                            ImplantRand%dm_SpaceDist_Implant,                         &
!                                                                            ImplantRand%dm_SizeDist_Implant,                          &
!                                                                            this%NACUT(1),                                            &
!                                                                            this%NACUT(2),                                            &
!                                                                            Dev_Boxes%dm_SEVirtualIndexBox,                           &
!                                                                            this%dm_ImplantInfo_DevPart%Dev_CompositWeight)
!            end if
!
!        END ASSOCIATE
!
!        return
!    end subroutine FillVirtualBoundary_GPU_Depth_Gauss
!
!    !**********************************************
!    attributes(global) subroutine Kernel_ImplantClusters_Depth_Gauss(TotalAllocateNC,            &
!                                                                    NewAllocateNCEachBox,        &
!                                                                    Dev_Clusters,                &
!                                                                    Dev_TypesEntities,           &
!                                                                    Dev_SingleAtomsDivideArrays, &
!                                                                    Nseeds,                      &
!                                                                    Dev_GrainSeeds,              &
!                                                                    Dev_RandArray_SpaceDist,     &
!                                                                    Dev_RandArray_SizeDist,      &
!                                                                    LNACUT,                      &
!                                                                    RNACUT,                      &
!                                                                    Dev_SEVirtualIndexBox,       &
!                                                                    Dev_CompositWeight)
!        implicit none
!        !---Dummy Vars---
!        integer, value::TotalAllocateNC
!        integer, value::NewAllocateNCEachBox
!        type(ACluster), device::Dev_Clusters(:)
!        type(DiffusorTypeEntity),device::Dev_TypesEntities(:)
!        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! If the two dimension array would be delivered to attributes(device), the first dimension must be known
!        integer,value::Nseeds
!        type(GrainSeed),device::Dev_GrainSeeds(:)
!        real(kind=KMCDF),device::Dev_RandArray_SpaceDist(:)
!        real(kind=KMCDF),device::Dev_RandArray_SizeDist(:)
!        real(kind=KMCDF),value::LNACUT
!        real(kind=KMCDF),value::RNACUT
!        integer, device::Dev_SEVirtualIndexBox(:,:)
!        real(kind=KMCDF),device::Dev_CompositWeight(:)
!        !---Local Vars---
!        integer::tid
!        integer::bid
!        integer::cid
!        integer::IBox
!        integer::cid0
!        integer::ICTRUE
!        integer::I
!        real(kind=KMCDF)::POS(3)
!        integer::IElement
!        type(DiffusorValue)::TheDiffusorValue
!        real(kind=KMCDF)::randSize
!        real(kind=KMCDF)::randDepth
!        !---Body---
!        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
!        bid = (blockidx%y - 1)*griddim%x + blockidx%x
!        cid = (bid - 1)*blockdim%x*blockdim%y + tid
!
!        IBox = (cid - 1)/NewAllocateNCEachBox + 1
!        cid0 = (IBox - 1)*NewAllocateNCEachBox + 1
!
!        if(cid .LE. TotalAllocateNC) then
!            ICTRUE = Dev_SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1 + (cid - cid0)
!
!            POS(1) = Dev_RandArray_SpaceDist(cid)*dm_BOXSIZE(1)+dm_BOXBOUNDARY(1,1)
!            POS(2) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC)*dm_BOXSIZE(2)+dm_BOXBOUNDARY(2,1)
!
!            randDepth = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*2)
!            if(randDepth .GT. dm_BOXBOUNDARY(3,2)) then
!                randDepth = 2*dm_BOXBOUNDARY(3,2) - randDepth - (dm_BOXBOUNDARY(3,2)-dm_BOXBOUNDARY(3,1))*floor((randDepth-dm_BOXBOUNDARY(3,2))/(dm_BOXBOUNDARY(3,2)-dm_BOXBOUNDARY(3,1)))
!            else if(randDepth .LT. dm_BOXBOUNDARY(3,1)) then
!                randDepth = 2*dm_BOXBOUNDARY(3,1) - randDepth - (dm_BOXBOUNDARY(3,2)-dm_BOXBOUNDARY(3,1))*floor((dm_BOXBOUNDARY(3,1)-randDepth)/(dm_BOXBOUNDARY(3,2)-dm_BOXBOUNDARY(3,1)))
!            end if
!
!            POS(3) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*2)*dm_BOXSIZE(3) + dm_BOXBOUNDARY(3,1)
!            !Initialize the position of clusters
!            Dev_Clusters(ICTRUE)%m_POS = POS
!
!            !Give the cluster an type(layer) ID for the convenience of visualization
!            Dev_Clusters(ICTRUE)%m_Layer = 1
!
!            !*** Initialize the size of the clusters
!            randSize = Dev_RandArray_SizeDist(cid)
!            if(randSize .GT. RNACUT) then
!                randSize = 2*RNACUT - randSize - (RNACUT-LNACUT)*floor((randSize-RNACUT)/(RNACUT-LNACUT))
!            else if(randSize .LT. LNACUT) then
!                randSize = 2*LNACUT - randSize - (RNACUT-LNACUT)*floor((LNACUT - randSize)/(RNACUT-LNACUT))
!            end if
!            DO IElement = 1,p_ATOMS_GROUPS_NUMBER
!                Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_NA = floor(randSize*Dev_CompositWeight(IElement)+0.5D0)
!                Dev_Clusters(ICTRUE)%m_Atoms(IElement)%m_ID = IElement
!            END DO
!
!            Dev_Clusters(ICTRUE)%m_GrainID(1) = GrainBelongsTo_Dev(Nseeds,Dev_GrainSeeds,POS)
!
!            Dev_Clusters(ICTRUE)%m_Statu = p_ACTIVEFREE_STATU
!
!            call Dev_GetValueFromDiffusorsMap(Dev_Clusters(ICTRUE),Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheDiffusorValue)
!
!            !-- In Current application, the simple implant distribution is only considered in free matrix, if you want to init the clusters in GB---
!            !---you should init the distribution by external file---
!            select case(TheDiffusorValue%ECRValueType_Free)
!                case(p_ECR_ByValue)
!                    Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_Free
!                case(p_ECR_ByBCluster)
!                    Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
!            end select
!
!            select case(TheDiffusorValue%DiffusorValueType_Free)
!                case(p_DiffuseCoefficient_ByValue)
!                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
!                case(p_DiffuseCoefficient_ByArrhenius)
!                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/dm_TKB)
!                case(p_DiffuseCoefficient_ByBCluster)
!                    ! Here we adopt a model that D=D0*(1/R)**Gama
!                    Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_FREESURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
!            end select
!
!        end if
!
!        return
!    end subroutine Kernel_ImplantClusters_Depth_Gauss
!
!    !*************************************************************
!    subroutine FillVirtualBoundary_GPU_FromFile(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ImplantSection)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        integer,intent(in)::NewAllocateNCEachBox
!        !---Local Vars---
!        integer::MultiBox
!        integer::TotalAllocateNC
!        type(dim3)::blocks
!        type(dim3)::threads
!        integer::NB
!        integer::NBX,NBY
!        Integer::BX
!        integer::BY
!        integer::err
!        !---Body---
!        ASSOCIATE(ImplantRand=>Dev_MigCoaleGVars%dm_MigCoale_RandDev)
!
!            if(NewAllocateNCEachBox .GT. 0) then
!
!                MultiBox = Host_SimuCtrlParam%MultiBox
!
!                TotalAllocateNC = MultiBox*NewAllocateNCEachBox
!
!                NB = (TotalAllocateNC - 1)/p_BLOCKSIZE + 1
!                NBX  = min(NB,p_BLOCKDIMX)
!                NBY = (NB - 1)/NBX + 1
!
!                !*** to determine the block size
!                BX = p_BLOCKSIZE
!                BY = 1
!                !*** to determine the dimension of blocks
!
!                blocks  = dim3(NBX, NBY, 1)
!                threads = dim3(BX,  BY,  1)
!
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_Layer,ImplantRand%dm_SpaceDist_Implant(1:TotalAllocateNC),TotalAllocateNC)
!                err = curandGenerateUniformDouble(ImplantRand%m_ranGen_ClustersSpaceDist_X,ImplantRand%dm_SpaceDist_Implant(TotalAllocateNC+1:2*TotalAllocateNC),TotalAllocateNC)
!
!                call Kernel_ImplantClusters_FromFile<<<blocks,threads>>>(TotalAllocateNC,                                         &
!                                                                        NewAllocateNCEachBox,                                     &
!                                                                        Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,                 &
!                                                                        Dev_Boxes%dm_DiffusorTypesMap%Dev_TypesEntities,          &
!                                                                        Dev_Boxes%dm_DiffusorTypesMap%Dev_SingleAtomsDivideArrays,&
!                                                                        Host_Boxes%m_GrainBoundary%GrainNum,                      &
!                                                                        Dev_Boxes%dm_GrainBoundary%dm_GrainSeeds,                 &
!                                                                        ImplantRand%dm_SpaceDist_Implant,                         &
!                                                                        Dev_Boxes%dm_SEVirtualIndexBox,                           &
!                                                                        this%dm_ImplantInfo_DevPart%Dev_LayerThick,               &
!                                                                        this%dm_ImplantInfo_DevPart%Dev_ClustersSample,           &
!                                                                        this%dm_ImplantInfo_DevPart%Dev_ClustersSampleRate)
!            end if
!
!        END ASSOCIATE
!
!        return
!    end subroutine FillVirtualBoundary_GPU_FromFile
!
!
!    !**********************************************
!    attributes(global) subroutine Kernel_ImplantClusters_FromFile(TotalAllocateNC,             &
!                                                                  NewAllocateNCEachBox,        &
!                                                                  Dev_Clusters,                &
!                                                                  Dev_TypesEntities,                &
!                                                                  Dev_SingleAtomsDivideArrays, &
!                                                                  Nseeds,                      &
!                                                                  Dev_GrainSeeds,              &
!                                                                  Dev_RandArray_SpaceDist,     &
!                                                                  Dev_SEVirtualIndexBox,       &
!                                                                  Dev_LayerThick,              &
!                                                                  Dev_ClustersSample,          &
!                                                                  Dev_ClustersSampleRate)
!        implicit none
!        !---Dummy Vars---
!        integer, value::TotalAllocateNC
!        integer, value::NewAllocateNCEachBox
!        type(ACluster), device::Dev_Clusters(:)
!        type(DiffusorTypeEntity),device::Dev_TypesEntities(:)
!        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! If the two dimension array would be delivered to attributes(device), the first dimension must be known
!        integer,value::Nseeds
!        type(GrainSeed),device::Dev_GrainSeeds(:)
!        real(kind=KMCDF),device::Dev_RandArray_SpaceDist(:)
!        integer, device::Dev_SEVirtualIndexBox(:,:)
!        real(kind=KMCDF),device::Dev_LayerThick(:)
!        type(ACluster),device::Dev_ClustersSample(:,:)
!        real(kind=KMCDF),device::Dev_ClustersSampleRate(:,:)
!        !---Local Vars---
!        integer::tid
!        integer::bid
!        integer::cid
!        integer::IBox
!        integer::cid0
!        integer::ICTRUE
!        real(kind=KMCDF)::POS(3)
!        integer::NLayer
!        integer::MaxGroups
!        integer::ILayer
!        integer::IGroup
!        logical::exitFlag
!        real(kind=KMCDF)::tempRand
!        real(kind=KMCDF)::GroupRateTemp
!        type(DiffusorValue)::TheDiffusorValue
!        !---Body---
!        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
!        bid = (blockidx%y - 1)*griddim%x + blockidx%x
!        cid = (bid - 1)*blockdim%x*blockdim%y + tid
!
!        IBox = (cid - 1)/NewAllocateNCEachBox + 1
!        cid0 = (IBox - 1)*NewAllocateNCEachBox + 1
!
!        if(cid .LE. TotalAllocateNC) then
!            NLayer = size(Dev_ClustersSampleRate,dim=1)
!            MaxGroups = size(Dev_ClustersSampleRate,dim=2)
!
!            tempRand = Dev_RandArray_SpaceDist(cid)
!
!            GroupRateTemp = 0.D0
!            exitFlag = .false.
!            DO ILayer = 1,NLayer
!
!                if(exitFlag .eq. .true.) then
!                    exit
!                end if
!
!                DO IGroup = 1,MaxGroups
!                    GroupRateTemp = GroupRateTemp + Dev_ClustersSampleRate(ILayer,IGroup)
!                    if(GroupRateTemp .GE. tempRand) then
!
!                        ICTRUE = Dev_SEVirtualIndexBox(IBox,2) - NewAllocateNCEachBox + 1 + (cid - cid0)
!
!                        Dev_Clusters(ICTRUE)%m_Atoms = Dev_ClustersSample(ILayer,IGroup)%m_Atoms
!
!                        !Initialize the position of clusters
!                        POS(1) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC)*dm_BOXSIZE(1)+dm_BOXBOUNDARY(1,1)
!                        POS(2) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*2)*dm_BOXSIZE(2)+dm_BOXBOUNDARY(2,1)
!                        POS(3) = Dev_RandArray_SpaceDist(cid + TotalAllocateNC*3)*Dev_LayerThick(ILayer) + sum(Dev_LayerThick(1:ILayer-1),dim=1) + dm_BOXBOUNDARY(3,1)
!                        Dev_Clusters(ICTRUE)%m_POS = POS
!
!                        !Give the cluster an type(layer) ID for the convenience of visualization
!                        Dev_Clusters(ICTRUE)%m_Layer = ILayer
!
!                        Dev_Clusters(ICTRUE)%m_Statu = Dev_ClustersSample(ILayer,IGroup)%m_Statu
!
!                        call Dev_GetValueFromDiffusorsMap(Dev_Clusters(ICTRUE),Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheDiffusorValue)
!
!                        if(Dev_Clusters(ICTRUE)%m_Statu .eq. p_ACTIVEFREE_STATU) then
!
!                            Dev_Clusters(ICTRUE)%m_GrainID(1) = GrainBelongsTo_Dev(Nseeds,Dev_GrainSeeds,POS)
!
!                            select case(TheDiffusorValue%ECRValueType_Free)
!                                case(p_ECR_ByValue)
!                                    Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_Free
!                                case(p_ECR_ByBCluster)
!                                    Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
!                            end select
!
!                            select case(TheDiffusorValue%DiffusorValueType_Free)
!                                case(p_DiffuseCoefficient_ByValue)
!                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_Free_Value
!                                case(p_DiffuseCoefficient_ByArrhenius)
!                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_Free*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_Free/dm_TKB)
!                                case(p_DiffuseCoefficient_ByBCluster)
!                                    ! Here we adopt a model that D=D0*(1/R)**Gama
!                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_FREESURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
!                            end select
!
!                        else if(Dev_Clusters(ICTRUE)%m_Statu .eq. p_ACTIVEINGB_STATU) then
!
!                            Dev_Clusters(ICTRUE)%m_GrainID = Dev_ClustersSample(ILayer,IGroup)%m_GrainID
!
!                            select case(TheDiffusorValue%ECRValueType_InGB)
!                                case(p_ECR_ByValue)
!                                    Dev_Clusters(ICTRUE)%m_RAD = TheDiffusorValue%ECR_InGB
!                                case(p_ECR_ByBCluster)
!                                    Dev_Clusters(ICTRUE)%m_RAD = DSQRT(sum(Dev_Clusters(ICTRUE)%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA,dim=1)/dm_RNFACTOR)
!                            end select
!
!                            select case(TheDiffusorValue%DiffusorValueType_InGB)
!                                case(p_DiffuseCoefficient_ByValue)
!                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%DiffuseCoefficient_InGB_Value
!                                case(p_DiffuseCoefficient_ByArrhenius)
!                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = TheDiffusorValue%PreFactor_InGB*exp(-C_EV2ERG*TheDiffusorValue%ActEnergy_InGB/dm_TKB)
!                                case(p_DiffuseCoefficient_ByBCluster)
!                                    ! Here we adopt a model that D=D0*(1/R)**Gama
!                                    Dev_Clusters(ICTRUE)%m_DiffCoeff = dm_GBSURDIFPRE*(Dev_Clusters(ICTRUE)%m_RAD**(-p_GAMMA))
!                            end select
!                        end if
!
!                        exitFlag = .true.
!                        exit
!
!                    end if
!                END DO
!            END DO
!
!        end if
!
!        return
!    end subroutine Kernel_ImplantClusters_FromFile
!
!    !*************************************************************
!    subroutine FillVirtualBoundary_GPU_FromExteFunc(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,Dev_MigCoaleGVars,NewAllocateNCEachBox)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ImplantSection)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        type(MigCoale_GVarsDev)::Dev_MigCoaleGVars
!        integer,intent(in)::NewAllocateNCEachBox
!        !---Body---
!
!        return
!    end subroutine FillVirtualBoundary_GPU_FromExteFunc


end module MIGCOALE_IMPLANTATION
