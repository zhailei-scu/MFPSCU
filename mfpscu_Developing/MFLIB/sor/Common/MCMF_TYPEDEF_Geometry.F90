module MCMF_TYPEDEF_GEOMETRY
    use MCMF_CONSTANTS
    use MCMF_TYPEDEF_SIMULATIONCTRLPARAM
    implicit none

    integer, parameter::p_GBIniConfig_Simple = 0
    integer, parameter::p_GBIniConfig_SpecialDistFromFile = 1
    integer, parameter::p_GBIniConfig_SpecialDistFromExteFunc = 2

    integer, parameter::p_GBInitSimple_BySeedCtl = 0
    integer, parameter::p_GBInitSimple_ByGVolumCtl = 1

    type,public::GrainSeed
        real(kind=KMCDF),dimension(3)::m_POS

        contains

        procedure,non_overridable,private,pass::CopySeedsFormOther
        Generic::Assignment(=)=>CopySeedsFormOther
    end type

    type,public::GrainSeedList
        type(GrainSeed)::Seed
        type(GrainSeedList),pointer::next=>null()

        integer::m_count = 0

        contains
        procedure,non_overridable,private,pass::AppendOne
        procedure,non_overridable,private,pass::AppendOneSeed
        procedure,non_overridable,private,pass::ConvertToArray
        procedure,non_overridable,private,nopass::CleanGrainSeedList
        procedure,non_overridable,private,pass::CopyGrainSeedLisFromOther
        Generic::ASSIGNMENT(=)=>CopyGrainSeedLisFromOther
        Final::DestroyGrainSeedList
    end type GrainSeedList

    type,public::GrainBoundary
        integer::GBInitType = p_GBIniConfig_Simple
        !---init by simple way---
        integer::GBInitSimple_Strategy = p_GBInitSimple_BySeedCtl
        real(kind=KMCDF)::Cutoff(2) = 0.D0                                  ! min and max distance cut-off
        integer::GrainNum = 0
        real(kind=KMCDF)::SeedsDistINI = 0.D0
        real(kind=KMCDF)::SeedsDistSD = 0.D0
        real(kind=KMCDF)::GVolumINI = 0.D0
        real(kind=KMCDF)::GVolumSD = 0.D0
        !---Init by external file----
        character*256::GBCfgFileName = ''
        !---Init by external function---

        !---Seeds---
        type(GrainSeed),dimension(:),allocatable::GrainSeeds

        contains
        procedure,non_overridable,public,pass::ConstructGrainBoundary
        procedure,non_overridable,public,pass::ConstructGrainBoundary_Simple
        procedure,non_overridable,public,pass::ConstructGrainBoundary_Simple_ByGSeedCtl
        procedure,non_overridable,public,pass::ConstructGrainBoundary_Simple_ByGVolumCtl
        procedure,non_overridable,public,pass::ConstructGrainBoundary_SpecialDistFromFile
        procedure,non_overridable,public,pass::ConstructGrainBoundary_SpecialDistFromExteFunc
        procedure,non_overridable,public,pass::RescaleGrainBoundary
        procedure,non_overridable,public,pass::GrainBelongsTo
        procedure,non_overridable,public,pass::Clean_Grainboundary
        procedure,non_overridable,public,pass::CopyGrainBoundaryFromOther
        Generic::ASSIGNMENT(=)=>CopyGrainBoundaryFromOther
        Final::CleanGrainboundary
    end type GrainBoundary

    private::CopySeedsFormOther
    private::AppendOne
    private::AppendOneSeed
    private::ConvertToArray
    private::CopyGrainSeedLisFromOther
    private::CleanGrainSeedList
    private::DestroyGrainSeedList
    private::ConstructGrainBoundary
    private::ConstructGrainBoundary_Simple
    private::ConstructGrainBoundary_Simple_ByGSeedCtl
    private::ConstructGrainBoundary_Simple_ByGVolumCtl
    private::ConstructGrainBoundary_SpecialDistFromFile
    private::ConstructGrainBoundary_SpecialDistFromExteFunc
    private::RescaleGrainBoundary
    private::GrainBelongsTo
    private::CopyGrainBoundaryFromOther
    private::Clean_Grainboundary
    private::CleanGrainboundary

    contains

    !******************************************
    function Get_MemoryConsuming_GrainSeed() result(TheSize)
        implicit none
        !---Dummy Vars---
        integer::TheSize
        !---Local Vars---
        type(GrainSeed)::OneGrainSeed
        !---Body---
        TheSize = sizeof(OneGrainSeed)
        return
    end function

    !******************************************
    subroutine CopySeedsFormOther(this,theOtherOne)
        implicit none
        !---Dummy Vars---
        CLASS(GrainSeed),intent(out)::this
        type(GrainSeed),intent(in)::theOtherOne
        !---Body---

        this%m_POS = theOtherOne%m_POS

        return
    end subroutine

    !******************************************
    subroutine AppendOne(this,newOne)
        implicit none
        !---Dummy Vars---
        CLASS(GrainSeedList),target::this
        type(GrainSeedList)::newOne
        !---Local Vars---
        type(GrainSeedList),pointer::Pcursor=>null()
        type(GrainSeedList),pointer::Cursor=>null()
        !---Body---

        if(this%m_count .eq. 0) then
            this%Seed = newOne%Seed
            this%next=>null()
        else
            Pcursor=>this
            Cursor=>Pcursor%next

            DO While(associated(Cursor))
                Pcursor=>Pcursor%next
                Cursor=>Pcursor%next
            END DO

            allocate(Cursor)
            Cursor%Seed = newOne%Seed
            Cursor%next=>null()

            Pcursor%next=>Cursor

            Nullify(Pcursor)
            Nullify(Cursor)

        end if

        this%m_count = this%m_count + 1

        return
    end subroutine

    !******************************************
    subroutine AppendOneSeed(this,newSeed)
        implicit none
        !---Dummy Vars---
        CLASS(GrainSeedList)::this
        type(GrainSeed)::newSeed
        !---Local Vars---
        type(GrainSeedList)::tempGrainSeedList
        !---Body---

        tempGrainSeedList%Seed = newSeed
        tempGrainSeedList%next=>null()

        call this%AppendOne(tempGrainSeedList)

        return
    end subroutine AppendOneSeed

    !*****************************************
    subroutine ConvertToArray(this,theArray)
        implicit none
        !---Dummy Vars---
        CLASS(GrainSeedList),target::this
        type(GrainSeed),dimension(:),allocatable::theArray
        !---Local Vars---
        type(GrainSeedList),pointer::Cursor=>null()
        integer::I
        !---Body---
        if(this%m_count .ne. size(theArray)) then
            write(*,*) "MCPSCUERROR: the convert from Grainseed List to Grainseeds array failed!"
            write(*,*) "MCPSCUERROR: becasue the size is not same!"
            write(*,*) "The size of grainseeds list is:",this%m_count
            write(*,*) "The size of grainseeds array is:",size(theArray)
            pause
            stop
        end if

        if(this%m_count .GT. 0) then
            Cursor=>this
            I = 1
            DO While(associated(Cursor))
                theArray(I) = Cursor%Seed
                Cursor=>Cursor%next
                I = I + 1
            END DO
        end if

        return
    end subroutine

    !*****************************************
    subroutine CopyGrainSeedLisFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(GrainSeedList),target,intent(out)::this
        type(GrainSeedList),target,intent(in)::other
        !---Local Vars---
        type(GrainSeedList),pointer::cursorThis=>null()
        type(GrainSeedList),pointer::cursorOther=>null()
        !---Body---
        cursorThis=>this
        call CleanGrainSeedList(cursorThis)

        cursorOther=>other
        if(other%m_count .GT. 0) then
            DO While(associated(cursorOther))
                call this%AppendOneSeed(cursorOther%Seed)

                cursorOther=>cursorOther%next
            END DO
        end if

        Nullify(cursorThis)
        cursorThis=>null()

        Nullify(cursorOther)
        cursorOther=>null()
        return
    end subroutine

    !******************************************
    subroutine CleanGrainSeedList(list)
        implicit none
        !---Dummy Vars---
        type(GrainSeedList),pointer::list
        !---Local Vars---
        type(GrainSeedList),pointer::PCursor=>null()
        type(GrainSeedList),pointer::cursor=>null()
        !---Body---
        if(associated(list)) then

            if(list%m_count .GT. 0) then
                PCursor=>list
                cursor=>list%next
                DO While(associated(cursor))
                    Nullify(PCursor%next)
                    PCursor%next=>null()
                    deallocate(PCursor)
                    PCursor=>cursor
                    cursor=>cursor%next
                END DO

                Nullify(cursor)
                cursor=>null()

                Nullify(PCursor)
                PCursor=>null()

            end if
        end if

        return
    end subroutine

    !******************************************
    subroutine DestroyGrainSeedList(this)
        implicit none
        !---Dummy Vars--
        type(GrainSeedList),target::this
        !---Local Vars---
        CLASS(GrainSeedList),pointer::ptr=>null()
        !---Body---
        ptr=>this

        call CleanGrainSeedList(ptr)

        Nullify(ptr)
        ptr=>null()

        return
    end subroutine


    !******************************************
    subroutine ConstructGrainBoundary(this,BOXBOUNDARY,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary)::this
        real(kind=KMCDF),intent(in)::BOXBOUNDARY(3,2)
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Body---
        select case(this%GBInitType)
            case(p_GBIniConfig_Simple)
                call this%ConstructGrainBoundary_Simple(BOXBOUNDARY)
            case(p_GBIniConfig_SpecialDistFromFile)
                call this%ConstructGrainBoundary_SpecialDistFromFile(Host_SimuCtrlParam)
            case(p_GBIniConfig_SpecialDistFromExteFunc)
                call this%ConstructGrainBoundary_SpecialDistFromExteFunc()
            case default
                write(*,*) "MCPSCUERROR: Unkonw way to construct the grain boundary in the boxes.",this%GBInitType
                pause
                stop
        end select

        return
    end subroutine ConstructGrainBoundary

    !******************************************
    subroutine ConstructGrainBoundary_Simple(this,BOXBOUNDARY)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary)::this
        real(kind=KMCDF),intent(in)::BOXBOUNDARY(3,2)
        !---Body---
        select case(this%GBInitSimple_Strategy)
            case(p_GBInitSimple_BySeedCtl)
                call this%ConstructGrainBoundary_Simple_ByGSeedCtl(BOXBOUNDARY)
            case(p_GBInitSimple_ByGVolumCtl)
                call this%ConstructGrainBoundary_Simple_ByGVolumCtl(BOXBOUNDARY)
            case default
                write(*,*) "MCPSCUERROR: Unkonw way to construct the grain boundary in the boxes.",this%GBInitType
                pause
                stop
        end select

        return
    end subroutine ConstructGrainBoundary_Simple

    !******************************************
    subroutine ConstructGrainBoundary_Simple_ByGSeedCtl(this,BOXBOUNDARY)
        use RAND32_MODULE
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary)::this
        real(kind=KMCDF),intent(in)::BOXBOUNDARY(3,2)
        !---Local Vars---
        real(kind=KMCDF)::BOXSIZE(3)
        real(kind=KMCDF)::BOXVOLUM
        real(kind=KMCDF)::POS(3)
        integer::N
        integer::MaxSeedsNum
        type(GrainSeedList),pointer::TheGrainSeedList=>null()
        type(GrainSeedList),pointer::cursor=>null()
        real(kind=KMCDF)::SEP(3)
        real(kind=KMCDF)::Dist2
        real(kind=KMCDF)::CutMin2
        logical::finded
        type(GrainSeed)::tempSeed
        integer::I
        integer::K
        !---Body---
        N = 0
        DO I = 1,3
            BOXSIZE(I) = BOXBOUNDARY(I,2) - BOXBOUNDARY(I,1)
        END DO

        BOXVOLUM = product(BOXSIZE)

        CutMin2 = this%Cutoff(1)**2

        MaxSeedsNum = ceiling(BOXVOLUM/((SQRT(2.D0)/12.D0)*this%Cutoff(1)**3))

        if(this%GrainNum .LT. 0 .or. this%GrainNum .GT. MaxSeedsNum) then       ! the grains is unsetted
            this%GrainNum = MaxSeedsNum
        end if

        allocate(TheGrainSeedList)

        DO K = 1,this%GrainNum

            finded = .false.

            DO While(finded .eq. .false.)
                cursor=>TheGrainSeedList
                DO I = 1,3
                    POS(I) = DRAND32()*BOXSIZE(I) + BOXBOUNDARY(I,1)
                END DO

                DO While(associated(cursor))

                    SEP = POS - cursor%Seed%m_POS

                    Dist2 = SEP(1)*SEP(1) + SEP(2)*SEP(2) + SEP(2)*SEP(2)

                    if(Dist2 .GT. CutMin2) then
                        finded = .true.
                        exit
                    end if

                    cursor=>cursor%next

                END DO

            END DO

            tempSeed%m_POS = POS

            call TheGrainSeedList%AppendOneSeed(tempSeed)

        END DO

        this%GrainNum = TheGrainSeedList%m_count

        if(this%GrainNum .GT. 0) then
            allocate(this%GrainSeeds(this%GrainNum))
        end if

        call TheGrainSeedList%ConvertToArray(this%GrainSeeds)

        call CleanGrainSeedList(TheGrainSeedList)

        return
    end subroutine ConstructGrainBoundary_Simple_ByGSeedCtl

    !******************************************
    subroutine ConstructGrainBoundary_Simple_ByGVolumCtl(this,BOXBOUNDARY)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary)::this
        real(kind=KMCDF),intent(in)::BOXBOUNDARY(3,2)
        !---Body---

    end subroutine ConstructGrainBoundary_Simple_ByGVolumCtl

    !******************************************
    subroutine ConstructGrainBoundary_SpecialDistFromFile(this,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary)::this
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        !---Local Vars---
        character*256::FilePath
        character*256::STR
        integer::LINE
        integer::hFile
        integer::ISeed
        integer::N
        character*20::STRTMP(10)
        integer::I
        !---Body---

        LINE = 0

        this%GrainNum = 0

        FilePath = INQUIREFILE(this%GBCfgFileName,Host_SimuCtrlParam%InputFilePath)

        hFile = OpenExistedFile(FilePath)



        DO While(.not. GETINPUTSTRLINE_New(hFile,STR,LINE,"!") )
            call RemoveComments(STR,"!")

            if(LENTRIM(adjustl(STR)) .LE. 0) then
                cycle
            end if

            this%GrainNum = this%GrainNum + 1
        END DO

        if(this%GrainNum .LE. 0) then
            write(*,*) "MCPSCUERROR: The grain number is less than 0"
            pause
            stop
        end if

        allocate(this%GrainSeeds(this%GrainNum))


        Rewind(hFile)
        LINE = 0

        ISeed = 0

        DO While(.not. GETINPUTSTRLINE_New(hFile,STR,LINE,"!") )
            call RemoveComments(STR,"!")

            if(LENTRIM(adjustl(STR)) .LE. 0) then
                cycle
            end if

            ISeed = ISeed  +  1

            call EXTRACT_NUMB(STR,3,N,STRTMP)

            if(N .LT. 3) then
                write(*,*) "MCPSCUERROR: The seed should be constructed by x, y and z, however, some dimension is lake at LINE: ",LINE
                write(*,*) STR
                pause
                stop
            end if

            DO I = 1,3
                this%GrainSeeds(ISeed)%m_POS(I) = DRSTR(STRTMP(I))
            END DO

        END DO


        return
    end subroutine ConstructGrainBoundary_SpecialDistFromFile

    !******************************************
    subroutine ConstructGrainBoundary_SpecialDistFromExteFunc(this)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary)::this
        !---Body---
        return
    end subroutine ConstructGrainBoundary_SpecialDistFromExteFunc

    !*******************************************
    subroutine RescaleGrainBoundary(this,DUPXYZ)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary)::this
        integer,intent(in)::DUPXYZ(3)
        !---Local Vars---
        integer::GrainSeedsNum
        type(GrainBoundary)::temp
        integer::dumpNum
        integer::IG
        integer::I,J,K
        integer::IP
        !---Body---
        GrainSeedsNum = this%GrainNum

        if(GrainSeedsNum .GT. 0) then

          dumpNum = product(DUPXYZ+1)

          allocate(temp%GrainSeeds(dumpNum*GrainSeedsNum))

          IP = 1
          DO IG = 1,GrainSeedsNum

            DO I = 0,DUPXYZ(1)
                DO J = 0,DUPXYZ(2)
                    DO K = 0,DUPXYZ(3)
                        ! the assignment(=) has been overrided
                        temp%GrainSeeds(IP) = this%GrainSeeds(IG)
                        IP = IP + 1
                    END DO
                END DO

            END DO

          END DO

          this%GrainNum = dumpNum*GrainSeedsNum

          deallocate(this%GrainSeeds)

          allocate(this%GrainSeeds(dumpNum*GrainSeedsNum))

          this%GrainSeeds = temp%GrainSeeds

          call temp%Clean_Grainboundary()

        end if

        return
    end subroutine

    !**********************************
    subroutine CopyGrainBoundaryFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary),intent(out)::this
        type(GrainBoundary),intent(in)::other
        !---Local Vars---
        !---Body---
        call this%Clean_Grainboundary()

        this%GBInitSimple_Strategy = other%GBInitSimple_Strategy
        this%Cutoff = other%Cutoff
        this%GrainNum = other%GrainNum
        this%SeedsDistINI = other%SeedsDistINI
        this%SeedsDistSD = other%SeedsDistSD
        this%GVolumINI = other%GVolumINI
        this%GVolumSD = other%GVolumSD

        this%GBCfgFileName = other%GBCfgFileName

        if(allocated(this%GrainSeeds)) then
            deallocate(this%GrainSeeds)
        end if
        if(size(other%GrainSeeds) .GT. 0) then
            allocate(this%GrainSeeds(size(other%GrainSeeds)))
        end if
        this%GrainSeeds = other%GrainSeeds

        return
    end subroutine

    !***********************************
    subroutine Clean_Grainboundary(this)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary)::this
        !---Body---
        if(allocated(this%GrainSeeds)) then
            deallocate(this%GrainSeeds)
        end if

        this%GBInitSimple_Strategy = p_GBInitSimple_BySeedCtl
        this%Cutoff = 0.D0
        this%GrainNum = 0
        this%SeedsDistINI = 0.D0
        this%SeedsDistSD = 0.D0
        this%GVolumINI = 0.D0
        this%GVolumSD = 0.D0

        this%GBCfgFileName = ''

        return
    end subroutine

    !***********************************
    subroutine CleanGrainboundary(this)
        !---Dummy Vars---
        type(GrainBoundary)::this
        !---Body---

        call this%Clean_Grainboundary()

        return
    end subroutine


end module MCMF_TYPEDEF_GEOMETRY
