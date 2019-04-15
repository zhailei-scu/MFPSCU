module MFLIB_UTILITIES

  use MiniUtilities,only:LENTRIM
  use MFLIB_CONSTANTS

  implicit none

  interface ResizeArray
    MODULE PROCEDURE ResizeArrayi_OneDim
    MODULE PROCEDURE ResizeArrayr_OneDim
    MODULE PROCEDURE ResizeArrayd_OneDim

    MODULE PROCEDURE ResizeArrayi_TwoDim
    MODULE PROCEDURE ResizeArrayr_TwoDim
    MODULE PROCEDURE ResizeArrayd_TwoDim

  end interface ResizeArray

  !-----------------------
  interface AllocateArray_Host
    MODULE PROCEDURE AllocateOneDimi_Host
    MODULE PROCEDURE AllocateOneDimr_Host
    MODULE PROCEDURE AllocateOneDimd_Host


    MODULE PROCEDURE AllocateTwoDimi_Host
    MODULE PROCEDURE AllocateTwoDimr_Host
    MODULE PROCEDURE AllocateTwoDimd_Host

    MODULE PROCEDURE AllocateThreeDimi_Host
    MODULE PROCEDURE AllocateThreeDimr_Host
    MODULE PROCEDURE AllocateThreeDimd_Host

  end interface AllocateArray_Host

  !------------------
  interface DeAllocateArray_Host
    MODULE PROCEDURE DeAllocateOneDimi_Host
    MODULE PROCEDURE DeAllocateOneDimr_Host
    MODULE PROCEDURE DeAllocateOneDimd_Host

    MODULE PROCEDURE DeAllocateTwoDimi_Host
    MODULE PROCEDURE DeAllocateTwoDimr_Host
    MODULE PROCEDURE DeAllocateTwoDimd_Host

    MODULE PROCEDURE DeAllocateThreeDimi_Host
    MODULE PROCEDURE DeAllocateThreeDimr_Host
    MODULE PROCEDURE DeAllocateThreeDimd_Host

  end interface DeAllocateArray_Host

  contains
  !************************************************************
  function INQUIREFILE(fileName,parentPath) result(truePath)
    implicit none
    !---Dummy Vars---
    character*(*),intent(in)::fileName
    character*(*),intent(in),optional::parentPath
    character*256::truePath
    !---Local Vars---
    logical::exits
    !---Body---
    exits = .false.

    truePath = adjustl(trim(fileName))
    INQUIRE(FILE=truePath,EXIST=exits)

    if(present(parentPath) .AND. .not. exits) then
        if(LENTRIM(adjustl(parentPath)) .LE. 0) then
            truePath = adjustl(trim(fileName))
        else
            truePath = adjustl(trim(parentPath))//adjustl(trim(fileName))
        end if

        INQUIRE(FILE=truePath,EXIST=exits)
    end if

    if(.not. exits) then
        write(*,*) "MCPSCU ERROR: The file do not exit :",fileName
        pause
        stop
    end if

    return
  end function INQUIREFILE

  !************************************************************
  function AvailableIOUnit() result(FileUnit)
    implicit none
    !---Dummy Vars---
    integer,intent(out)::FileUnit
    !---Local Vars---
    logical::opened
    !---Body---
    opened = .true.

    DO FileUnit = 10,99
       INQUIRE(UNIT=FileUnit,OPENED=opened)
       if(.not. opened) exit
    END DO

    return
  end function AvailableIOUnit

  !*********************************************************
  function IsStrEqual(SubjectStr,ObjectStr) result(theResult)
    implicit none
    !---Dummy Vars---
    character*(*),intent(in)::SubjectStr
    character*(*),intent(in)::ObjectStr
    logical,intent(out)::theResult
    !---Local Vars---
    character(len=256)::tempSubjectStr
    character(len=256)::tempObjectStr
    character,dimension(:),allocatable::tempSubjectStrAllo
    character,dimension(:),allocatable::tempObjectStrAllo
    integer::SubjectStrlength
    integer::ObjectStrlength
    logical::SubjectuseArray
    logical::ObjectuseArray
    integer::I
    !---Body---

    theResult = .true.

    SubjectuseArray = .false.

    ObjectuseArray = .false.

    SubjectStrlength = LENTRIM(adjustl(SubjectStr))

    ObjectStrlength = LENTRIM(adjustl(ObjectStr))

    if(SubjectStrlength .ne. ObjectStrlength) then
        theResult = .false.
        return
    end if

    if(SubjectStrlength .GT. len(tempSubjectStr)) then
        allocate(tempSubjectStrAllo(SubjectStrlength))
        tempSubjectStrAllo = trim(adjustl(SubjectStr))
        SubjectuseArray = .true.
    end if

    if(ObjectStrlength .GT. len(tempObjectStr)) then
        allocate(tempObjectStrAllo(ObjectStrlength))
        tempObjectStrAllo = trim(adjustl(ObjectStr))
        ObjectuseArray = .true.
    end if

    tempSubjectStr = adjustl(trim(SubjectStr))
    tempObjectStr = adjustl(trim(ObjectStr))

    if(SubjectuseArray .AND. (.not. ObjectuseArray)) then
        DO I = 1,SubjectStrlength
            if(tempSubjectStrAllo(I) .ne. tempObjectStr(I:I)) then
                theResult = .false.
                exit
            end if
        END DO

    else if((.not. SubjectuseArray) .AND. ObjectuseArray) then
        DO I = 1,SubjectStrlength
            if(tempSubjectStr(I:I) .ne. tempObjectStrAllo(I)) then
                theResult = .false.
                exit
            end if
        END DO

    else if(SubjectuseArray .AND. ObjectuseArray) then
        DO I = 1,SubjectStrlength
            if(tempSubjectStrAllo(I) .ne. tempObjectStrAllo(I)) then
                theResult = .false.
                exit
            end if
        END DO
    else
        if(tempSubjectStr(1:SubjectStrlength) .ne. tempObjectStr(1:ObjectStrlength)) then
            theResult = .false.
        end if
    end if

    return
  end function IsStrEqual

!*********************************************************
  subroutine separateStrByString(STR,separateSTR,SegmentArray,segNumber)
    !   ***Purpose: to sperate the string by the givin seperate character
    !          STR: the string
    !    separateSTR: the seperate string
    ! SegmentArray: the result array
    !    segNumber: the number of segments
    implicit none
    !---Dummy Vars---
    character*(*), intent(in)::STR
    character*(*), intent(in)::separateSTR
    integer::segNumber
    character*(*),dimension(:)::SegmentArray

    !---Local Vars---
    integer::IPosBefore,IPosAfter
    integer::I,LastI
    integer::Length,tempLength
    integer::ArraySize

    !---Body---

    SegmentArray = ""

    ArraySize = size(SegmentArray,DIM=1)

    if(ArraySize .LE. 0) then
        return
    end if

    Length = LENTRIM(adjustl(STR))

    segNumber = 0
    IPosBefore = 1
    IPosAfter = 1
    DO I = 1,Length

      IPosAfter = INDEX(STR(IPosBefore:Length),separateSTR)

      tempLength = IPosAfter - 1

      if(tempLength .GT. 0) then

        if(segNumber .GT. ArraySize) then
            write(*,*) "MCPSCU ERROR: The size of segment array is not enough. Process would be stop.",ArraySize
            write(*,*) STR
            pause
            stop
        end if

        segNumber = segNumber + 1
        SegmentArray(segNumber)(1:tempLength) = STR(IPosBefore:IPosBefore+tempLength - 1)

      else if(tempLength .LT. 0) then

        if(IPosBefore .LE. Length) then
            segNumber = segNumber + 1
            SegmentArray(segNumber)(1:Length-IPosBefore+1) = STR(IPosBefore:Length)
        end if
        exit
      end if

      IPosBefore = IPosBefore + tempLength + 1

    END DO

    return
  end subroutine separateStrByString

  !*************************************************************
  subroutine resolveLongFileName(FileLongName,FilePath,FileShortName)
    !***     Purpose: to resolve the File Long Name(include full Path) to a pure Path and a pure file name
    !   FileLongName: the File Long Name(include full Path)
    !       FilePath: the pure file path
    !  FileShortName: the pure file name
    implicit none
    !---Dummy Vars---
    character*(*), intent(in)::FileLongName
    character*256::FilePath
    character*256::FileShortName
    !---Local Vars---
    character*256::tempLongName
    integer::I,IFind,Length
    character*2::CharacterSet1
    character*2::CharacterSet2
    character*2::CharacterSet3
    !---Body---
    CharacterSet1 = "\\"
    CharacterSet2 = "\/"
    CharacterSet3 = "//"

    tempLongName = adjustl(trim(FileLongName))

    Length = LENTRIM(tempLongName)
    if(Length .LE. 0) then
        return
    end if

    IFind = 0

    IFind = max(SCAN(trim(tempLongName),CharacterSet1,.true.), &
                SCAN(trim(tempLongName),CharacterSet2,.true.), &
                SCAN(trim(tempLongName),CharacterSet3,.true.))

    if(IFind .LE. 0) then
        FilePath = ""
        FileShortName = tempLongName
    else
        FilePath = tempLongName(1:IFind)
        FileShortName = tempLongName(IFind+1:Length)
    end if

    return
  end subroutine resolveLongFileName

  !****************************************
  subroutine resolveExePrefixName(ExeFullName,ExePrefixName)
    implicit none
    !---Dummy Vars---
    character*(*),intent(in)::ExeFullName
    character*(*),intent(inout)::ExePrefixName
    !---Local Vars---
    character(len=128),dimension(10)::seperatedStrsArray
    integer::seperatedNum
    integer::I
    !---Body---

    seperatedStrsArray = ''

    call separateStrByString(ExeFullName,".",seperatedStrsArray,seperatedNum)

    ExePrefixName = ''

    DO I = 1,seperatedNum
        ExePrefixName = ExePrefixName(1:LENTRIM(adjustl(ExePrefixName)))//seperatedStrsArray(I)(1:LENTRIM(adjustl(seperatedStrsArray(I))))
    END DO

    return
  end subroutine resolveExePrefixName

  !*************************************************************
  function OpenExistedFile(fileName,thePosition) result(fileUnit)
    implicit none
    !---Dummy Vars---
    character*(*),intent(in)::fileName
    character*(*),optional,intent(in)::thePosition
    integer,intent(out)::fileUnit
    !---Local Vars---
    integer::ISTAT
    character*256::openInfo
    !---Body---
    fileUnit = AvailableIOUnit()

    if(present(thePosition)) then
        open(Unit=fileUnit,File=fileName,STATUS="old",POSITION=thePosition(1:LENTRIM(thePosition)),iostat=ISTAT,IOMSG=openInfo)
    else
        open(Unit=fileUnit,File=fileName,STATUS="old",iostat=ISTAT,IOMSG=openInfo)
    end if

    if(ISTAT .ne. 0) then
        write(*,*) "MCPSCUERROR: open file failed:",fileName
        write(*,*) openInfo
        pause
        close(fileUnit)
        stop
    end if

    return
  end function OpenExistedFile

  !*************************************************************
  function CreateNewFile(fileName,thePosition) result(fileUnit)
    implicit none
    !---Dummy Vars---
    character*(*),intent(in)::fileName
    character*(*),optional,intent(in)::thePosition
    integer,intent(out)::fileUnit
    !---Local Vars---
    integer::ISTAT
    character*256::openInfo
    !---Body---
    fileUnit = AvailableIOUnit()

    if(present(thePosition)) then
        open(Unit=fileUnit,File=fileName,STATUS="replace",POSITION=thePosition(1:LENTRIM(thePosition)),iostat=ISTAT,IOMSG=openInfo)
    else
        open(Unit=fileUnit,File=fileName,STATUS="replace",iostat=ISTAT,IOMSG=openInfo)
    end if

    if(ISTAT .ne. 0) then
        write(*,*) "MCPSCUERROR: create file failed:",fileName
        write(*,*) openInfo
        pause
        close(fileUnit)
        stop
    end if

    return
  end function CreateNewFile

  !*************************************************************
  function CreateOrOpenExistedFile(fileName,thePosition) result(fileUnit)
    implicit none
    !---Dummy Vars---
    character*(*),intent(in)::fileName
    character*(*),optional,intent(in)::thePosition
    integer,intent(out)::fileUnit
    !---Local Vars---
    integer::ISTAT
    character*256::openInfo
    logical::exits
    !---Body---
    exits = .false.

    INQUIRE(FILE=fileName(1:LENTRIM(adjustl(fileName))),EXIST=exits)

    fileUnit = AvailableIOUnit()

    if(exits) then

        if(present(thePosition)) then
            fileUnit = OpenExistedFile(fileName(1:LENTRIM(fileName)),thePosition)
        else
            fileUnit = OpenExistedFile(fileName(1:LENTRIM(fileName)))
        end if
    else
        if(present(thePosition)) then
            fileUnit = CreateNewFile(fileName(1:LENTRIM(fileName)),thePosition)
        else
            fileUnit = CreateNewFile(fileName(1:LENTRIM(fileName)))
        end if

    end if

    return
  end function CreateOrOpenExistedFile

  !*****************************************************************
  function CreateDataFolder(distPath,parentPath) result(resultPath)
    implicit none
    !---Dummy Vars---
    character*(*)::distPath
    character*(*)::parentPath
    character*256::resultPath
    !---Local Vars---
    character*256::truePath
    logical::exits
    character(len=64),dimension(20)::seperatedStrsArray
    integer::segNumber
    character*256::resolvedParentPath
    integer::I
    logical::AbsolutePath
    !---Body---
    exits = .false.

    AbsolutePath = .false.

    truePath = adjustl(trim(distPath))

    INQUIRE(FILE=truePath,EXIST=exits)

    if(.not. exits) then

        if(LENTRIM(adjustl(parentPath)) .GT. 0) then
            truePath = adjustl(trim(parentPath))//FolderSpe//adjustl(trim(distPath))
        end if

        if(truePath(1:1) .eq. FolderSpe) then
            AbsolutePath = .true.
        end if

        seperatedStrsArray = ''
        call separateStrByString(truePath,FolderSpe,seperatedStrsArray,segNumber)

        truePath = ''

        DO I = 1,segNumber

            if(LENTRIM(adjustl(truePath)) .GT. 0) then
                truePath = truePath(1:LENTRIM(adjustl(truePath)))//FolderSpe//seperatedStrsArray(I)(1:LENTRIM(adjustl(seperatedStrsArray(I))))
            else
                if(AbsolutePath .eq. .true.) then
                    truePath = FolderSpe//seperatedStrsArray(I)(1:LENTRIM(adjustl(seperatedStrsArray(I))))
                end if

            end if

            INQUIRE(FILE=adjustl(trim(truePath)),EXIST=exits)

            if(.not. exits) then
                write(*,*) "MCPSCU: making data directory: ",truePath
                call SYSTEM("mkdir "//truePath(1:LENTRIM(truePath)))
            end if

        END DO

    end if

    resultPath = adjustl(trim(truePath))

    return
  end function CreateDataFolder


  !*************************************************************
  subroutine AllocateOneDimi_Host(Array,Length,Name)
    implicit none
    !---Dummy Vars---
    integer,dimension(:),allocatable::Array
    integer,intent(in)::Length
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateOneDimi_Host(Array,Name)

    if(Length .GT. 0) then
        allocate(Array(Length),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateOneDimi_Host

  !*************************************************************
  subroutine AllocateOneDimr_Host(Array,Length,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCSF),dimension(:),allocatable::Array
    integer,intent(in)::Length
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateOneDimr_Host(Array,Name)

    if(Length .GT. 0) then
        allocate(Array(Length),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateOneDimr_Host

  !*************************************************************
  subroutine AllocateOneDimd_Host(Array,Length,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCDF),dimension(:),allocatable::Array
    integer,intent(in)::Length
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateOneDimd_Host(Array,Name)

    if(Length .GT. 0) then
        allocate(Array(Length),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateOneDimd_Host

  !*************************************************************
  subroutine AllocateTwoDimi_Host(Array,LengthX,LengthY,Name)
    implicit none
    !---Dummy Vars---
    integer,dimension(:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateTwoDimi_Host(Array,Name)

    if(LengthX .GT. 0 .AND. LengthY .GT. 0) then
        allocate(Array(LengthX,LengthY),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateTwoDimi_Host

  !*************************************************************
  subroutine AllocateTwoDimr_Host(Array,LengthX,LengthY,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCSF),dimension(:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateTwoDimr_Host(Array,Name)

    if(LengthX .GT. 0 .AND. LengthY .GT. 0) then
        allocate(Array(LengthX,LengthY),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateTwoDimr_Host

  !*************************************************************
  subroutine AllocateTwoDimd_Host(Array,LengthX,LengthY,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCDF),dimension(:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateTwoDimd_Host(Array,Name)

    if(LengthX .GT. 0 .AND. LengthY .GT. 0) then
        allocate(Array(LengthX,LengthY),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateTwoDimd_Host


  !*************************************************************
  subroutine AllocateThreeDimi_Host(Array,LengthX,LengthY,LengthZ,Name)
    implicit none
    !---Dummy Vars---
    integer,dimension(:,:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    integer,intent(in)::LengthZ
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateThreeDimi_Host(Array,Name)

    if(LengthX .GT. 0 .AND. LengthY .GT. 0 .AND. LengthZ .GT. 0) then
        allocate(Array(LengthX,LengthY,LengthZ),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateThreeDimi_Host

  !*************************************************************
  subroutine AllocateThreeDimr_Host(Array,LengthX,LengthY,LengthZ,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCSF),dimension(:,:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    integer,intent(in)::LengthZ
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateThreeDimr_Host(Array,Name)

    if(LengthX .GT. 0 .AND. LengthY .GT. 0 .AND. LengthZ .GT. 0) then
        allocate(Array(LengthX,LengthY,LengthZ),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateThreeDimr_Host

  !*************************************************************
  subroutine AllocateThreeDimd_Host(Array,LengthX,LengthY,LengthZ,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCDF),dimension(:,:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    integer,intent(in)::LengthZ
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateThreeDimd_Host(Array,Name)

    if(LengthX .GT. 0 .AND. LengthY .GT. 0 .AND. LengthZ .GT. 0) then
        allocate(Array(LengthX,LengthY,LengthZ),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateThreeDimd_Host


  !*************************************************************
  subroutine DeAllocateOneDimi_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    integer,dimension(:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateOneDimi_Host


  !*************************************************************
  subroutine DeAllocateOneDimr_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCSF),dimension(:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateOneDimr_Host

  !*************************************************************
  subroutine DeAllocateOneDimd_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCDF),dimension(:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateOneDimd_Host

  !*************************************************************
  subroutine DeAllocateTwoDimi_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    integer,dimension(:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateTwoDimi_Host

  !*************************************************************
  subroutine DeAllocateTwoDimr_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCSF),dimension(:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateTwoDimr_Host

  !*************************************************************
  subroutine DeAllocateTwoDimd_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCDF),dimension(:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateTwoDimd_Host


  !*************************************************************
  subroutine DeAllocateThreeDimi_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    integer,dimension(:,:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateThreeDimi_Host

  !*************************************************************
  subroutine DeAllocateThreeDimr_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCSF),dimension(:,:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateThreeDimr_Host

  !*************************************************************
  subroutine DeAllocateThreeDimd_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KMCDF),dimension(:,:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateThreeDimd_Host


  !*****************************************
  subroutine ResizeArrayi_OneDim(TheArray,NewSize)
    implicit none
    !---Dummy Vars---
    integer, dimension(:), allocatable,intent(inout)::TheArray
    integer,intent(in)::NewSize
    !---Local Vars---
    integer::OldSize(1)
    integer, dimension(:), allocatable::tempArray
    integer::istat
    !---Body----
    OldSize = shape(TheArray)

    if(OldSize(1) .ne. NewSize) then

        ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
        !       the tempArray need not need to be allocated, the compiler would
        !       allocate the array automatic that is marked as "allocatable" based on the
        !       assigned array's size
        call AllocateArray_Host(tempArray,OldSize(1),"tempArray")

        tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

        call AllocateArray_Host(TheArray,NewSize,"TheArray")

        TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[0])

        call DeAllocateArray_Host(tempArray,"tempArray")

    end if


    return
  end subroutine ResizeArrayi_OneDim

  !*****************************************
  subroutine ResizeArrayr_OneDim(TheArray,NewSize)
    implicit none
    !---Dummy Vars---
    real(kind=KMCSF), dimension(:), allocatable,intent(inout)::TheArray
    integer,intent(in)::NewSize
    !---Local Vars---
    integer::OldSize(1)
    real(kind=KMCSF), dimension(:), allocatable::tempArray
    integer::istat
    real(kind=KMCSF)::zero_Single
    !---Body----
    zero_Single = 0.D0
    OldSize = shape(TheArray)

    if(OldSize(1) .ne. NewSize) then

        ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
        !       the tempArray need not need to be allocated, the compiler would
        !       allocate the array automatic that is marked as "allocatable" based on the
        !       assigned array's size
        call AllocateArray_Host(tempArray,OldSize(1),"tempArray")

        tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

        call AllocateArray_Host(TheArray,NewSize,"TheArray")

        TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[zero_Single])

        call DeAllocateArray_Host(tempArray,"tempArray")

    end if


    return
  end subroutine ResizeArrayr_OneDim

  !*****************************************
  subroutine ResizeArrayd_OneDim(TheArray,NewSize)
    implicit none
    !---Dummy Vars---
    real(kind=KMCDF), dimension(:), allocatable,intent(inout)::TheArray
    integer,intent(in)::NewSize
    !---Local Vars---
    integer::OldSize(1)
    real(kind=KMCDF), dimension(:), allocatable::tempArray
    integer::istat
    real(kind=KMCDF)::zero_Double
    !---Body----
    zero_Double = 0.D0
    OldSize = shape(TheArray)

    if(OldSize(1) .ne. NewSize) then

        ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
        !       the tempArray need not need to be allocated, the compiler would
        !       allocate the array automatic that is marked as "allocatable" based on the
        !       assigned array's size
        call AllocateArray_Host(tempArray,OldSize(1),"tempArray")

        tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

        call AllocateArray_Host(TheArray,NewSize,"TheArray")

        TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[zero_Double])

        call DeAllocateArray_Host(tempArray,"tempArray")

    end if

    return
  end subroutine ResizeArrayd_OneDim

  !**********************************************************
  subroutine ResizeArrayi_TwoDim(TheArray,NewX,NewY)
    implicit none
    !---Dummy Vars---
    integer, dimension(:,:), allocatable::TheArray
    integer::NewX
    integer::NewY
    !---Local Vars---
    integer::oldShape(2)
    integer, dimension(:,:), allocatable::tempArray
    integer::istat
    !---Body---

    oldShape = shape(TheArray)

    if(oldShape(1) .ne. NewX .and. oldShape(2) .ne. NewY) then


        if(NewX .LE. oldShape(1)) then

            if(NewY .LE. oldShape(2)) then

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateArray_Host(tempArray,NewX,NewY,"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(1:NewX,1:NewY)],SHAPE=[NewX,NewY])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY])
            else
                ! note: Here, we use the "reshape" to copy the array. The main reason for this is that we need to consider the
                !       case that compiler option "-Mallocatable=03" is used, for this compiler options, base on our test, if we
                !       try to copy the array1 to array2 by the form of array2 = array1 and the array2 is marked as "allocatable",
                !       if these two array own same shape, the array2 would be reshaped to the size saming with array1, for instant,
                !       if array1 with size of 1xs1 and array with size of 1xs2 and s1 .ne. s2, in this case, after the operation,
                !       the array2 would be changed to the size of 1xs1 !!!!!
                !       Thus, We apply the "reshape" here, reshape([source1,source2,...],[dim1,dim2,dim3...]) to copy the value and
                !       ensure that array2 would keep the size of 1xs2 after the copying operation.
                !       If s1< s2, we use the form of array2 = reshape([array1],[s2]) here,infact ,there are tow steps are
                !       implicated in this opertion:
                !           1.The reshape operation would construct return an "tempArray" with size 1xs2, and the arrays elements
                !             from "tempArray(1)" to "tempArray(s1)" would be filled by array1
                !           2, array2 = "tempArray" would be executed,
                !       So, based on this startegy, we can keep the size for array2 and copy value from array1 at the same time.
                !       (If s1> s2, we should use array2 = reshape([array1(1:s2)],[s2]) to avoid overflow)
                !
                !       Additional, if s1< s2, the pgfortran compiler would do two things:
                !            1: copy the contents with the size of s1 from to array1 to array2
                !            2: the reminded contents with the sizeof (s2 -s1 +1) in array2
                !               would be assigned to some value automaticlly(the reminded
                !               contents would be filled based on the value in array1, for
                !               instant, if the value in array1 is 1,2,3...., the remineded
                !               contents would be assigned to s1+1,s1+2,....,s2)
                !       In fact, we hope that the value in array2 would be array(1:s1),0(s1+1:s2) ,which means the remined value
                !       shoud not be assigned by some values, thus we should do the additional operation: array2(s1+1:s2) = 0

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateArray_Host(tempArray,NewX,oldShape(2),"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(NewX,oldShape(2))],SHAPE=[NewX,oldShape(2)])


                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[0],ORDER=[1,2])

            end if

        else

            if(NewY .LE. oldShape(2)) then

                call AllocateArray_Host(tempArray,NewY,oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[NewY,oldShape(1)],PAD=[0],ORDER=[2,1])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[0],ORDER=[2,1])

            else

                call AllocateArray_Host(tempArray,oldShape(2),oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[oldShape(2),oldShape(1)],PAD=[0],ORDER=[2,1])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[reshape(SOURCE=[tempArray],&
                                                   SHAPE=[NewX,oldShape(2)],&
                                                   PAD=[0],&
                                                   ORDER=[2,1])],&
                                   SHAPE=[NewX,NewY],&
                                   PAD=[0],&
                                   ORDER=[1,2])

            end if

        end if

        call DeAllocateArray_Host(tempArray,"tempArray")

    end if

    return
  end subroutine ResizeArrayi_TwoDim

  !**********************************************************
  subroutine ResizeArrayr_TwoDim(TheArray,NewX,NewY)
    implicit none
    !---Dummy Vars---
    real(kind=KMCSF), dimension(:,:), allocatable::TheArray
    integer::NewX
    integer::NewY
    !---Local Vars---
    integer::oldShape(2)
    real(kind=KMCSF), dimension(:,:), allocatable::tempArray
    real(kind=KMCSF)::zero_Single
    integer::istat
    !---Body---

    zero_Single = 0.D0

    oldShape = shape(TheArray)

    if(oldShape(1) .ne. NewX .and. oldShape(2) .ne. NewY) then


        if(NewX .LE. oldShape(1)) then

            if(NewY .LE. oldShape(2)) then

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateArray_Host(tempArray,NewX,NewY,"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(1:NewX,1:NewY)],SHAPE=[NewX,NewY])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY])
            else
                ! note: Here, we use the "reshape" to copy the array. The main reason for this is that we need to consider the
                !       case that compiler option "-Mallocatable=03" is used, for this compiler options, base on our test, if we
                !       try to copy the array1 to array2 by the form of array2 = array1 and the array2 is marked as "allocatable",
                !       if these two array own same shape, the array2 would be reshaped to the size saming with array1, for instant,
                !       if array1 with size of 1xs1 and array with size of 1xs2 and s1 .ne. s2, in this case, after the operation,
                !       the array2 would be changed to the size of 1xs1 !!!!!
                !       Thus, We apply the "reshape" here, reshape([source1,source2,...],[dim1,dim2,dim3...]) to copy the value and
                !       ensure that array2 would keep the size of 1xs2 after the copying operation.
                !       If s1< s2, we use the form of array2 = reshape([array1],[s2]) here,infact ,there are tow steps are
                !       implicated in this opertion:
                !           1.The reshape operation would construct return an "tempArray" with size 1xs2, and the arrays elements
                !             from "tempArray(1)" to "tempArray(s1)" would be filled by array1
                !           2, array2 = "tempArray" would be executed,
                !       So, based on this startegy, we can keep the size for array2 and copy value from array1 at the same time.
                !       (If s1> s2, we should use array2 = reshape([array1(1:s2)],[s2]) to avoid overflow)
                !
                !       Additional, if s1< s2, the pgfortran compiler would do two things:
                !            1: copy the contents with the size of s1 from to array1 to array2
                !            2: the reminded contents with the sizeof (s2 -s1 +1) in array2
                !               would be assigned to some value automaticlly(the reminded
                !               contents would be filled based on the value in array1, for
                !               instant, if the value in array1 is 1,2,3...., the remineded
                !               contents would be assigned to s1+1,s1+2,....,s2)
                !       In fact, we hope that the value in array2 would be array(1:s1),0(s1+1:s2) ,which means the remined value
                !       shoud not be assigned by some values, thus we should do the additional operation: array2(s1+1:s2) = 0

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateArray_Host(tempArray,NewX,oldShape(2),"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(NewX,oldShape(2))],SHAPE=[NewX,oldShape(2)])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[zero_Single],ORDER=[1,2])

            end if

        else

            if(NewY .LE. oldShape(2)) then

                call AllocateArray_Host(tempArray,NewY,oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[NewY,oldShape(1)],PAD=[zero_Single],ORDER=[2,1])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[zero_Single],ORDER=[2,1])

            else
                call AllocateArray_Host(tempArray,oldShape(2),oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[oldShape(2),oldShape(1)],PAD=[zero_Single],ORDER=[2,1])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[reshape(SOURCE=[tempArray],&
                                                   SHAPE=[NewX,oldShape(2)],&
                                                   PAD=[zero_Single],&
                                                   ORDER=[2,1])],&
                                   SHAPE=[NewX,NewY],&
                                   PAD=[zero_Single],&
                                   ORDER=[1,2])

            end if

        end if

        call DeAllocateArray_Host(tempArray,"tempArray")

    end if

    return
  end subroutine ResizeArrayr_TwoDim


  !**********************************************************
  subroutine ResizeArrayd_TwoDim(TheArray,NewX,NewY)
    implicit none
    !---Dummy Vars---
    real(kind=KMCDF), dimension(:,:), allocatable::TheArray
    integer::NewX
    integer::NewY
    !---Local Vars---
    integer::oldShape(2)
    real(kind=KMCDF), dimension(:,:), allocatable::tempArray
    real(kind=KMCDF)::zero_Double
    integer::istat
    !---Body---

    zero_Double = 0.D0

    oldShape = shape(TheArray)

    if(oldShape(1) .ne. NewX .and. oldShape(2) .ne. NewY) then


        if(NewX .LE. oldShape(1)) then

            if(NewY .LE. oldShape(2)) then

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateArray_Host(tempArray,NewX,NewY,"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(1:NewX,1:NewY)],SHAPE=[NewX,NewY])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY])
            else
                ! note: Here, we use the "reshape" to copy the array. The main reason for this is that we need to consider the
                !       case that compiler option "-Mallocatable=03" is used, for this compiler options, base on our test, if we
                !       try to copy the array1 to array2 by the form of array2 = array1 and the array2 is marked as "allocatable",
                !       if these two array own same shape, the array2 would be reshaped to the size saming with array1, for instant,
                !       if array1 with size of 1xs1 and array with size of 1xs2 and s1 .ne. s2, in this case, after the operation,
                !       the array2 would be changed to the size of 1xs1 !!!!!
                !       Thus, We apply the "reshape" here, reshape([source1,source2,...],[dim1,dim2,dim3...]) to copy the value and
                !       ensure that array2 would keep the size of 1xs2 after the copying operation.
                !       If s1< s2, we use the form of array2 = reshape([array1],[s2]) here,infact ,there are tow steps are
                !       implicated in this opertion:
                !           1.The reshape operation would construct return an "tempArray" with size 1xs2, and the arrays elements
                !             from "tempArray(1)" to "tempArray(s1)" would be filled by array1
                !           2, array2 = "tempArray" would be executed,
                !       So, based on this startegy, we can keep the size for array2 and copy value from array1 at the same time.
                !       (If s1> s2, we should use array2 = reshape([array1(1:s2)],[s2]) to avoid overflow)
                !
                !       Additional, if s1< s2, the pgfortran compiler would do two things:
                !            1: copy the contents with the size of s1 from to array1 to array2
                !            2: the reminded contents with the sizeof (s2 -s1 +1) in array2
                !               would be assigned to some value automaticlly(the reminded
                !               contents would be filled based on the value in array1, for
                !               instant, if the value in array1 is 1,2,3...., the remineded
                !               contents would be assigned to s1+1,s1+2,....,s2)
                !       In fact, we hope that the value in array2 would be array(1:s1),0(s1+1:s2) ,which means the remined value
                !       shoud not be assigned by some values, thus we should do the additional operation: array2(s1+1:s2) = 0

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateArray_Host(tempArray,NewX,oldShape(2),"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(NewX,oldShape(2))],SHAPE=[NewX,oldShape(2)])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[zero_Double],ORDER=[1,2])

            end if

        else

            if(NewY .LE. oldShape(2)) then

                call AllocateArray_Host(tempArray,NewY,oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[NewY,oldShape(1)],PAD=[zero_Double],ORDER=[2,1])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[zero_Double],ORDER=[2,1])

            else

                call AllocateArray_Host(tempArray,oldShape(2),oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[oldShape(2),oldShape(1)],PAD=[zero_Double],ORDER=[2,1])

                call AllocateArray_Host(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[reshape(SOURCE=[tempArray],&
                                                   SHAPE=[NewX,oldShape(2)],&
                                                   PAD=[zero_Double],&
                                                   ORDER=[2,1])],&
                                   SHAPE=[NewX,NewY],&
                                   PAD=[zero_Double],&
                                   ORDER=[1,2])

            end if

        end if

        call DeAllocateArray_Host(tempArray,"tempArray")

    end if

    return
  end subroutine ResizeArrayd_TwoDim


    !*******************************************
    subroutine DumplicateArrayi_OneDim(TheArray,DumplicateNum)
        implicit none
        !---Dummy Vars---
        integer, dimension(:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DumplicateNum
        !---Local Vars---
        integer::OldSize(1)
        integer, dimension(:), allocatable::tempArray
        integer::NewSize
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        NewSize = DumplicateNum*OldSize(1)

        if(DumplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateArray_Host(tempArray,OldSize(1),"tempArray")

            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

            call AllocateArray_Host(TheArray,NewSize,"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[tempArray])

            call DeAllocateArray_Host(tempArray,"tempArray")

        end if

        return

    end subroutine DumplicateArrayi_OneDim

    !*******************************************
    subroutine DumplicateArrayr_OneDim(TheArray,DumplicateNum)
        implicit none
        !---Dummy Vars---
        real(kind=KMCSF), dimension(:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DumplicateNum
        !---Local Vars---
        integer::OldSize(1)
        real(kind=KMCSF), dimension(:), allocatable::tempArray
        integer::NewSize
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        NewSize = DumplicateNum*OldSize(1)

        if(DumplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateArray_Host(tempArray,OldSize(1),"tempArray")

            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

            call AllocateArray_Host(TheArray,NewSize,"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[tempArray])

            call DeAllocateArray_Host(tempArray,"tempArray")

        end if

        return

    end subroutine DumplicateArrayr_OneDim

    !*******************************************
    subroutine DumplicateArrayd_OneDim(TheArray,DumplicateNum)
        implicit none
        !---Dummy Vars---
        real(kind=KMCDF), dimension(:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DumplicateNum
        !---Local Vars---
        integer::OldSize(1)
        real(kind=KMCDF), dimension(:), allocatable::tempArray
        integer::NewSize
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        NewSize = DumplicateNum*OldSize(1)

        if(DumplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateArray_Host(tempArray,OldSize(1),"tempArray")

            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

            call AllocateArray_Host(TheArray,NewSize,"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[tempArray])

            call DeAllocateArray_Host(tempArray,"tempArray")

        end if

        return

    end subroutine DumplicateArrayd_OneDim

    !*******************************************
    subroutine DumplicateArrayi_TwoDim(TheArray,DumplicateNum)
        implicit none
        !---Dummy Vars---
        integer, dimension(:,:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DumplicateNum
        !---Local Vars---
        integer::OldSize(2)
        integer, dimension(:,:), allocatable::tempArray
        integer::NewSize(2)
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        NewSize(1) = DumplicateNum*OldSize(1)
        NewSize(2) = OldSize(2)

        if(DumplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateArray_Host(tempArray,OldSize(2),OldSize(1),"tempArray")

            !---Restore and reverse
            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(2),OldSize(1)],PAD=[0],ORDER=[2,1])

            call AllocateArray_Host(TheArray,NewSize(1),NewSize(2),"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize(1),NewSize(2)],PAD=[tempArray],ORDER=[2,1])

            call DeAllocateArray_Host(tempArray,"tempArray")

        end if

        return

    end subroutine DumplicateArrayi_TwoDim

    !*******************************************
    subroutine DumplicateArrayr_TwoDim(TheArray,DumplicateNum)
        implicit none
        !---Dummy Vars---
        real(kind=KMCSF), dimension(:,:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DumplicateNum
        !---Local Vars---
        integer::OldSize(2)
        real(kind=KMCSF), dimension(:,:), allocatable::tempArray
        real(kind=KMCSF)::zero_Single
        integer::NewSize(2)
        integer::istat
        !---Body----

        zero_Single = 0.D0

        OldSize = shape(TheArray)

        NewSize(1) = DumplicateNum*OldSize(1)
        NewSize(2) = OldSize(2)

        if(DumplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateArray_Host(tempArray,OldSize(2),OldSize(1),"tempArray")

            !---Restore and reverse
            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(2),OldSize(1)],PAD=[zero_Single],ORDER=[2,1])

            call AllocateArray_Host(TheArray,NewSize(1),NewSize(2),"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize(1),NewSize(2)],PAD=[tempArray],ORDER=[2,1])

            call DeAllocateArray_Host(tempArray,"tempArray")

        end if

        return

    end subroutine DumplicateArrayr_TwoDim


    !*******************************************
    subroutine DumplicateArrayd_TwoDim(TheArray,DumplicateNum)
        implicit none
        !---Dummy Vars---
        real(kind=KMCDF), dimension(:,:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DumplicateNum
        !---Local Vars---
        integer::OldSize(2)
        real(kind=KMCDF), dimension(:,:), allocatable::tempArray
        real(kind=KMCDF)::zero_Double
        integer::NewSize(2)
        integer::istat
        !---Body----

        zero_Double = 0.D0

        OldSize = shape(TheArray)

        NewSize(1) = DumplicateNum*OldSize(1)
        NewSize(2) = OldSize(2)

        if(DumplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateArray_Host(tempArray,OldSize(2),OldSize(1),"tempArray")

            !---Restore and reverse
            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(2),OldSize(1)],PAD=[zero_Double],ORDER=[2,1])

            call AllocateArray_Host(TheArray,NewSize(1),NewSize(2),"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize(1),NewSize(2)],PAD=[tempArray],ORDER=[2,1])

            call DeAllocateArray_Host(tempArray,"tempArray")

        end if

        return

    end subroutine DumplicateArrayd_TwoDim

    function BinarySearch_GE(InputNum,TheArray,LeftBound,RightBound) result(ResultIndex)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::InputNum
        integer,dimension(:)::TheArray
        integer,intent(in)::LeftBound
        integer,intent(in)::RightBound
        integer,intent(out)::ResultIndex
        !---Local Vars---
        integer::ILeft,IRight,IMiddle
        !---Body---
        ResultIndex = -1

        ILeft = LeftBound
        IRight = RightBound

        DO While(ILeft .LE. IRight)

            IMiddle = (ILeft+IRight)/2

            if(InputNum .LE. TheArray(IMiddle)) then
                if(IMiddle .eq. LeftBound .or. InputNum .GT. TheArray(IMiddle-1) ) then
                    ResultIndex = IMiddle
                    exit
                else
                    IRight = IMiddle - 1
                    cycle
                end if
            else
                ILeft = IMiddle + 1
            end if

        END DO

        return
    end function BinarySearch_GE

    !*****************************************************
    function BinarySearch_EQ(InputNum,TheArray,LeftBound,RightBound) result(ResultIndex)
        implicit none
        !---Dummy Vars---
        integer,intent(in)::InputNum
        integer,dimension(:)::TheArray
        integer,intent(in)::LeftBound
        integer,intent(in)::RightBound
        integer,intent(out)::ResultIndex
        !---Local Vars---
        integer::ILeft,IRight,IMiddle
        !---Body---
        ResultIndex = -1

        ILeft = LeftBound
        IRight = RightBound

        DO While(ILeft .LE. IRight)
            IMiddle = (ILeft+IRight)/2

            if(InputNum .LT. TheArray(IMiddle)) then
                IRight = IMiddle - 1
            else if(InputNum .GT. TheArray(IMiddle)) then
                ILeft = IMiddle + 1
            else
                ResultIndex = IMiddle
                exit
            end if

        END DO

        return
    end function BinarySearch_EQ

end module MFLIB_UTILITIES
