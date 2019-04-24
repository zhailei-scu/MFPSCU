module MCLIB_UTILITIES
  USE MCLIB_CONSTANTS
  USE MCLIB_TYPEDEF_ACLUSTER
  use MCLIB_TYPEDEF_ATOMSLIST
  USE MCLIB_UTILITIES_FORMER
  use MiniUtilities,only:LENTRIM,ISTR
  #ifdef MC_PROFILING
  USE MCLIB_TimeProfile
  #endif

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
    MODULE PROCEDURE AllocateOneDimACluster_Host


    MODULE PROCEDURE AllocateTwoDimi_Host
    MODULE PROCEDURE AllocateTwoDimr_Host
    MODULE PROCEDURE AllocateTwoDimd_Host
    MODULE PROCEDURE AllocateTwoDimACluster_Host

    MODULE PROCEDURE AllocateThreeDimi_Host
    MODULE PROCEDURE AllocateThreeDimr_Host
    MODULE PROCEDURE AllocateThreeDimd_Host

  end interface AllocateArray_Host

  !------------------
  interface DeAllocateArray_Host
    MODULE PROCEDURE DeAllocateOneDimi_Host
    MODULE PROCEDURE DeAllocateOneDimr_Host
    MODULE PROCEDURE DeAllocateOneDimd_Host
    MODULE PROCEDURE DeAllocateOneDimACluster_Host

    MODULE PROCEDURE DeAllocateTwoDimi_Host
    MODULE PROCEDURE DeAllocateTwoDimr_Host
    MODULE PROCEDURE DeAllocateTwoDimd_Host
    MODULE PROCEDURE DeAllocateTwoDimACluster_Host

    MODULE PROCEDURE DeAllocateThreeDimi_Host
    MODULE PROCEDURE DeAllocateThreeDimr_Host
    MODULE PROCEDURE DeAllocateThreeDimd_Host

  end interface DeAllocateArray_Host

  contains

  !*************************************************
  function ResolveSymbol2AtomsSetRange(symbol,BasicAtomsList) result(TheAtomsSetRange)
        implicit none
        !---Dummy Vars---
        character(*)::symbol
        type(AtomsList),intent(in)::BasicAtomsList
        type(AtomsSetRange)::TheAtomsSetRange
        !---Local Vars---
        character(len=20)::ElementsStrs(10)
        character(len=20)::symbolANDNumRangeStr(2)
        character(len=10)::NumRangeStr(2)
        integer::ElemtsGroup
        integer::SepNum
        integer::ElementIndex
        integer::I
        !---Body---

        ElementsStrs = ''

        call TheAtomsSetRange%ReleaseSetsRange()

        call separateStrByString(symbol,p_ElementsTypeSpe,ElementsStrs,ElemtsGroup)

        if(ElemtsGroup .GT. p_ATOMS_GROUPS_NUMBER) then
            write(*,*) "MCPSCUERROR: The kinds of atoms in cluster "//adjustl(trim(symbol))//" is ",ElemtsGroup
            write(*,*) "MCPSCUERROR: Which has been greater than defined max atoms kinds: ",p_ATOMS_GROUPS_NUMBER
            pause
            stop
        end if

        DO I = 1,p_ATOMS_GROUPS_NUMBER
            TheAtomsSetRange%m_SetsRange(I)%m_ID = I
        END DO

        DO I = 1,ElemtsGroup

            symbolANDNumRangeStr = ""

            NumRangeStr = ""

            call separateStrByString(ElementsStrs(I),p_ElementsNumSpe,symbolANDNumRangeStr,SepNum)
            if(SepNum .NE. 2) then
                write(*,*) "MCPSCUERROR: The Element "//ElementsStrs(I)//" define is not correct"
                write(*,*) "In cluster :",symbol
                write(*,*) symbolANDNumRangeStr
                pause
                stop
            end if

            ElementIndex = BasicAtomsList%FindIndexBySymbol(symbolANDNumRangeStr(1))

            if(ElementIndex .LE. 0) then
                write(*,*) "MCPSCUERROR: The element symbol is not defineded: ",symbolANDNumRangeStr(1)
                write(*,*) "In cluster: ",symbol
                pause
                stop
            end if

            if(TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_From .GT. 0) then
                write(*,*) "MCPSCUERROR: Cannot define two ranges for one same element in one cluster symbol: ",symbol
                pause
                stop
            end if

            if(ElementIndex .ne. TheAtomsSetRange%m_SetsRange(ElementIndex)%m_ID) then
                write(*,*) "The pre-putted element index is not true"
                write(*,*) "For cluster: ",symbol
                write(*,*) "In position: ",ElementIndex
                write(*,*) "The pre-putted element index is: ",TheAtomsSetRange%m_SetsRange(ElementIndex)%m_ID
                pause
                stop
            end if

            call separateStrByString(symbolANDNumRangeStr(2),p_NumRangeSpe,NumRangeStr,SepNum)

            if(SepNum .eq. 1) then
                TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_From = ISTR(NumRangeStr(1))
                TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_To = ISTR(NumRangeStr(1))
            else if(SepNum .LE. 0) then
                write(*,*) "MCPSCUERROR: you must special the atoms compents number for cluster: ",symbol
                pause
                stop
            else if(SepNum .GE. 3) then
                write(*,*) "MCPSCUERROR: only up and down limits can be accepted, you have specialied too much limit for cluster :",symbol
                pause
                stop
            else if(SepNum .eq. 2) then
                if(IsStrEqual(NumRangeStr(2),p_InfStr)) then
                    TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_To = 1.D32
                    TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_From = ISTR(NumRangeStr(1))
                else if(IsStrEqual(NumRangeStr(1),p_InfStr)) then
                    TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_To = 1.D32
                    TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_From = ISTR(NumRangeStr(2))
                else
                    TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_From = min(ISTR(NumRangeStr(1)),ISTR(NumRangeStr(2)))
                    TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_To = max(ISTR(NumRangeStr(1)),ISTR(NumRangeStr(2)))
                end if
            end if

            if(TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_From .LE. 0 .AND. TheAtomsSetRange%m_SetsRange(ElementIndex)%m_NA_To .GT. 0) then
                write(*,*) "MCPSCUERROR: The atoms number down limit cannot less than 1 in cluster: ",symbol
                write(*,*) NumRangeStr
                pause
                stop
            end if

        END DO

        return
    end function ResolveSymbol2AtomsSetRange


  !****************************************************************
  real(kind=KMCDF) function Calc_RCUT_Old(MultiBox,CUTREGIONEXTEND,BOXVOLUM,NCAct,RMAX)
    !---Purpose: To calculate the radium while searching nerighbor list
    !---Dummy Vars---
    integer,intent(in)::MultiBox
    real(kind=KMCDF),intent(in)::CUTREGIONEXTEND
    real(kind=KMCDF),intent(in)::BOXVOLUM
    integer,intent(in)::NCAct
    real(kind=KMCDF),intent(in)::RMAX
    !---Local Vars---
    real(kind=KMCDF)::SEP

    #ifdef MC_PROFILING
    call Time_Start(T_Calc_RCUT_Old_Start)
    #endif
    !---Body---
    SEP = MultiBox*BOXVOLUM/dble(NCAct)
    SEP = SEP**(0.33333333333333D0)
    SEP = SEP+2*RMAX
    Calc_RCUT_Old = CUTREGIONEXTEND*SEP

    #ifdef MC_PROFILING
    call Time_Accumulate(T_Calc_RCUT_Old_Start,T_Calc_RCUT_Old)
    #endif
    return
  end function Calc_RCUT_Old

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
  subroutine AllocateOneDimACluster_Host(Array,Length,Name)
    implicit none
    !---Dummy Vars---
    type(ACluster),dimension(:),allocatable::Array
    integer,intent(in)::Length
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateOneDimACluster_Host(Array,Name)

    if(Length .GT. 0) then
        allocate(Array(Length),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateOneDimACluster_Host

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
  subroutine AllocateTwoDimACluster_Host(Array,LengthX,LengthY,Name)
    implicit none
    !---Dummy Vars---
    type(ACluster),dimension(:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateTwoDimACluster_Host(Array,Name)

    if(LengthX .GT. 0 .AND. LengthY .GT. 0) then
        allocate(Array(LengthX,LengthY),STAT=istat)
        if(istat /=0) then
            write(*,*) "MCPSCUERROR: The Array :",Name,"allocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine AllocateTwoDimACluster_Host

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
  subroutine DeAllocateOneDimACluster_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    type(ACluster),dimension(:),allocatable::Array
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
  end subroutine DeAllocateOneDimACluster_Host

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
  subroutine DeAllocateTwoDimACluster_Host(Array,Name)
    implicit none
    !---Dummy Vars---
    type(ACluster),dimension(:,:),allocatable::Array
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
  end subroutine DeAllocateTwoDimACluster_Host


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

  !*******************************************************
  subroutine ResizeClustersArray_OneDim(TheArray,NewSize)
        implicit none
        !---Dummy Vars---
        type(ACluster), dimension(:), allocatable,intent(inout)::TheArray
        integer,intent(in)::NewSize
        !---Local Vars---
        integer::OldSize(1)
        type(ACluster), dimension(:), allocatable::tempArray
        type(ACluster)::zero_Cluster
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        if(OldSize(1) .ne. NewSize) then
            zero_Cluster%m_POS = 0
            zero_Cluster%m_Statu = p_Empty
            zero_Cluster%m_Layer = 1
            zero_Cluster%m_RAD = 0
            zero_Cluster%m_DiffCoeff = 0.D0
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateArray_Host(tempArray,OldSize(1),"tempArray")

            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

            call AllocateArray_Host(TheArray,NewSize,"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[zero_Cluster])

            call DeAllocateArray_Host(tempArray,"tempArray")

        end if

        return
    end subroutine ResizeClustersArray_OneDim

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
    subroutine DumplicateClustersArray_OneDim(TheArray,DumplicateNum)
        implicit none
        !---Dummy Vars---
        type(ACluster), dimension(:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DumplicateNum
        !---Local Vars---
        integer::OldSize(1)
        type(ACluster), dimension(:), allocatable::tempArray
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

    end subroutine DumplicateClustersArray_OneDim

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


end module MCLIB_UTILITIES
