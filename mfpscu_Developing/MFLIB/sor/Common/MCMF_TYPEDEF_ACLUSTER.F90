module MCMF_TYPEDEF_ACLUSTER
    #ifdef MC_PROFILING
    USE MCMF_TimeProfile
    #endif
    USE MCMF_CONSTANTS
    implicit none

    TYPE,PUBLIC::Single_AtomsSet
        integer::m_ID = 0
        integer::m_NA = 0
    END TYPE Single_AtomsSet

    TYPE,PUBLIC::Single_AtomsSetRange

        integer::m_ID = 0

        integer::m_NA_From = 0

        integer::m_NA_To = 0

    END TYPE Single_AtomsSetRange

    TYPE,PUBLIC::AtomsSetRange
        type(Single_AtomsSetRange),dimension(p_ATOMS_GROUPS_NUMBER)::m_SetsRange

        contains
        procedure,public,pass,non_overridable::ReleaseSetsRange
        procedure,private,pass,non_overridable::PermutationAtomsSetRange2ClusterList
        procedure,public,pass,non_overridable::AtomsSetRange2ClusterList
    END TYPE AtomsSetRange


    TYPE,PUBLIC::ACluster
         type(Single_AtomsSet),dimension(p_ATOMS_GROUPS_NUMBER)::m_Atoms
         real(kind=KMCDF),dimension(3)::m_POS
         integer::m_Layer=1
         real(kind=KMCDF)::m_RAD=0
         integer::m_Statu = p_Empty
         integer::m_GrainID(2) = 0
         real(kind=KMCDF)::m_DiffCoeff = 0.D0

         contains
         procedure,non_overridable,pass,public::CopyClusterFromOther
         procedure,non_overridable,pass,public::Clean_Cluster
         Generic::Assignment(=)=>CopyClusterFromOther
         !*********Important note: The PGI CUDA Fortran not support the Final symbol, the final symbol would cause the
         !*********Compiler error************
         !Final::CleanCluster
    END TYPE ACluster

    TYPE,PUBLIC::AClusterList
        type(ACluster)::TheCluster
        integer,private::ListCount = 0
        type(AClusterList),pointer::next=>null()

        contains
        procedure,public,pass,non_overridable::AppendOneCluster
        procedure,public,pass,non_overridable::AppendOtherClusterList
        procedure,public,pass,non_overridable::GetList_Count=>GetClustersList_Count
        procedure,public,pass,non_overridable::CopyClustersListFromOther
        procedure,public,pass,non_overridable::Clean_ClusterList
        Generic::Assignment(=)=>CopyClustersListFromOther
        Final::CleanClusterList
    END TYPE

    private::CopyClusterFromOther
    private::Clean_Cluster
    private::CleanCluster
    private::ReleaseSetsRange
    private::PermutationAtomsSetRange2ClusterList
    private::AtomsSetRange2ClusterList
    private::AppendOneCluster
    private::AppendOtherClusterList
    private::GetClustersList_Count
    private::CopyClustersListFromOther
    private::Clean_ClusterList
    private::CleanClusterList

    contains

    !*************For AtomsSetRange******************
    subroutine ReleaseSetsRange(this)
        implicit none
        !---Dummy Vars---
        Class(AtomsSetRange)::this
        !---Local Vars---
        integer::I
        !---Body---

        DO I = 1,size(this%m_SetsRange)
            this%m_SetsRange(I)%m_ID = 0
            this%m_SetsRange(I)%m_NA_From = 0
            this%m_SetsRange(I)%m_NA_To = 0
        END DO

        return
    end subroutine ReleaseSetsRange

    !*************For ACluster***********************

    subroutine CopyClusterFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(ACluster),intent(out)::this
        CLASS(ACluster),intent(in)::other
        !---Local Vars---
        integer::IElement
        !---Body---
        DO IElement = 1,p_ATOMS_GROUPS_NUMBER
            this%m_Atoms(IElement) = other%m_Atoms(IElement)
        END DO

        this%m_POS = other%m_POS

        this%m_RAD = other%m_RAD

        this%m_Layer = other%m_Layer

        this%m_Statu = other%m_Statu

        this%m_GrainID = other%m_GrainID

        this%m_DiffCoeff = other%m_DiffCoeff

        return
    end subroutine CopyClusterFromOther


    subroutine Clean_Cluster(this)
        implicit none
        !---Dummy Vars---
        CLASS(ACluster)::this
        !---Local Vars---
        integer::I
        !---Body---

        DO I = 1,p_ATOMS_GROUPS_NUMBER
            this%m_Atoms(I)%m_ID = 0
            this%m_Atoms(I)%m_NA = 0
        END DO

        this%m_POS = 0
        this%m_Layer = 1
        this%m_RAD = 0
        this%m_Statu = p_Empty
        this%m_GrainID(2) = 0
        this%m_DiffCoeff = 0.D0

        return
    end subroutine

    subroutine CleanCluster(this)
        implicit none
        !---Dummy Vars---
        type(ACluster)::this
        !---Local Vars---
        call this%Clean_Cluster()

        return
    end subroutine


    integer function Get_MemoryConsuming_ClusterType()
        implicit none
        type(ACluster)::aCluster

        Get_MemoryConsuming_ClusterType = sizeof(aCluster)
        return
    end function Get_MemoryConsuming_ClusterType

    !**************************************
    function AtomsSetRange2ClusterList(this,SingleAtomsDivideArrays) result(List)
        !---Dummy Vars---
        CLASS(AtomsSetRange)::this
        integer,dimension(:,:),allocatable::SingleAtomsDivideArrays
        type(AClusterList)::List
        !---Local Vars---
        type(ACluster)::Cluster
        !---body---

        call List%Clean_ClusterList()

        call this%PermutationAtomsSetRange2ClusterList(List,Cluster,SingleAtomsDivideArrays,1,1)

        return
    end function

    !**************************************************
    recursive subroutine PermutationAtomsSetRange2ClusterList(this,List,Cluster,SingleAtomsDivideArrays,TheLevel,choosen)
        !---Dummy Vars---
        CLASS(AtomsSetRange)::this
        type(AClusterList)::List
        type(ACluster)::Cluster
        integer,dimension(:,:),allocatable::SingleAtomsDivideArrays
        integer::TheLevel
        integer::choosen
        !---Local Vars---
        integer::IGroup
        !---Body---

        if(TheLevel .GT. 1) then

            Cluster%m_Atoms(TheLevel-1)%m_NA = SingleAtomsDivideArrays(TheLevel-1,choosen)

            Cluster%m_Atoms(TheLevel-1)%m_ID = this%m_SetsRange(TheLevel-1)%m_ID
        end if

        if(TheLevel .eq. (size(this%m_SetsRange)+1)) then
            call List%AppendOneCluster(Cluster)
            return
        else
            DO IGroup = 1,size(SingleAtomsDivideArrays(TheLevel,:))
                if(this%m_SetsRange(TheLevel)%m_NA_From .GT. SingleAtomsDivideArrays(TheLevel,IGroup)) then
                    cycle
                end if

                if(this%m_SetsRange(TheLevel)%m_NA_To .LT. SingleAtomsDivideArrays(TheLevel,IGroup)) then
                    if(IGroup .GT. 1) then
                        exit
                    end if
                end if

                call this%PermutationAtomsSetRange2ClusterList(List,Cluster,SingleAtomsDivideArrays,TheLevel+1,IGroup)

                if(this%m_SetsRange(TheLevel)%m_NA_From .eq. this%m_SetsRange(TheLevel)%m_NA_To) then
                    exit
                end if
            END DO

        end if

        return
    end subroutine

    !***************************************
    subroutine AppendOneCluster(this,newOne)
        implicit none
        !---Dummy Vars---
        CLASS(AClusterList),target::this
        type(ACluster)::newOne
        !---Local Vars---
        type(AClusterList),pointer::cursor=>null(),cursorP=>null()
        !---Body---

        cursorP=>this

        if(.not. associated(cursorP)) then
            write(*,*) "MCPSCUERROR: you need to init the AAClusterList first!"
            pause
            stop
        end if

        if(this%GetList_Count() .LE. 0) then
            this%ListCount = 1
            this%TheCluster = newOne
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
            cursor%next=>null()
            ! The assignment(=) had been overrided
            cursor%TheCluster = newOne
            cursorP%next=>cursor
        end if

        Nullify(cursorP)
        cursorP=>null()
        Nullify(cursor)
        cursor=>null()
        return
    end subroutine AppendOneCluster

    !**************************************
    subroutine AppendOtherClusterList(this,OtherList)
        implicit none
        !---Dummy Vars---
        CLASS(AClusterList),target::this
        type(AClusterList),target::OtherList
        !---Local Vars---
        type(AClusterList),pointer::cursorThis=>null()
        type(AClusterList),pointer::cursorOther=>null()
        !---Body---
        cursorThis=>this

        if(.not. associated(cursorThis)) then
            write(*,*) "MCPSCUERROR: you need to init the AAClusterList first!"
            pause
            stop
        end if

        cursorOther=>OtherList

        if(.not. associated(cursorOther)) then
            return
        end if

        if(cursorOther%GetList_Count() .LE. 0) then
            return
        end if

        DO While(associated(cursorOther))
            call this%AppendOneCluster(cursorOther%TheCluster)
            cursorOther=>cursorOther%next
        END DO

        return
    end subroutine

    !**************************************
    integer function GetClustersList_Count(this)
        implicit none
        !---Dummy Vars---
        CLASS(AClusterList),target::this
        !---Local Vars---
        type(AClusterList),pointer::cursor=>null()
        !---Body---

        cursor=>this

        if(.not. associated(cursor)) then
            write(*,*) "MCPSCUERROR: you need to init the AAClusterList first!"
            pause
            stop
        end if

        GetClustersList_Count = this%ListCount

        return
    end function

    !**************************************
    subroutine CopyClustersListFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(AClusterList),intent(out),target::this
        CLASS(AClusterList),intent(in),target::other
        !---Local Vars---
        type(AClusterList),pointer::thisCursor=>null()
        type(AClusterList),pointer::otherCursor=>null()
        type(AClusterList),pointer::thisCursorP=>null()
        type(AClusterList),pointer::otherCursorP=>null()
        !---Body---

        thisCursorP=>this
        if(.not. associated(thisCursorP)) then
            write(*,*) "MCPSCUERROR: You must allocate the list first !"
            pause
            stop
        end if

        call this%Clean_ClusterList()

        otherCursorP=>other
        if(.not. associated(otherCursorP)) then
            return
        end if

        if(otherCursorP%GetList_Count() .LE. 0) then
            return
        end if

        ! The assignment(=) had been override
        thisCursorP%TheCluster = otherCursorP%TheCluster

        this%ListCount = this%ListCount + 1

        thisCursor=>thisCursorP%next
        otherCursor=>otherCursorP%next
        DO While(associated(otherCursor))

            allocate(thisCursor)
            ! The assignment(=) had been override
            thisCursor%TheCluster = otherCursor%TheCluster

            this%ListCount = this%ListCount + 1

            thisCursorP%next=>thisCursor

            thisCursorP=>thisCursor
            otherCursorP=>otherCursor

            otherCursor=>otherCursor%next
            thisCursor=>thisCursor%next
        END DO

        Nullify(thisCursor)
        thisCursor=>null()
        Nullify(thisCursorP)
        thisCursorP=>null()
        Nullify(otherCursor)
        otherCursor=>null()
        Nullify(otherCursorP)
        otherCursorP=>null()
        return
    end subroutine

    !**************************************
    subroutine Clean_ClusterList(this)
        implicit none
        !---Dummy Vars---
        CLASS(AClusterList),target::this
        !---Local Vars---
        type(AClusterList),pointer::cursor=>null()
        type(AClusterList),pointer::next=>null()
        !---Body---

        cursor=>this

        if(.not. associated(cursor)) then
            return
        end if

        cursor=>this%next

        call this%TheCluster%Clean_Cluster()

        DO While(associated(cursor))
            next=>cursor%next
            call Clean_Cluster(cursor%TheCluster)
            Nullify(cursor)
            deallocate(cursor)
            cursor=>next
        END DO

        this%next=>null()

        this%ListCount = 0

        Nullify(cursor)
        Nullify(next)
        cursor=>null()
        next=>null()

        return
    end subroutine Clean_ClusterList

    !************************************
    subroutine CleanClusterList(this)
        implicit none
        !---Dummy Vars---
        type(AClusterList)::this
        !---Body---

        call this%Clean_ClusterList()

        return
    end subroutine CleanClusterList

end module MCMF_TYPEDEF_ACLUSTER
