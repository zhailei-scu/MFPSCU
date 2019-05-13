module MCMF_TYPEDEF_ReactionsDefine_GPU
    use cudafor
    use MCMF_TYPEDEF_REACTIONSVALUE
    use MCMF_UTILITIES_GPU
    implicit none

    type,public::Dev_ReactionsMap

        type(ReactionEntity),device,dimension(:),allocatable::Dev_RecordsEntities

        integer,device,dimension(:,:),allocatable::Dev_SingleAtomsDivideArrays

        contains
        procedure,public,non_overridable,pass::Init=>InitReactionsMap_Dev
        procedure,public,non_overridable,pass::copyFromHost=>copyReactionsMapFromHost
        procedure,public,non_overridable,pass::Clean=>Clean_ReactionMap_Dev
        Final::CleanReactionMap_Dev
    end type Dev_ReactionsMap

    integer,private,constant::dm_MaxDivideGroups_SingleElement_Reactions
    integer,private,constant::dm_MapBitLength_Reactions
    integer,private,constant::dm_MapLength_Reactions

    !---Constructor---
    private::InitReactionsMap_Dev
    private::copyReactionsMapFromHost
    private::Clean_ReactionMap_Dev
    private::CleanReactionMap_Dev

    contains

    !*********************************
    subroutine InitReactionsMap_Dev(this,Host_ReactionsMap)
        implicit none
        !---Dummy Vars---
        CLASS(Dev_ReactionsMap)::this
        type(ReactionsMap)::Host_ReactionsMap
        !---Local Vars---

        !---Body---
        if(allocated(this%Dev_RecordsEntities)) then
            deallocate(this%Dev_RecordsEntities)
        end if
        allocate(this%Dev_RecordsEntities(Host_ReactionsMap%MapLength))

        call DeAllocateArray_GPU(this%Dev_SingleAtomsDivideArrays,"Dev_SingleAtomsDivideArrays")
        call AllocateArray_GPU(this%Dev_SingleAtomsDivideArrays,p_ATOMS_GROUPS_NUMBER,Host_ReactionsMap%MaxDivideGroups_SingleElement,"Dev_SingleAtomsDivideArrays")

        return
    end subroutine


    !**********************************
    subroutine copyReactionsMapFromHost(this,Host_ReactionsMap)
        implicit none
        !---Dummy Vars---
        CLASS(Dev_ReactionsMap)::this
        type(ReactionsMap)::Host_ReactionsMap
        !---Local Vars---
        integer::err
        type(c_ptr)::hp
        type(c_devptr)::dp
        !---Body---

        dm_MapLength_Reactions = Host_ReactionsMap%MapLength

        dm_MapBitLength_Reactions = Host_ReactionsMap%MapBitLength

        dm_MaxDivideGroups_SingleElement_Reactions = Host_ReactionsMap%MaxDivideGroups_SingleElement

        hp = c_loc(Host_ReactionsMap%RecordsEntities)
        dp = c_devloc(this%Dev_RecordsEntities)

        err = cudaMemcpy(dp,hp,sizeof(Host_ReactionsMap%RecordsEntities))

        this%Dev_SingleAtomsDivideArrays = Host_ReactionsMap%SingleAtomsDivideArrays

        return
    end subroutine copyReactionsMapFromHost

    !**********************************
    subroutine Clean_ReactionMap_Dev(this)
        implicit none
        !---Dummy Vars---
        CLASS(Dev_ReactionsMap)::this
        !---Body---

        if(allocated(this%Dev_RecordsEntities)) then
            deallocate(this%Dev_RecordsEntities)
        end if

        call DeAllocateArray_GPU(this%Dev_SingleAtomsDivideArrays,"Dev_SingleAtomsDivideArrays")

        return
    end subroutine

    !*********************************
    subroutine CleanReactionMap_Dev(this)
        implicit none
        !---Dummy Vars---
        TYPE(Dev_ReactionsMap)::this
        !---Body---
        call this%Clean()

        return
    end subroutine

    !**********************************
    attributes(device) subroutine Dev_GetValueFromReactionsMap(KeySubject,KeyObject,Dev_RecordsEntities,Dev_SingleAtomsDivideArrays,TheValue)
        implicit none
        !---Dummy Vars---
        type(ACluster)::KeySubject
        type(ACluster)::KeyObject
        type(ReactionEntity),device::Dev_RecordsEntities(*) ! When the nollvm compiler option is used, the attributes(device) dummy vars array should write as (*) for one dimension,cannot be (:)
        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! When the nollvm compiler option is used, the attributes(device) dummy vars array should write as (x,*) for two dimension, cannot be (:,:)
        type(ReactionValue)::TheValue
        !---Local Vars---
        integer(kind=KMCLINT)::SubjectCode,ObjectCode
        integer(kind=KMCLINT)::reSparedCode
        integer(kind=KMCLINT)::IndexFor
        integer(kind=KMCLINT)::NextIndex
        !---Body---
        call Dev_GetCode(KeySubject%m_Atoms,Dev_SingleAtomsDivideArrays,SubjectCode)
        call Dev_GetCode(KeyObject%m_Atoms,Dev_SingleAtomsDivideArrays,ObjectCode)

        call Dev_Hash(SubjectCode,ObjectCode,reSparedCode)

        call Dev_GetIndexFor(reSparedCode,IndexFor)

        DO While(IndexFor .GT. 0)

            if(Dev_RecordsEntities(IndexFor)%SubjectCode .eq. SubjectCode .AND. Dev_RecordsEntities(IndexFor)%ObjectCode .eq. ObjectCode) then
                TheValue = Dev_RecordsEntities(IndexFor)%TheValue
                exit
            end if

            !---here, we consider that the reaction is symmetrical which means if A can react with B , vice versa
            if(Dev_RecordsEntities(IndexFor)%SubjectCode .eq. ObjectCode .AND. Dev_RecordsEntities(IndexFor)%ObjectCode .eq. SubjectCode) then
                TheValue = Dev_RecordsEntities(IndexFor)%TheValue
                exit
            end if

            IndexFor = Dev_RecordsEntities(IndexFor)%NextIndex
        END DO


        return
    end subroutine Dev_GetValueFromReactionsMap

    !**********************************
    attributes(device) subroutine Dev_GetCode(Atoms,Dev_SingleAtomsDivideArrays,Code)
        implicit none
        !---Dummy Vars---
        type(Single_AtomsSet)::Atoms(p_ATOMS_GROUPS_NUMBER)
        integer::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*)  ! When the nollvm compiler option is used, the attributes(device) dummy vars array should write as (x,*) for two dimension, cannot be (:,:)
        integer(kind=KMCLINT)::Code
        !---Local Vars---
        integer::I
        integer::J
        !---Body---
        Code = 0

        DO I = 1,p_ATOMS_GROUPS_NUMBER

            J = BinarySearch_GE_DEV(Atoms(I)%m_NA,p_ATOMS_GROUPS_NUMBER,Dev_SingleAtomsDivideArrays,I,1,dm_MaxDivideGroups_SingleElement_Reactions)

            Code = ISHFT(Code,dm_MapBitLength_Reactions) + J

        END DO

        return
    end subroutine Dev_GetCode

    !**********************************
    attributes(device) subroutine Dev_Hash(SubjectCode,ObjectCode,reSparedCode)
        implicit none
        ! Purpose: to spare the code to be more uniform
        !---Dummy Vars---
        integer(kind=KMCLINT)::SubjectCode
        integer(kind=KMCLINT)::ObjectCode
        integer(kind=KMCLINT)::reSparedCode
        !---Local Vars---
        integer(kind=KMCLINT)::TempCode
        !---Body---
        reSparedCode = IOR(SubjectCode,ObjectCode)

        TempCode = reSparedCode
        TempCode = ISHFT(TempCode,-dm_MapBitLength_Reactions)

        DO While(TempCode .GT. 0)

            reSparedCode = IOR(reSparedCode,IBITS(TempCode,0,dm_MapBitLength_Reactions-1))

            TempCode = ISHFT(TempCode,-dm_MapBitLength_Reactions)

        END DO

        return
    end subroutine Dev_Hash

    !************************************
    attributes(device) subroutine Dev_GetIndexFor(Code,IndexFor)
        implicit none
        !---Dummy Vars---
        integer(kind=KMCLINT)::Code
        integer(kind=KMCLINT)::IndexFor
        !---Body---
        IndexFor = IAND(Code,dm_MapLength_Reactions)

        return
    end subroutine Dev_GetIndexFor


end module MCMF_TYPEDEF_ReactionsDefine_GPU
