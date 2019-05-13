module MCMF_TYPEDEF_DiffusorsDefine_GPU
    use cudafor
    use MCMF_TYPEDEF_DiffusorsValue
    use MCMF_UTILITIES_GPU
    implicit none

    !---Params---

    type,public::Dev_DiffusorTypesMap

        type(DiffusorTypeEntity),device,dimension(:),allocatable::Dev_TypesEntities

        integer,device,dimension(:,:),allocatable::Dev_SingleAtomsDivideArrays

        contains
        procedure,public,non_overridable,pass::Init=>InitDiffuosrMap_Dev
        procedure,public,non_overridable,pass::copyFromHost=>copyDiffusorTypesMapFromHost
        procedure,public,non_overridable,pass::Clean=>Clean_DiffuosrMap_Dev
        Final::CleanDiffuosrMap_Dev
    end type Dev_DiffusorTypesMap

    integer,private,constant::dm_MaxDivideGroups_SingleElement_Diffusors
    integer,private,constant::dm_MapBitLength_Diffusors
    integer,private,constant::dm_MapLength_Diffusors

    !---Constructor---
    private::InitDiffuosrMap_Dev
    private::copyDiffusorTypesMapFromHost
    private::Clean_DiffuosrMap_Dev
    private::CleanDiffuosrMap_Dev

    contains

    !*********************************
    subroutine InitDiffuosrMap_Dev(this,Host_DiffusorTypesMap)
        implicit none
        !---Dummy Vars---
        CLASS(Dev_DiffusorTypesMap)::this
        type(DiffusorTypesMap)::Host_DiffusorTypesMap
        !---Local Vars---

        !---Body---
        if(allocated(this%Dev_TypesEntities)) then
            deallocate(this%Dev_TypesEntities)
        end if
        allocate(this%Dev_TypesEntities(Host_DiffusorTypesMap%MapLength))

        call DeAllocateArray_GPU(this%Dev_SingleAtomsDivideArrays,"Dev_SingleAtomsDivideArrays")
        call AllocateArray_GPU(this%Dev_SingleAtomsDivideArrays,p_ATOMS_GROUPS_NUMBER,Host_DiffusorTypesMap%MaxDivideGroups_SingleElement,"Dev_SingleAtomsDivideArrays")

        return
    end subroutine


    !**********************************
    subroutine copyDiffusorTypesMapFromHost(this,Host_DiffusorTypesMap)
        implicit none
        !---Dummy Vars---
        CLASS(Dev_DiffusorTypesMap)::this
        type(DiffusorTypesMap)::Host_DiffusorTypesMap
        !---Local Vars---
        integer::err
        type(c_ptr)::hp_Clusters
        type(c_devptr)::dp_Clusters
        !---Body---

        dm_MapLength_Diffusors = Host_DiffusorTypesMap%MapLength

        dm_MapBitLength_Diffusors = Host_DiffusorTypesMap%MapBitLength

        dm_MaxDivideGroups_SingleElement_Diffusors = Host_DiffusorTypesMap%MaxDivideGroups_SingleElement

        hp_Clusters = c_loc(Host_DiffusorTypesMap%TypesEntities)
        dp_Clusters = c_devloc(this%Dev_TypesEntities)

        err = cudaMemcpy(dp_Clusters,hp_Clusters,sizeof(Host_DiffusorTypesMap%TypesEntities))

        this%Dev_SingleAtomsDivideArrays = Host_DiffusorTypesMap%SingleAtomsDivideArrays


        return
    end subroutine copyDiffusorTypesMapFromHost

    !**********************************
    subroutine Clean_DiffuosrMap_Dev(this)
        implicit none
        !---Dummy Vars---
        CLASS(Dev_DiffusorTypesMap)::this
        !---Body---

        if(allocated(this%Dev_TypesEntities)) then
            deallocate(this%Dev_TypesEntities)
        end if

        call DeAllocateArray_GPU(this%Dev_SingleAtomsDivideArrays,"Dev_SingleAtomsDivideArrays")

        return
    end subroutine

    !*********************************
    subroutine CleanDiffuosrMap_Dev(this)
        implicit none
        !---Dummy Vars---
        TYPE(Dev_DiffusorTypesMap)::this
        !---Body---
        call this%Clean()

        return
    end subroutine

    !**********************************
    attributes(device) subroutine Dev_GetValueFromDiffusorsMap(Key,Dev_TypesEntities,Dev_SingleAtomsDivideArrays,TheValue)
        implicit none
        !---Dummy Vars---
        type(ACluster)::Key
        type(DiffusorTypeEntity),device::Dev_TypesEntities(*) ! When the nollvm compiler option is used, the attributes(device) dummy vars array should write as (*) for one dimension,cannot be (:)
        integer,device::Dev_SingleAtomsDivideArrays(p_ATOMS_GROUPS_NUMBER,*) ! When the nollvm compiler option is used, the attributes(device) dummy vars array should write as (x,*) for two dimension, cannot be (:,:)
        type(DiffusorValue)::TheValue
        !---Local Vars---
        integer(kind=KMCLINT)::Code
        integer(kind=KMCLINT)::reSparedCode
        integer(kind=KMCLINT)::IndexFor
        integer(kind=KMCLINT)::NextIndex
        !---Body---
        call Dev_GetCode(Key%m_Atoms,Dev_SingleAtomsDivideArrays,Code)

        call Dev_Hash(Code,reSparedCode)

        call Dev_GetIndexFor(reSparedCode,IndexFor)

        DO While(IndexFor .GT. 0)

            if(Dev_TypesEntities(IndexFor)%Code .eq. Code) then
                TheValue = Dev_TypesEntities(IndexFor)%TheValue
                exit
            end if

            IndexFor = Dev_TypesEntities(IndexFor)%NextIndex
        END DO


        return
    end subroutine Dev_GetValueFromDiffusorsMap

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

            J = BinarySearch_GE_DEV(Atoms(I)%m_NA,p_ATOMS_GROUPS_NUMBER,Dev_SingleAtomsDivideArrays,I,1,dm_MaxDivideGroups_SingleElement_Diffusors)

            Code = ISHFT(Code,dm_MapBitLength_Diffusors) + J

        END DO

        return
    end subroutine Dev_GetCode

    !**********************************
    attributes(device) subroutine Dev_Hash(Code,reSparedCode)
        implicit none
        ! Purpose: to spare the code to be more uniform
        !---Dummy Vars---
        integer(kind=KMCLINT)::Code
        integer(kind=KMCLINT)::reSparedCode
        !---Local Vars---
        integer(kind=KMCLINT)::TempCode
        !---Body---
        reSparedCode = Code
        TempCode = Code

        TempCode = ISHFT(TempCode,-dm_MapBitLength_Diffusors)

        DO While(TempCode .GT. 0)

            reSparedCode = IOR(reSparedCode,IBITS(TempCode,0,dm_MapBitLength_Diffusors-1))

            TempCode = ISHFT(TempCode,-dm_MapBitLength_Diffusors)

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
        IndexFor = IAND(Code,dm_MapLength_Diffusors)

        return
    end subroutine Dev_GetIndexFor


end module MCMF_TYPEDEF_DiffusorsDefine_GPU
