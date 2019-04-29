module MFLIB_TYPEDEF_ClustersInfo_CPU

    USE MCMF_CONSTANTS
    USE MCMF_UTILITIES
    USE MCMF_TYPEDEF_ACLUSTER

    implicit none

    type,public::ClustersInfo_CPU
        !---The concentrates distribution---
        real(kind=KMCDF),dimension(:,:),allocatable::Concentrate

        !---The clusters kinds define---
        type(ACluster),dimension(:),allocatable::ClustersKindArray

        contains

        procedure,NON_OVERRIDABLE,pass,public::AllocateClustersInfo_CPU=>Allocate_ClustersInfo_CPU
        procedure,NON_OVERRIDABLE,pass,public::GetClustersInfo_ArraySize=>Get_ClustersInfo_ArraySize
        procedure,NON_OVERRIDABLE,pass,private::Copy_ClustersInfo_CPU
        procedure,NON_OVERRIDABLE,pass,public::Clean=>Clean_ClustersInfo_CPU
        !---overload some operator---
        GENERIC::ASSIGNMENT(=)=>Copy_ClustersInfo_CPU

        !---DeConstructor function---
        FINAL::CleanClustersInfo_CPU

    end type ClustersInfo_CPU

    private::Allocate_ClustersInfo_CPU
    private::Get_ClustersInfo_ArraySize
    private::Copy_ClustersInfo_CPU
    private::Clean_ClustersInfo_CPU
    private::CleanClustersInfo_CPU

    contains

    !******************************************************************
    subroutine Allocate_ClustersInfo_CPU(this,ClustersKinds,NNodes)
        implicit none
        !------Dummy Vars-------
        CLASS(ClustersInfo_CPU)::this
        integer,intent(in)::ClustersKinds
        integer,intent(in)::NNodes
        !------Local Vars------
        integer::istat
        !---------Body---------
        !---Cluster In host
        if(ClustersKinds .GT. 0) then

            call AllocateArray_Host(this%ClustersKindArray,ClustersKinds,"ClustersKindArray")

            call AllocateArray_Host(this%Concentrate,ClustersKinds,NNodes,"Concentrate")
        end if
        return
    end subroutine Allocate_ClustersInfo_CPU

    !******************************************************************
    subroutine Clean_ClustersInfo_CPU(this)
        implicit none
        !------Dummy Vars-------
        CLASS(ClustersInfo_CPU)::this
        !---------Body---------
        !---Cluster In host
        call DeAllocateArray_Host(this%ClustersKindArray,"ClustersKindArray")

        call DeAllocateArray_Host(this%Concentrate,"Concentrate")

        return
    end subroutine Clean_ClustersInfo_CPU

    !******************************************************************
    subroutine CleanClustersInfo_CPU(this)
        implicit none
        !------Dummy Vars-------
        type(ClustersInfo_CPU)::this
        !---------Body---------
        call this%Clean()

        return
    end subroutine CleanClustersInfo_CPU

    !*****************************************************************
    subroutine Get_ClustersInfo_ArraySize(this,ClustersKinds,NNodes)
        implicit none
        !---Dummy Vars---
        CLASS(ClustersInfo_CPU)::this
        integer,intent(out)::ClustersKinds
        integer,intent(out)::NNodes
        !---Local Vars---
        integer::ClustersKindsNum
        !---Body---
        ClustersKindsNum = size(this%ClustersKindArray)

        if(ClustersKindsNum .ne. size(this%Concentrate,dim=1)) then
            write(*,*) "MFPSCUERROR: The clusters kinds number for ClustersKindArray array is not same with concentrate."
            pause
            stop
        end if

        NNodes = size(this%Concentrate,dim=2)

        return
    end subroutine Get_ClustersInfo_ArraySize

    !*****************************************************************
    subroutine Copy_ClustersInfo_CPU(Dist_Info,Source_Info)
        implicit none
        !---Dummy Vars---
        CLASS(ClustersInfo_CPU),intent(out)::Dist_Info
        CLASS(ClustersInfo_CPU),intent(in)::Source_Info
        !---Local Vars---
        integer::ClustersKinds
        integer::NNodes
        !---Body----
        call Dist_Info%Clean()

        call Get_ClustersInfo_ArraySize(Source_Info,ClustersKinds,NNodes)

        if(ClustersKinds .GT. 0) then
            Dist_Info%ClustersKindArray = reshape(SOURCE=[Source_Info%ClustersKindArray],SHAPE=[ClustersKinds])


            if(NNodes .GT. 0) then
                Dist_Info%Concentrate = reshape(SOURCE=[Source_Info%Concentrate],SHAPE=[ClustersKinds,NNodes])
            end if
        end if

        return
    end subroutine Copy_ClustersInfo_CPU

end module MFLIB_TYPEDEF_ClustersInfo_CPU
