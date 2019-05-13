module MFLIB_TYPEDEF_ClustersInfo_GPU
!    USE cudafor
!    USE MFLIB_TYPEDEF_ClustersInfo_CPU
!    USE MCMF_Utilities_GPU
!    implicit none
!
!    type,public::ClustersInfo_GPU
!        !---Cluster In Devices
!        type(Acluster),device,dimension(:),allocatable::dm_Clusters
!        !---NeighborList-Info in Device---
!        integer, device, dimension(:,:), allocatable::dm_INDI
!        integer, device, dimension(:), allocatable::dm_KVOIS
!
!        !---Merging Table---
!        integer, device, dimension(:,:), allocatable::dm_MergeINDI
!        integer, device, dimension(:), allocatable::dm_MergeKVOIS
!        !---Active status---
!        integer, device, dimension(:), allocatable::dm_ActiveStatus
!        !---Active Index---
!        integer, device, dimension(:), allocatable::dm_ActiveIndex
!
!        !---Record Some status
!        integer::NLUpdateCount_Dev = 0
!
!        contains
!
!        procedure,NON_OVERRIDABLE,pass,public::AllocateClustersInfo_GPU=>Allocate_ClustersInfo_GPU
!        procedure,NON_OVERRIDABLE,pass,public::ReleaseClustersInfo_GPU=>Release_ClustersInfo_GPU
!        procedure,non_overridable,pass,public::CopyInFromHost=>CopyIn_FromHost
!        procedure,non_overridable,pass,public::CopyOutToHost=>CopyOut_ToHost
!        procedure,non_overridable,pass,public::GetMemConsumOneClusterInfo=>Get_MemoryConsuming_OneClusterInfo_GPU
!        procedure,non_overridable,public,pass::GetNLUpdateCount_Dev=>Get_NLUpdateCount_Dev
!        procedure,non_overridable,public,pass::SetNLUpdateCount_Dev=>Set_NLUpdateCount_Dev
!        procedure,non_overridable,public,pass::IncreaseOneNLUpdateCount_Dev=>IncreaseOne_NLUpdateCount_Dev
!
!        !---De-Constructor function
!        FINAL::CleanClusterInfo_GPU
!
!    end type ClustersInfo_GPU
!
!    private::Allocate_ClustersInfo_GPU
!    private::Release_ClustersInfo_GPU
!    private::CopyIn_FromHost
!    private::CopyOut_ToHost
!    private::Get_MemoryConsuming_OneClusterInfo_GPU
!    private::Get_NLUpdateCount_Dev
!    private::Set_NLUpdateCount_Dev
!    private::IncreaseOne_NLUpdateCount_Dev
!    private::CleanClusterInfo_GPU
!
!    contains
!
!    !******************************************************************
!    subroutine CopyIn_FromHost(this,Host_ClustersInfo,NSIZE,IfCpyNL)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ClustersInfo_GPU)::this
!        type(ClustersInfo_CPU)::Host_ClustersInfo
!        integer,intent(in)::NSIZE
!        logical::IfCpyNL        ! determine if it is need to copy neighbor-list
!        !---Body---
!
!        call copyInClustersSync(Host_ClustersInfo%m_Clusters,this%dm_Clusters,NSIZE)
!
!        call copyInOneDimSync(Host_ClustersInfo%m_ActiveIndex,this%dm_ActiveIndex,NSIZE)
!
!        !***In must case, the neighbor-list is calculated by GPU,so we may not need to copy neighbor-list from host
!        if(IfCpyNL .eq. .true.) then
!            call copyNeighborListFromCPUTOGPU(NSIZE,this%dm_INDI,this%dm_KVOIS,Host_ClustersInfo%m_list)
!            call this%SetNLUpdateCount_Dev(Host_ClustersInfo%m_list%GetNLUpdateCount_Host())
!        end if
!
!
!        return
!    end subroutine CopyIn_FromHost
!
!    !******************************************************************
!    subroutine CopyOut_ToHost(this,Host_ClustersInfo,NSIZE,IfCpyNL)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ClustersInfo_GPU)::this
!        type(ClustersInfo_CPU)::Host_ClustersInfo
!        integer,intent(in)::NSIZE
!        logical::IfCpyNL        ! determine if it is need to copy neighbor-list
!        !---Body---
!
!        call copyOutClustersSync(Host_ClustersInfo%m_Clusters,this%dm_Clusters,NSIZE)
!
!        call copyOutOneDimSync(Host_ClustersInfo%m_ActiveIndex,this%dm_ActiveIndex,NSIZE)
!
!        !***In must case, the neighbor-list is calculated by GPU,so we may not need to copy neighbor-list from host
!        if(IfCpyNL .eq. .true.) then
!            call copyNeighborListFromGPUTOCPU(NSIZE,this%dm_INDI,this%dm_KVOIS,Host_ClustersInfo%m_list)
!            call Host_ClustersInfo%m_list%SetNLUpdateCount_Host(this%GetNLUpdateCount_Dev())
!        end if
!
!
!        return
!    end subroutine CopyOut_ToHost
!
!    !******************************************************************
!    subroutine CleanClusterInfo_GPU(this)
!        implicit none
!        type(ClustersInfo_GPU) this
!
!        !----Body----
!        !---Cluster In Devices
!        call DeAllocateArray_GPU(this%dm_Clusters,"dm_Clusters")
!
!        !---NeighborList-Info in Device
!        call DeAllocateArray_GPU(this%dm_INDI,"dm_INDI")
!
!        call DeAllocateArray_GPU(this%dm_KVOIS,"dm_KVOIS")
!
!        !---Merging Table
!        call DeAllocateArray_GPU(this%dm_MergeINDI,"dm_MergeINDI")
!
!        call DeAllocateArray_GPU(this%dm_MergeKVOIS,"dm_MergeKVOIS")
!
!        !---Active status---
!        call DeAllocateArray_GPU(this%dm_ActiveStatus,"dm_ActiveStatus")
!
!        !---Active Index---
!        call DeAllocateArray_GPU(this%dm_ActiveIndex,"dm_ActiveIndex")
!
!        call this%SetNLUpdateCount_Dev(0)
!
!        return
!    end subroutine
!
!    !******************************************************************
!    subroutine Allocate_ClustersInfo_GPU(this,AllocSize,NeighborhoodsNum)
!        implicit none
!        !------Dummy Vars-------
!        CLASS(ClustersInfo_GPU)::this
!        integer,intent(in)::AllocSize
!        integer,intent(in)::NeighborhoodsNum
!
!        !------Local Vars------
!        integer::istat
!        !---------Body---------
!        if(AllocSize .GT. 0) then
!            !---Cluster In Devices
!            call AllocateArray_GPU(this%dm_Clusters,AllocSize,"dm_Clusters")
!
!            !---NeighborList-Info in Device
!            call AllocateArray_GPU(this%dm_INDI,AllocSize,NeighborhoodsNum,"dm_INDI")
!
!            call AllocateArray_GPU(this%dm_KVOIS,AllocSize,"dm_KVOIS")
!
!            !---Merging Table
!            call AllocateArray_GPU(this%dm_MergeINDI,AllocSize,NeighborhoodsNum,"dm_MergeINDI")
!
!            call AllocateArray_GPU(this%dm_MergeKVOIS,AllocSize,"dm_MergeKVOIS")
!
!            !---Active status---
!            call AllocateArray_GPU(this%dm_ActiveStatus,AllocSize,"dm_ActiveStatus")
!
!            !---Active Index---
!            call AllocateArray_GPU(this%dm_ActiveIndex,AllocSize,"dm_ActiveIndex")
!
!            call this%SetNLUpdateCount_Dev(0)
!
!        end if
!
!        return
!    end subroutine Allocate_ClustersInfo_GPU
!
!    !******************************************************************
!    subroutine Release_ClustersInfo_GPU(this)
!        implicit none
!        !------Dummy Vars-------
!        CLASS(ClustersInfo_GPU)::this
!        !------Local Vars------
!        integer::istat
!        !---------Body---------
!
!        !---Cluster In Devices
!        call DeAllocateArray_GPU(this%dm_Clusters,"dm_Clusters")
!
!        !---NeighborList-Info in Device
!        call DeAllocateArray_GPU(this%dm_INDI,"dm_INDI")
!
!        call DeAllocateArray_GPU(this%dm_KVOIS,"dm_KVOIS")
!
!        !---Merging Table
!        call DeAllocateArray_GPU(this%dm_MergeINDI,"dm_MergeINDI")
!
!        call DeAllocateArray_GPU(this%dm_MergeKVOIS,"dm_MergeKVOIS")
!
!        !---Active status---
!        call DeAllocateArray_GPU(this%dm_ActiveStatus,"dm_ActiveStatus")
!
!        !---Active Index---
!        call DeAllocateArray_GPU(this%dm_ActiveIndex,"dm_ActiveIndex")
!
!        call this%SetNLUpdateCount_Dev(0)
!
!        return
!    end subroutine Release_ClustersInfo_GPU
!
!    !******************************************************************
!    integer function Get_MemoryConsuming_OneClusterInfo_GPU(this,NeighborhoodsNum)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ClustersInfo_GPU)::this
!        integer,intent(in)::NeighborhoodsNum
!        !---Local Vars---
!        type(ClustersInfo_GPU)::oneClusterInfo_GPU
!        !---Body---
!        call oneClusterInfo_GPU%AllocateClustersInfo_GPU(1,NeighborhoodsNum)
!
!        !Get_MemoryConsuming_OneClusterInfo_GPU = sizeof(oneClusterInfo_GPU)
!
!        Get_MemoryConsuming_OneClusterInfo_GPU = sizeof(oneClusterInfo_GPU%dm_Clusters) + &
!                                                 sizeof(oneClusterInfo_GPU%dm_INDI) + &
!                                                 sizeof(oneClusterInfo_GPU%dm_KVOIS) + &
!                                                 sizeof(oneClusterInfo_GPU%dm_MergeINDI) + &
!                                                 sizeof(oneClusterInfo_GPU%dm_MergeKVOIS) + &
!                                                 sizeof(oneClusterInfo_GPU%dm_ActiveStatus) + &
!                                                 sizeof(oneClusterInfo_GPU%dm_ActiveIndex)
!
!
!
!        call CleanClusterInfo_GPU(oneClusterInfo_GPU)
!
!        return
!    end function Get_MemoryConsuming_OneClusterInfo_GPU
!
!    !**************************************************************************************
!    function Get_NLUpdateCount_Dev(this) result(TheCount)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ClustersInfo_GPU)::this
!        integer,intent(out)::TheCount
!        !---Body---
!
!        TheCount = this%NLUpdateCount_Dev
!        return
!    end function Get_NLUpdateCount_Dev
!
!    !**************************************************************************************
!    subroutine Set_NLUpdateCount_Dev(this,TheCount)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ClustersInfo_GPU)::this
!        integer,intent(in)::TheCount
!        !---Body---
!
!        this%NLUpdateCount_Dev = TheCount
!        return
!    end subroutine Set_NLUpdateCount_Dev
!
!    !**************************************************************************************
!    subroutine IncreaseOne_NLUpdateCount_Dev(this)
!        implicit none
!        !---Dummy Vars---
!        CLASS(ClustersInfo_GPU)::this
!        !---Body---
!        this%NLUpdateCount_Dev = this%NLUpdateCount_Dev + 1
!        return
!    end subroutine IncreaseOne_NLUpdateCount_Dev
!
end module MFLIB_TYPEDEF_ClustersInfo_GPU
