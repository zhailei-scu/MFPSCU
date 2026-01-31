!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MFLIB_TYPEDEF_SIMULATIONBOXARRAY_GPU
!  use MCLIB_GLOBAL
!  use MCLIB_TYPEDEF_ClustersInfo_GPU
!  use MCLIB_TYPEDEF_SIMULATIONBOXARRAY
!  use MCLIB_TYPEDEF_GEOMETRY_GPU
!  use MCLIB_TYPEDEF_DiffusorsDefine_GPU
!  use MCLIB_TYPEDEF_ReactionsDefine_GPU
!  implicit none
!
!  integer,private,device,dimension(:,:),allocatable::dm_CountsNCArray
!  integer,private,dimension(:,:),allocatable::m_CountsNCArray
!  integer(kind=KMCLINT),private,device,dimension(:,:),allocatable::dm_CountsNAArray
!  integer(kind=KMCLINT),private,dimension(:,:),allocatable::m_CountsNAArray
!
!  type,public::SimulationBoxes_GPU
!
!    type(Dev_DiffusorTypesMap)::dm_DiffusorTypesMap
!
!    type(Dev_ReactionsMap)::dm_ReactionsMap
!
!    type(ClustersInfo_GPU)::dm_ClusterInfo_GPU
!
!    !************Info about Geometry**********
!    type(GrainBoundary_Dev)::dm_GrainBoundary
!
!    integer, device, dimension(:,:), allocatable::dm_SEActIndexBox                ! the start and end active clusters Index in each Box
!                                                                                  ! (this array is got while each box is re-scanned and
!                                                                                  ! dm_SEActIndexBox(1,1)=1,
!                                                                                  ! dm_SEActIndexBox(1,2)=dm_SEActIndexBox(1,1) + m_Boxes%Box(1)%SNC(p_ActiveFree_statu) + m_Boxes%Box(1)%SNC(p_ActiveINGB_statu)- 1,
!                                                                                  ! dm_SEActIndexBox(2,1)=dm_SEActIndexBox(1,2) + 1,
!                                                                                  ! dm_SEActIndexBox(2,2)=dm_SEActIndexBox(2,1) + m_Boxes%Box(2)%SNC(p_ActiveFree_statu) + m_Boxes%Box(2)%SNC(p_ActiveINGB_statu) - 1)
!
!    integer, device, dimension(:,:), allocatable::dm_SEUsedIndexBox               ! the actual simulated start and end index cluster in each Box
!
!    integer, device, dimension(:,:), allocatable::dm_SEVirtualIndexBox             ! while there are clusters are implanted to the box, to improve the calculation efficiency, we would pre-allocate a bigger block of
!                                                                                  ! array (bigger than required) to leave some free memory space for each box, and based on the implantation type, we would put lots of
!                                                                                  ! clusters in these free memory one-time. So we need not put clusters in each step, in each step, we need only to move the end cluster index
!                                                                                  ! of each box. The dm_SEVirtualIndexBox is the start and end index for the bigger memory space for each box. The dm_SEUsedIndexBox is the actual start
!                                                                                  ! and end index for each box. By other words, we need only move the dm_SEUsedIndexBox(:,2) for each box in each step.
!
!    integer, device, dimension(:,:), allocatable::dm_SEExpdIndexBox              ! In implantation situation, as the total clusters in each box is changing in each box, it is necessary to update neighbor-list in each step,
!                                                                                  ! obviously, it is pretty inefficiency. To improve the calculation efficiency, as described for dm_SEVirtualIndexBox, a block of bigger memory space
!                                                                                  ! is pre-allocated for each box, we can let the neighbor-list calculation to cover all of this block, but it seems like have over-calculation,
!                                                                                  ! because some clusters are not used. To be a compromise strategy, we use dm_SEExpdIndexBox to indicate the clusters index range to be calculated
!                                                                                  ! the dm_SEExpdIndexBox(IBox,1) = dm_SEUsedIndexBox(IBox,1) = dm_SEVirtualIndexBox(IBox,1),
!                                                                                  !  = dm_SEUsedIndexBox(IBox,2) < dm_SEExpdIndexBox(IBox,2) < dm_SEVirtualIndexBox(IBox,2)
!
!
!    integer, device, dimension(:,:), allocatable::dm_SEAddedClustersBoxes         ! the start and end index for added clusters in each Box
!
!
!    contains
!    procedure,non_overridable,public,pass::InitSimulationBoxes_Dev=>Init_SimulationBoxes_Dev
!    procedure,non_overridable,public,pass::InitBoxesInfo_GPU=>Init_BoxesInfo_GPU
!    procedure,non_overridable,public,pass::CleanBoxesInfo_GPU
!    procedure,non_overridable,public,pass::CopyInBoxesInfoFromHost=>CopyIn_BoxesInfoFromHost
!    procedure,non_overridable,public,pass::CopyInBoxesArrayFromHost=>CopyInBoxesArray_FromHost
!    procedure,non_overridable,public,pass::RescaleBoxes_GPUToCPU=>Rescale_Boxes_GPUToCPU
!    procedure,non_overridable,public,pass::GetBoxesBasicStatistic_AllStatu_GPU
!    procedure,NON_OVERRIDABLE,pass,public::ExpandClustersInfor_CPUToGPU_BoxByBox=>Expand_ClustersInfor_CPUToGPU_BoxByBox
!    procedure,NON_OVERRIDABLE,pass,public::ExpandClustersInfor_CPUToGPU_EqualNum=>Expand_ClustersInfor_CPUToGPU_EqualNum
!    procedure,NON_OVERRIDABLE,pass,public::ExpandClustersInfo_GPUToCPU_BoxByBox=>Expand_ClustersInfor_GPUToCPU_BoxByBox
!    procedure,NON_OVERRIDABLE,pass,public::ExpandClustersInfo_GPUToCPU_EqualNum=>Expand_ClustersInfor_GPUToCPU_EqualNum
!    procedure,NON_OVERRIDABLE,pass,public::SweepUnActiveMemory_GPUToCPU=>Sweep_UnActiveMemory_GPUToCPU
!    procedure,NON_OVERRIDABLE,pass,public::SweepUnActiveMemory_CPUToGPU=>Sweep_UnActiveMemory_CPUToGPU
!    Final::CleanSimulationBoxes_GPU
!
!  end type SimulationBoxes_GPU
!
!  private::Init_SimulationBoxes_Dev
!  private::Init_BoxesInfo_GPU
!  private::CleanBoxesInfo_GPU
!  private::CopyIn_BoxesInfoFromHost
!  private::CopyInBoxesArray_FromHost
!  private::Rescale_Boxes_GPUToCPU
!  private::GetBoxesBasicStatistic_AllStatu_GPU
!  private::Expand_ClustersInfor_CPUToGPU_BoxByBox
!  private::Expand_ClustersInfor_CPUToGPU_EqualNum
!  private::Expand_ClustersInfor_GPUToCPU_BoxByBox
!  private::Expand_ClustersInfor_GPUToCPU_EqualNum
!  private::Sweep_UnActiveMemory_GPUToCPU
!  private::Sweep_UnActiveMemory_CPUToGPU
!  public::CleanSimulationBoxes_GPU
!
!  contains
!
!  !**********************************************
!  subroutine Init_SimulationBoxes_Dev(this,Host_Boxes,Host_SimuCtrlParams)
!    implicit none
!    !---Dummy vars--
!    CLASS(SimulationBoxes_GPU)::this
!    type(SimulationBoxes)::Host_Boxes
!    type(SimulationCtrlParam)::Host_SimuCtrlParams
!    !---Local Vars---
!    integer::MultiBox
!    integer::TotalSize
!    !---Body--
!
!    MultiBox = Host_SimuCtrlParams%MultiBox
!
!    if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!        TotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!    else
!        TotalSize = 0
!    end if
!
!    call this%dm_ClusterInfo_GPU%ReleaseClustersInfo_GPU()
!
!    call this%dm_GrainBoundary%Clean_GrainBoundary_Dev()
!
!    call this%dm_GrainBoundary%InitGrainBoundary_Dev(Host_Boxes%m_GrainBoundary%GrainNum)
!
!    call this%CleanBoxesInfo_GPU()
!
!    call this%InitBoxesInfo_GPU(Host_SimuCtrlParams)
!
!    if(Host_SimuCtrlParams%FreeDiffusion .eq. .false.) then
!        call this%dm_ClusterInfo_GPU%AllocateClustersInfo_GPU(TotalSize,Host_SimuCtrlParams%MAXNEIGHBORNUM)
!    else
!        call this%dm_ClusterInfo_GPU%AllocateClustersInfo_GPU(TotalSize,0)
!    end if
!
!    call this%dm_DiffusorTypesMap%Clean()
!
!    call this%dm_DiffusorTypesMap%Init(Host_Boxes%m_DiffusorTypesMap)
!
!    call this%dm_ReactionsMap%Clean()
!
!    call this%dm_ReactionsMap%Init(Host_Boxes%m_ReactionsMap)
!
!    return
!  end subroutine Init_SimulationBoxes_Dev
!
!  !**********************************************
!  subroutine CopyInBoxesArray_FromHost(this,Host_Boxes,NSIZE,IfCpyNL)
!    CLASS(SimulationBoxes_GPU)::this
!    type(SimulationBoxes)::Host_Boxes
!    integer,intent(in)::NSIZE
!    logical::IfCpyNL
!    !---Body---
!
!    call this%dm_ClusterInfo_GPU%CopyInFromHost(Host_Boxes%m_ClustersInfo_CPU,NSIZE,IfCpyNL)
!
!    call this%dm_GrainBoundary%CopyInGrainBoundaryFromHost(Host_Boxes%m_GrainBoundary)
!
!    call this%CopyInBoxesInfoFromHost(Host_Boxes)
!
!    call this%dm_DiffusorTypesMap%copyFromHost(Host_Boxes%m_DiffusorTypesMap)
!
!    call this%dm_ReactionsMap%copyFromHost(Host_Boxes%m_ReactionsMap)
!
!    call copyInBoxParamsConstant(Host_Boxes%BOXBOUNDARY,Host_Boxes%BOXSIZE,Host_Boxes%HBOXSIZE)
!
!    return
!  end subroutine CopyInBoxesArray_FromHost
!
!  !**********************************************
!  subroutine Init_BoxesInfo_GPU(this,Host_SimuCtrlParams)
!    implicit none
!    !---Dummy Vars---
!    CLASS(SimulationBoxes_GPU)::this
!    type(SimulationCtrlParam)::Host_SimuCtrlParams
!    !---Local Vars---
!    integer::MULTIBOX
!    integer::istat
!    !---Body---
!
!    MULTIBOX = Host_SimuCtrlParams%MultiBox
!
!    call AllocateArray_GPU(this%dm_SEActIndexBox,MULTIBOX,2,"dm_SEActIndexBox")
!
!    call AllocateArray_GPU(this%dm_SEUsedIndexBox,MULTIBOX,2,"dm_SEUsedIndexBox")
!
!    call AllocateArray_GPU(this%dm_SEVirtualIndexBox,MULTIBOX,2,"dm_SEVirtualIndexBox")
!
!    call AllocateArray_GPU(this%dm_SEExpdIndexBox,MULTIBOX,2,"dm_SEExpdIndexBox")
!
!    call AllocateArray_GPU(this%dm_SEAddedClustersBoxes,MULTIBOX,2,"dm_SEAddedClustersBoxes")
!
!    return
!  end subroutine Init_BoxesInfo_GPU
!
!  !**********************************************
!  subroutine CleanBoxesInfo_GPU(this)
!    implicit none
!    !---Dummy Vars---
!    CLASS(SimulationBoxes_GPU) this
!    !---Local Vars---
!    integer::istat
!    !---Body---
!
!    call DeAllocateArray_GPU(this%dm_SEActIndexBox,"dm_SEActIndexBox")
!
!    call DeAllocateArray_GPU(this%dm_SEUsedIndexBox,"dm_SEUsedIndexBox")
!
!    call DeAllocateArray_GPU(this%dm_SEVirtualIndexBox,"dm_SEVirtualIndexBox")
!
!    call DeAllocateArray_GPU(this%dm_SEExpdIndexBox,"dm_SEExpdIndexBox")
!
!    call DeAllocateArray_GPU(this%dm_SEAddedClustersBoxes,"dm_SEAddedClustersBoxes")
!
!    return
!  end subroutine CleanBoxesInfo_GPU
!
!  !**************************************************
!  subroutine CopyIn_BoxesInfoFromHost(this,Host_Boxes)
!    implicit none
!    !---Dummy Vars---
!    CLASS(SimulationBoxes_GPU)::this
!    type(SimulationBoxes)::Host_Boxes
!    !---Local Vars---
!
!    this%dm_SEActIndexBox = Host_Boxes%m_BoxesInfo%SEActIndexBox
!    this%dm_SEUsedIndexBox = Host_Boxes%m_BoxesInfo%SEUsedIndexBox
!    this%dm_SEVirtualIndexBox = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox
!    this%dm_SEExpdIndexBox = Host_Boxes%m_BoxesInfo%SEExpdIndexBox
!    this%dm_SEAddedClustersBoxes = Host_Boxes%m_BoxesInfo%SEAddedClustersBoxes
!
!
!    return
!  end subroutine CopyIn_BoxesInfoFromHost
!
!
!    !*****************************************************************
!    subroutine Expand_ClustersInfor_CPUToGPU_BoxByBox(this,Host_Boxes,Host_SimuCtrlParams,ExpandNCNum)
!        implicit none
!        !-----Dummy Vars-------
!        CLASS(SimulationBoxes_GPU)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParams
!        integer,intent(in),dimension(:),allocatable::ExpandNCNum
!        !---Local Vars---
!        integer::NewTotalSize
!        integer::MultiBox
!        integer::NeighborsNum
!        !---Body-----------
!
!        MultiBox = Host_SimuCtrlParams%MultiBox
!
!        NeighborsNum = Host_SimuCtrlParams%MAXNEIGHBORNUM
!
!        call Host_Boxes%ExpandClustersInfor_CPU(Host_SimuCtrlParams,ExpandNCNum)
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            NewTotalSize = 0
!        end if
!
!        call CleanSimulationBoxes_GPU(this)
!
!        call this%InitSimulationBoxes_Dev(Host_Boxes,Host_SimuCtrlParams)
!
!        call this%CopyInBoxesArrayFromHost(Host_Boxes,NewTotalSize,IfCpyNL=.false.)
!
!        return
!    end subroutine Expand_ClustersInfor_CPUToGPU_BoxByBox
!
!
!    !*****************************************************************
!    subroutine Expand_ClustersInfor_CPUToGPU_EqualNum(this,Host_Boxes,Host_SimuCtrlParams,ExpandNCNum)
!        implicit none
!        !-----Dummy Vars-------
!        CLASS(SimulationBoxes_GPU)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParams
!        integer,intent(in)::ExpandNCNum
!        !---Local Vars---
!        integer::NewTotalSize
!        integer::MultiBox
!        integer::NeighborsNum
!        !---Body-----------
!
!        MultiBox = Host_SimuCtrlParams%MultiBox
!
!        NeighborsNum = Host_SimuCtrlParams%MAXNEIGHBORNUM
!
!        call Host_Boxes%ExpandClustersInfor_CPU(Host_SimuCtrlParams,ExpandNCNum)
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            NewTotalSize = 0
!        end if
!
!        call CleanSimulationBoxes_GPU(this)
!
!        call this%InitSimulationBoxes_Dev(Host_Boxes,Host_SimuCtrlParams)
!
!        call this%CopyInBoxesArrayFromHost(Host_Boxes,NewTotalSize,IfCpyNL=.false.)
!
!        return
!    end subroutine Expand_ClustersInfor_CPUToGPU_EqualNum
!
!
!    !*****************************************************************
!    subroutine Expand_ClustersInfor_GPUToCPU_BoxByBox(this,Host_Boxes,Host_SimuCtrlParams,ExpandNCNum)
!        implicit none
!        !-----Dummy Vars-------
!        CLASS(SimulationBoxes_GPU)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParams
!        integer,intent(in),dimension(:),allocatable::ExpandNCNum
!        !---Local Vars---
!        integer::MultiBox
!        integer::NeighborsNum
!        integer::OldTotalSize
!        integer::NewTotalSize
!        !---Body-----------
!
!        MultiBox = Host_SimuCtrlParams%MultiBox
!
!        NeighborsNum = Host_SimuCtrlParams%MAXNEIGHBORNUM
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            OldTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            OldTotalSize = 0
!        end if
!
!        call this%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,OldTotalSize,IfCpyNl=.false.)
!
!        call Host_Boxes%ExpandClustersInfor_CPU(Host_SimuCtrlParams,ExpandNCNum)
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            NewTotalSize = 0
!        end if
!
!        call CleanSimulationBoxes_GPU(this)
!
!        call this%InitSimulationBoxes_Dev(Host_Boxes,Host_SimuCtrlParams)
!
!        call this%CopyInBoxesArrayFromHost(Host_Boxes,NewTotalSize,IfCpyNL=.false.)
!
!        return
!    end subroutine Expand_ClustersInfor_GPUToCPU_BoxByBox
!
!    !*****************************************************************
!    subroutine Expand_ClustersInfor_GPUToCPU_EqualNum(this,Host_Boxes,Host_SimuCtrlParams,ExpandNCNum)
!        implicit none
!        !-----Dummy Vars-------
!        CLASS(SimulationBoxes_GPU)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParams
!        integer,intent(in)::ExpandNCNum
!        !---Local Vars---
!        integer::MultiBox
!        integer::NeighborsNum
!        integer::OldTotalSize
!        integer::NewTotalSize
!        !---Body-----------
!
!        MultiBox = Host_SimuCtrlParams%MultiBox
!
!        NeighborsNum = Host_SimuCtrlParams%MAXNEIGHBORNUM
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            OldTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            OldTotalSize = 0
!        end if
!
!        call this%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,OldTotalSize,IfCpyNl=.false.)
!
!        call Host_Boxes%ExpandClustersInfor_CPU(Host_SimuCtrlParams,ExpandNCNum)
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            NewTotalSize = 0
!        end if
!
!        write(*,*) "NewTotalSize",NewTotalSize
!
!        call CleanSimulationBoxes_GPU(this)
!
!        call this%InitSimulationBoxes_Dev(Host_Boxes,Host_SimuCtrlParams)
!
!        call this%CopyInBoxesArrayFromHost(Host_Boxes,NewTotalSize,IfCpyNL=.false.)
!
!        return
!    end subroutine Expand_ClustersInfor_GPUToCPU_EqualNum
!
!    !************************************************************
!    subroutine Rescale_Boxes_GPUToCPU(this,Host_Boxes,Host_SimuCtrlParams,DUPXYZ)
!        implicit none
!        !---Dummy Vars---
!        CLASS(SimulationBoxes_GPU)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParams
!        integer,intent(in)::DUPXYZ(3)
!        !---Local Vars---
!        integer::OldTotalSize
!        integer::MultiBox
!        integer::NewTotalSize
!        !---Body---
!        MultiBox = Host_SimuCtrlParams%MultiBox
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            OldTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            OldTotalSize = 0
!        end if
!
!        call this%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,OldTotalSize,IfCpyNl=.false.)
!
!        call Host_Boxes%RescaleBoxes_CPU(Host_SimuCtrlParams,DUPXYZ)
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            NewTotalSize = 0
!        end if
!
!        call CleanSimulationBoxes_GPU(this)
!
!        call this%InitSimulationBoxes_Dev(Host_Boxes,Host_SimuCtrlParams)
!
!        call this%CopyInBoxesArrayFromHost(Host_Boxes,NewTotalSize,IfCpyNL=.false.)
!
!        return
!    end subroutine Rescale_Boxes_GPUToCPU
!
!    !************************************************************
!    subroutine Sweep_UnActiveMemory_GPUToCPU(this,Host_Boxes,Host_SimuCtrlParams)
!        !---Dummy Vars---
!        CLASS(SimulationBoxes_GPU)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParams
!        !---Local Vars---
!        integer::MultiBox
!        integer::OldTotalSize
!        integer::NewTotalSize
!        !---Body---
!
!        MultiBox = Host_SimuCtrlParams%MultiBox
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            OldTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            OldTotalSize = 0
!        end if
!
!        call this%dm_ClusterInfo_GPU%CopyOutToHost(Host_Boxes%m_ClustersInfo_CPU,OldTotalSize,IfCpyNl=.false.)
!
!        call Host_Boxes%SweepUnActiveMemory_CPU(Host_SimuCtrlParams)
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            NewTotalSize = 0
!        end if
!
!        call CleanSimulationBoxes_GPU(this)
!
!        call this%InitSimulationBoxes_Dev(Host_Boxes,Host_SimuCtrlParams)
!
!        call this%CopyInBoxesArrayFromHost(Host_Boxes,NewTotalSize,IfCpyNL=.false.)
!
!        return
!    end subroutine Sweep_UnActiveMemory_GPUToCPU
!
!    !************************************************************
!    subroutine Sweep_UnActiveMemory_CPUToGPU(this,Host_Boxes,Host_SimuCtrlParams)
!        !---Dummy Vars---
!        CLASS(SimulationBoxes_GPU)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParams
!        !---Local Vars---
!        integer::MultiBox
!        integer::NewTotalSize
!        !---Body---
!
!        MultiBox = Host_SimuCtrlParams%MultiBox
!
!        call Host_Boxes%SweepUnActiveMemory_CPU(Host_SimuCtrlParams)
!
!        if(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) .GT. 0) then
!            NewTotalSize = Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(MultiBox,2) - Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(1,1) + 1
!        else
!            NewTotalSize = 0
!        end if
!
!        call CleanSimulationBoxes_GPU(this)
!
!        call this%InitSimulationBoxes_Dev(Host_Boxes,Host_SimuCtrlParams)
!
!        call this%CopyInBoxesArrayFromHost(Host_Boxes,NewTotalSize,IfCpyNL=.false.)
!
!        return
!    end subroutine Sweep_UnActiveMemory_CPUToGPU
!
!    !******************************************
!    subroutine GetBoxesBasicStatistic_AllStatu_GPU(this,Host_Boxes,Host_SimuCtrlParam)
!        implicit none
!        !---Dummy Vars---
!        CLASS(SimulationBoxes_GPU)::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        !---Local Vars---
!        integer::MultiBox
!        integer::IBox
!        integer::IBFROM
!        integer::IBTO
!        integer::NCCount
!        integer::BlockNumEachBox
!        integer::BlockNumEachBoxExpd
!        integer::BlockNumEachBoxVirtual
!        integer::NB
!        integer::NBExpd
!        integer::NBVirtual
!        integer::BX,BY
!        type(dim3)::blocks
!        type(dim3)::threads
!        integer::NBAllocate
!        integer::err
!        integer::IStatu
!        !---Body---
!        MultiBox = Host_SimuCtrlParam%MultiBox
!
!        if(maxval(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(:,2)-Host_Boxes%m_BoxesInfo%SEUsedIndexBox(:,1)) .LT. 0) then
!            return
!        end if
!
!        BlockNumEachBox = (maxval(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(:,2)-Host_Boxes%m_BoxesInfo%SEUsedIndexBox(:,1)))/p_Reduce_BLOCKSIZE + 1
!        BlockNumEachBox = (BlockNumEachBox -1)/2 + 1
!        NB = BlockNumEachBox*MultiBox
!        BX = p_Reduce_BLOCKSIZE
!        BY = 1
!        blocks = dim3(NB,1,1)
!        threads = dim3(BX,BY,1)
!
!        BlockNumEachBoxExpd = (maxval(Host_Boxes%m_BoxesInfo%SEExpdIndexBox(:,2)-Host_Boxes%m_BoxesInfo%SEExpdIndexBox(:,1)))/p_Reduce_BLOCKSIZE + 1
!        BlockNumEachBoxExpd = (BlockNumEachBoxExpd -1)/2 + 1
!        NBExpd = BlockNumEachBoxExpd*MultiBox
!
!        BlockNumEachBoxVirtual = (maxval(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(:,2)-Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(:,1)))/p_Reduce_BLOCKSIZE + 1
!        BlockNumEachBoxVirtual = (BlockNumEachBoxVirtual -1)/2 + 1
!        NBVirtual = BlockNumEachBoxVirtual*MultiBox
!
!        !---Here, the array size is domain by the virtual or expd situation, for used situation,
!        !   The array size maybe a little bigger than used size
!        !   The main purpose to let what happen is that we do not want to allocate a suitable memory size
!        !   for used situation in each step, because for implant, the used size is keeping changing and we should not adjustment below
!        !   memory each step, so we use the virtual or expd size, which mean, while GetBoxesMigCoaleStat_Expd_GPU or
!        !   GetBoxesMigCoaleStat_Virtual_GPU is used, we had get a bigger block of memory size that can ensure the usage
!        !   for next N steps and need not to adjustment memory size.
!        NBAllocate = max(NB,NBExpd,NBVirtual)
!
!        if(.not. allocated(m_CountsNCArray)) then
!            allocate(m_CountsNCArray(p_NUMBER_OF_STATU,NBAllocate))
!        else if(NBAllocate .GT. size(m_CountsNCArray,dim=2)) then
!            deallocate(m_CountsNCArray)
!            allocate(m_CountsNCArray(p_NUMBER_OF_STATU,NBAllocate))
!        end if
!
!        if(.not. allocated(dm_CountsNCArray)) then
!            allocate(dm_CountsNCArray(p_NUMBER_OF_STATU,NBAllocate))
!        else if(NBAllocate .GT. size(dm_CountsNCArray,dim=2)) then
!            deallocate(dm_CountsNCArray)
!            allocate(dm_CountsNCArray(p_NUMBER_OF_STATU,NBAllocate))
!        end if
!
!        if(.not. allocated(m_CountsNAArray)) then
!            allocate(m_CountsNAArray(p_NUMBER_OF_STATU,NBAllocate))
!        else if(NBAllocate .GT. size(m_CountsNAArray,dim=2)) then
!            deallocate(m_CountsNAArray)
!            allocate(m_CountsNAArray(p_NUMBER_OF_STATU,NBAllocate))
!        end if
!
!        if(.not. allocated(dm_CountsNAArray)) then
!            allocate(dm_CountsNAArray(p_NUMBER_OF_STATU,NBAllocate))
!        else if(NBAllocate .GT. size(dm_CountsNAArray,dim=2)) then
!            deallocate(dm_CountsNAArray)
!            allocate(dm_CountsNAArray(p_NUMBER_OF_STATU,NBAllocate))
!        end if
!
!        call Kernel_StatisticBasic3<<<blocks,threads>>>(BlockNumEachBox,                      &
!                                                        this%dm_ClusterInfo_GPU%dm_Clusters,  &
!                                                        this%dm_SEUsedIndexBox,               &
!                                                        dm_CountsNCArray,                     &
!                                                        dm_CountsNAArray)
!
!        m_CountsNCArray = dm_CountsNCArray
!        m_CountsNAArray = dm_CountsNAArray
!
!        DO IStatu = 1,p_NUMBER_OF_STATU
!
!            DO IBox = 1,MultiBox
!                IBFROM = (IBox -1)*BlockNumEachBox + 1
!                IBTO = IBox*BlockNumEachBox
!
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(IStatu) = sum(m_CountsNCArray(IStatu,IBFROM:IBTO))
!                Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NA(IStatu) = sum(m_CountsNAArray(IStatu,IBFROM:IBTO))
!            END DO
!
!            Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(IStatu) = sum(m_CountsNCArray(IStatu,1:NB))
!            Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NA(IStatu) = sum(m_CountsNAArray(IStatu,1:NB))
!
!        END DO
!
!        return
!    end subroutine GetBoxesBasicStatistic_AllStatu_GPU
!
!    !*******************************************
!    attributes(global) subroutine Kernel_StatisticBasic3(BlockNumEachBox,DevArray,Dev_SEIndexBox,       &
!                                                         ResultCountsNCArray,ResultCountsNAArray)
!        !use libm
!        implicit none
!        !---Dummy Vars---
!        integer,value::BlockNumEachBox
!        type(ACluster),device::DevArray(:)
!        integer,device::Dev_SEIndexBox(:,:)
!        integer,device::ResultCountsNCArray(p_NUMBER_OF_STATU,*)
!        integer(kind=KMCLINT),device::ResultCountsNAArray(p_NUMBER_OF_STATU,*)
!        !---Local Vars---
!        integer::tid
!        integer::bid
!        integer::cid
!        integer::bid0
!        integer::IBox
!        integer::scid
!        integer::ecid
!        integer::IC
!        integer,shared::Share_CountsNCOneStatu(p_Reduce_BLOCKSIZE)
!        integer(kind=KMCLINT),shared::Share_CountsNAOneStatu(p_Reduce_BLOCKSIZE)
!        integer::I
!        integer::IAtomsGroup
!        integer::IStatu
!        !---Body---
!        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
!        bid = (blockidx%y  - 1)*griddim%x  + blockidx%x
!        cid = (bid -1)*p_Reduce_BLOCKSIZE + tid
!
!        IBox = (bid - 1)/BlockNumEachBox + 1
!
!        bid0 = (IBox - 1)*BlockNumEachBox + 1
!
!        scid = Dev_SEIndexBox(IBox,1)
!
!        ecid = Dev_SEIndexBox(IBox,2)
!
!        IC = scid + (bid - bid0)*p_Reduce_BLOCKSIZE*2 + tid - 1
!
!        DO IStatu = 1,p_NUMBER_OF_STATU
!
!            Share_CountsNCOneStatu(tid) = 0
!            Share_CountsNAOneStatu(tid) = 0
!
!            if(IC .LE. ecid) then
!                if(DevArray(IC)%m_Statu .eq. IStatu) then
!
!                    DO IAtomsGroup = 1,p_ATOMS_GROUPS_NUMBER
!                        Share_CountsNAOneStatu(tid) = Share_CountsNAOneStatu(tid) + DevArray(IC)%m_Atoms(IAtomsGroup)%m_NA
!                    END DO
!                    Share_CountsNCOneStatu(tid) = 1
!                end if
!
!            end if
!
!            if((IC + p_Reduce_BLOCKSIZE) .LE. ecid) then
!
!                if(DevArray(IC + p_Reduce_BLOCKSIZE)%m_Statu .eq. IStatu) then
!
!                    DO IAtomsGroup = 1,p_ATOMS_GROUPS_NUMBER
!                        Share_CountsNAOneStatu(tid) = Share_CountsNAOneStatu(tid) + DevArray(IC + p_Reduce_BLOCKSIZE)%m_Atoms(IAtomsGroup)%m_NA
!                    END DO
!
!                    Share_CountsNCOneStatu(tid) = Share_CountsNCOneStatu(tid) + 1
!                end if
!
!            end if
!
!            call syncthreads()
!
!            I = p_Reduce_BLOCKSIZE/2
!            DO While(I .GT. 0)
!
!                if(tid .LE. I) then
!                    Share_CountsNCOneStatu(tid) = Share_CountsNCOneStatu(tid) + Share_CountsNCOneStatu(tid+I)
!                    Share_CountsNAOneStatu(tid) = Share_CountsNAOneStatu(tid) + Share_CountsNAOneStatu(tid+I)
!                end if
!
!                call syncthreads()
!
!                I = I/2
!            END DO
!
!            if(tid .eq. 1) then
!                ResultCountsNCArray(IStatu,bid)  = Share_CountsNCOneStatu(1)
!                ResultCountsNAArray(IStatu,bid) = Share_CountsNAOneStatu(1)
!            end if
!
!            call syncthreads()
!        END DO
!
!        return
!    end subroutine Kernel_StatisticBasic3
!
!!        !************************************************************
!!    subroutine GetBoxesBasicStatistic_AllStatu_GPU(this,Host_Boxes,Host_SimuCtrlParam)
!!      implicit none
!!      !---Dummy Vars---
!!      CLASS(SimulationBoxes_GPU)::this
!!      type(SimulationBoxes)::Host_Boxes
!!      type(SimulationCtrlParam)::Host_SimuCtrlParam
!!      !---Local Vars---
!!      integer::Statu
!!      !---Body---
!!      DO Statu = 1,p_NUMBER_OF_STATU
!!        call StatisticBasicOneStatu_Way3_1(this,Statu,Host_Boxes,Host_SimuCtrlParam)
!!      END DO
!!
!!      return
!!    end subroutine GetBoxesBasicStatistic_AllStatu_GPU
!!
!!    !******************************************
!!    subroutine StatisticBasicOneStatu_Way3_1(Dev_Boxes,Statu,Host_Boxes,Host_SimuCtrlParam)
!!        implicit none
!!        !---Dummy Vars---
!!        type(SimulationBoxes_GPU)::Dev_Boxes
!!        integer,intent(in)::Statu
!!        type(SimulationBoxes)::Host_Boxes
!!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!!        !---Local Vars---
!!        integer::MultiBox
!!        integer::IBox
!!        integer::IBFROM
!!        integer::IBTO
!!        integer::NCCount
!!        integer::BlockNumEachBox
!!        integer::BlockNumEachBoxExpd
!!        integer::BlockNumEachBoxVirtual
!!        integer::NB
!!        integer::NBExpd
!!        integer::NBVirtual
!!        integer::BX,BY
!!        type(dim3)::blocks
!!        type(dim3)::threads
!!        integer::NBAllocate
!!        integer::err
!!        !---Body---
!!        MultiBox = Host_SimuCtrlParam%MultiBox
!!
!!        if(maxval(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(:,2)-Host_Boxes%m_BoxesInfo%SEUsedIndexBox(:,1)) .LT. 0) then
!!            return
!!        end if
!!
!!        BlockNumEachBox = (maxval(Host_Boxes%m_BoxesInfo%SEUsedIndexBox(:,2)-Host_Boxes%m_BoxesInfo%SEUsedIndexBox(:,1)))/p_Reduce_BLOCKSIZE + 1
!!        BlockNumEachBox = (BlockNumEachBox -1)/2 + 1
!!        NB = BlockNumEachBox*MultiBox
!!        BX = p_Reduce_BLOCKSIZE
!!        BY = 1
!!        blocks = dim3(NB,1,1)
!!        threads = dim3(BX,BY,1)
!!
!!        BlockNumEachBoxExpd = (maxval(Host_Boxes%m_BoxesInfo%SEExpdIndexBox(:,2)-Host_Boxes%m_BoxesInfo%SEExpdIndexBox(:,1)))/p_Reduce_BLOCKSIZE + 1
!!        BlockNumEachBoxExpd = (BlockNumEachBoxExpd -1)/2 + 1
!!        NBExpd = BlockNumEachBoxExpd*MultiBox
!!
!!        BlockNumEachBoxVirtual = (maxval(Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(:,2)-Host_Boxes%m_BoxesInfo%SEVirtualIndexBox(:,1)))/p_Reduce_BLOCKSIZE + 1
!!        BlockNumEachBoxVirtual = (BlockNumEachBoxVirtual -1)/2 + 1
!!        NBVirtual = BlockNumEachBoxVirtual*MultiBox
!!
!!        !---Here, the array size is domain by the virtual or expd situation, for used situation,
!!        !   The array size maybe a little bigger than used size
!!        !   The main purpose to let what happen is that we do not want to allocate a suitable memory size
!!        !   for used situation in each step, because for implant, the used size is keeping changing and we should not adjustment below
!!        !   memory each step, so we use the virtual or expd size, which mean, while GetBoxesMigCoaleStat_Expd_GPU or
!!        !   GetBoxesMigCoaleStat_Virtual_GPU is used, we had get a bigger block of memory size that can ensure the usage
!!        !   for next N steps and need not to adjustment memory size.
!!        NBAllocate = max(NB,NBExpd,NBVirtual)
!!
!!        if(.not. allocated(m_CountsNCOneStatuArray)) then
!!            allocate(m_CountsNCOneStatuArray(NBAllocate))
!!        else if(NBAllocate .GT. size(m_CountsNCOneStatuArray)) then
!!            deallocate(m_CountsNCOneStatuArray)
!!            allocate(m_CountsNCOneStatuArray(NBAllocate))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNCOneStatuArray)) then
!!            allocate(dm_CountsNCOneStatuArray(NBAllocate))
!!        else if(NBAllocate .GT. size(dm_CountsNCOneStatuArray)) then
!!            deallocate(dm_CountsNCOneStatuArray)
!!            allocate(dm_CountsNCOneStatuArray(NBAllocate))
!!        end if
!!
!!        if(.not. allocated(m_CountsNAOneStatuArray)) then
!!            allocate(m_CountsNAOneStatuArray(NBAllocate))
!!        else if(NBAllocate .GT. size(m_CountsNAOneStatuArray)) then
!!            deallocate(m_CountsNAOneStatuArray)
!!            allocate(m_CountsNAOneStatuArray(NBAllocate))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNAOneStatuArray)) then
!!            allocate(dm_CountsNAOneStatuArray(NBAllocate))
!!        else if(NBAllocate .GT. size(dm_CountsNAOneStatuArray)) then
!!            deallocate(dm_CountsNAOneStatuArray)
!!            allocate(dm_CountsNAOneStatuArray(NBAllocate))
!!        end if
!!
!!        call Kernel_StatisticBasicOneStatu3<<<blocks,threads>>>(Statu,                                        &
!!                                                                   BlockNumEachBox,                           &
!!                                                                   Dev_Boxes%dm_ClusterInfo_GPU%dm_Clusters,  &
!!                                                                   Dev_Boxes%dm_SEUsedIndexBox,               &
!!                                                                   dm_CountsNCOneStatuArray,                  &
!!                                                                   dm_CountsNAOneStatuArray)
!!
!!        m_CountsNCOneStatuArray = dm_CountsNCOneStatuArray
!!        m_CountsNAOneStatuArray = dm_CountsNAOneStatuArray
!!
!!        DO IBox = 1,MultiBox
!!            IBFROM = (IBox -1)*BlockNumEachBox + 1
!!            IBTO = IBox*BlockNumEachBox
!!
!!            Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC(Statu) = sum(m_CountsNCOneStatuArray(IBFROM:IBTO))
!!            Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NA(Statu) = sum(m_CountsNAOneStatuArray(IBFROM:IBTO))
!!        END DO
!!
!!        Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(Statu) = sum(m_CountsNCOneStatuArray(1:NB))
!!        Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NA(Statu) = sum(m_CountsNAOneStatuArray(1:NB))
!!
!!        return
!!    end subroutine StatisticBasicOneStatu_Way3_1
!!
!!    !*******************************************
!!    attributes(global) subroutine Kernel_StatisticBasicOneStatu3(Statu,BlockNumEachBox,DevArray,Dev_SEIndexBox,       &
!!                                                                 ResultCountsNCOneStatuArray,ResultCountsNAOneStatuArray)
!!        !use libm
!!        implicit none
!!        !---Dummy Vars---
!!        integer,value::Statu
!!        integer,value::BlockNumEachBox
!!        type(ACluster),device::DevArray(:)
!!        integer,device::Dev_SEIndexBox(:,:)
!!        integer,device::ResultCountsNCOneStatuArray(:)
!!        integer(kind=KMCLINT),device::ResultCountsNAOneStatuArray(:)
!!        !---Local Vars---
!!        integer::tid
!!        integer::bid
!!        integer::cid
!!        integer::bid0
!!        integer::IBox
!!        integer::scid
!!        integer::ecid
!!        integer::IC
!!        integer,shared::Share_CountsNCOneStatu(p_Reduce_BLOCKSIZE)
!!        integer(kind=KMCLINT),shared::Share_CountsNAOneStatu(p_Reduce_BLOCKSIZE)
!!        integer::I
!!        integer::IAtomsGroup
!!        !---Body---
!!        tid = (threadidx%y - 1)*blockdim%x + threadidx%x
!!        bid = (blockidx%y  - 1)*griddim%x  + blockidx%x
!!        cid = (bid -1)*p_Reduce_BLOCKSIZE + tid
!!
!!        IBox = (bid - 1)/BlockNumEachBox + 1
!!
!!        bid0 = (IBox - 1)*BlockNumEachBox + 1
!!
!!        scid = Dev_SEIndexBox(IBox,1)
!!
!!        ecid = Dev_SEIndexBox(IBox,2)
!!
!!        IC = scid + (bid - bid0)*p_Reduce_BLOCKSIZE*2 + tid - 1
!!
!!        Share_CountsNCOneStatu(tid) = 0
!!        Share_CountsNAOneStatu(tid) = 0
!!
!!        if(IC .LE. ecid) then
!!            if(DevArray(IC)%m_Statu .eq. Statu) then
!!
!!                DO IAtomsGroup = 1,p_ATOMS_GROUPS_NUMBER
!!                    Share_CountsNAOneStatu(tid) = Share_CountsNAOneStatu(tid) + DevArray(IC)%m_Atoms(IAtomsGroup)%m_NA
!!                END DO
!!                Share_CountsNCOneStatu(tid) = 1
!!            end if
!!
!!        end if
!!
!!        if((IC + p_Reduce_BLOCKSIZE) .LE. ecid) then
!!
!!            if(DevArray(IC + p_Reduce_BLOCKSIZE)%m_Statu .eq. Statu) then
!!
!!                DO IAtomsGroup = 1,p_ATOMS_GROUPS_NUMBER
!!                    Share_CountsNAOneStatu(tid) = Share_CountsNAOneStatu(tid) + DevArray(IC + p_Reduce_BLOCKSIZE)%m_Atoms(IAtomsGroup)%m_NA
!!                END DO
!!
!!                Share_CountsNCOneStatu(tid) = Share_CountsNCOneStatu(tid) + 1
!!            end if
!!
!!        end if
!!
!!        call syncthreads()
!!
!!        I = p_Reduce_BLOCKSIZE/2
!!        DO While(I .GT. 0)
!!
!!            if(tid .LE. I) then
!!                Share_CountsNCOneStatu(tid) = Share_CountsNCOneStatu(tid) + Share_CountsNCOneStatu(tid+I)
!!                Share_CountsNAOneStatu(tid) = Share_CountsNAOneStatu(tid) + Share_CountsNAOneStatu(tid+I)
!!            end if
!!
!!            call syncthreads()
!!
!!            I = I/2
!!        END DO
!!
!!        if(tid .eq. 1) then
!!            ResultCountsNCOneStatuArray(bid)  = Share_CountsNCOneStatu(1)
!!            ResultCountsNAOneStatuArray(bid) = Share_CountsNAOneStatu(1)
!!        end if
!!
!!        return
!!    end subroutine Kernel_StatisticBasicOneStatu3
!
!!    !******************************************
!!    subroutine StatisticBasicOneStatu_Way3(Statu,N0,DevArray,CountsNCOneStatu,CountsNAOneStatu)
!!        implicit none
!!        !---Dummy Vars---
!!        integer,intent(in)::Statu
!!        integer,intent(in)::N0
!!        type(ACluster),device::DevArray(:)
!!        integer::CountsNCOneStatu
!!        integer::CountsNAOneStatu
!!        !---Local Vars---
!!        integer::N
!!        integer::NB
!!        integer::BX,BY
!!        integer::NBX,NBY
!!        type(dim3)::blocks
!!        type(dim3)::threads
!!        integer::ITYPE
!!        !---Body---
!!        N = N0
!!        NB = ((N - 1)/p_Reduce_BLOCKSIZE)/2 + 1
!!        NBX = min(NB,p_BLOCKDIMX)
!!        NBY = (NB - 1)/NBX + 1
!!        BX = p_Reduce_BLOCKSIZE
!!        BY = 1
!!        blocks = dim3(NBX,NBY,1)
!!        threads = dim3(BX,BY,1)
!!
!!        if(.not. allocated(m_CountsNCOneStatuArray)) then
!!            allocate(m_CountsNCOneStatuArray(NB))
!!        else if(NB .ne. size(m_CountsNCOneStatuArray)) then
!!            deallocate(m_CountsNCOneStatuArray)
!!            allocate(m_CountsNCOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNCOneStatuArray)) then
!!            allocate(dm_CountsNCOneStatuArray(NB))
!!        else if(NB .ne. size(dm_CountsNCOneStatuArray)) then
!!            deallocate(dm_CountsNCOneStatuArray)
!!            allocate(dm_CountsNCOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(m_CountsNAOneStatuArray)) then
!!            allocate(m_CountsNAOneStatuArray(NB))
!!        else if(NB .ne. size(m_CountsNAOneStatuArray)) then
!!            deallocate(m_CountsNAOneStatuArray)
!!            allocate(m_CountsNAOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNAOneStatuArray)) then
!!            allocate(dm_CountsNAOneStatuArray(NB))
!!        else if(NB .ne. size(dm_CountsNAOneStatuArray)) then
!!            deallocate(dm_CountsNAOneStatuArray)
!!            allocate(dm_CountsNAOneStatuArray(NB))
!!        end if
!!
!!
!!        call Kernel_StatisticBasicOneStatu3<<<blocks,threads>>>(Statu,                      &
!!                                                                   N,                       &
!!                                                                   DevArray,                &
!!                                                                   dm_CountsNCOneStatuArray,&
!!                                                                   dm_CountsNAOneStatuArray)
!!
!!        DO While(N .GT. 1)
!!            N = NB
!!            NB = ((N - 1)/p_Reduce_BLOCKSIZE)/2 + 1
!!            NBX = min(NB,p_BLOCKDIMX)
!!            NBY = (NB - 1)/NBX + 1
!!            BX = p_Reduce_BLOCKSIZE
!!            BY = 1
!!            blocks = dim3(NBX,NBY,1)
!!            threads = dim3(BX,BY,1)
!!            call Kernel_StatisticCount3<<<blocks,threads>>>(N,                        &
!!                                                            dm_CountsNCOneStatuArray, &
!!                                                            dm_CountsNAOneStatuArray)
!!        END DO
!!
!!        m_CountsNCOneStatuArray = dm_CountsNCOneStatuArray
!!        CountsNCOneStatu = m_CountsNCOneStatuArray(1)
!!        m_CountsNAOneStatuArray = dm_CountsNAOneStatuArray
!!        CountsNAOneStatu = m_CountsNAOneStatuArray(1)
!!
!!        return
!!    end subroutine StatisticBasicOneStatu_Way3
!!
!!    !******************************************
!!    subroutine StatisticBasicOneStatu_Way3_1(Statu,N0,DevArray,CountsNCOneStatu,CountsNAOneStatu)
!!        implicit none
!!        !---Dummy Vars---
!!        integer,intent(in)::Statu
!!        integer,intent(in)::N0
!!        type(ACluster),device::DevArray(:)
!!        integer::CountsNCOneStatu
!!        integer::CountsNAOneStatu
!!        !---Local Vars---
!!        integer::N
!!        integer::NB
!!        integer::BX,BY
!!        integer::NBX,NBY
!!        type(dim3)::blocks
!!        type(dim3)::threads
!!        integer::ITYPE
!!        integer::err
!!        !---Body---
!!        N = N0
!!        NB = ((N - 1)/p_Reduce_BLOCKSIZE)/2 + 1
!!        NBX = min(NB,p_BLOCKDIMX)
!!        NBY = (NB - 1)/NBX + 1
!!        BX = p_Reduce_BLOCKSIZE
!!        BY = 1
!!        blocks = dim3(NBX,NBY,1)
!!        threads = dim3(BX,BY,1)
!!
!!
!!        if(.not. allocated(m_CountsNCOneStatuArray)) then
!!            allocate(m_CountsNCOneStatuArray(NB))
!!        else if(NB .ne. size(m_CountsNCOneStatuArray)) then
!!            deallocate(m_CountsNCOneStatuArray)
!!            allocate(m_CountsNCOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNCOneStatuArray)) then
!!            allocate(dm_CountsNCOneStatuArray(NB))
!!        else if(NB .ne. size(dm_CountsNCOneStatuArray)) then
!!            deallocate(dm_CountsNCOneStatuArray)
!!            allocate(dm_CountsNCOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(m_CountsNAOneStatuArray)) then
!!            allocate(m_CountsNAOneStatuArray(NB))
!!        else if(NB .ne. size(m_CountsNAOneStatuArray)) then
!!            deallocate(m_CountsNAOneStatuArray)
!!            allocate(m_CountsNAOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNAOneStatuArray)) then
!!            allocate(dm_CountsNAOneStatuArray(NB))
!!        else if(NB .ne. size(dm_CountsNAOneStatuArray)) then
!!            deallocate(dm_CountsNAOneStatuArray)
!!            allocate(dm_CountsNAOneStatuArray(NB))
!!        end if
!!
!!        call Kernel_StatisticBasicOneStatu3<<<blocks,threads>>>(Statu,                 &
!!                                                                   N,                       &
!!                                                                   DevArray,                &
!!                                                                   dm_CountsNCOneStatuArray,&
!!                                                                   dm_CountsNAOneStatuArray)
!!
!!        m_CountsNCOneStatuArray = dm_CountsNCOneStatuArray
!!        m_CountsNAOneStatuArray = dm_CountsNAOneStatuArray
!!
!!        CountsNCOneStatu = sum(m_CountsNCOneStatuArray)
!!        CountsNAOneStatu = sum(m_CountsNAOneStatuArray)
!!
!!        return
!!    end subroutine StatisticBasicOneStatu_Way3_1
!!
!!
!!    !*******************************************
!!    attributes(global) subroutine Kernel_StatisticBasicOneStatu3(Statu,N,DevArray,&
!!                                                                      ResultCountsNCOneStatuArray, &
!!                                                                      ResultCountsNAOneStatuArray)
!!        !use libm
!!        implicit none
!!        !---Dummy Vars---
!!        integer,value::Statu
!!        integer,value::N
!!        type(ACluster),device::DevArray(N)
!!        integer,device::ResultCountsNCOneStatuArray(:)
!!        integer,device::ResultCountsNAOneStatuArray(:)
!!        !---Local Vars---
!!        integer::IT,IB,IC
!!        integer,shared::Share_CountsNCOneStatu(p_Reduce_BLOCKSIZE)
!!        integer,shared::Share_CountsNAOneStatu(p_Reduce_BLOCKSIZE)
!!        integer::I
!!        integer::IAtomsGroup
!!        !---Body---
!!        IT = (threadidx%y - 1)*blockdim%x + threadidx%x
!!        IB = (blockidx%y - 1)*griddim%x + blockidx%x
!!        IC = (IB - 1)*(blockdim%x*blockdim%y)*2 + IT
!!
!!        Share_CountsNCOneStatu(IT) = 0
!!        Share_CountsNAOneStatu(IT) = 0
!!
!!        if(IC .LE. N) then
!!            if(DevArray(IC)%m_Statu .eq. Statu) then
!!
!!                DO IAtomsGroup = 1,p_ATOMS_GROUPS_NUMBER
!!                    Share_CountsNAOneStatu(IT) = Share_CountsNAOneStatu(IT) + DevArray(IC)%m_Atoms(IAtomsGroup)%m_NA
!!                END DO
!!                Share_CountsNCOneStatu(IT) = 1
!!            end if
!!
!!        end if
!!
!!        if((IC + p_Reduce_BLOCKSIZE) .LE. N) then
!!
!!            if(DevArray(IC + p_Reduce_BLOCKSIZE)%m_Statu .eq. Statu) then
!!
!!                DO IAtomsGroup = 1,p_ATOMS_GROUPS_NUMBER
!!                    Share_CountsNAOneStatu(IT) = Share_CountsNAOneStatu(IT) + DevArray(IC + p_Reduce_BLOCKSIZE)%m_Atoms(IAtomsGroup)%m_NA
!!                END DO
!!
!!                Share_CountsNCOneStatu(IT) = Share_CountsNCOneStatu(IT) + 1
!!            end if
!!
!!        end if
!!
!!        call syncthreads()
!!
!!        I = p_Reduce_BLOCKSIZE/2
!!        DO While(I .GT. 0)
!!
!!            if(IT .LE. I) then
!!                Share_CountsNCOneStatu(IT) = Share_CountsNCOneStatu(IT) + Share_CountsNCOneStatu(IT+I)
!!                Share_CountsNAOneStatu(IT) = Share_CountsNAOneStatu(IT) + Share_CountsNAOneStatu(IT+I)
!!            end if
!!
!!            call syncthreads()
!!
!!            I = I/2
!!        END DO
!!
!!
!!        if(IT .eq. 1) then
!!            ResultCountsNCOneStatuArray(IB) = Share_CountsNCOneStatu(1)
!!            ResultCountsNAOneStatuArray(IB) = Share_CountsNAOneStatu(1)
!!        end if
!!
!!        return
!!    end subroutine Kernel_StatisticBasicOneStatu3
!!
!!    !*******************************************
!!    attributes(global) subroutine Kernel_StatisticCount3(N,CountsNCOneStatuArray,CountsNAOneStatuArray)
!!        !use libm
!!        implicit none
!!        !---Dummy Vars---
!!        integer,value::N
!!        integer,device::CountsNCOneStatuArray(:)
!!        integer,device::CountsNAOneStatuArray(:)
!!        !---Local Vars---
!!        integer::IT,IB,IC
!!        integer,shared::Share_CountsNCOneStatu(p_Reduce_BLOCKSIZE)
!!        integer,shared::Share_CountsNAOneStatu(p_Reduce_BLOCKSIZE)
!!        integer::I
!!        !---Body---
!!        IT = (threadidx%y - 1)*blockdim%x + threadidx%x
!!        IB = (blockidx%y - 1)*griddim%x + blockidx%x
!!        IC = (IB - 1)*(blockdim%x*blockdim%y)*2 + IT
!!
!!        Share_CountsNCOneStatu(IT) = 0
!!        Share_CountsNAOneStatu(IT) = 0
!!        if(IC .LE. N) then
!!            Share_CountsNCOneStatu(IT) = CountsNCOneStatuArray(IC)
!!            Share_CountsNAOneStatu(IT) = CountsNAOneStatuArray(IC)
!!        end if
!!
!!        if((IC + p_Reduce_BLOCKSIZE) .LE. N) then
!!            Share_CountsNCOneStatu(IT) = Share_CountsNCOneStatu(IT) + CountsNCOneStatuArray(IC + p_Reduce_BLOCKSIZE)
!!            Share_CountsNAOneStatu(IT) = Share_CountsNAOneStatu(IT) + CountsNAOneStatuArray(IC + p_Reduce_BLOCKSIZE)
!!        end if
!!
!!        call syncthreads()
!!
!!        I = p_Reduce_BLOCKSIZE/2
!!
!!        DO While(I .GT. 0)
!!
!!            if(IT .LE. I) then
!!
!!                Share_CountsNCOneStatu(IT) = Share_CountsNCOneStatu(IT) + Share_CountsNCOneStatu(IT+I)
!!
!!                Share_CountsNAOneStatu(IT) = Share_CountsNAOneStatu(IT) + Share_CountsNAOneStatu(IT+I)
!!            end if
!!
!!            call syncthreads()
!!
!!            I = I/2
!!        END DO
!!
!!        if(IT .eq. 1) then
!!            CountsNCOneStatuArray(IB)  = Share_CountsNCOneStatu(1)
!!            Share_CountsNAOneStatu(IB) = Share_CountsNAOneStatu(1)
!!        end if
!!
!!        return
!!    end subroutine Kernel_StatisticCount3
!!
!!    !******************************************
!!    subroutine StatisticBasicOneStatu_Way2(Statu,N0,DevArray,CountsNCOneStatu,CountsNAOneStatu)
!!        implicit none
!!        !---Dummy Vars---
!!        integer,intent(in)::Statu
!!        integer,intent(in)::N0
!!        type(ACluster),device::DevArray(:)
!!        integer::CountsNCOneStatu
!!        integer::CountsNAOneStatu
!!        !---Local Vars---
!!        integer::N
!!        integer::NB
!!        integer::BX,BY
!!        integer::NBX,NBY
!!        type(dim3)::blocks
!!        type(dim3)::threads
!!        integer::ITYPE
!!        !---Body---
!!        N = N0
!!        NB = (N - 1)/p_Reduce_BLOCKSIZE + 1
!!        NBX = min(NB,p_BLOCKDIMX)
!!        NBY = (NB - 1)/NBX + 1
!!        BX = p_Reduce_BLOCKSIZE
!!        BY = 1
!!        blocks = dim3(NBX,NBY,1)
!!        threads = dim3(BX,BY,1)
!!
!!        if(.not. allocated(m_CountsNCOneStatuArray)) then
!!            allocate(m_CountsNCOneStatuArray(NB))
!!        else if(NB .ne. size(m_CountsNCOneStatuArray)) then
!!            deallocate(m_CountsNCOneStatuArray)
!!            allocate(m_CountsNCOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNCOneStatuArray)) then
!!            allocate(dm_CountsNCOneStatuArray(NB))
!!        else if(NB .ne. size(dm_CountsNCOneStatuArray)) then
!!            deallocate(dm_CountsNCOneStatuArray)
!!            allocate(dm_CountsNCOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(m_CountsNAOneStatuArray)) then
!!            allocate(m_CountsNAOneStatuArray(NB))
!!        else if(NB .ne. size(m_CountsNAOneStatuArray)) then
!!            deallocate(m_CountsNAOneStatuArray)
!!            allocate(m_CountsNAOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNAOneStatuArray)) then
!!            allocate(dm_CountsNAOneStatuArray(NB))
!!        else if(NB .ne. size(dm_CountsNAOneStatuArray)) then
!!            deallocate(dm_CountsNAOneStatuArray)
!!            allocate(dm_CountsNAOneStatuArray(NB))
!!        end if
!!
!!        call Kernel_StatisticBasicOneStatu2<<<blocks,threads>>>(Statu,                      &
!!                                                                   N,                       &
!!                                                                   DevArray,                &
!!                                                                   dm_CountsNCOneStatuArray,&
!!                                                                   dm_CountsNAOneStatuArray)
!!
!!        DO While(N .GT. 1)
!!            N = NB
!!            NB = (N - 1)/p_Reduce_BLOCKSIZE + 1
!!            NBX = min(NB,p_BLOCKDIMX)
!!            NBY = (NB - 1)/NBX + 1
!!            BX = p_Reduce_BLOCKSIZE
!!            BY = 1
!!            blocks = dim3(NBX,NBY,1)
!!            threads = dim3(BX,BY,1)
!!            call Kernel_StatisticCount2<<<blocks,threads>>>(N,                        &
!!                                                            dm_CountsNCOneStatuArray, &
!!                                                            dm_CountsNAOneStatuArray)
!!        END DO
!!
!!        m_CountsNCOneStatuArray = dm_CountsNCOneStatuArray
!!        CountsNCOneStatu = m_CountsNCOneStatuArray(1)
!!        m_CountsNAOneStatuArray = dm_CountsNAOneStatuArray
!!        CountsNAOneStatu = m_CountsNAOneStatuArray(1)
!!        return
!!    end subroutine StatisticBasicOneStatu_Way2
!!
!!    !******************************************
!!    subroutine StatisticBasicOneStatu_Way2_1(Statu,N0,DevArray,CountsNCOneStatu,CountsNAOneStatu)
!!        implicit none
!!        !---Dummy Vars---
!!        integer,intent(in)::Statu
!!        integer,intent(in)::N0
!!        type(ACluster),device::DevArray(:)
!!        integer::CountsNCOneStatu
!!        integer::CountsNAOneStatu
!!        !---Local Vars---
!!        integer::N
!!        integer::NB
!!        integer::BX,BY
!!        integer::NBX,NBY
!!        type(dim3)::blocks
!!        type(dim3)::threads
!!        integer::ITYPE
!!        !---Body---
!!        N = N0
!!        NB = (N - 1)/p_Reduce_BLOCKSIZE + 1
!!        NBX = min(NB,p_BLOCKDIMX)
!!        NBY = (NB - 1)/NBX + 1
!!        BX = p_Reduce_BLOCKSIZE
!!        BY = 1
!!        blocks = dim3(NBX,NBY,1)
!!        threads = dim3(BX,BY,1)
!!
!!        if(.not. allocated(m_CountsNCOneStatuArray)) then
!!            allocate(m_CountsNCOneStatuArray(NB))
!!        else if(NB .ne. size(m_CountsNCOneStatuArray)) then
!!            deallocate(m_CountsNCOneStatuArray)
!!            allocate(m_CountsNCOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNCOneStatuArray)) then
!!            allocate(dm_CountsNCOneStatuArray(NB))
!!        else if(NB .ne. size(dm_CountsNCOneStatuArray)) then
!!            deallocate(dm_CountsNCOneStatuArray)
!!            allocate(dm_CountsNCOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(m_CountsNAOneStatuArray)) then
!!            allocate(m_CountsNAOneStatuArray(NB))
!!        else if(NB .ne. size(m_CountsNAOneStatuArray)) then
!!            deallocate(m_CountsNAOneStatuArray)
!!            allocate(m_CountsNAOneStatuArray(NB))
!!        end if
!!
!!        if(.not. allocated(dm_CountsNAOneStatuArray)) then
!!            allocate(dm_CountsNAOneStatuArray(NB))
!!        else if(NB .ne. size(dm_CountsNAOneStatuArray)) then
!!            deallocate(dm_CountsNAOneStatuArray)
!!            allocate(dm_CountsNAOneStatuArray(NB))
!!        end if
!!
!!
!!        call Kernel_StatisticBasicOneStatu2<<<blocks,threads>>>(Statu,                      &
!!                                                                   N,                       &
!!                                                                   DevArray,                &
!!                                                                   dm_CountsNCOneStatuArray,&
!!                                                                   dm_CountsNAOneStatuArray)
!!
!!        m_CountsNCOneStatuArray = dm_CountsNCOneStatuArray
!!        m_CountsNAOneStatuArray = dm_CountsNAOneStatuArray
!!
!!        CountsNCOneStatu = sum(m_CountsNCOneStatuArray)
!!        CountsNAOneStatu = sum(m_CountsNAOneStatuArray)
!!
!!        return
!!    end subroutine StatisticBasicOneStatu_Way2_1
!!
!!
!!    !*******************************************
!!    attributes(global) subroutine Kernel_StatisticBasicOneStatu2(Statu,N,DevArray,&
!!                                                                      ResultCountsNCOneStatuArray, &
!!                                                                      ResultCountsNAOneStatuArray)
!!        !use libm
!!        implicit none
!!        !---Dummy Vars---
!!        integer,value::Statu
!!        integer,value::N
!!        type(ACluster),device::DevArray(N)
!!        integer,device::ResultCountsNCOneStatuArray(:)
!!        integer,device::ResultCountsNAOneStatuArray(:)
!!        !---Local Vars---
!!        integer::IT,IB,IC
!!        integer,shared::Share_CountsNCOneStatu(p_Reduce_BLOCKSIZE)
!!        integer,shared::Share_CountsNAOneStatu(p_Reduce_BLOCKSIZE)
!!        integer::I
!!        integer::IAtomsGroup
!!        !---Body---
!!        IT = (threadidx%y - 1)*blockdim%x + threadidx%x
!!        IB = (blockidx%y - 1)*griddim%x + blockidx%x
!!        IC = (IB - 1)*(blockdim%x*blockdim%y)*2 + IT
!!
!!        Share_CountsNCOneStatu(IT) = 0
!!        Share_CountsNAOneStatu(IT) = 0
!!
!!        if(IC .LE. N) then
!!            if(DevArray(IC)%m_Statu .eq. Statu) then
!!
!!                DO IAtomsGroup = 1,p_ATOMS_GROUPS_NUMBER
!!                    Share_CountsNAOneStatu(IT) = Share_CountsNAOneStatu(IT) + DevArray(IC)%m_Atoms(IAtomsGroup)%m_NA
!!                END DO
!!
!!                Share_CountsNCOneStatu(IT) = 1
!!            end if
!!
!!        end if
!!
!!        call syncthreads()
!!
!!        I = p_Reduce_BLOCKSIZE/2
!!        DO While(I .GT. 0)
!!
!!            if(IT .LE. I) then
!!                Share_CountsNCOneStatu(IT) = Share_CountsNCOneStatu(IT) + Share_CountsNCOneStatu(IT+I)
!!                Share_CountsNAOneStatu(IT) = Share_CountsNAOneStatu(IT) + Share_CountsNAOneStatu(IT+I)
!!            end if
!!
!!            call syncthreads()
!!
!!            I = I/2
!!        END DO
!!
!!        if(IT .eq. 1) then
!!            ResultCountsNCOneStatuArray(IB) = Share_CountsNCOneStatu(1)
!!            ResultCountsNAOneStatuArray(IB) = Share_CountsNAOneStatu(1)
!!        end if
!!
!!        return
!!    end subroutine Kernel_StatisticBasicOneStatu2
!!
!!    !*******************************************
!!    attributes(global) subroutine Kernel_StatisticCount2(N,CountsNCOneStatuArray,CountsNAOneStatuArray)
!!        !use libm
!!        implicit none
!!        !---Dummy Vars---
!!        integer,value::N
!!        integer,device::CountsNCOneStatuArray(:)
!!        integer,device::CountsNAOneStatuArray(:)
!!        !---Local Vars---
!!        integer::IT,IB,IC
!!        integer,shared::Share_CountsNCOneStatu(p_Reduce_BLOCKSIZE)
!!        integer,shared::Share_CountsNAOneStatu(p_Reduce_BLOCKSIZE)
!!        integer::I
!!        !---Body---
!!        IT = (threadidx%y - 1)*blockdim%x + threadidx%x
!!        IB = (blockidx%y - 1)*griddim%x + blockidx%x
!!        IC = (IB - 1)*(blockdim%x*blockdim%y) + IT
!!
!!        Share_CountsNCOneStatu(IT) = 0
!!        Share_CountsNAOneStatu(IT) = 0
!!        if(IC .LE. N) then
!!            Share_CountsNCOneStatu(IT) = CountsNCOneStatuArray(IC)
!!            Share_CountsNAOneStatu(IT) = CountsNAOneStatuArray(IC)
!!        end if
!!
!!        call syncthreads()
!!
!!        I = p_Reduce_BLOCKSIZE/2
!!        DO While(I .GT. 0)
!!
!!            if(IT .LE. I) then
!!
!!                Share_CountsNCOneStatu(IT) = Share_CountsNCOneStatu(IT) + Share_CountsNCOneStatu(IT+I)
!!                Share_CountsNAOneStatu(IT) = Share_CountsNAOneStatu(IT) + Share_CountsNAOneStatu(IT+I)
!!            end if
!!
!!            call syncthreads()
!!
!!            I = I/2
!!        END DO
!!
!!
!!        if(IT .eq. 1) then
!!            CountsNCOneStatuArray(IB) = Share_CountsNCOneStatu(1)
!!            CountsNAOneStatuArray(IB) = Share_CountsNAOneStatu(1)
!!        end if
!!
!!        return
!!    end subroutine Kernel_StatisticCount2
!
!    !**********************************************************
!    subroutine CleanSimulationBoxes_GPU(this)
!        !---Dummy Vars---
!        type(SimulationBoxes_GPU)::this
!        !---Body---
!        call this%dm_ClusterInfo_GPU%ReleaseClustersInfo_GPU()
!
!        call this%CleanBoxesInfo_GPU()
!
!        return
!    end subroutine
!
end module MFLIB_TYPEDEF_SIMULATIONBOXARRAY_GPU
