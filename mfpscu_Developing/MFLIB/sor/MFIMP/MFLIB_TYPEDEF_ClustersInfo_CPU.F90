!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
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
        procedure,NON_OVERRIDABLE,pass,public::GetKindIndexByCluster=>Get_KindIndexByCluster
        procedure,NON_OVERRIDABLE,pass,private::Copy_ClustersInfo_CPU
        procedure,NON_OVERRIDABLE,pass,public::Clean=>Clean_ClustersInfo_CPU
        !---overload some operator---
        GENERIC::ASSIGNMENT(=)=>Copy_ClustersInfo_CPU

        !---DeConstructor function---
        FINAL::CleanClustersInfo_CPU

    end type ClustersInfo_CPU

    private::Allocate_ClustersInfo_CPU
    private::Get_ClustersInfo_ArraySize
    private::Get_KindIndexByCluster
    private::Copy_ClustersInfo_CPU
    private::Clean_ClustersInfo_CPU
    private::CleanClustersInfo_CPU

    contains

    !******************************************************************
    subroutine Allocate_ClustersInfo_CPU(this,NNodes,ClustersKinds)
        implicit none
        !------Dummy Vars-------
        CLASS(ClustersInfo_CPU)::this
        integer,intent(in)::NNodes
        integer,intent(in)::ClustersKinds
        !------Local Vars------
        integer::istat
        !---------Body---------
        !---Cluster In host
        if(ClustersKinds .GT. 0) then

            call AllocateArray_Host(this%ClustersKindArray,ClustersKinds,"ClustersKindArray")

            call AllocateArray_Host(this%Concentrate,NNodes,ClustersKinds,"Concentrate")

            this%Concentrate = 0.D0
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
    subroutine Get_ClustersInfo_ArraySize(this,NNodes,ClustersKinds)
        implicit none
        !---Dummy Vars---
        CLASS(ClustersInfo_CPU)::this
        integer,intent(out)::NNodes
        integer,intent(out)::ClustersKinds
        !---Local Vars---
        integer::ClustersKindsNum
        !---Body---
        ClustersKindsNum = size(this%ClustersKindArray)

        if(ClustersKindsNum .ne. size(this%Concentrate,dim=2)) then
            write(*,*) "MFPSCUERROR: The clusters kinds number for ClustersKindArray array is not same with concentrate."
            pause
            stop
        end if

        NNodes = size(this%Concentrate,dim=1)

        ClustersKinds = ClustersKindsNum

        return
    end subroutine Get_ClustersInfo_ArraySize

    !*****************************************************************
    function Get_KindIndexByCluster(this,TheCluster) result(TheIndex)
        implicit none
        !---Dummy Vars---
        CLASS(ClustersInfo_CPU),intent(in)::this
        type(ACluster),intent(in)::TheCluster
        integer::TheIndex
        !---Local Vars---
        integer::KindNum
        integer::IKind
        integer::IElement
        !---Body---
        KindNum= size(this%ClustersKindArray)

        DO IKind = 1,KindNum
            DO IElement = 1,p_NUMBER_OF_STATU
                if(this%ClustersKindArray(IKind)%m_Atoms(IElement)%m_ID .ne. TheCluster%m_Atoms(IElement)%m_ID .OR. &
                   this%ClustersKindArray(IKind)%m_Atoms(IElement)%m_NA .ne. TheCluster%m_Atoms(IElement)%m_NA) then
                    TheIndex = 0
                    exit
                else
                    TheIndex = IKind
                end if
            END DO

            if(TheIndex .ne. 0) then
                exit
            end if

        END DO

        if(TheIndex .eq. 0) then
            write(*,*) "MFPSCUERROR: The cluster kind is not defined."
            write(*,*) TheCluster%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA
            pause
            stop
        end if

        return
    end function Get_KindIndexByCluster

    !*****************************************************************
    subroutine Copy_ClustersInfo_CPU(Dist_Info,Source_Info)
        implicit none
        !---Dummy Vars---
        CLASS(ClustersInfo_CPU),intent(out)::Dist_Info
        CLASS(ClustersInfo_CPU),intent(in)::Source_Info
        !---Local Vars---
        integer::NNodes
        integer::ClustersKinds
        !---Body----
        call Dist_Info%Clean()

        call Get_ClustersInfo_ArraySize(Source_Info,NNodes,ClustersKinds)

        if(ClustersKinds .GT. 0) then
            Dist_Info%ClustersKindArray = reshape(SOURCE=[Source_Info%ClustersKindArray],SHAPE=[ClustersKinds])


            if(NNodes .GT. 0) then
                Dist_Info%Concentrate = reshape(SOURCE=[Source_Info%Concentrate],SHAPE=[NNodes,ClustersKinds])
            end if
        end if

        return
    end subroutine Copy_ClustersInfo_CPU

end module MFLIB_TYPEDEF_ClustersInfo_CPU
