module NUCLEATION_SPACEDIST_GPU
    use cudafor
    use MFLIB_GLOBAL
    use MIGCOALE_IMPLANTATION
    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
    use MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
    use MIGCOALE_ADDONDATA_HOST
    implicit none
    !---Local Vars---
    character*256::fileName
    !---Body---
    real(kind=KMCDF),dimension(:),allocatable::MatA
    real(kind=KMCDF),dimension(:),allocatable::MatB
    real(kind=KMCDF),dimension(:),allocatable::MatC
    real(kind=KMCDF),dimension(:),allocatable::MatD
    real(kind=KMCDF),dimension(:),allocatable::MatW
    real(kind=KMCDF),dimension(:),allocatable::MatH

    integer::Dumplicate = 1

    integer::LastOutPutConifgIndex = 0

    integer::m_StatisticFile

    !--------
    type(ACluster),device,dimension(:),allocatable::dm_ClustersKindArray
    real(kind=KMCDF),device,dimension(:,:),allocatable::dm_Concentrate
    real(kind=KMCDF),device,dimension(:),allocatable::dm_NodeSpace

    contains

    subroutine InitSimu_SpaceDist_GPU(Host_Boxes,Host_SimuCtrlParam)
        implicit none
        !---Dummy Vars---
        type(SimulationBoxes)::Host_Boxes
        type(SimulationCtrlParam)::Host_SimuCtrlParam
        integer::I
        integer::NNodes
        integer::CKind
        !---Body---

        NNodes = Host_Boxes%NNodes

        CKind = Host_Boxes%CKind

        if(allocated(dm_ClustersKindArray)) then
            deallocate(dm_ClustersKindArray)
        end if
        allocate(dm_ClustersKindArray(CKind))
        dm_ClustersKindArray = Host_Boxes%m_ClustersInfo_CPU%ClustersKindArray

        if(allocated(dm_Concentrate)) then
            deallocate(dm_Concentrate)
        end if
        allocate(dm_Concentrate(NNodes,CKind))
        dm_Concentrate = Host_Boxes%m_ClustersInfo_CPU%Concentrate

        if(allocated(dm_NodeSpace)) then
            deallocate(dm_NodeSpace)
        end if
        allocate(dm_NodeSpace(NNodes))
        dm_NodeSpace = Host_Boxes%NodeSpace

        if(allocated(MatA)) then
            deallocate(MatA)
        end if
        allocate(MatA(NNodes))

        if(allocated(MatB)) then
            deallocate(MatB)
        end if
        allocate(MatB(NNodes))

        if(allocated(MatC)) then
            deallocate(MatC)
        end if
        allocate(MatC(NNodes))

        if(allocated(MatD)) then
            deallocate(MatD)
        end if
        allocate(MatD(NNodes))

        if(allocated(MatW)) then
            deallocate(MatW)
        end if
        allocate(MatW(NNodes))

        if(allocated(MatH)) then
            deallocate(MatH)
        end if
        allocate(MatH(NNodes))


        call AvailableIOUnit(m_StatisticFile)

        fileName = Host_SimuCtrlParam%OutFilePath(1:len_trim(Host_SimuCtrlParam%OutFilePath))//FolderSpe//"Statistic.dat"

        open(Unit=m_StatisticFile,file=fileName(1:len_trim(fileName)))

        write(m_StatisticFile,FMT="(15(A15,1x))") "Step","Time","TStep","NPOWER0Ave","NPOWER1DIV2Ave","NPOWER1Ave","NPOWER3DIV2Ave","N1","N2","N3","Rave(nm)"

        return
    end subroutine InitSimu_SpaceDist_GPU




end module NUCLEATION_SPACEDIST_GPU
