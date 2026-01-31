!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MIGCOALE_STATISTIC_CPU
  use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
  use MIGCOALE_TYPEDEF_STATISTICINFO
  use MIGCOALE_ADDONDATA_HOST
  contains

  !*************************************************************
!  subroutine GetBoxesMigCoaleStat_CPU(Host_Boxes,Host_SimuCtrlParam,TheMigCoaleStatisticInfo)
!    !***                 Purpose: To get average size of clusters, maxma size, and so on, at current time(for all boxs)
!    !                 Host_Boxes: the boxes information in host
!    implicit none
!    !---Dummy Vars---
!    type(SimulationBoxes)::Host_Boxes
!    type(SimulationCtrlParam)::Host_SimuCtrlParam
!    type(MigCoaleStatisticInfo_Used)::TheMigCoaleStatisticInfo
!    !---Local Vars---
!    integer::MultiBox
!    integer::IBox
!    integer::IStatu
!    real(kind=KMCDF)::RAVA(p_NUMBER_OF_STATU)
!    integer::NActCountTemp(p_NUMBER_OF_STATU)
!    integer::NCCount(p_NUMBER_OF_STATU)
!    !---Body---
!
!    call Host_Boxes%GetBoxesBasicStatistic_AllStatu_CPU(Host_SimuCtrlParam)
!
!    NCCount = 0
!
!    RAVA = 0.D0
!
!    ASSOCIATE(TMigStatInfo=>TheMigCoaleStatisticInfo%statistic_IntegralBox)
!
!    TMigStatInfo%RMAX = -1.D32
!    TMigStatInfo%RMIN = 1.D32
!    TMigStatInfo%RAVA = 0.D0
!    TMigStatInfo%DiffusorValueMax = -1.D32
!
!    MultiBox = Host_SimuCtrlParam%MultiBox
!
!    DO IBox = 1, MultiBox
!
!        NActCountTemp = 0
!
!        ASSOCIATE(SMigStatInfo=>TheMigCoaleStatisticInfo%statistic_SingleBoxes(IBox))
!            call GetOneBoxMigCoaleStat_Used_CPU(IBox,Host_Boxes,Host_SimuCtrlParam,SMigStatInfo,NActCountTemp)
!
!            DO IStatu = 1, p_NUMBER_OF_STATU
!                if(TMigStatInfo%RMAX(IStatu) .LT. SMigStatInfo%RMAX(IStatu)) then
!                    TMigStatInfo%RMAX(IStatu) = SMigStatInfo%RMAX(IStatu)
!                    TMigStatInfo%ICMAX(IStatu) = SMigStatInfo%ICMAX(IStatu)
!                end if
!
!                if(TMigStatInfo%DiffusorValueMax(IStatu) .LT. SMigStatInfo%DiffusorValueMax(IStatu)) then
!                    TMigStatInfo%DiffusorValueMax(IStatu) = SMigStatInfo%DiffusorValueMax(IStatu)
!                end if
!
!                if(TMigStatInfo%RMIN(IStatu) .GT. SMigStatInfo%RMIN(IStatu) .AND. NActCountTemp(IStatu) .GT. 0) then
!                    TMigStatInfo%RMIN(IStatu) = SMigStatInfo%RMIN(IStatu)
!                end if
!
!                RAVA(IStatu) = RAVA(IStatu) + SMigStatInfo%RAVA(IStatu)*NActCountTemp(IStatu)
!            END DO
!
!            NCCount = NCCount + NActCountTemp
!
!        END ASSOCIATE
!
!    END DO
!
!    DO IStatu = 1, p_NUMBER_OF_STATU
!      if(NCCount(IStatu) .GT. 0) then
!        TMigStatInfo%RAVA(IStatu) = RAVA(IStatu)/NCCount(IStatu)
!      end if
!
!    END DO
!
!
!    END ASSOCIATE
!
!    return
!  end subroutine GetBoxesMigCoaleStat_CPU
!
!  !**************************************************************
!  subroutine GetOneBoxMigCoaleStat_CPU(IBox, Host_Boxes,Host_SimuCtrlParam,TheMigCoaleStatisticOneBox,NActCount)
!    !   ***   Purpose        : To get average size of clusters, maxma size, and so on, at current time and this single box
!    !                    IBox: the index of simulation box
!    !              Host_Boxes: the boxes information in host
!    implicit none
!    !---Dummy Vars---
!    integer, intent(in)::IBox
!    type(SimulationBoxes)::Host_Boxes
!    type(SimulationCtrlParam)::Host_SimuCtrlParam
!    type(MigCoaleStatisticOneBox)::TheMigCoaleStatisticOneBox
!    integer,optional::NActCount(p_NUMBER_OF_STATU)
!    !---Local Vars---
!    integer::IC, ICFROM, ICTO
!    integer::IStatu
!    integer::ActiveFlag,DisappearFlag
!    integer::NActCountTemp(p_NUMBER_OF_STATU)
!    !---Body---
!    TheMigCoaleStatisticOneBox%RMAX = -1.D32
!    TheMigCoaleStatisticOneBox%RMIN = 1.D32
!    TheMigCoaleStatisticOneBox%RAVA  = 0.D0
!    TheMigCoaleStatisticOneBox%DiffusorValueMax = -1.D32
!
!    ICFROM = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,1)
!    ICTO   = Host_Boxes%m_BoxesInfo%SEUsedIndexBox(IBox,2)
!
!    ActiveFlag = ICFROM
!    DisappearFlag = ICTO
!
!    if(ICTO .GT. 0) then
!      DO IC = ICFROM, ICTO
!        IStatu = Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_Statu
!
!        if(IStatu .eq. p_Empty) then
!            cycle
!        end if
!
!        if(IStatu .eq. p_ACTIVEFREE_STATU .or. IStatu .eq. p_ACTIVEINGB_STATU) then
!            ! The index for active clusters
!            Host_Boxes%m_ClustersInfo_CPU%m_ActiveIndex(ActiveFlag) = IC
!            ActiveFlag = ActiveFlag + 1
!        else
!            ! The index for unactive clusters
!            Host_Boxes%m_ClustersInfo_CPU%m_ActiveIndex(DisappearFlag) = IC
!            DisappearFlag = DisappearFlag - 1
!        end if
!
!        TheMigCoaleStatisticOneBox%RAVA(IStatu) = TheMigCoaleStatisticOneBox%RAVA(IStatu) + Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD
!
!        if(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff .GT. TheMigCoaleStatisticOneBox%DiffusorValueMax(IStatu)) then
!            TheMigCoaleStatisticOneBox%DiffusorValueMax(IStatu) = Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_DiffCoeff
!        end if
!
!        if(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD .GT. TheMigCoaleStatisticOneBox%RMAX(IStatu)) then
!            TheMigCoaleStatisticOneBox%RMAX(IStatu) = Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD
!            TheMigCoaleStatisticOneBox%ICMAX(IStatu) = IC
!        end if
!
!        if(Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD .LT. TheMigCoaleStatisticOneBox%RMIN(IStatu)) then
!            TheMigCoaleStatisticOneBox%RMIN(IStatu) = Host_Boxes%m_ClustersInfo_CPU%m_Clusters(IC)%m_RAD
!        end if
!
!      END DO
!    end if
!
!    ! move the unactive clusters indexes array position and inverse direction
!    Host_Boxes%m_ClustersInfo_CPU%m_ActiveIndex(ActiveFlag:ActiveFlag + ICTO - DisappearFlag - 1) = Host_Boxes%m_ClustersInfo_CPU%m_ActiveIndex(ICTO:DisappearFlag+1:-1)
!
!    ! The index for empty clusters space
!    !FORALL(IC=DisappearFlag:ActiveFlag:-1)
!    FORALL(IC=ActiveFlag + ICTO - DisappearFlag:ICTO)
!        Host_Boxes%m_ClustersInfo_CPU%m_ActiveIndex(IC) = IC
!    END FORALL
!
!    ! Do some check
!
!    if((DisappearFlag + 1 - ActiveFlag) .LT. 0) then
!        write(*,*)  "MFPSCUERROR: The statistic error. The active and unactive clusters index would overlap."
!        pause
!        stop
!    end if
!
!    NActCountTemp = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Single(IBox)%NC
!
!    DO IStatu = 1,p_NUMBER_OF_STATU
!        if(NActCountTemp(IStatu) .GT. 0) then
!            TheMigCoaleStatisticOneBox%RAVA(IStatu) = TheMigCoaleStatisticOneBox%RAVA(IStatu)/NActCountTemp(IStatu)
!        end if
!    END DO
!
!    if(present(NActCount)) then
!        NActCount = NActCountTemp
!    end if
!
!    return
!  end subroutine GetOneBoxMigCoaleStat_CPU


end module MIGCOALE_STATISTIC_CPU
