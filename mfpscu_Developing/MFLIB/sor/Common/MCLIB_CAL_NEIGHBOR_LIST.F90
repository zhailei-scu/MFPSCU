module MCLIB_CAL_NEIGHBOR_LIST

  USE MCLIB_CONSTANTS
  USE MCLIB_TYPEDEF_SIMULATIONBOXARRAY
  USE MCLIB_TYPEDEF_SIMULATIONCTRLPARAM
  USE MCLIB_TYPEDEF_ACLUSTER
  USE MCLIB_TYPEDEF_NEIGHBOR_LIST
  #ifdef MC_PROFILING
  USE MCLIB_TimeProfile
  #endif
  implicit none

  contains

  !**********************************************************************************
  subroutine Cal_Neighbor_List_CPU(Host_Boxes,Host_SimuCtrlParam,Record,IfDirectly,RMAX)
    implicit none
    !---Dummy Vars---
    type(SimulationBoxes), intent(in)::Host_Boxes
    type(SimulationCtrlParam),intent(in)::Host_SimuCtrlParam
    CLASS(SimulationRecord)::Record
    logical,optional::IfDirectly
    real(kind=KMCDF),optional::RMAX
    !---Local Vars---
    logical::SureToUpdateNL
    integer::NAct
    !---Body---
    if(Host_SimuCtrlParam%FreeDiffusion .eq. .true.) then
        return
    end if

    SureToUpdateNL = .false.

    NAct = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEFREE_STATU) + &
           Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEINGB_STATU)

    if(present(IfDirectly)) then
        if(IfDirectly .eq. .true.) then
            SureToUpdateNL = .true.
        end if
    end if

    if(SureToUpdateNL .eq. .false.) then
        select case(Host_SimuCtrlParam%NEIGHBORUPDATESTRATEGY)
            case(mp_NEIGHBORUPDATEBYNCREMIND)
                if(dble(NAct)/dble(Record%GetLastUpdateNLNC0()) .LE. Host_SimuCtrlParam%NEIGHBORUPDATE) then
                        SureToUpdateNL = .true.
                end if
            case(mp_NEIGHBORUPDATEBYSTEP)
                if((Record%GetSimuSteps() - Record%GetLastUpdateNLTime()) .GE. Host_SimuCtrlParam%NEIGHBORUPDATE) then
                    SureToUpdateNL = .true.
                end if
        end select
    end if

    if(SureToUpdateNL .eq. .true.) then
        select case(Host_SimuCtrlParam%NEIGHBORCALWAY)
            case(mp_CalcNeighborList_NNEAREST)
                call Cal_Neighbore_Table_NNearest_CPU(Host_Boxes,Host_SimuCtrlParam)
            case(mp_CalcNeighborList_RCUT)
                if(.not. present(RMAX)) then
                    write(*,*) "MCPSCUERROR: Must special the max radius if the cut-off neighbor-list is used."
                    pause
                    stop
                end if
                call Cal_Neighbore_Table_Old_CPU(Host_Boxes,Host_SimuCtrlParam,RMAX)
            case default
                write(*,*) "MCPSCUERROR: Unknown strategy to update neighbor-list."
                pause
                stop
        end select

        call Record%SetLastUpdateNLNC0(NAct)
        call Record%SetLastUpdateNLTime((dble(Record%GetSimuSteps())))

        call Host_Boxes%m_ClustersInfo_CPU%m_list%IncreaseOneNLUpdateCount_Host()
    end if

    return
  end subroutine Cal_Neighbor_List_CPU

  ! **************************************************************
  subroutine Cal_Neighbore_Table_Old_CPU(Host_Boxes, Host_SimuCtrlParam,RMAX)
    !***             PURPOSE: to update the neighbore list of atoms
    !             Host_Boxes: the simulation boxs
    !      Host_SimuCtrlParam: the control parameters
    !
    implicit none
    !---Dummy vars---
    type(SimulationBoxes)::Host_Boxes
    type(SimulationCtrlParam),intent(in)::Host_SimuCtrlParam
    real(kind=KMCDF),intent(in)::RMAX
    !---Local Vars---
    real(kind=KMCDF)::RCUT
    integer::MultiBox
    integer::IC, ICFROM, ICTO, J, K, IW
    integer::MAXNEIGHBORNUM
    integer::IBox
    real(kind=KMCDF)::SEP(3), R2, RCUT2
    integer, dimension(:), allocatable::INDI
    real(kind=KMCDF)::tempBOXSIZE(3)
    real(kind=KMCDF)::tempHBOXSIZE(3)
    integer::tempPERIOD(3)
    integer::NAct

    #ifdef MC_PROFILING
    call Time_Start(T_Cal_Neighbore_Table_Old_CPU_Start)
    N_Invoke_Cal_Neighbor_CPU = N_Invoke_Cal_Neighbor_CPU + 1
    #endif
    !---Body---

    MAXNEIGHBORNUM = Host_SimuCtrlParam%MAXNEIGHBORNUM

    NAct = Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEFREE_STATU) + Host_Boxes%m_BoxesBasicStatistic%BoxesStatis_Integral%NC(p_ACTIVEINGB_STATU)

    RCUT = Calc_RCUT_Old(Host_SimuCtrlParam%MultiBox,Host_SimuCtrlParam%CUTREGIONEXTEND,Host_Boxes%BOXVOLUM,NAct,RMAX)

    RCUT2 = RCUT*RCUT

    call AllocateArray_Host(INDI,MAXNEIGHBORNUM,"INDI")

    MultiBox = Host_SimuCtrlParam%MultiBox

    tempBOXSIZE = Host_Boxes%BOXSIZE

    tempHBOXSIZE = Host_Boxes%HBOXSIZE

    tempPERIOD = Host_SimuCtrlParam%PERIOD


    Associate(Clusters=>Host_Boxes%m_ClustersInfo_CPU%m_Clusters)

      DO IBox = 1, MultiBox

        ICFROM = Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,1)

        ICTO   = Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2)

        DO IC=ICFROM, ICTO
          IW=0
          if(Clusters(IC)%m_Statu .ne. p_ACTIVEFREE_STATU .AND. Clusters(IC)%m_Statu .ne. p_ACTIVEINGB_STATU) cycle

          DO  J=IC+1,ICTO
             if(Clusters(J)%m_Statu .ne. p_ACTIVEFREE_STATU .AND. Clusters(J)%m_Statu .ne. p_ACTIVEINGB_STATU) cycle

             SEP = Clusters(IC)%m_POS - Clusters(J)%m_POS
             DO K=1, 3
                if(DABS(SEP(K)) .GT. tempHBOXSIZE(K) .AND. tempPERIOD(K)) then
                   SEP(K) = SEP(K)-DSIGN(tempBOXSIZE(K),SEP(K))
                end if
             END DO

             R2=SUM(SEP*SEP)
             if(R2.GE.RCUT2) cycle

             IW=IW+1

             if(IW .GT. MAXNEIGHBORNUM) then
                write(*,*) "MCPSCU Error: Number of neighbores larger than permitted value:",MAXNEIGHBORNUM
                write(*,*) "Process to be stoped"
                stop
             end if

             INDI(IW)=J
          END DO

          Host_Boxes%m_ClustersInfo_CPU%m_list%m_KVOIS(IC)=IW
          if(IW .GT. 0) then
             Host_Boxes%m_ClustersInfo_CPU%m_list%m_INDI(IC,1:IW) = INDI(1:IW)
          end if

        END DO

      END DO

    END Associate

    call Host_Boxes%m_ClustersInfo_CPU%m_list%IncreaseOneNLUpdateCount_Host()

    call DeAllocateArray_Host(INDI,"INDI")

    #ifdef MC_PROFILING
    call Time_Accumulate(T_Cal_Neighbore_Table_Old_CPU_Start,T_Cal_Neighbore_Table_Old_CPU)
    #endif

    return
  end subroutine Cal_Neighbore_Table_Old_CPU

  ! **************************************************************
  subroutine Cal_Neighbore_Table_NNearest_CPU(Host_Boxes,Host_SimuCtrlParam)
    !***PURPOSE: to update the neighbore list of atoms(use NNearest)
    !   INPUT:            Host_Boxes       , the simulation boxs
    !            Host_SimuCtrlParam         , the control parameters
    !
    implicit none
    !---Dummy Vars---
    type(SimulationBoxes)::Host_Boxes
    type(SimulationCtrlParam),intent(in)::Host_SimuCtrlParam
    !---Local Vars---
    integer::MultiBox
    integer::NNearest
    integer::IBox
    integer::IC, ICFROM, ICTO, J, K
    integer::NN
    real(kind=KMCDF)::SEP(3), R2
    integer::Nearest_maxDistIndex
    real(kind=KMCDF)::Nearest_maxDIST2,temp_DIST2
    integer, dimension(:),allocatable::Nearest_INDI
    real(kind=KMCDF),dimension(:),allocatable::Nearest_DIST2
    real(kind=KMCDF)::tempBOXSIZE(3)
    real(kind=KMCDF)::tempHBOXSIZE(3)
    integer::tempPERIOD(3)
    !---Body---
    #ifdef MC_PROFILING
    call Time_Start(T_Cal_Neighbore_Table_NNearest_CPU_Start)
    N_Invoke_Cal_Neighbor_CPU = N_Invoke_Cal_Neighbor_CPU + 1
    #endif

    MultiBox = Host_SimuCtrlParam%MultiBox

    NNearest = Host_SimuCtrlParam%MAXNEIGHBORNUM

    tempBOXSIZE = Host_Boxes%BOXSIZE

    tempHBOXSIZE = Host_Boxes%HBOXSIZE

    tempPERIOD = Host_SimuCtrlParam%PERIOD

    call AllocateArray_Host(Nearest_INDI,NNearest,"Nearest_INDI")

    call AllocateArray_Host(Nearest_DIST2,NNearest,"Nearest_DIST2")

    Associate(Clusters=>Host_Boxes%m_ClustersInfo_CPU%m_Clusters)

    !---Body---
    DO IBox = 1, MultiBox

      ICFROM = Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox, 1)

      ICTO   = Host_Boxes%m_BoxesInfo%SEExpdIndexBox(IBox,2)

      DO IC = ICFROM, ICTO
        if(Clusters(IC)%m_Statu .ne. p_ACTIVEFREE_STATU .AND. Clusters(IC)%m_Statu .ne. p_ACTIVEINGB_STATU) cycle

         NN = 0
         Nearest_maxDIST2 = -1.D0
         DO  J=ICFROM, ICTO   !mutual
             if(IC .eq. J .or. (Clusters(J)%m_Statu .ne. p_ACTIVEFREE_STATU .AND. Clusters(J)%m_Statu .ne. p_ACTIVEINGB_STATU)) cycle

             SEP = Clusters(IC)%m_POS - Clusters(J)%m_POS
             DO K=1, 3
                IF(DABS(SEP(K)) .GT. tempHBOXSIZE(K) .AND. tempPERIOD(K)) THEN
                   SEP(K) = SEP(K)-DSIGN(tempBOXSIZE(K),SEP(K))
                END IF
             END DO

             R2=SUM(SEP*SEP)

             if(NN .LT. NNearest) then

                  NN = NN + 1
                  Nearest_INDI(NN) = J
                  Nearest_DIST2(NN) = R2

                  if(Nearest_maxDIST2 .LT. R2) then

                    Nearest_maxDIST2 = R2
                    Nearest_maxDISTIndex = NN
                  end if

             else if(Nearest_maxDIST2 .GT. R2) then

                  Nearest_INDI(Nearest_maxDISTIndex) = J
                  Nearest_DIST2(Nearest_maxDISTIndex) = R2
                  Nearest_maxDIST2 = -1.0
                  DO K = 1,NNearest
                     temp_DIST2 = Nearest_DIST2(K)
                     if(temp_DIST2 .GT. Nearest_maxDIST2) then

                       Nearest_maxDIST2 = temp_DIST2
                       Nearest_maxDISTIndex = K
                     end if
                  END DO
             end if

         END DO

         Host_Boxes%m_ClustersInfo_CPU%m_list%m_KVOIS(IC)=NN

         Host_Boxes%m_ClustersInfo_CPU%m_list%m_INDI(IC,1:NN) = Nearest_INDI(1:NN)
      END DO

    END DO

    END Associate

    call Host_Boxes%m_ClustersInfo_CPU%m_list%IncreaseOneNLUpdateCount_Host()

    call DeAllocateArray_Host(Nearest_INDI,"Nearest_INDI")

    call DeAllocateArray_Host(Nearest_DIST2,"Nearest_DIST2")

    #ifdef MC_PROFILING
    call Time_Accumulate(T_Cal_Neighbore_Table_NNearest_CPU_Start,T_Cal_Neighbore_Table_NNearest_CPU)
    #endif

    RETURN
  end subroutine Cal_Neighbore_Table_NNearest_CPU

end module MCLIB_CAL_NEIGHBOR_LIST
