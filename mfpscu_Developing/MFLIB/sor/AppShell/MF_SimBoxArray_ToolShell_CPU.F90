module MF_SimBoxArray_ToolShell
!    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
!    use MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
!
!    implicit none
!    type(SimulationBoxes)::m_SimBoxes
!    type(SimulationCtrlParam)::m_CtrlParam
!
!    contains
!
!    !*********************************************************
!    subroutine Main_ANALYSIS(NMPI,processid,INICONFIGPROC)
!        implicit none
!        !---Dummy Vars---
!        integer,intent(in)::NMPI
!        integer,intent(in)::processid
!        optional::INICONFIGPROC
!        external::INICONFIGPROC
!        interface
!            subroutine INICONFIGPROC()
!                implicit none
!            end subroutine INICONFIGPROC
!        end interface
!        !---Local Vars---
!!        character*256::ARG
!!        character*256::filePath
!!        integer::err
!!        integer::arg_Num
!!        integer::start_Index_Dev = 0
!!        integer::num_use_Device = 1
!!        integer::TestLoops
!!        integer::ILoop
!!        !-----------Body--------------
!!
!!        arg_Num = COMMAND_ARGUMENT_COUNT()
!!
!!        if(arg_NUM .GE. 1) THEN
!!
!!            call GET_COMMAND_ARGUMENT(0,ARG)
!!
!!            call GET_COMMAND_ARGUMENT(1,ARG)
!!            Read(ARG,fmt="(A256)") filePath
!!        end if
!!
!!        !*********Create/Open log file********************
!!        call OpenLogFile(m_hFILELOG)
!!
!!        !********Load Global vars from input file**************
!!        call Initialize_Global_Variables(m_CtrlParam,m_SimBoxes)
!!
!!        !*******Init the simulation boxes*****************
!!        call m_SimBoxes%InitSimulationBox(m_CtrlParam)
!!
!!        call Print_Global_Variables(6,m_CtrlParam,m_SimBoxes)
!
!        return
!    end subroutine Main_ANALYSIS
!
end module MF_SimBoxArray_ToolShell
