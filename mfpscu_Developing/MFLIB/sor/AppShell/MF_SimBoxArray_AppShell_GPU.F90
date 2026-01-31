!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MF_SimBoxArray_AppShell_GPU
    use MFLIB_GLOBAL_GPU
    use MF_MethodClass_Factory_GPU

    type(SimulationBoxes)::m_SimBoxes
    type(SimulationCtrlParam)::m_CtrlParam
    type(MFMethodClassGPU)::m_MethodClass

    contains

    !*********************************************************
    subroutine AppShell_Main_GPU(NMPI,processid,INICONFIGPROC)
        use RAND32SEEDLIB_MODULE,only:GetSeed_RAND32SEEDLIB
        use RAND32_MODULE,only:DRAND32_PUTSEED
        implicit none
        !---Dummy Vars---
        integer,intent(in)::NMPI
        integer,intent(in)::processid
        optional::INICONFIGPROC
        external::INICONFIGPROC
        interface
            subroutine INICONFIGPROC()
                implicit none
            end subroutine INICONFIGPROC
        end interface
        !---Local Vars---
        character*256::ARG
        character*256::filePath
        integer::err
        integer::arg_Num
        integer::start_Index_Dev = 0
        integer::num_use_Device = 1
        integer::ISEED0,ISEED(2)

        character(len=256)::ExePath
        character(len=256)::path
        character(len=256)::ExeName
        character(len=256)::ExePrefixName
        logical::exits
        integer::ISTAT
        integer::TestLoops
        integer::ILoop
        !-----------Body--------------

        arg_Num = COMMAND_ARGUMENT_COUNT()

        if(arg_NUM .GE. 1) THEN

            call GET_COMMAND_ARGUMENT(0,ARG)

            call GET_COMMAND_ARGUMENT(1,ARG)
            Read(ARG,fmt="(A256)") filePath

            if(arg_NUM .GE. 2) THEN
                call GET_COMMAND_ARGUMENT(2,ARG)
                Read(ARG,*) start_Index_Dev
            end if

            if(arg_NUM .GE. 3) THEN
                call GET_COMMAND_ARGUMENT(3,ARG)
                Read(ARG,*) num_use_Device
            end if
        end if

        !*********Init device setting*********************
        call Init_Device_Setting(start_Index_Dev,num_use_Device)

        !*********Create/Open log file********************
        call OpenLogFile(m_hFILELOG)

        !********Load Global vars from input file**************
        call Initialize_Global_Variables(m_CtrlParam,m_SimBoxes)

        !*******Init the simulation boxes*****************
        call m_SimBoxes%InitSimulationBox(m_CtrlParam)

        !********Init the simulation methods*******************
        call m_MethodClass%Register_Method_Class(m_AppType,m_SimBoxes,m_CtrlParam)

        ISEED0 = m_CtrlParam%RANDSEED(1)
        call GetSeed_RAND32SEEDLIB(ISEED0,ISEED(1),ISEED(2))
        ISEED0 = ISEED0 + processid - 1
        call GetSeed_RAND32SEEDLIB(ISEED0,ISEED(1),ISEED(2))
        call DRAND32_PUTSEED(ISEED)

        call Print_Global_Variables(6,m_CtrlParam,m_SimBoxes)

        if(m_CtrlParam%INDEPBOX) then
            TestLoops = m_CtrlParam%TOTALBOX/m_CtrlParam%MultiBox
        else
            TestLoops = 1
        end if

        !call m_OneStepProcudureList%AppendOne()

        DO ILoop = 1,TestLoops
            call m_MethodClass%ForOneTest(m_SimBoxes,m_CtrlParam,ILoop)
        END DO

        return
    end subroutine AppShell_Main_GPU

end module MF_SimBoxArray_AppShell_GPU
