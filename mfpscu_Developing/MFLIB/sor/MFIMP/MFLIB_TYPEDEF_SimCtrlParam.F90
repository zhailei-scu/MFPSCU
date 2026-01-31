!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
  !*** Description: This module is created for the control parameters in simulation
  !
  !    HISTORY: Created by Zhai Lei in May, 2018
  !
  !    Reference: MDPSCU module MD_TYPEDEF_SIMULATIONCTRLPARAM
  use MCMF_CONSTANTS
  use MCMF_UTILITIES
  use MSM_TYPEDEF_InputPaser
  use MiniUtilities, only:EXTRACT_NUMB,GETINPUTSTRLINE,GETKEYWORD,UPCASE,DRSTR,ISTR

  implicit none

  character(len=7), parameter, private::m_CTLSTARTFLAG = "&MFCTLF"

  type,public::SimulationCtrlParam

     !***Run status
     integer::RESTARTAT = 0                                             !=0,   start a new session
                                                                        !=-1,  restart from the end point of the previous calculation
                                                                        !=xxx(>0), restart from a xxx time step of the previouss calculation
     !***Information about simulation boxs of jobs
     integer::MultiBox = 1                                              ! the number of simulation boxs in one run
     integer::TOTALBOX = 1                                              ! the total number of indentical boxes, this variable would be used in GPU computing
     integer::INDEPBOX = 1                                              ! if all the box are independent

     !***Information for random number
     integer(kind=KMCDF)::RANDSEED(2) = (/43434, 54454532/)             ! the inputed random seed

     !***PERIOD boundary************
     integer::PERIOD(3) = (/1,1,1/)                                     ! determine if PERIOD condition used
     integer::BDCTYPE(3,2) = p_Neumann_BDC


     !***Informaiton about temperature
     real(kind=KMCDF)::TEMP = 300.D0                                    ! temperature
     real(kind=KMCDF)::TKB = 300D0*C_KB                                 ! the kinetic energy

     !****Parameters for mean filed control***
     real(kind=KMCDF)::DumplicateFactor = 1.D-6
     real(kind=KMCDF)::MaxReactChangeRate = 0.05
     real(kind=KMCDF)::MaxDiffuseChangeRate = 0.2

     !***Implantation section*********
     integer::ImplantSectID = 0                                         ! the implantation section index

     !***Information about Time
     integer::TermTFlag = mp_TermTimeFlag_ByRealTime                    ! = 0 stans for by steps,flag = 1 by time(s)
     real::TermTValue = 3000                                            ! terminate time

     integer::UPDATETSTEPSTRATEGY = mp_SelfAdjustlStep_NearestSep       ! flag = 0 the time step is determined by average distance of nearest cluster
                                                                        ! flag = 1 is by fixed time-step
                                                                        ! flag = 2 the time step is determined by volume average distance and suppose the clusters distribute uniform in the box
     real::FixedTimeStepValue = 1                                       ! Fixed step time length for mp_FixedTimeStep strategy
     real::EnlageTStepScale = 0.01                                      ! adjustl time-step enlarge-scale mp_SelfAdjustlStep_NearestSep strategy and mp_SelfAdjustlStep_AveSep

     integer::TUpdateStatisFlag = mp_UpdateStatisFlag_ByIntervalSteps   ! flag = 0 for output each interval steps,flag = 1 for output each interval time(s)
     real::TUpdateStatisValue = 10                                      ! the time interval to update statistic

     integer::OutPutConfFlag = mp_OutTimeFlag_ByIntervalSteps           ! flag = 0 for output each interval steps,flag = 1 for output each interval time(s), flag = 2 by magnification
     real::OutPutConfValue = 100                                        ! output configuration interval value

     integer::OutPutSCFlag = mp_OutTimeFlag_ByIntervalSteps             ! flag = 0 for output each interval steps,flag = 1 for output each interval time(s), flag = 2 by magnification
     real::OutPutSCValue_IntegralBox = 100                              ! output statistic configuration interval value(for integral box)
     real::OutPutSCValue_EachBox = 100                                  ! output statistic configuration interval value(for integral box)

     integer::OutPutFuncSFlag = mp_OutTimeFlag_ByIntervalSteps          ! flag = 0 for output each interval steps,flag = 1 for output each interval time(s), flag = 2 by magnification
     real::OutPutFuncSValue = 100                                       ! output function statistic interval value

     integer::OutPutSwapFlag = mp_OutTimeFlag_ByIntervalSteps           ! flag = 0 for output each interval steps,flag = 1 for output each interval time(s), flag = 2 by magnification
     real::OutPutSwapValue = 100                                        ! output swap value interval value

     !*****restored statments for special problem*************
     type(Statementlist),pointer::AddOnData=>null()

     !***File path and name restore*****
     character*256::InputFilePath = ""
     character*256::InputFileshortName = ""
     character*256::OutFilePath = ""

     !*********Parameters for analysis processes*************
     integer::STARTJOB = -1
     integer::ENDJOB = -1
     integer::JOBSTEP = 0
     integer::STARTTSECTION = -1
     integer::ENDTSECTION = -1
     integer::TSECTIONSTEP = 0
     integer::STARTCFG = -1
     integer::ENDCFG = -1
     integer::CFGSTEP = 0
     integer::STARTBOX = -1
     integer::ENDBOX = -1
     integer::BOXSTEP = 0

     !---Determine whether the reaction are considered
     logical::FreeDiffusion = .false.

     !********************************************
     type(SimulationCtrlParam),pointer::next=>null()

     contains
     procedure,non_overridable,pass,public::AppendOne_SimulationCtrlParam
     procedure,non_overridable,pass,private::CopyFromOther
     procedure,non_overridable,nopass,private::CleanSimulationCtrlParam
     procedure,non_overridable,pass,public::DefaultValue_CtrlParam
     procedure,non_overridable,pass,public::Load_Ctrl_Parameters
     procedure,non_overridable,pass,public::Print_CtrlParameters
     procedure,non_overridable,pass,private::Load_Ctrl_CommParameter
     procedure,non_overridable,pass,private::Load_Ctrl_AnalyParameter
     procedure,non_overridable,pass,private::Load_Ctrl_SectionParameter
     procedure,non_overridable,pass,private::Load_Ctrl_Temperature
     procedure,non_overridable,pass,private::Load_Ctrl_Boundary
     procedure,non_overridable,pass,private::Load_Ctrl_Implant
     procedure,non_overridable,pass,private::Load_Ctrl_TimeStep
     procedure,non_overridable,pass,private::Load_Ctrl_MeanFiledRate
     procedure,non_overridable,pass,private::Load_AddOnDataStatments
     Generic::Assignment(=)=>CopyFromOther
     Final::Clean

  end type

  private::AppendOne_SimulationCtrlParam
  private::CopyFromOther
  private::CleanOneSimulationCtrlParam
  private::CleanSimulationCtrlParam
  private::DefaultValue_CtrlParam
  private::Load_Ctrl_Parameters
  private::Print_CtrlParameters
  private::Load_Ctrl_CommParameter
  private::Load_Ctrl_AnalyParameter
  private::Load_Ctrl_SectionParameter
  private::Load_Ctrl_Temperature
  private::Load_Ctrl_Boundary
  private::Load_Ctrl_Implant
  private::Load_Ctrl_TimeStep
  private::Load_Ctrl_MeanFiledRate
  private::Load_AddOnDataStatments
  private::Clean


  contains

  !***************************************
  subroutine AppendOne_SimulationCtrlParam(this,newOne)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam),target::this
    type(SimulationCtrlParam)::newOne
    !---Local Vars---
    type(SimulationCtrlParam),pointer::cursor=>null(),cursorP=>null()
    !---Body---
    cursor=>this%next
    cursorP=>this
    DO while(associated(cursor))
       cursor=>cursor%next
       cursorP=>cursorP%next
    END DO

    allocate(cursor)
    ! The assignment(=) had been overrided
    cursor = newOne
    cursorP%next=>cursor
    return
  end subroutine AppendOne_SimulationCtrlParam

  !****************************************************************
  subroutine CopyFromOther(this,otherOne)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam),intent(out)::this
    type(SimulationCtrlParam), intent(in)::otherOne
    !---Body----

    this%RESTARTAT = otherOne%RESTARTAT

    !***Information about simulation boxs of jobs
    this%MultiBox = otherOne%MultiBox
    this%TOTALBOX = otherOne%TOTALBOX
    this%INDEPBOX = otherOne%INDEPBOX

    !***Information for random number
    this%RANDSEED = otherOne%RANDSEED

     !***PERIOD boundary************
     this%PERIOD = otherOne%PERIOD
     this%BDCTYPE = otherOne%BDCTYPE

    !***Informaiton about temperature
    this%TEMP = otherOne%TEMP
    this%TKB = otherOne%TKB

    !---Parameters for mean filed control---
    this%DumplicateFactor = otherOne%DumplicateFactor
    this%MaxReactChangeRate = otherOne%MaxReactChangeRate
    this%MaxDiffuseChangeRate = otherOne%MaxDiffuseChangeRate

    !***Implantation sectin*********
    this%ImplantSectID = otherOne%ImplantSectID

    !***Information about Time
    this%TermTFlag = otherOne%TermTFlag
    this%TermTValue = otherOne%TermTValue

    this%UPDATETSTEPSTRATEGY = otherOne%UPDATETSTEPSTRATEGY
    this%FixedTimeStepValue = otherOne%FixedTimeStepValue
    this%EnlageTStepScale = otherOne%EnlageTStepScale

    this%TUpdateStatisFlag = otherOne%TUpdateStatisFlag
    this%TUpdateStatisValue = otherOne%TUpdateStatisValue

    this%OutPutConfFlag = otherOne%OutPutConfFlag
    this%OutPutConfValue = otherOne%OutPutConfValue

    this%OutPutSCFlag = otherOne%OutPutSCFlag
    this%OutPutSCValue_IntegralBox = otherOne%OutPutSCValue_IntegralBox
    this%OutPutSCValue_EachBox = otherOne%OutPutSCValue_EachBox

    this%OutPutFuncSFlag = otherOne%OutPutFuncSFlag
    this%OutPutFuncSValue = otherOne%OutPutFuncSValue

    this%OutPutSwapFlag = otherOne%OutPutSwapFlag
    this%OutPutSwapValue = otherOne%OutPutSwapValue

    !***File path and name restore*****
    this%InputFilePath = otherOne%InputFilePath
    this%InputFileshortName = otherOne%InputFileshortName
    this%OutFilePath = otherOne%OutFilePath

    !---Determine whether the reaction are considered
    this%FreeDiffusion = otherOne%FreeDiffusion

    !*********Parameters for analysis processes*************
    this%STARTJOB = otherOne%STARTJOB
    this%ENDJOB = otherOne%ENDJOB
    this%JOBSTEP = otherOne%JOBSTEP
    this%STARTTSECTION = otherOne%STARTTSECTION
    this%ENDTSECTION = otherOne%ENDTSECTION
    this%TSECTIONSTEP = otherOne%TSECTIONSTEP
    this%STARTCFG = otherOne%STARTCFG
    this%ENDCFG = otherOne%ENDCFG
    this%CFGSTEP = otherOne%CFGSTEP
    this%STARTBOX = otherOne%STARTBOX
    this%ENDBOX = otherOne%ENDBOX
    this%BOXSTEP = otherOne%BOXSTEP


    call Copy_StatementList(otherOne%AddOnData,this%AddOnData)

    return
  end subroutine CopyFromOther

  !****************************************************************
  subroutine CleanOneSimulationCtrlParam(this)
    implicit none
    !---Dummy Vars---
    type(SimulationCtrlParam)::this
    !---Body---
    !***Run status
     this%RESTARTAT = 0

    !***Information about simulation boxs of jobs
     this%MultiBox = 1
     this%TOTALBOX = 1
     this%INDEPBOX = 1

     !***Information for random number
     this%RANDSEED = (/43434, 54454532/)

     !***PERIOD boundary************
     this%PERIOD = (/1,1,1/)
     this%BDCTYPE = p_Neumann_BDC

     !***Informaiton about temperature
     this%TEMP = 300.D0
     this%TKB = 300D0*C_KB

     !---Parameters for mean filed control---
     this%DumplicateFactor = 1.D-6
     this%MaxReactChangeRate = 0.05
     this%MaxDiffuseChangeRate = 0.2

     !***Implantation section*********
     this%ImplantSectID = 0

     !***Information about Time
     this%TermTFlag = mp_TermTimeFlag_ByRealTime
     this%TermTValue = 3000

     this%UPDATETSTEPSTRATEGY = mp_SelfAdjustlStep_NearestSep
     this%FixedTimeStepValue = 1
     this%EnlageTStepScale = 0.01

     this%TUpdateStatisFlag = mp_UpdateStatisFlag_ByIntervalSteps
     this%TUpdateStatisValue = 10

     this%OutPutConfFlag = mp_OutTimeFlag_ByIntervalSteps
     this%OutPutConfValue = 100

     this%OutPutSCFlag = mp_OutTimeFlag_ByIntervalSteps
     this%OutPutSCValue_IntegralBox = 100
     this%OutPutSCValue_EachBox = 100

     this%OutPutFuncSFlag = mp_OutTimeFlag_ByIntervalSteps
     this%OutPutFuncSValue = 100

     this%OutPutSwapFlag = mp_OutTimeFlag_ByIntervalSteps
     this%OutPutSwapValue = 100

     !***File path and name restore*****
     this%InputFilePath = ""
     this%InputFileshortName = ""
     this%OutFilePath = ""

     !---Determine whether the reaction are considered
     this%FreeDiffusion = .false.

     !*********Parameters for analysis processes*************
     this%STARTJOB = -1
     this%ENDJOB = -1
     this%JOBSTEP = 0
     this%STARTTSECTION = -1
     this%ENDTSECTION = -1
     this%TSECTIONSTEP = 0
     this%STARTCFG = -1
     this%ENDCFG = -1
     this%CFGSTEP = 0
     this%STARTBOX = -1
     this%ENDBOX = -1
     this%BOXSTEP = 0

     call Release_StatementList(this%AddOnData)
     this%AddOnData=>null()

  end subroutine CleanOneSimulationCtrlParam


  !****************************************************************
  subroutine CleanSimulationCtrlParam(this)
    implicit none
    !---Dummy Vars---
    type(SimulationCtrlParam),pointer::this
    !---Local Vars--
    type(SimulationCtrlParam),pointer::cursorP=>null()
    type(SimulationCtrlParam),pointer::cursor=>null()
    !---Body---

    if(.not. associated(this)) then
        return
    end if

    cursorP=>this
    cursor=>this%next

    if(associated(cursor)) then
        DO While(associated(cursor))
            Nullify(cursorP%next)
            cursorP%next=>null()
            call CleanOneSimulationCtrlParam(cursorP)
            deallocate(cursorP)
            cursorP=>cursor
            cursor=>cursor%next
        END DO

        Nullify(cursor)
        cursor=>null()

    else
        call CleanOneSimulationCtrlParam(cursorP)
    end if

    Nullify(cursorP)
    cursorP=>null()

     return
  end subroutine

  !****************************************************************
  subroutine Clean(this)
    implicit none
    !---Dummy Vars---
    type(SimulationCtrlParam),target::this
    !---Local Vars---
    type(SimulationCtrlParam),pointer::ptr=>null()
    !---Body---
    ptr=>this

    call CleanSimulationCtrlParam(ptr)

    Nullify(ptr)
    ptr=>null()

    return
  end subroutine Clean


  !****************************************************************
  subroutine DefaultValue_CtrlParam(this)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    !---Body---
    !***Run status
     this%RESTARTAT = 0

    !***Information about simulation boxs of jobs
     this%MultiBox = 1
     this%TOTALBOX = 1
     this%INDEPBOX = 1

     !***Information for random number
     this%RANDSEED = (/43434, 54454532/)

     !***PERIOD boundary************
     this%PERIOD = (/1,1,1/)
     this%BDCTYPE = p_Neumann_BDC

     !***Informaiton about temperature
     this%TEMP = 300.D0
     this%TKB = 300D0*C_KB

     !---Parameters for mean filed control---
     this%DumplicateFactor = 1.D-6
     this%MaxReactChangeRate = 0.05
     this%MaxDiffuseChangeRate = 0.2

     !***Implantation section*********
     this%ImplantSectID = 0

     !***Information about Time
     this%TermTFlag = mp_TermTimeFlag_ByRealTime
     this%TermTValue = 3000

     this%UPDATETSTEPSTRATEGY = mp_SelfAdjustlStep_NearestSep
     this%FixedTimeStepValue = 1
     this%EnlageTStepScale = 0.01

     this%TUpdateStatisFlag = mp_UpdateStatisFlag_ByIntervalSteps
     this%TUpdateStatisValue = 10

     this%OutPutConfFlag = mp_OutTimeFlag_ByIntervalSteps
     this%OutPutConfValue = 100

     this%OutPutSCFlag = mp_OutTimeFlag_ByIntervalSteps
     this%OutPutSCValue_IntegralBox = 100
     this%OutPutSCValue_EachBox = 100

     this%OutPutFuncSFlag = mp_OutTimeFlag_ByIntervalSteps
     this%OutPutFuncSValue = 100

     this%OutPutSwapFlag = mp_OutTimeFlag_ByIntervalSteps
     this%OutPutSwapValue = 100

     !***File path and name restore*****
     this%InputFilePath = ""
     this%InputFileshortName = ""
     this%OutFilePath = ""


     !---Determine whether the reaction are considered
     this%FreeDiffusion = .false.

     !*********Parameters for analysis processes*************
     this%STARTJOB = -1
     this%ENDJOB = -1
     this%JOBSTEP = 0
     this%STARTTSECTION = -1
     this%ENDTSECTION = -1
     this%TSECTIONSTEP = 0
     this%STARTCFG = -1
     this%ENDCFG = -1
     this%CFGSTEP = 0
     this%STARTBOX = -1
     this%ENDBOX = -1
     this%BOXSTEP = 0

     call Release_StatementList(this%AddOnData)
     this%AddOnData=>null()

     this%next=>null()

    return
  end subroutine DefaultValue_CtrlParam

  !******************************************************
  subroutine Load_Ctrl_Parameters(this,hFile)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam),target::this
    integer, intent(in)::hFile
    !---Local Vars---
    integer::LINE
    character*256::STR
    character*32::KEYWORD
    type(SimulationCtrlParam),pointer::cursor=>null()
    type(SimulationCtrlParam)::tempCtrlParam
    integer::ISect
    !---Body---

    call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
    call RemoveComments(STR,"!")
    STR = adjustl(STR)
    call GETKEYWORD("&", STR, KEYWORD)
    call UPCASE(KEYWORD)
    if(KEYWORD(1:LENTRIM(KEYWORD)) .ne. m_CTLSTARTFLAG) then
      write(*,*) "MFPSCUERROR: The Start Flag of simulation Control Parameters is Illegal: ",KEYWORD(1:LENTRIM(KEYWORD))
      pause
      stop
    end if

    ISect = 0
    DO While(.TRUE.)
      call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
      call RemoveComments(STR,"!")
      STR = adjustl(STR)
      call GETKEYWORD("&", STR, KEYWORD)
      call UPCASE(KEYWORD)

      select case(KEYWORD(1:LENTRIM(KEYWORD)))
        case("&ENDMFCTLF")
          exit

        case("&COMMSUBCTL")
          call this%Load_Ctrl_CommParameter(hFile,*100)
          cursor=>this
          DO While(.true.)
            if(.not. associated(cursor)) then
                exit
            end if
            cursor%MultiBox = this%MultiBox
            cursor%TOTALBOX = this%TOTALBOX
            cursor%INDEPBOX = this%INDEPBOX
            cursor%RANDSEED = this%RANDSEED

            cursor=>cursor%next
          END DO

        case("&ANALYSUBCTL")
            call this%Load_Ctrl_AnalyParameter(hFile,*100)
            cursor=>this
            DO While(.true.)
                if(.not. associated(cursor)) then
                    exit
                end if
                cursor%STARTJOB = this%STARTJOB
                cursor%ENDJOB = this%ENDJOB
                cursor%JOBSTEP = this%JOBSTEP
                cursor%STARTTSECTION = this%STARTTSECTION
                cursor%ENDTSECTION = this%ENDTSECTION
                cursor%TSECTIONSTEP = this%TSECTIONSTEP
                cursor%STARTCFG = this%STARTCFG
                cursor%ENDCFG = this%ENDCFG
                cursor%CFGSTEP = this%CFGSTEP
                cursor%STARTBOX = this%STARTBOX
                cursor%ENDBOX = this%ENDBOX
                cursor%BOXSTEP = this%BOXSTEP


                cursor=>cursor%next
            END DO

        case("&SECTSUBCTL")

          if(ISect .eq. 0) then
            call this%Load_Ctrl_SectionParameter(hFile,*100)
          else
            tempCtrlParam = this
            tempCtrlParam%next=>null()
            call tempCtrlParam%Load_Ctrl_SectionParameter(hFile,*100)
            call this%AppendOne_SimulationCtrlParam(tempCtrlParam)
          end if
          ISect = ISect + 1

        case default
          write(*,*) "MCPSCU ERROR: The Illegal Flag: ",KEYWORD(1:LENTRIM(KEYWORD))
          write(*,*) "Please Check Control File at Line: ",LINE
          stop
      end select

    END DO

    return
    !-----------------------------------------------------
    100 write(*,*) "MCPSCU ERROR: Failer to read Simulation Control Parameters."
        write(*,*) "The process would be stop."
        stop
  end subroutine Load_Ctrl_Parameters

  !*****************************************
  subroutine Load_Ctrl_CommParameter(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    integer, intent(in)::hFile
    !---Local Vars---
    character*256::STR
    character*256::KEYWORD
    character*32::STRNUMB(10)
    integer::N
    integer::I
    integer::LINE
    !---Body---

    DO While(.TRUE.)
      call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
      call RemoveComments(STR,"!")
      STR = adjustl(STR)
      call GETKEYWORD("&", STR, KEYWORD)
      call UPCASE(KEYWORD)

      select case(KEYWORD(1:LENTRIM(KEYWORD)))
        case("&ENDSUBCTL")
          exit
        case("&BOX")
           call EXTRACT_NUMB(STR,3,N,STRNUMB)

           if(N .LT. 3) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for MultiBox Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&BOX box in one test = , total num = , independent ='."
             stop
           else
             this%MultiBox = ISTR(STRNUMB(1))
             this%TOTALBOX = ISTR(STRNUMB(2))
             this%INDEPBOX = ISTR(STRNUMB(3))

             if(this%MultiBox .LE. 0) then
               write(*,*) "MCPSCU ERROR: The number of box can not less than 0 .",this%MultiBox
               stop
             end if

             if(this%TOTALBOX .LT. this%MultiBox) then
               write(*,*) "MCPSCU ERROR: The number of total box should great than box in one test. "
               write(*,*) "At control file line: ",LINE
               write(*,*) "box in one test =",this%MultiBox," total num =",this%TOTALBOX
               pause
               stop
             end if
           end if
        case("&RANDSEED")

            call EXTRACT_NUMB(STR,2,N,STRNUMB)

            if(N .LT. size(this%RANDSEED)) then
                write(*,*) "MFPSCUERROR: The random seeds number shoud be at lease: ",size(this%RANDSEED)
                pause
                stop
            end if

            DO I = 1,size(this%RANDSEED)
                this%RANDSEED(I) = ISTR(STRNUMB(I))
            END DO

        case default
          write(*,*) "MCPSCU ERROR: The Illegal Flag: ",KEYWORD(1:LENTRIM(KEYWORD)),LINE
          write(*,*) "Please Check Control File at Line: ",LINE
          stop
      end select
    END DO

    return

    100 return 1
  end subroutine Load_Ctrl_CommParameter

  !*****************************************
  subroutine Load_Ctrl_AnalyParameter(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    integer, intent(in)::hFile
    !---Local Vars---
    character*256::STR
    character*256::KEYWORD
    character*32::STRNUMB(10)
    integer::N
    integer::I
    integer::LINE
    !---Body---

    DO While(.TRUE.)
      call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
      call RemoveComments(STR,"!")
      STR = adjustl(STR)
      call GETKEYWORD("&", STR, KEYWORD)
      call UPCASE(KEYWORD)

      select case(KEYWORD(1:LENTRIM(KEYWORD)))
        case("&ENDSUBCTL")
          exit
        case("&JOBSEL")
           call EXTRACT_NUMB(STR,3,N,STRNUMB)

           if(N .GT. 0) then
                this%STARTJOB = ISTR(STRNUMB(1))
           end if

           if(N .GE. 2) then
                this%ENDJOB = ISTR(STRNUMB(2))
           end if

           if(N .GE. 3) then
                this%JOBSTEP = ISTR(STRNUMB(3))
           end if

        case("&CFGSEL")
           call EXTRACT_NUMB(STR,3,N,STRNUMB)

           if(N .GT. 0) then
                this%STARTCFG = ISTR(STRNUMB(1))
           end if

           if(N .GE. 2) then
                this%ENDCFG = ISTR(STRNUMB(2))
           end if

           if(N .GE. 3) then
                this%CFGSTEP = ISTR(STRNUMB(3))
           end if

        case("&BOXSEL")
           call EXTRACT_NUMB(STR,3,N,STRNUMB)

           if(N .GT. 0) then
                this%STARTBOX = ISTR(STRNUMB(1))
           end if

           if(N .GE. 2) then
                this%ENDBOX = ISTR(STRNUMB(2))
           end if

           if(N .GE. 3) then
                this%BOXSTEP = ISTR(STRNUMB(3))
           end if

        case("&TSECTIONSEL")
           call EXTRACT_NUMB(STR,3,N,STRNUMB)

           if(N .GT. 0) then
                this%STARTTSECTION = ISTR(STRNUMB(1))
           end if

           if(N .GE. 2) then
                this%ENDTSECTION = ISTR(STRNUMB(2))
           end if

           if(N .GE. 3) then
                this%TSECTIONSTEP = ISTR(STRNUMB(3))
           end if


        case default
          write(*,*) "MCPSCU ERROR: The Illegal Flag: ",KEYWORD(1:LENTRIM(KEYWORD)),LINE
          write(*,*) "Please Check Control File at Line: ",LINE
          stop
      end select
    END DO

    return

    100 return 1
  end subroutine Load_Ctrl_AnalyParameter

  !*****************************************
  subroutine Load_Ctrl_SectionParameter(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    integer, intent(in)::hFile
    !---Local Vars---
    character*256::STR
    integer::LINE
    character*32::KEYWORD
    !---Body---

    DO while(.true.)
        call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&", STR, KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case default
                write(*,*) "MFPSCUERROR: Illegl symbol : ",KEYWORD,LINE
                pause
                stop
            case("&ENDSUBCTL")
                exit
            case("&TEMPSUBCTL")
                call this%Load_Ctrl_Temperature(hFile,*100)
            case("&BOUNDSUBCTL")
                call this%Load_Ctrl_Boundary(hFile,*100)
            case("&RATIOSUBCTL")
                call this%Load_Ctrl_MeanFiledRate(hFile,*100)
            case("&IMPLANTSUBCTL")
                call this%Load_Ctrl_Implant(hFile,*100)
            case("&TIMESUBCTL")
                call this%Load_Ctrl_TimeStep(hFile,*100)
            case("&ADDONDATA")
                call this%Load_AddOnDataStatments(hFile,*100)
        end select
    END DO

    return

    100 return 1
  end subroutine Load_Ctrl_SectionParameter

  !*****************************************
  subroutine Load_Ctrl_Temperature(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    integer, intent(in)::hFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRNUMB(10)

    DO While(.TRUE.)
      call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
      call RemoveComments(STR,"!")
      STR = adjustl(STR)
      call GETKEYWORD("&", STR, KEYWORD)
      call UPCASE(KEYWORD)

      select case(KEYWORD(1:LENTRIM(KEYWORD)))
        case("&ENDSUBCTL")
          exit
        case("&TEMPERATURE")
           call EXTRACT_NUMB(STR,1,N,STRNUMB)

           if(N .LT. 1) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for TEMPERATURE Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&TEMPERATURE SYSTEM SIMULATION TEMPERATURE = (K)'."
             stop
           else
             this%TEMP = DRSTR(STRNUMB(1))
             this%TKB = this%TEMP*C_KB

             if(this%TEMP .LE. 0) then
               write(*,*) "MCPSCU ERROR: The absoloute temperature cannot less than 0 .",this%TEMP
               stop
             end if

           end if

        case default
          write(*,*) "MCPSCU ERROR: The Illegal Flag: ",KEYWORD(1:LENTRIM(KEYWORD))
          write(*,*) "Please Check Control File at Line: ",LINE
          stop
      end select
    END DO

    return

    100 return 1
  end subroutine Load_Ctrl_Temperature


  !*****************************************
  subroutine Load_Ctrl_Boundary(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    integer, intent(in)::hFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRNUMB(10)
    integer::TotalBDC
    integer::I
    integer::J
    !---Body---
    DO while(.true.)
        call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case default
                write(*,*) "MFPSCUERROR: Illegl flag: ",KEYWORD,LINE
                pause
                stop
            case("&ENDSUBCTL")
                exit
            case("&PERIDIC")
                call EXTRACT_NUMB(STR,3,N,STRNUMB)

                if(N .LT. 3) then
                    write(*,*) "MFPSCUERROR: please special the boundary condition in three direction."
                    write(*,*) "At control file line: ",LINE
                    write(*,*) "Should be '&PERIDIC If use periodic boundary condition: X = , Y = , Z = '."
                    pause
                    stop
                end if
                this%PERIOD(1) = ISTR(STRNUMB(1))
                this%PERIOD(2) = ISTR(STRNUMB(2))
                this%PERIOD(3) = ISTR(STRNUMB(3))
            case("&BDTYPE")
                TotalBDC = 0

                call EXTRACT_SUBSTR(STR,6,N,STRNUMB)

                DO I = 1,3
                    if(this%PERIOD(I) .LE. 0) then
                        TotalBDC = TotalBDC + 2
                    end if
                END DO

                if(N .ne. TotalBDC) then
                    write(*,*) "MFPSCUERROR: You had open ",TotalBDC/2," dimension"
                    write(*,*) "However, you specialed ",N, " surface boundary condition"
                    write(*,*) "You should special ",TotalBDC," surface boundary condition"
                    pause
                    stop
                end if

                TotalBDC = 0
                DO I = 1,3
                    if(this%PERIOD(I) .LE. 0) then
                        Do J = 1,2
                            TotalBDC = TotalBDC + 1
                            if(IsStrEqual(STRNUMB(TotalBDC),p_CDirichlet_BDC)) then
                                this%BDCTYPE(I,J) = p_Dirichlet_BDC
                            else if(IsStrEqual(STRNUMB(TotalBDC),p_CNeumann_BDC)) then
                                this%BDCTYPE(I,J) = p_Neumann_BDC
                            else
                                write(*,*) "MFPSCUERROR: unknown boundary condition: ",STRNUMB(TotalBDC)
                                pause
                                stop
                            end if
                        END DO
                    end if
                END DO
        end select
    END DO

    return

    100 return 1
  end subroutine


  subroutine Load_Ctrl_MeanFiledRate(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    integer, intent(in)::hFile
    !---Local Vars---
    integer::LINE
    character*256::STR
    character*32::KEYWORD
    character*20::STRTMP(10)
    integer::N
    !---Body---

    DO While(.true.)
        call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&MAXREACTCHANGERATE")
                call EXTRACT_NUMB(STR,1,N,STRTMP)

                if(N .LT. 1) then
                    write(*,*) "MFPSCUERROR: To few parameters for mean field max change rate induced by reaction setting"
                    write(*,*) "At line: ",LINE
                    write(*,*) STR
                    pause
                    stop
                end if
                this%MaxReactChangeRate = DRSTR(STRTMP(1))

                if(this%MaxReactChangeRate .GT. 1 .or. this%MaxReactChangeRate .LT. 0.D0) then
                    write(*,*) "MFPSCUERROR: The max change rate induced by reaction should between 0 and 1."
                    write(*,*) this%MaxReactChangeRate
                    pause
                    stop
                end if

            case("&MAXDIFFUSECHANGERATE")
                call EXTRACT_NUMB(STR,1,N,STRTMP)

                if(N .LT. 1) then
                    write(*,*) "MFPSCUERROR: To few parameters for mean field max change rate induced by diffusion setting"
                    write(*,*) "At line: ",LINE
                    write(*,*) STR
                    pause
                    stop
                end if
                this%MaxDiffuseChangeRate = DRSTR(STRTMP(1))

                if(this%MaxDiffuseChangeRate .GT. 1 .or. this%MaxDiffuseChangeRate .LT. 0.D0) then
                    write(*,*) "MFPSCUERROR: The max change rate induced by diffusion should between 0 and 1."
                    write(*,*) this%MaxDiffuseChangeRate
                    pause
                    stop
                end if

            case("&DUMPLICATEFACTOR")
                call EXTRACT_NUMB(STR,1,N,STRTMP)

                if(N .LT. 1) then
                    write(*,*) "MFPSCUERROR: To few parameters for mean field clusters range dumplicate setting"
                    write(*,*) "At line: ",LINE
                    write(*,*) STR
                    pause
                    stop
                end if
                this%DumplicateFactor = DRSTR(STRTMP(1))
            case default
                write(*,*) "MCPSCU ERROR: The Illegal Flag: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "Please Check Control File at Line: ",LINE
                pause
                stop
        end select

    END DO

    return
    100 return 1
  end subroutine

  !*****************************************
  subroutine Load_Ctrl_Implant(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    integer, intent(in)::hFile
    !---Local Vars---
    character*256::STR
    character*32::KEYWORD
    character*20::SUBNUM(10)
    integer::LINE
    integer::N
    !---Body---
    DO While(.true.)
        call GETINPUTSTRLINE(hFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&IMPLANTID")
                call EXTRACT_NUMB(STR,1,N,SUBNUM)
                if(N .LT. 1) then
                    write(*,*) "MFPSCUERROR: Too few parameters for the Keyword: ",KEYWORD
                    write(*,*) "You should special: '&IMPLANTID The used implantation id in current time section = '"
                    pause
                    stop
                end if
                this%ImplantSectID = ISTR(SUBNUM(1))
            case default
                write(*,*) "MFPSCUERROR: The flag is illegal: ",KEYWORD
                pause
                stop
        end select
    END DO

    return

    100 return 1
  end subroutine Load_Ctrl_Implant

  !*****************************************
  subroutine Load_Ctrl_TimeStep(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam)::this
    integer, intent(in)::hFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRNUMB(10)

    DO While(.TRUE.)
      call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
      call RemoveComments(STR,"!")
      STR = adjustl(STR)
      call GETKEYWORD("&", STR, KEYWORD)
      call UPCASE(KEYWORD)

      select case(KEYWORD(1:LENTRIM(KEYWORD)))
        case("&ENDSUBCTL")
          exit
        case("&TERMINATE")
           call EXTRACT_NUMB(STR,2,N,STRNUMB)

           if(N .LT. 2) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for TERMINATE Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be ' &TERMINATE Maxma simulation flag =, the time ='."
             pause
             stop
           else
             this%TermTFlag = ISTR(STRNUMB(1))
             if(this%TermTFlag .ne. mp_TermTimeFlag_ByStep .AND. this%TermTFlag .ne. mp_TermTimeFlag_ByRealTime) then
               write(*,*) "MCPSCU ERROR: The TERMINATE flag cannot is not defined.",this%TermTFlag
               pause
               stop
             end if

             this%TermTValue = DRSTR(STRNUMB(2))
             if(this%TermTValue .LT. 0) then
               write(*,*) "MCPSCU ERROR: The TERMINATE value cannot less than 0.",this%TermTValue
               pause
               stop
             end if

           end if

        case("&TSTEPSTRATEGY")
           call EXTRACT_NUMB(STR,2,N,STRNUMB)

           if(N .LT. 2) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for TSTEPSTRATEGY Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&TSTEPSTRATEGY The update time-step strategy = , the correspond value = '."
             pause
             stop
           else
             this%UPDATETSTEPSTRATEGY = ISTR(STRNUMB(1))

             select case(this%UPDATETSTEPSTRATEGY)
                case(mp_SelfAdjustlStep_NearestSep)
                    this%EnlageTStepScale = DRSTR(STRNUMB(2))
                    if(this%EnlageTStepScale .LT. 0) then
                        write(*,*) "MCPSCU ERROR: The time-step-enlarge-value cannot less than 0.",this%EnlageTStepScale
                        pause
                        stop
                    end if

                case(mp_FixedTimeStep)
                    this%FixedTimeStepValue = DRSTR(STRNUMB(2))
                    if(this%FixedTimeStepValue .LT. 0) then
                        write(*,*) "MCPSCU ERROR: The fixed time-step-value cannot less than 0.",this%FixedTimeStepValue
                        pause
                        stop
                    end if

                case(mp_SelfAdjustlStep_AveSep)
                    this%EnlageTStepScale = DRSTR(STRNUMB(2))
                    if(this%EnlageTStepScale .LT. 0) then
                        write(*,*) "MCPSCU ERROR: The time-step-enlarge-value cannot less than 0.",this%EnlageTStepScale
                        pause
                        stop
                    end if
                case default
                    write(*,*) "MCPSCU ERROR: The TSTEPSTRATEGY flag cannot is not defined.",this%UPDATETSTEPSTRATEGY
                    pause
                    stop
             end select
           end if

        case("&UPDATESTATISTIC")
           call EXTRACT_NUMB(STR,2,N,STRNUMB)

           if(N .LT. 2) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for UPDATESTATISTIC Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&UPDATESTATISTIC  Use fixd step flag = , the correspond value =       '."
             pause
             stop
           else
             this%TUpdateStatisFlag = ISTR(STRNUMB(1))
             if(this%TUpdateStatisFlag .ne. mp_UpdateStatisFlag_ByIntervalSteps .AND. this%TUpdateStatisFlag .ne. mp_UpdateStatisFlag_ByIntervalRealTime) then
               write(*,*) "MCPSCU ERROR: The TUpdateStatisFlag flag cannot is not defined.",this%TUpdateStatisFlag
               pause
               stop
             end if

             this%TUpdateStatisValue = DRSTR(STRNUMB(2))
             if(this%TUpdateStatisValue .LT. 0) then
               write(*,*) "MCPSCU ERROR: The TUpdateStatis value cannot less than 0.",this%TUpdateStatisValue
               pause
               stop
             end if

           end if

        case("&OUTPUT_CONF")
           call EXTRACT_NUMB(STR,2,N,STRNUMB)

           if(N .LT. 2) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for OUTPUT_CONF Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&OUTPUT_CONF Output instant configuration flag =, the interval = '."
             pause
             stop
           else
             this%OutPutConfFlag = ISTR(STRNUMB(1))
             if(this%OutPutConfFlag .ne. mp_OutTimeFlag_ByIntervalSteps .AND. &
                this%OutPutConfFlag .ne. mp_OutTimeFlag_ByIntervalRealTime .AND. &
                this%OutPutConfFlag .ne. mp_OutTimeFlag_ByIntervalTimeMagnification) then
               write(*,*) "MCPSCU ERROR: The OUTPUT_CONF flag cannot is not defined.",this%OutPutConfFlag
               pause
               stop
             end if

             this%OutPutConfValue = DRSTR(STRNUMB(2))
             if(this%OutPutConfValue .LT. 0) then
               write(*,*) "MCPSCU ERROR: The OUTPUT_CONF value cannot less than 0.",this%OutPutConfValue
               pause
               stop
             end if

           end if

        case("&OUTPUT_SC")
           call EXTRACT_NUMB(STR,3,N,STRNUMB)

           if(N .LT. 3) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for OUTPUT_SC Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&OUTPUT_SC Output instant configuration flag =, the interval for integral box =, the interval for each box =  '."
             pause
             stop
           else
             this%OutPutSCFlag = ISTR(STRNUMB(1))
             if(this%OutPutSCFlag .ne. mp_OutTimeFlag_ByIntervalSteps .AND. &
                this%OutPutSCFlag .ne. mp_OutTimeFlag_ByIntervalRealTime .AND. &
                this%OutPutSCFlag .ne. mp_OutTimeFlag_ByIntervalTimeMagnification) then
               write(*,*) "MCPSCU ERROR: The OUTPUT_SC flag cannot is not defined.",this%OutPutSCFlag
               pause
               stop
             end if

             this%OutPutSCValue_IntegralBox = DRSTR(STRNUMB(2))
             if(this%OutPutSCValue_IntegralBox .LT. 0) then
               write(*,*) "MCPSCU ERROR: The OUTPUT_SC value cannot less than 0 for integral box.",this%OutPutSCValue_IntegralBox
               pause
               stop
             end if

            this%OutPutSCValue_EachBox = DRSTR(STRNUMB(3))
             if(this%OutPutSCValue_EachBox .LT. 0) then
               write(*,*) "MCPSCU ERROR: The OUTPUT_SC value cannot less than 0 for each box.",this%OutPutSCValue_EachBox
               pause
               stop
             end if

           end if

        case("&OUTPUT_FUNCS")
           call EXTRACT_NUMB(STR,2,N,STRNUMB)

           if(N .LT. 2) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for OUTPUT_FUNCS Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&OUTPUT_FUNCS Output instant configuration flag =, the interval = '."
             pause
             stop
           else
             this%OutPutFuncSFlag = ISTR(STRNUMB(1))
             if(this%OutPutFuncSFlag .ne. mp_OutTimeFlag_ByIntervalSteps .AND. &
                this%OutPutFuncSFlag .ne. mp_OutTimeFlag_ByIntervalRealTime .AND. &
                this%OutPutFuncSFlag .ne. mp_OutTimeFlag_ByIntervalTimeMagnification) then
               write(*,*) "MCPSCU ERROR: The OUTPUT_FUNCS flag cannot is not defined.",this%OutPutFuncSFlag
               pause
               stop
             end if

             this%OutPutFuncSValue = DRSTR(STRNUMB(2))
             if(this%OutPutFuncSValue .LT. 0) then
               write(*,*) "MCPSCU ERROR: The OUTPUT_FUNCS value cannot less than 0.",this%OutPutFuncSValue
               pause
               stop
             end if

           end if

        case("&SAVE")
           call EXTRACT_NUMB(STR,2,N,STRNUMB)

           if(N .LT. 2) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for SAVE Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&SAVE Output instant configuration flag =, the interval = '."
             pause
             stop
           else
             this%OutPutSwapFlag = ISTR(STRNUMB(1))
             if(this%OutPutSwapFlag .ne. mp_OutTimeFlag_ByIntervalSteps .AND.  &
                this%OutPutSwapFlag .ne. mp_OutTimeFlag_ByIntervalRealTime .AND. &
                this%OutPutSwapFlag .ne. mp_OutTimeFlag_ByIntervalTimeMagnification) then
               write(*,*) "MCPSCU ERROR: The SAVE flag cannot is not defined.",this%OutPutSwapFlag
               pause
               stop
             end if

             this%OutPutSwapValue = DRSTR(STRNUMB(2))
             if(this%OutPutSwapValue .LT. 0) then
               write(*,*) "MCPSCU ERROR: The SAVE value cannot less than 0.",this%OutPutSwapValue
               pause
               stop
             end if

           end if

        case default
          write(*,*) "MCPSCU ERROR: The Illegal Flag: ",KEYWORD(1:LENTRIM(KEYWORD))
          write(*,*) "Please Check Control File at Line: ",LINE
          pause
          stop
      end select
    END DO

    return

    100 return 1
  end subroutine Load_Ctrl_TimeStep

  !********************************************
  subroutine Load_AddOnDataStatments(this,hFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam),target::this
    integer,intent(in)::hFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    !---Body---
    DO While(.true.)
        call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&", STR, KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case default
                call Add_StatementList(this%AddOnData, STR, LINE)
        end select

    END DO

    return
    100 return 1
  end subroutine

  !********************************************
  subroutine Print_CtrlParameters(this,hFile)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationCtrlParam),target::this
    integer,intent(in)::hFile
    !---Local Vars---
    type(SimulationCtrlParam),pointer::cursor=>null()
    integer::ISect
    integer::I
    integer::J
    !---Body---

    write(hFile,*) "!****************Control file information***********************"
    write(hFile,fmt="('!',A70,'!',2x,I10)") "Box in one test =",this%MultiBox
    write(hFile,fmt="('!',A70,'!',2x,I10)") "Total boxes num =",this%TOTALBOX
    write(hFile,fmt="('!',A70,'!',2x,I10)") "Boxes independent =",this%INDEPBOX

    write(hFile,fmt="('!',A70,'!',2x,3I10)") "Random Seeds =",this%RANDSEED


    cursor=>this

    ISect = 1
    DO while(associated(cursor))
        write(hFile,fmt="('*******************SubSection #: ',I10)") ISect

        write(hFile,fmt="('!',A70,'!',2x,1PE10.4)") "SYSTEM SIMULATION TEMPERATURE =",cursor%TEMP

        write(hFile,fmt="('!',A70,'!',2x,3I10)") "PERIDIC condition =",cursor%PERIOD

        write(hFile,fmt="('!',A70)") "----The boundary condition----"
        DO I = 1,3
            DO J = 1,2
                if(cursor%BDCTYPE(I,J) .eq. p_Neumann_BDC) then
                    write(hFile,fmt="('!',A70)") p_CNeumann_BDC
                else if(cursor%BDCTYPE(I,J) .eq. p_Dirichlet_BDC) then
                    write(hFile,fmt="('!',A70)") p_CDirichlet_BDC
                else
                    write(*,*) "MFPSCUERROR: Unknown boundary condition ",cursor%BDCTYPE(I,J)
                    pause
                    stop
                end if
            END DO
        END DO

        !***Information about Implantation******************
        write(hFile,fmt="('!',A70,'!',2x,I10)") "The implantation section index is :", cursor%ImplantSectID

        !***Information about Time
        write(hFile,fmt="('!',A70,'!',2x,I10,2x,1PE10.4)") "Maxma simulation flag = , the time =",cursor%TermTFlag,cursor%TermTValue

        select case(this%UPDATETSTEPSTRATEGY)
            case(mp_SelfAdjustlStep_NearestSep)
                write(hFile,fmt="('!',A70,'!',2x,I10,2x,1PE10.4)") "Use Time-update step strategy =, the correspond value =", &
                                                                    cursor%UPDATETSTEPSTRATEGY,cursor%EnlageTStepScale
            case(mp_FixedTimeStep)
                write(hFile,fmt="('!',A70,'!',2x,I10,2x,1PE10.4)") "Use Time-update step strategy =, the correspond value =", &
                                                                    cursor%UPDATETSTEPSTRATEGY,cursor%FixedTimeStepValue
            case(mp_SelfAdjustlStep_AveSep)
                write(hFile,fmt="('!',A70,'!',2x,I10,2x,1PE10.4)") "Use Time-update step strategy =, the correspond value =", &
                                                                    cursor%UPDATETSTEPSTRATEGY,cursor%EnlageTStepScale
        end select

        write(hFile,fmt="('!',A70,'!',2x,I10,2x,1PE10.4)") "The update statistic frequency flag =, the correspond value = ",cursor%TUpdateStatisFlag,cursor%TUpdateStatisValue

        write(hFile,fmt="('!',A70,'!',2x,I10,2x,1PE10.4)") "Output instant configuration flag = , the interval =",cursor%OutPutConfFlag,cursor%OutPutConfValue

        write(hFile,fmt="('!',A70,'!',2x,I10,2(2x,1PE10.4))") "Output instant size statistic information flag =, the interval for integral box =, the interval for each box =",           &
                                                               cursor%OutPutSCFlag,cursor%OutPutSCValue_IntegralBox,cursor%OutPutSCValue_EachBox

        write(hFile,fmt="('!',A70,'!',2x,I10,2x,1PE10.4)") "Output instant function statistic information flag =, the interval =",cursor%OutPutFuncSFlag,cursor%OutPutFuncSValue

        write(hFile,fmt="('!',A70,'!',2x,I10,2x,1PE10.4)") "Output instant information for restart flag =,the interval =",cursor%OutPutSwapFlag,cursor%OutPutSwapValue

        write(hFile,fmt="('*******************END SubSection #: ',I10)") ISect

        cursor=>cursor%next
        ISect = ISect + 1

    END DO


    return
  end subroutine

end module MFLIB_TYPEDEF_SIMULATIONCTRLPARAM

