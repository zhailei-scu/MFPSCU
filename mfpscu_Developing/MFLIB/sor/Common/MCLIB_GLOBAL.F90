module MCLIB_GLOBAL
!  *****************************************************************
!
!  THIS MODULE DEFINE THE CONSTANTS AND FUNCTIONS USED IN KMC SIMUATION
!
!  ********************************************************************

   USE MCLIB_CONSTANTS
   USE MCLIB_UTILITIES
   USE MCLIB_TYPEDEF_SIMULATIONCTRLPARAM
   USE MCLIB_TYPEDEF_SIMULATIONBOXARRAY
   implicit none

   character*32::m_AppType = "MIGCOALE_CLUSTER_GPU"     ! the type of application

   integer::m_hFILELOG = 0

   contains

  !***********************************************
  subroutine Initialize_Global_Variables(CtrlParam,SimBoxes)
    !*** Purpose: Load the initialize parameters for simulation
    use MiniUtilities, only:EXTRACT_NUMB,EXTRACT_SUBSTR,GETINPUTSTRLINE, ISTR, DRSTR
    implicit none
    !---Dummy Vars---
    type(SimulationCtrlParam)::CtrlParam
    type(SimulationBoxes)::SimBoxes
    !---Local Vars---
    character*256::sampleFilePath,STR
    character*256::shortFileName
    integer::fileUnit,length,ISTAT
    integer::LINE
    character*32::KEYWORD
    character*256::ctlFile,boxFile,initFile,impFile,outPath
    character*256::STRTMP(5)
    integer::N
    !---Body---
    call CtrlParam%DefaultValue_CtrlParam()

    if(COMMAND_ARGUMENT_COUNT() .LT. 1) then
       write(*,*) "MCPSCU ERROR: Need to output the sample File ."
       pause
       stop
    end if

    call GET_COMMAND_ARGUMENT(1,sampleFilePath,length,STATUS=ISTAT)
    Read(sampleFilePath,fmt="(A256)") sampleFilePath
    if(ISTAT .LT. 0) then
       write(*,*) "MCPSCU: Faile to read the command arguments of sample file."
       stop
    end if

    LINE = 0

    fileUnit = OpenExistedFile(sampleFilePath)

    call resolveLongFileName(sampleFilePath,CtrlParam%InputFilePath,CtrlParam%InputFileShortName)

    call GETINPUTSTRLINE(fileUnit,STR,LINE,'!',*100)
    call RemoveComments(STR,'!')
    STR = adjustl(STR)
    call GETKEYWORD("&", STR, KEYWORD)
    call UPCASE(KEYWORD)

    select case(KEYWORD(1:LENTRIM(KEYWORD)))
        case("&START_MIGCOALE_CLUSTER_GPU")
          CtrlParam%RESTARTAT = 0
          m_AppType = KEYWORD(LENTRIM("&START_")+1:LENTRIM(KEYWORD))
        case("&RESTART_MIGCOALE_CLUSTER_GPU")
          CtrlParam%RESTARTAT = 1
          m_AppType = KEYWORD(LENTRIM("&RESTART_")+1:LENTRIM(KEYWORD))
        case default
          write(*,*) "MCPSCUERROR: Illegal flag in sample file: ",KEYWORD(1:LENTRIM(KEYWORD))
          pause
          close(fileUnit)
          stop
    end select

    DO While(.TRUE.)
      call GETINPUTSTRLINE(fileUnit,STR,LINE,'!',*100)
      call RemoveComments(STR,'!')
      STR = adjustl(STR)
      call GETKEYWORD("&", STR, KEYWORD)
      call UPCASE(KEYWORD)

      select case(KEYWORD(1:LENTRIM(KEYWORD)))
          case("&END")
            exit
          case default
            write(*,*) "MCPSCU ERROR: Illegal flag in sample file: ",KEYWORD(1:LENTRIM(KEYWORD))
            pause
            close(fileUnit)
            stop
          case("&CTLF")
            call EXTRACT_SUBSTR(STR,1,N,STRTMP)
            if(N .LT. 1) then
                write(*,*) "MCPSCUERROR: You must special the control file name or path"
                write(*,*) STR
                write(*,*) "At line: ",LINE
                pause
                stop
            end if
            ctlFile = INQUIREFILE(adjustl(trim(STRTMP(1))),CtrlParam%InputFilePath)

          case("&BOXF")
            call EXTRACT_SUBSTR(STR,1,N,STRTMP)
            if(N .LT. 1) then
                write(*,*) "MCPSCUERROR: You must special the box file name or path"
                write(*,*) STR
                write(*,*) "At line: ",LINE
                pause
                stop
            end if
            boxFile = INQUIREFILE(adjustl(trim(STRTMP(1))),CtrlParam%InputFilePath)

          case("&INIF")
            call EXTRACT_SUBSTR(STR,1,N,STRTMP)
            if(N .LT. 1) then
                write(*,*) "MCPSCUERROR: You must special the Initialization file name or path"
                write(*,*) STR
                write(*,*) "At line: ",LINE
                pause
                stop
            end if
            initFile = INQUIREFILE(adjustl(trim(STRTMP(1))),CtrlParam%InputFilePath)
            SimBoxes%IniConfig = adjustl(trim(initFile))

          case("&IMPF")
            call EXTRACT_SUBSTR(STR,1,N,STRTMP)
            if(N .LT. 1) then
                write(*,*) "MCPSCUERROR: You must special the Implantation file name or path"
                write(*,*) STR
                write(*,*) "At line: ",LINE
                pause
                stop
            end if
            impFile = INQUIREFILE(adjustl(trim(STRTMP(1))),CtrlParam%InputFilePath)
            SimBoxes%ImpFile = adjustl(trim(impFile))

          case("&COUT")
            call EXTRACT_SUBSTR(STR,1,N,STRTMP)
            outPath = STRTMP(1)
            outPath = adjustl(outPath)
            if(LENTRIM(adjustl(CtrlParam%InputFilePath)) .GT. 0) then
                CtrlParam%OutFilePath = CreateDataFolder(adjustl(trim(CtrlParam%InputFilePath))//FolderSpe//adjustl(trim(outPath)))
            else
                CtrlParam%OutFilePath = CreateDataFolder(adjustl(trim(outPath)))
            end if
      end select
    END DO

    close(fileUnit)

    fileUnit = OpenExistedFile(ctlFile)

    call CtrlParam%Load_Ctrl_Parameters(fileUnit)

    close(fileUnit)

    fileUnit = OpenExistedFile(boxFile)

    call SimBoxes%LoadParameter_SimulationBoxes(fileUnit)

    close(fileUnit)

    call CheckSimulationParamters(CtrlParam,SimBoxes)

    return
    !------------------------------------
    100 write(*,*) "MCPSCUERROR: Failed to read line: ",LINE,"in file: ",sampleFilePath
    return
  end subroutine Initialize_Global_Variables

  !**********************************************
  subroutine OpenLogFile(hFILELOG)
    implicit none
    !---Dummy Vars---
    integer,intent(inout)::hFILELOG
    !---Local Vars---
    character(len=256)::ARG
    character(len=256)::ExePath
    character(len=256)::SampleFilePath
    character(len=256)::path
    character(len=256)::ExeName
    character(len=256)::ExePrefixName
    character(len=256)::fileName
    integer::length
    integer::ISTAT
    logical::exits
    !---Body---
    call GET_COMMAND_ARGUMENT(0,ARG)

    Read(ARG,fmt="(A256)") ExePath

    call resolveLongFileName(ExePath,path,ExeName)

    call resolveExePrefixName(ExeName,ExePrefixName)

    call GET_COMMAND_ARGUMENT(1,ARG,length,STATUS=ISTAT)

    Read(ARG,fmt="(A256)") SampleFilePath
    if(ISTAT .LT. 0) then
       write(*,*) "MCPSCU: Faile to read the command arguments of sample file."
       pause
       stop
    end if

    call resolveLongFileName(SampleFilePath,path,fileName)

    if(LENTRIM(adjustl(path)) .LE. 0) then
        hFILELOG = CreateOrOpenExistedFile(ExePrefixName(1:LENTRIM(ExePrefixName))//".log","APPEND")
    else
        hFILELOG = CreateOrOpenExistedFile(path(1:LENTRIM(path))//FolderSpe//ExePrefixName(1:LENTRIM(ExePrefixName))//".log","APPEND")
    end if

    return
  end subroutine OpenLogFile

  !***********************************************
  subroutine Print_Global_Variables(hFile,CtrlParam,SimBoxes)
    implicit none
    !---Dummy Vars---
    integer,intent(in)::hFile
    type(SimulationCtrlParam)::CtrlParam
    type(SimulationBoxes)::SimBoxes

    !---Body---
    call CtrlParam%Print_CtrlParameters(hFile)

    call SimBoxes%Print_Parameter_SimulationBoxes(hFile)
    return
  end subroutine


  !***********************************************
  subroutine CheckSimulationParamters(CtrlParam,SimBoxes)
    implicit none
    !---Dummy Vars---
    type(SimulationCtrlParam),target::CtrlParam
    type(SimulationBoxes)::SimBoxes
    !---Local Vars---
    integer::ISection
    type(SimulationCtrlParam),pointer::PCtrlParam=>null()
    !---Body---
    ISection = 1

    PCtrlParam=>CtrlParam

    DO while(associated(PCtrlParam))
      !---We can exclude some time section where the reaction need not to be considered and thus the neighbor-lists can not required, then, the---
      !---calculation would be accelerated---

      if(associated(SimBoxes%ReadReactionProp_List)) then
         PCtrlParam%FreeDiffusion = SimBoxes%ReadReactionProp_List%WhetherFreeDiffusion(PCtrlParam%TKB)
      else
         PCtrlParam%FreeDiffusion = .true.
      end if

      if(PCtrlParam%FreeDiffusion .eq. .true.) then
         PCtrlParam%UPDATETSTEPSTRATEGY = mp_SelfAdjustlStep_AveSep

         write(*,*) "***********************************************************************************************"
         write(*,*) "MCPSCU Info: The Diffusion would not be considered in time section: ",ISection
         write(*,*) "MCPSCU Info: The time step strategy would be changed to by average separation."
         write(*,*) "***********************************************************************************************"
      end if

      PCtrlParam=>PCtrlParam%next

      ISection = ISection + 1
    END DO

    return
  end subroutine CheckSimulationParamters

end module MCLIB_GLOBAL
