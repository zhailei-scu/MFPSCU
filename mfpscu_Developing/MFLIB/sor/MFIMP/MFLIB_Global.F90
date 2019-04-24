module MFLIB_GLOBAL
    USE MCMF_CONSTANTS
    USE MCMF_UTILITIES
    USE MCMF_TYPEDEF_SIMULATIONCTRLPARAM
    USE MCMF_TYPEDEF_SIMULATIONBOXARRAY
    implicit none

    character*32::m_AppType = "MIGCOALE_CLUSTER_GPU"     ! the type of application

    integer::m_hFILELOG = 0

    character*7,private::MFCTRL = "&MFCTRL"

    real(kind=KMCDF)::m_Temperature = 1000                ! (K)

    real(kind=KMCDF)::m_MatrixAtomVolum = 1.5855D-23      ! (cm^3)

    real(kind=KMCDF)::m_TKB = 1000*CP_KB         ! ERG

    real(kind=KMCDF)::m_MaxChangeRate = 0.05

    real(kind=KMCDF)::m_ConCentrat0 = 0.D0 ! 5.D19        ! 1/(cm^3)

    real(kind=KMCDF)::m_TargetConCentrat = 0.D0           ! 1/(cm^3)

    real(kind=KMCDF)::m_Flux =  1.D26                    ! 1/(cm^3*s)

    integer::m_BKind = 200

    real(kind=KMCDF)::m_DumplicateFactor = 1.D-6

    real(kind=KMCDF)::m_TermTime = 3000               ! (s)

    integer::m_NNodes = 10

    real(kind=KMCDF)::m_NodeSpace = 1.D-5             ! (cm)

    real(kind=KMCDF)::m_SURFE = 1700                   ! (ERG)/cm^2

    real(kind=KMCDF)::m_RNFACTOR

    real(kind=KMCDF)::m_DIFCOESPRE(3)     ! (cm^2)/s
    real(kind=KMCDF)::m_DIFCOESES(3)      ! (ev)
    real(kind=KMCDF)::m_DIFCOES(3)        ! (cm^2)/s
    real(kind=KMCDF)::m_SURDIFPRE         ! (cm^2)/s

    character*256::SampleFileShortName = ""

    character*256::InputFilePath = ""

    character*256::OutFilePath = ""

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



    !**************************************************
    subroutine Load_SimulationParams()
        use MiniUtilities, only:EXTRACT_NUMB,EXTRACT_SUBSTR,GETINPUTSTRLINE, GETKEYWORD, UPCASE, ISTR, DRSTR
        implicit none
        !---Local Vars---
        character*256::args
        character*256::sampleFilePath
        integer::hFile
        character*256::STR
        character*32::KEYWORD
        integer::LINE
        character*32::STRNUMB(10)
        integer::N
        real(kind=KMCDF)::ATOMV0
        integer::I
        !---Body---
        LINE = 0

        if(COMMAND_ARGUMENT_COUNT() .LT. 1) then
            write(*,*) "MFPSCU ERROR: Need to input the sample File ."
            pause
            stop
        end if

        call GET_COMMAND_ARGUMENT(1,args)

        read(args,fmt="(A256)") sampleFilePath

        call resolveLongFileName(sampleFilePath,InputFilePath,SampleFileShortName)

        hFile = OpenExistedFile(sampleFilePath)

        call GETINPUTSTRLINE(hFile,STR,LINE,'!',*100)
        STR = adjustl(STR)
        call GETKEYWORD("&", STR, KEYWORD)
        call UPCASE(KEYWORD)

        if(KEYWORD(1:len_trim(KEYWORD)) .ne. MFCTRL) then
            write(*,*) "MFPSCU ERROR: the symbol is unknown"
            write(*,*) "At line: ",LINE
            write(*,*) KEYWORD
            pause
            stop
        end if

        DO While(.true.)
            call GETINPUTSTRLINE(hFile,STR,LINE,'!',*100)
            STR = adjustl(STR)
            call GETKEYWORD("&", STR, KEYWORD)
            call UPCASE(KEYWORD)

            select case(KEYWORD(1:len_trim(KEYWORD)))
                case("&ENDMFCTRL")
                    exit

                case("&TEMPERATURE")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &TEMPERATURE"
                        write(*,*) "You shoud splecial: '&TEMPERATURE The simulation temperature =  (K)' "
                        pause
                        stop
                    end if
                    m_Temperature = DRSTR(STRNUMB(1))
                    m_TKB = m_Temperature*CP_KB
                    write(*,*) "m_Temperature",m_Temperature
                    write(*,*) "m_TKB",m_TKB

                case("&MATRIXVOLUM")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &MATRIXVOLUM"
                        write(*,*) "You shoud splecial: '&MATRIXVOLUM The matrix atoms radius = (nm)' "
                        pause
                        stop
                    end if
                    m_MatrixAtomVolum = DRSTR(STRNUMB(1))*(C_NM2CM**3)
                    write(*,*) "m_MatrixAtomVolum",m_MatrixAtomVolum

                case("&MAXCHANGERATE")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &MAXCHANGERATE"
                        write(*,*) "You shoud splecial: '&MAXCHANGERATE The max change rate in each step is ' "
                        pause
                        stop
                    end if
                    m_MaxChangeRate = DRSTR(STRNUMB(1))
                    write(*,*) "m_MaxChangeRate",m_MaxChangeRate

                case("&STARTCONCENTRATE")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &STARTCONCENTRATE"
                        write(*,*) "You shoud splecial: '&STARTCONCENTRATE The start atoms concentrate is ' "
                        pause
                        stop
                    end if
                    m_ConCentrat0 = DRSTR(STRNUMB(1))
                    write(*,*) "m_ConCentrat0",m_ConCentrat0

                case("&TARGETCONCENTRATE")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &TARGETCONCENTRATE"
                        write(*,*) "You shoud splecial: '&TARGETCONCENTRATE The target atoms concentrate is ' "
                        pause
                        stop
                    end if
                    m_TargetConCentrat = DRSTR(STRNUMB(1))
                    write(*,*) "m_TargetConCentrat",m_TargetConCentrat

                case("&IMPLANTRATE")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &IMPLANTRATE"
                        write(*,*) "You shoud splecial: '&IMPLANTRATE The implant flux (1/(cm^2*s)) = ' "
                        pause
                        stop
                    end if
                    m_Flux = DRSTR(STRNUMB(1))
                    write(*,*) "m_Flux",m_Flux

                case("&MAXNUCLEATION")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &MAXNUCLEATION"
                        write(*,*) "You shoud splecial: '&MAXNUCLEATION The max atoms number for one bubble =' "
                        pause
                        stop
                    end if
                    m_BKind = ISTR(STRNUMB(1))
                    write(*,*) "m_BKind",m_BKind

                case("&TERMTIME")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &TERMTIME"
                        write(*,*) "You shoud splecial: '&TERMTIME The termated time = ' "
                        pause
                        stop
                    end if
                    m_TermTime = DRSTR(STRNUMB(1))
                    write(*,*) "m_TermTime",m_TermTime

                case("&DUMPLICATEFACTOR")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &DUMPLICATEFACTOR"
                        write(*,*) "You shoud splecial: '&DUMPLICATEFACTOR The factor for dumplicate the size is ' "
                        pause
                        stop
                    end if
                    m_DumplicateFactor = DRSTR(STRNUMB(1))
                    write(*,*) "m_DumplicateFactor",m_DumplicateFactor

                case("&NNODES")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &NNODES"
                        write(*,*) "You shoud splecial: '&NNODES The spatical nodes number = ' "
                        pause
                        stop
                    end if
                    m_NNodes = ISTR(STRNUMB(1))
                    write(*,*) "m_NNodes",m_NNodes

                case("&NODESPACE")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &NODESPACE"
                        write(*,*) "You shoud splecial: '&NODESPACE The space for each node = (nm)' "
                        pause
                        stop
                    end if
                    m_NodeSpace = DRSTR(STRNUMB(1))*C_NM2CM
                    write(*,*) "m_NodeSpace",m_NodeSpace

                case("&SURENG")
                    call EXTRACT_NUMB(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &SURENG"
                        write(*,*) "You shoud splecial: '&SURENG THE SURFACE ENERGY OF A BUBBLE = ! (ERG/CM^2))' "
                        pause
                        stop
                    end if
                    m_SURFE = DRSTR(STRNUMB(1))
                    m_RNFACTOR = 8.D0*PI*m_SURFE/(3.D0*m_TKB)
                    write(*,*) "m_SURFE",m_SURFE
                    write(*,*) "m_RNFACTOR",m_RNFACTOR

                case("&SURDIF")
                    call EXTRACT_NUMB(STR,6,N,STRNUMB)
                    if(N .LT. 6) then
                        write(*,*) "MCPSCUERROR: Too few parameters for surface diffusion oarameters at line: ",LINE
                        write(*,*) STR
                        write(*,*) "You should special: &SURDIF  THE Surface Diffusion coefficiens, prefactor (cm^2/s) and ES(ev):"
                        pause
                        stop
                    end if
                    DO I=1,3
                        m_DIFCOESPRE(I) = DRSTR(STRNUMB(2*I-1))
                        m_DIFCOESES(I) =  DRSTR(STRNUMB(2*I))
                        m_DIFCOES(I) = m_DIFCOESPRE(I)*DEXP(-m_DIFCOESES(I)*C_EV2ERG/m_TKB)
                    END DO
                    m_SURDIFPRE = (3.D0/(2.D0*PI))*(m_MatrixAtomVolum**C_FOURBYTHREE)*m_DIFCOES(1)
                    write(*,*) "m_SURDIFPRE",m_SURDIFPRE

                case("&OUTDIR")
                    call EXTRACT_SUBSTR(STR,1,N,STRNUMB)
                    if(N .LT. 1) then
                        write(*,*) "MFPSCU ERROR: Too few paramters for &OUTDIR"
                        write(*,*) "You shoud splecial: '&OUTDIR The output directorary is ' "
                        pause
                        stop
                    end if
                    OutFilePath = STRNUMB(1)
                    OutFilePath = adjustl(OutFilePath)
                    OutFilePath = CreateDataFolder(OutFilePath,InputFilePath)

                case default
                    write(*,*) "MFPSCU ERROR: Unknown flag: ",KEYWORD
                    pause
                    stop
            end select

        END DO

        m_Flux = m_Flux/(m_NodeSpace*m_NNodes)

        close(hFile)

        return

        100 write(*,*) "Load MeanFiled sample file faild!"
            write(*,*) "At line: ", LINE
            pause
            stop
    end subroutine


end module MFLIB_GLOBAL
