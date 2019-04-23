module MFLIB_GLOBAL
    use MFLIB_CONSTANTS
    use MFLIB_UTILITIES
    implicit none

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
                        write(*,*) "You shoud splecial: '&IMPLANTRATE The implant atoms number per seconds unit volum = ' "
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

        close(hFile)

        return

        100 write(*,*) "Load MeanFiled sample file faild!"
            write(*,*) "At line: ", LINE
            pause
            stop
    end subroutine


end module MFLIB_GLOBAL
