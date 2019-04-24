module MCMF_TYPEDEF_SIMULATIONBOXARRAY
  use MCMF_CONSTANTS
  use MCMF_TYPEDEF_ATOMSLIST
  use MCMF_TYPEDEF_SIMULATIONCTRLPARAM
  use MCMF_TYPEDEF_USUAL
  use MCMF_TYPEDEF_DiffusorPropList
  use MCMF_TYPEDEF_ReactionPropList
  use MCMF_UTILITIES
  use MCMF_TYPEDEF_GEOMETRY
  use MiniUtilities, only:EXTRACT_NUMB,EXTRACT_SUBSTR,GETINPUTSTRLINE, GETKEYWORD, UPCASE, ISTR, DRSTR
  implicit none

  character(len=5), parameter, private::m_BOXSTARTFLAG = "&BOXF"


  type,public::SimulationBoxes

    !***********Diffusor list memory mapping*************
    type(DiffusorTypesMap)::m_DiffusorTypesMap

    !**********Reaction list memory mapping**************
    type(ReactionsMap)::m_ReactionsMap

    !************Info about Geometry**********
    type(GrainBoundary)::m_GrainBoundary

    !************Info for boxsize*************
    real(kind=KMCDF)::LatticeLength = 3.14D-8                          ! Lattice length (cm)
    real(kind=KMCDF)::BOXBOUNDARY(3,2) = 0                             ! siumulation box boundary, in unit of atomic radiua
    real(kind=KMCDF)::BOXSIZE(3) = 0                                   ! simulation box size
    real(kind=KMCDF)::HBOXSIZE(3) = 0                                  ! half box size
    real(kind=KMCDF)::BOXVOLUM = 0                                     ! voulme of the box

    !************Info for matrix*************
    type(ATOM)::MatrixAtom

    !***********Info for atoms***************
    type(AtomsList),pointer::Atoms_list=>null()

    !**************Init file****************
    character*256::IniConfig = ""

    !**********Implantation file************
    character*256::ImpFile = ""

    !***********Info for diffusor************
    type(ReadDiffusorPropList),pointer::ReadDiffusorProp_List=>null()

    !**********Info for reactions************
    type(ReadReactionPropList),pointer::ReadReactionProp_List=>null()

    contains

    procedure,non_overridable,public,pass::DefaultValueSimulationBoxes=>DefaultValue_SimulationBoxes
    procedure,non_overridable,public,pass::LoadParameter_SimulationBoxes=>Load_Parameter_SimulationBoxes
    procedure,non_overridable,public,pass::Print_Parameter_SimulationBoxes
    procedure,non_overridable,private,pass::Load_Box_Shape
    procedure,non_overridable,private,pass::Load_Box_AtomsDefine
    procedure,non_overridable,private,pass::Load_OneSecton_AtomDefine
    procedure,non_overridable,private,pass::Load_Box_Diffusors
    procedure,non_overridable,private,pass::LoadDiffusorsValue
    procedure,non_overridable,private,pass::LoadOneDiffusors
    procedure,non_overridable,private,pass::LoadDiffusorsValueFromScript
    procedure,non_overridable,private,pass::LoadReactions
    procedure,non_overridable,private,pass::LoadOneReaction
    procedure,non_overridable,private,pass::LoadReactionsFromScript
    procedure,non_overridable,private,pass::Load_Box_GrainBoundary
    procedure,non_overridable,private,pass::Load_GB_Simple
    procedure,non_overridable,private,pass::Load_GB_Simple_Distribution_ByGSeedCtl
    procedure,non_overridable,private,pass::Load_GB_Simple_Distribution_ByGVolumCtl
    procedure,non_overridable,private,pass::Load_GB_SpecialDistFromFile
    procedure,non_overridable,private,pass::Load_GB_SpecialDistFromExteFunc
    procedure,non_overridable,public,pass::InitSimulationBox=>Init_SimulationBox
    procedure,non_overridable,private,pass::CopySimulationBoxesFromOther
    Generic::Assignment(=)=>CopySimulationBoxesFromOther
    procedure,non_overridable,public,pass::Clean=>CleanSimulationBoxes
    Final::DestorySimulationBoxes

  end type SimulationBoxes

  private::DefaultValue_SimulationBoxes
  private::Load_Parameter_SimulationBoxes
  private::Print_Parameter_SimulationBoxes
  private::Load_Box_Shape
  private::Load_Box_AtomsDefine
  private::Load_OneSecton_AtomDefine
  private::Load_Box_Diffusors
  private::LoadDiffusorsValue
  private::LoadOneDiffusors
  private::LoadDiffusorsValueFromScript
  private::LoadReactions
  private::LoadOneReaction
  private::LoadReactionsFromScript
  private::Load_Box_GrainBoundary
  private::Load_GB_Simple
  private::Load_GB_Simple_Distribution_ByGSeedCtl
  private::Load_GB_Simple_Distribution_ByGVolumCtl
  private::Load_GB_SpecialDistFromFile
  private::Load_GB_SpecialDistFromExteFunc
  private::Init_SimulationBox
  private::CleanSimulationBoxes
  private::DestorySimulationBoxes
  private::CopySimulationBoxesFromOther

  contains

  !***************************************************************
  subroutine Init_SimulationBox(this,Host_SimuCtrlParam)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    type(SimulationCtrlParam)::Host_SimuCtrlParam
    !---Local Vars---
    !---Body---

    call this%m_GrainBoundary%ConstructGrainBoundary(this%BOXBOUNDARY,Host_SimuCtrlParam)

    return
  end subroutine

  !****************************************************************
  subroutine DefaultValue_SimulationBoxes(this)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this

    !***Box Shape
    this%BOXBOUNDARY = 0
    this%BOXSIZE  = 0
    this%HBOXSIZE = 0
    this%BOXVOLUM    = 0

    !***Peridic boundary
    this%IniConfig = ""

    this%ImpFile = ""

    return
  end subroutine DefaultValue_SimulationBoxes

  !*****************************************************************
  subroutine CopySimulationBoxesFromOther(this,Other)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes),intent(out)::this
    type(SimulationBoxes),intent(in)::Other
    !---Body---

    call DestorySimulationBoxes(this)

    ! The Assignment had been override
    this%m_DiffusorTypesMap = other%m_DiffusorTypesMap

    ! The Assignment had been override
    this%m_ReactionsMap = other%m_ReactionsMap

    ! The Assignment had been override
    this%m_GrainBoundary = Other%m_GrainBoundary

    this%LatticeLength = Other%LatticeLength
    this%BOXBOUNDARY = Other%BOXBOUNDARY
    this%BOXSIZE = Other%BOXSIZE
    this%HBOXSIZE = Other%HBOXSIZE
    this%BOXVOLUM = Other%BOXVOLUM

    ! The Assignment(=) had been override
    this%MatrixAtom = Other%MatrixAtom

    if(associated(this%Atoms_list)) then
        ! The Assignment(=) had been override
        this%Atoms_list = Other%Atoms_list
    end if

    this%IniConfig = Other%IniConfig

    this%ImpFile = Other%ImpFile

    ! The Assignment(=) had been override
    if(associated(this%ReadDiffusorProp_List)) then
        this%ReadDiffusorProp_List = Other%ReadDiffusorProp_List
    end if

    ! The Assignment(=) had been override
    if(associated(this%ReadReactionProp_List)) then
        this%ReadReactionProp_List = Other%ReadReactionProp_List
    end if

    return
  end subroutine

  !*****************************************************************
  subroutine Load_Parameter_SimulationBoxes(this,hBoxFile)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---

    integer::LINE
    character*256::STR
    character*32::KEYWORD
    !---Body---

    call GETINPUTSTRLINE(hBoxFile,STR, LINE, "!", *100)
    call RemoveComments(STR,"!")
    STR = adjustl(STR)
    call GETKEYWORD("&", STR, KEYWORD)
    call UPCASE(KEYWORD)
    if(KEYWORD(1:LENTRIM(KEYWORD)) .ne. m_BOXSTARTFLAG) then
      write(*,*) "MCPSCUERROR: The Start Flag of simulation box Parameters is Illegal: ",KEYWORD(1:LENTRIM(KEYWORD))
      pause
      stop
    end if

    DO While(.TRUE.)
      call GETINPUTSTRLINE(hBoxFile,STR, LINE, "!", *100)
      call RemoveComments(STR,"!")
      STR = adjustl(STR)
      call GETKEYWORD("&", STR, KEYWORD)
      call UPCASE(KEYWORD)

      select case(KEYWORD(1:LENTRIM(KEYWORD)))
        case("&ENDBOXF")
            exit
        case("&BOXSUBCTL")
          call this%Load_Box_Shape(hBoxFile,*100)
        case("&ATOMSUBCTL")
          call this%Load_Box_AtomsDefine(hBoxFile,*100)
        case("&DIFFUSORSUBCTL")
          call this%Load_Box_Diffusors(hBoxFile,*100)
        case("&GBSUBCTL")
          call this%Load_Box_GrainBoundary(hBoxFile,*100)
        case default
          write(*,*) "MCPSCU ERROR: The Illegal Flag: ",KEYWORD(1:LENTRIM(KEYWORD))
          write(*,*) "Please Check Box File at Line: ",LINE
          pause
          stop
      end select

    END DO

    return
    !-----------------------------------------------------
    100 write(*,*) "MCPSCU ERROR: Failer to read Simulation box Parameters."
        write(*,*) "The process would be stop."
        stop
  end subroutine Load_Parameter_SimulationBoxes

  !****************************************
  subroutine Print_Parameter_SimulationBoxes(this,hFile)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hFile
    !---Local Vars---
    type(AtomsList),pointer::atomsListCursor=>null()
    type(ReadDiffusorPropList),pointer::diffusorListCursor=>null()
    !---Body---

    write(hFile,*) "!****************Siumation Boxes file information***********************"
    write(hFile,fmt="('!',A70,'!',2x,1PE10.4)") "Simulation Box lattice length(cm) =",this%LatticeLength
    write(hFile,fmt="('!',A70,'!',2x,3(1PE10.4,2x))") "Simulation Box size(cm) =",this%BOXSIZE

    atomsListCursor=>this%Atoms_list

    if(atomsListCursor%Get_ListCount() .GT. 0) then
        DO While(associated(atomsListCursor))
            write(hFile,fmt="('!',A30,I10,2x,'!',A30,I10)") "------The atom inner index = ",atomsListCursor%m_Atom%m_ID,&
                                                            "------The number of atoms = ",atomsListCursor%m_AtomsNumber
            write(hFile,fmt="('!',A30,A10,2x,'!',A30,I10,2(2x,A30,1PE10.4))") "The atom symbol = ",atomsListCursor%m_Atom%m_Symbol, &
                                                                             "The element index = ",atomsListCursor%m_Atom%m_ElementIndex, &
                                                                             "The element mass = ",atomsListCursor%m_Atom%m_AtomMass, &
                                                                             "The atomic m_Volum(cm^3) = ",atomsListCursor%m_Atom%m_Volum
            atomsListCursor=>atomsListCursor%next
        END DO
    end if

    write(hFile,fmt="('!',A70,'!',2x,A10)") "Simulation Box Matrix symbol =",this%MatrixAtom%m_Symbol

    write(hFile,fmt="('!',A70,'!',2x,1PE10.4)") "Simulation Box Matrix atom volum(cm^3) =",this%MatrixAtom%m_Volum


    diffusorListCursor=>this%ReadDiffusorProp_List
    if(diffusorListCursor%GetList_Count() .GT. 0) then
        Do While(associated(diffusorListCursor))

            write(hFile,fmt="('!','The diffusor symbol = ',A20,2x,  &
                              '!','CoefficentsGenerate way in free matrix =',I1,2x, &
                              '!','DiffusionCiefficents value in free matrix =',1PE10.4,2x, &
                              '!','PreFactor in free matrix = ',1PE10.4,2x, &
                              '!','ActEnergy in free matrix = ',1PE10.4,2x, &
                              '!','ECR Generate way in free matrix = ',I1,2x, &
                              '!','ECR Value in free matrix = ',1PE10.4,&
                              '!','CoefficentsGenerate way in GB =',I1,2x, &
                              '!','DiffusionCiefficents value in GB =',1PE10.4,2x, &
                              '!','PreFactor in GB = ',1PE10.4,2x, &
                              '!','ActEnergy in GB = ',1PE10.4,2x, &
                              '!','ECR Generate way in GB = ',I1,2x, &
                              '!','ECR Value in GB = ',1PE10.4)")              diffusorListCursor%Diffusor%symbol, &
                                                                               diffusorListCursor%Diffusor%DiffusorValueType_Free, &
                                                                               diffusorListCursor%Diffusor%DiffuseCoefficient_Free_Value,  &
                                                                               diffusorListCursor%Diffusor%PreFactor_Free, &
                                                                               diffusorListCursor%Diffusor%ActEnergy_Free, &
                                                                               diffusorListCursor%Diffusor%ECRValueType_Free, &
                                                                               diffusorListCursor%Diffusor%ECR_Free,&
                                                                               diffusorListCursor%Diffusor%DiffusorValueType_InGB, &
                                                                               diffusorListCursor%Diffusor%DiffuseCoefficient_InGB_Value,  &
                                                                               diffusorListCursor%Diffusor%PreFactor_InGB, &
                                                                               diffusorListCursor%Diffusor%ActEnergy_InGB, &
                                                                               diffusorListCursor%Diffusor%ECRValueType_InGB, &
                                                                               diffusorListCursor%Diffusor%ECR_InGB
            diffusorListCursor=>diffusorListCursor%next
        End Do
    end if

    !---Check the diffusorList---

    write(*,*) "**************************************************************************************************"
    write(*,*) "*                                                                                                *"
    write(*,*) "***********************Start to Check The diffusors map*******************************************"
    write(*,*) "*                                                                                                *"
    write(*,*) "**************************************************************************************************"
    call this%ReadDiffusorProp_List%PrintOutCheckingResult(hFile,this%Atoms_list,this%m_DiffusorTypesMap)

    !---Check the ReactionList---
    write(*,*) "**************************************************************************************************"
    write(*,*) "*                                                                                                *"
    write(*,*) "***********************Start to Check The reactions map*******************************************"
    write(*,*) "*                                                                                                *"
    write(*,*) "**************************************************************************************************"
    call this%ReadReactionProp_List%PrintOutCheckingResult(hFile,this%Atoms_list,this%m_ReactionsMap)

    Nullify(atomsListCursor)
    atomsListCursor=>null()

    Nullify(diffusorListCursor)
    diffusorListCursor=>null()

    return
  end subroutine Print_Parameter_SimulationBoxes

  !*****************************************
  subroutine Load_Box_Shape(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    integer::LINE
    integer::I
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRNUMB(10)
    real(kind=KMCDF)::BOXSIZE(3)
    !---Body---

    DO While(.TRUE.)
      call GETINPUTSTRLINE(hBoxFile,STR, LINE, "!", *100)
      call RemoveComments(STR,"!")
      STR = adjustl(STR)
      call GETKEYWORD("&", STR, KEYWORD)
      call UPCASE(KEYWORD)

      select case(KEYWORD(1:LENTRIM(KEYWORD)))
        case("&ENDSUBCTL")
          exit
        case("&SIZE")
           call EXTRACT_NUMB(STR,3,N,STRNUMB)

           if(N .LT. 3) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for BOXSIZE Setting."
             write(*,*) "At control file line: ",LINE
             write(*,*) "Should be '&SIZE   bx(LU)= ,  by(LU) = , bz(LU) = '."
             stop
           else
             BOXSIZE(1) = DRSTR(STRNUMB(1))
             BOXSIZE(2) = DRSTR(STRNUMB(2))
             BOXSIZE(3) = DRSTR(STRNUMB(3))

             if(any(BOXSIZE .LT. 0)) then
               write(*,*) "MCPSCU ERROR: The value of BOXSIZE can not less than 0 .",this%BOXSIZE
               stop
             end if

           end if

        case("&LATT")
           call EXTRACT_NUMB(STR,1,N,STRNUMB)

           if(N .LT. 1) then
             write(*,*) "MCPSCU ERROR: Too Few Parameters for LATT Setting."
             write(*,*) "At box file line: ",LINE
             write(*,*) "Should be '&LATT latiice constant(nm) = '."
             stop
           else
             this%LatticeLength = DRSTR(STRNUMB(1))*C_NM2CM

           end if

        case default
          write(*,*) "MCPSCU ERROR: The Illegal Flag: ",KEYWORD(1:LENTRIM(KEYWORD))
          write(*,*) "Please Check box File at Line: ",LINE
          stop
      end select
    END DO

    this%BOXSIZE(1) = this%LatticeLength*BOXSIZE(1)
    this%BOXSIZE(2) = this%LatticeLength*BOXSIZE(2)
    this%BOXSIZE(3) = this%LatticeLength*BOXSIZE(3)

    DO I = 1,3
        this%HBOXSIZE(I) = 0.5*this%BOXSIZE(I)
        this%BOXBOUNDARY(I,1) = -0.5*this%BOXSIZE(I)
        this%BOXBOUNDARY(I,2) = 0.5*this%BOXSIZE(I)
    END DO

    this%BOXVOLUM = this%BOXSIZE(1)*this%BOXSIZE(2)*this%BOXSIZE(3)

    return

    100 return 1
  end subroutine Load_Box_Shape

  !*****************************************
  subroutine Load_Box_AtomsDefine(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRNUMB(10)
    !---Body---

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&GROUPSUBCTL")
                call this%Load_OneSecton_AtomDefine(hBoxFile,*100)
            case default
                write(*,*) "MCPSCUERROR: The illegal flag: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "Please check box file at Line: ",LINE
                pause
                stop
        end select


    END DO

    return

    100 return 1

  end subroutine Load_Box_AtomsDefine


  !*****************************************
  subroutine Load_OneSecton_AtomDefine(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    character*256::STR
    character*32::KEYWORD
    character*32::STRNUMB(10)
    type(ATOM)::tempAtom
    integer::AtomNumb
    integer::N
    integer::LINE
    logical::isMatrixAtom
    !---Body---

    AtomNumb = 0

    isMatrixAtom = .FALSE.

    call tempAtom%CleanAtom()

    if(.not. associated(this%Atoms_list)) then
        allocate(this%Atoms_list)
    end if

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case ("&ENDSUBCTL")
                exit
            case("&NATOM")
                call EXTRACT_NUMB(STR,1,N,STRNUMB)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR : Too few parameters for Atoms number"
                    write(*,*) "You should special: '&NATOM    the number of atoms in the group = 1' "
                    pause
                    stop
                else
                    AtomNumb = ISTR(STRNUMB(1))
                end if
            case("&ATOMP")
                call EXTRACT_SUBSTR(STR,1,N,STRNUMB)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR : You must define the Atom symbol"
                    write(*,*) STR
                    write(*,*) "You should special: '&ATOMP atomic symbol = (symbol), element index = ,  atomic mass= ' "
                    pause
                    stop
                end if

                tempAtom%m_Symbol =  trim(adjustl(STRNUMB(1)))

                call UPCASE(tempAtom%m_Symbol)

                call EXTRACT_NUMB(STR,2,N,STRNUMB)
                if(N .LT. 2) then
                    write(*,*) "MCPSCUERROR : Too few parameters for Atoms define"
                    write(*,*) STR
                    write(*,*) "You should special: '&ATOMP    atomic symbol = (symbol), element index = ,  atomic mass= ' "
                    pause
                    stop
                end if

                tempAtom%m_ElementIndex = ISTR(STRNUMB(1))
                tempAtom%m_AtomMass = DRSTR(STRNUMB(2))

            case("&ATOMVOLUM")
                call EXTRACT_NUMB(STR,1,N,STRNUMB)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: Too few parameters for matrix atom,you should specify the volum(in nm^3)"
                    write(*,*) STR
                    write(*,*) "You should specify : '&ATOMVOLUM   Volum of matrix atom (in nm^3) = ' "
                    pause
                    stop
                end if

                tempAtom%m_Volum = DRSTR(STRNUMB(1))*(C_NM2CM**3)
                isMatrixAtom = .true.

            case default
                write(*,*) "MCPSCUERROR: The Illegal flag: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "At box file Line: ",LINE
                pause
                stop
            case("&STAT")
                ! @todo (zhail#1#):
        end select
    END DO

    call this%Atoms_list%AppendOne(tempAtom,AtomNumb)

    if(isMatrixAtom .eq. .true.) then
        this%MatrixAtom = tempAtom
    end if

    return

    100 return 1
  end subroutine Load_OneSecton_AtomDefine



  !*********************************************
  subroutine Load_Box_Diffusors(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    character*256::STR
    character*32::KEYWORD
    integer::LINE
    !---Body---

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&DIFFUSORDEFSUBCTL")
                call this%LoadDiffusorsValue(hBoxFile,*100)
            case("&REACTDEFSUBCTL")
                call this%LoadReactions(hBoxFile,*100)
            case default
                write(*,*) "MCPSCUERROR: Illegal symbol:",KEYWORD
                write(*,*) "Please check box file at Line: ",LINE
                pause
                stop
        end select

    END DO

    return

    100 return 1
  end subroutine Load_Box_Diffusors


  !************************************************
  subroutine LoadDiffusorsValue(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    character*256::STR
    character*32::KEYWORD
    integer::LINE
    character*32::STRNUMB(10)
    type(ReadDiffusorPropList),pointer::cursor=>null()
    !---Body---
    allocate(this%ReadDiffusorProp_List)

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)
        SELECT CASE(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&DIFFUSOR")
                call this%LoadOneDiffusors(hBoxFile,*100)
            case("&FUNCSUBCTL")
                call this%LoadDiffusorsValueFromScript(hBoxFile,*100)
            case default
                write(*,*) "MCPSCUERROR: Illegal Keyword: ",KEYWORD
                write(*,*) "Please check box file at Line: ",LINE
                pause
                STOP

        END SELECT

    END DO

    cursor=>this%ReadDiffusorProp_List

    DO While(associated(cursor))
        call UPCASE(cursor%Diffusor%symbol)
        cursor=>cursor%next
    END DO

    Nullify(cursor)

    call this%ReadDiffusorProp_List%ConvertToDiffusorsTypesMap(this%Atoms_list,this%m_DiffusorTypesMap)

    return
    100 return 1
  end subroutine LoadDiffusorsValue

  !*******************************************
  subroutine LoadOneDiffusors(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    character*32::KEYWORD
    character*256::STR
    integer::LINE
    character*32::STRNUMB(10)
    type(ReadedDiffusorValue)::newDiffusor
    integer::N
    !---Body---
    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)
        SELECT CASE(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&SYMBOL")
                call EXTRACT_SUBSTR(STR,1,N,STRNUMB)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: You must specialize the diffusors symbol by 'Element1'#'number of Element1'@'Element2'#'number of Element2' "
                    pause
                    stop
                end if
                newDiffusor%symbol = trim(adjustl(STRNUMB(1)))

            case("&DIFFCOEFFVALUE_FREE")
                call EXTRACT_NUMB(STR,1,N,STRNUMB)

                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: You must special the diffusor value type in free matrix."
                    write(*,*) "At Line: ",LINE
                    pause
                    stop
                end if

                newDiffusor%DiffusorValueType_Free = ISTR(STRNUMB(1))

                if(newDiffusor%DiffusorValueType_Free .eq. p_DiffuseCoefficient_ByValue) then
                    call EXTRACT_NUMB(STR,2,N,STRNUMB)
                    if(N .LT. 2) then
                        write(*,*) "MCPSCUERROR: If you had used the by-diffusionValue strategy, you should give the diffusor value."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newDiffusor%DiffuseCoefficient_Free_Value = DRSTR(STRNUMB(2))
                else if(newDiffusor%DiffusorValueType_Free .eq. p_DiffuseCoefficient_ByArrhenius) then
                    call EXTRACT_NUMB(STR,3,N,STRNUMB)

                    if(N .LT. 3) then
                        write(*,*) "MCPSCUERROR: If you had used the by-Arrhenius strategy, you should give the prefacotr and active energy."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newDiffusor%PreFactor_Free = DRSTR(STRNUMB(2))
                    newDiffusor%ActEnergy_Free = DRSTR(STRNUMB(3))
                else if(newDiffusor%DiffusorValueType_Free .ne. p_DiffuseCoefficient_ByBCluster) then
                    write(*,*) "MCPSCUERROR: unknown diffusor value type :",newDiffusor%DiffusorValueType_Free
                    write(*,*) "At line: ",LINE
                    pause
                    stop
                end if

            case("&ECR_FREE")
                call EXTRACT_NUMB(STR,1,N,STRNUMB)

                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: You must special the ECR value type in free matrix."
                    write(*,*) "At Line: ",LINE
                    pause
                    stop
                end if

                newDiffusor%ECRValueType_Free = ISTR(STRNUMB(1))

                if(newDiffusor%ECRValueType_Free .eq. p_ECR_ByValue) then
                    call EXTRACT_NUMB(STR,2,N,STRNUMB)
                    if(N .LT. 2) then
                        write(*,*) "MCPSCUERROR: If you had used the by-ECRValue strategy, you should give the ECR value."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newDiffusor%ECR_Free = DRSTR(STRNUMB(2))*this%LatticeLength
                end if

            case("&DIFFCOEFFVALUE_INGB")
                call EXTRACT_NUMB(STR,1,N,STRNUMB)

                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: You must special the diffusor value type in GB."
                    write(*,*) "At Line: ",LINE
                    pause
                    stop
                end if

                newDiffusor%DiffusorValueType_InGB = ISTR(STRNUMB(1))

                if(newDiffusor%DiffusorValueType_InGB .eq. p_DiffuseCoefficient_ByValue) then
                    call EXTRACT_NUMB(STR,2,N,STRNUMB)
                    if(N .LT. 2) then
                        write(*,*) "MCPSCUERROR: If you had used the by-diffusionValue strategy, you should give the diffusor value."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newDiffusor%DiffuseCoefficient_InGB_Value = DRSTR(STRNUMB(2))
                else if(newDiffusor%DiffusorValueType_InGB .eq. p_DiffuseCoefficient_ByArrhenius) then
                    call EXTRACT_NUMB(STR,3,N,STRNUMB)

                    if(N .LT. 3) then
                        write(*,*) "MCPSCUERROR: If you had used the by-Arrhenius strategy, you should give the prefacotr and active energy."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newDiffusor%PreFactor_InGB = DRSTR(STRNUMB(2))
                    newDiffusor%ActEnergy_InGB = DRSTR(STRNUMB(3))
                else if(newDiffusor%DiffusorValueType_InGB .ne. p_DiffuseCoefficient_ByBCluster) then
                    write(*,*) "MCPSCUERROR: unknown diffusor value type :",newDiffusor%DiffusorValueType_InGB
                    write(*,*) "At line: ",LINE
                    pause
                    stop
                end if

            case("&ECR_INGB")
                call EXTRACT_NUMB(STR,1,N,STRNUMB)

                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: You must special the ECR value type in GB."
                    write(*,*) "At Line: ",LINE
                    pause
                    stop
                end if

                newDiffusor%ECRValueType_InGB = ISTR(STRNUMB(1))

                if(newDiffusor%ECRValueType_InGB .eq. p_ECR_ByValue) then
                    call EXTRACT_NUMB(STR,2,N,STRNUMB)
                    if(N .LT. 2) then
                        write(*,*) "MCPSCUERROR: If you had used the by-ECRValue strategy, you should give the ECR value."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newDiffusor%ECR_InGB = DRSTR(STRNUMB(2))*this%LatticeLength
                end if

            case default
                write(*,*) "MCPSCUERROR: The unknown symbol: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "Please check box file at Line: ",LINE
                pause
                stop
        END SELECT

    END DO

    call this%ReadDiffusorProp_List%AppendOne_ReadDiffusorPropList(newDiffusor)
    return

    100 return 1
  end subroutine LoadOneDiffusors

  !*******************************************
  subroutine LoadDiffusorsValueFromScript(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    character*256::STR
    character*32::KEYWORD
    character(kind=c_char,len=10000)::scriptStr
    integer::LINE
    !---Body---
    scriptStr = ''

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)
        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
        end select

        scriptStr = scriptStr(1:LEN_TRIM(scriptStr))//STR(1:LEN_TRIM(STR))//"\n"

    END DO

    scriptStr(LEN_TRIM(scriptStr):LEN_TRIM(scriptStr)+1) = CHAR(0)

    call ResloveDiffusorsValueFromCScript(scriptStr,this%ReadDiffusorProp_List)

    return
    100 return 1
  end subroutine LoadDiffusorsValueFromScript

  !************************************************
  subroutine LoadReactions(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    character*256::STR
    character*32::KEYWORD
    integer::LINE
    character*32::STRNUMB(10)
    type(ReadReactionPropList),pointer::cursor=>null()
    !---Body---
    allocate(this%ReadReactionProp_List)

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)
        SELECT CASE(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&REACTION")
                call this%LoadOneReaction(hBoxFile,*100)
            case("&FUNCSUBCTL")
                call this%LoadReactionsFromScript(hBoxFile,*100)
            case default
                write(*,*) "MCPSCUERROR: Illegal Keyword: ",KEYWORD
                write(*,*) "Please check box file at Line: ",LINE
                pause
                STOP

        END SELECT

    END DO

    cursor=>this%ReadReactionProp_List

    DO While(associated(cursor))
        call UPCASE(cursor%Reaction%SubjectSymbol)
        call UPCASE(cursor%Reaction%ObjectSymbol)
        cursor=>cursor%next
    END DO

    Nullify(cursor)

    call this%ReadReactionProp_List%ConvertToReactionsMap(this%Atoms_list,this%m_ReactionsMap)

    return
    100 return 1
  end subroutine LoadReactions

  !*******************************************
  subroutine LoadOneReaction(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    character*32::KEYWORD
    character*256::STR
    integer::LINE
    character*32::STRNUMB(10)
    type(ReadReactionPair)::newReactionPair
    integer::N
    !---Body---
    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)
        SELECT CASE(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&REACTPAIRS")
                call EXTRACT_SUBSTR(STR,2,N,STRNUMB)
                if(N .LT. 2) then
                    write(*,*) "MCPSCUERROR: You must specialize the &REACTPAIRS The reaction pairs by 'symbol of subject cluster', 'symbol of oubject cluster' "
                    pause
                    stop
                end if
                newReactionPair%SubjectSymbol = trim(adjustl(STRNUMB(1)))
                newReactionPair%ObjectSymbol = trim(adjustl(STRNUMB(2)))
            case("&REACTCOEFF")
                call EXTRACT_NUMB(STR,1,N,STRNUMB)

                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: You must special the reaction coefficients type."
                    write(*,*) "At Line: ",LINE
                    pause
                    stop
                end if

                newReactionPair%ReactionCoefficientType = ISTR(STRNUMB(1))

                if(newReactionPair%ReactionCoefficientType .eq. p_ReactionCoefficient_ByValue) then
                    call EXTRACT_NUMB(STR,2,N,STRNUMB)
                    if(N .LT. 2) then
                        write(*,*) "MCPSCUERROR: If you had used the reaction coefficients by-value strategy, you should give the corresponded value."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newReactionPair%ReactionCoefficient_Value = DRSTR(STRNUMB(2))
                else if(newReactionPair%ReactionCoefficientType .eq. p_ReactionCoefficient_ByArrhenius) then
                    call EXTRACT_NUMB(STR,3,N,STRNUMB)

                    if(N .LT. 3) then
                        write(*,*) "MCPSCUERROR: If you had used reaction coefficients by-Arrhenius strategy , you should give the prefacotr and active energy."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newReactionPair%PreFactor = DRSTR(STRNUMB(2))
                    newReactionPair%ActEnergy = DRSTR(STRNUMB(3))
                else
                    write(*,*) "MCPSCUERROR: unknown reaction coefficients type :",newReactionPair%ReactionCoefficientType
                    write(*,*) "At line: ",LINE
                    pause
                    stop
                end if

            case("&ECR")
                call EXTRACT_NUMB(STR,1,N,STRNUMB)

                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: You must special the ECR value type."
                    write(*,*) "At Line: ",LINE
                    pause
                    stop
                end if

                newReactionPair%ECRValueType = ISTR(STRNUMB(1))

                if(newReactionPair%ECRValueType .eq. p_ECR_ByValue) then
                    call EXTRACT_NUMB(STR,2,N,STRNUMB)
                    if(N .LT. 2) then
                        write(*,*) "MCPSCUERROR: If you had used the by-ECRValue strategy, you should give the ECR value (LU)."
                        write(*,*) "At Line: ",LINE
                        pause
                        stop
                    end if
                    newReactionPair%ECR = DRSTR(STRNUMB(2))*this%LatticeLength
                end if

            case default
                write(*,*) "MCPSCUERROR: The unknown symbol: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "Please check box file at Line: ",LINE
                pause
                stop
        END SELECT

    END DO

    call this%ReadReactionProp_List%AppendOne_ReadReactionPropList(newReactionPair)
    return

    100 return 1
  end subroutine LoadOneReaction

  !*******************************************
  subroutine LoadReactionsFromScript(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    character*256::STR
    character*32::KEYWORD
    character(kind=c_char,len=10000)::scriptStr
    integer::LINE
    !---Body---
    scriptStr = ''

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)
        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
        end select

        scriptStr = scriptStr(1:LEN_TRIM(scriptStr))//STR(1:LEN_TRIM(STR))//"\n"

    END DO

    scriptStr(LEN_TRIM(scriptStr):LEN_TRIM(scriptStr)+1) = CHAR(0)

    call ResloveReactionsFromCScript(scriptStr,this%ReadReactionProp_List)

    return
    100 return 1
  end subroutine LoadReactionsFromScript

  !**********************************************
  subroutine Load_Box_GrainBoundary(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRTEMP(10)
    !---Body---

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&SIMPLEDISTSUBCTL")
                this%m_GrainBoundary%GBInitType = p_GBIniConfig_Simple
                call this%Load_GB_Simple(hBoxFile,*100)

            case("&FILEDISTSUBCTL")
                this%m_GrainBoundary%GBInitType = p_GBIniConfig_SpecialDistFromFile
                call this%Load_GB_SpecialDistFromFile(hBoxFile,*100)

            case("&EXTFUNCDISTSUBCTL")
                this%m_GrainBoundary%GBInitType = p_GBIniConfig_SpecialDistFromExteFunc
                call this%Load_GB_SpecialDistFromExteFunc(hBoxFile,*100)

            case default
                write(*,*) "MCPSCUERROR: unKnown type to for grain boundary distribution!"
                write(*,*) KEYWORD
                pause
                stop
        end select

    END DO

    return
    100 return 1
  end subroutine Load_Box_GrainBoundary

  !*********************************************
  subroutine Load_GB_Simple(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRTEMP(10)
    !---Body---
    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&GRAINSNUMBER")
                call EXTRACT_NUMB(STR,1,N,STRTEMP)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: Too few parameters for the minial (cut-off) distance between seeds."
                    write(*,*) "At line: ",LINE
                    write(*,*) "You should special: &CUTOFF   The minial (cut-off) distance between seeds ="
                    pause
                    stop
                end if
                this%m_GrainBoundary%GrainNum = ISTR(STRTEMP(1))
            case("&BYSEEDSUBCTL")
                this%m_GrainBoundary%GBInitSimple_Strategy = p_GBInitSimple_BySeedCtl
                call this%Load_GB_Simple_Distribution_ByGSeedCtl(hBoxFile,*100)
            case("&BYGVOLUMSUBCTL")
                this%m_GrainBoundary%GBInitSimple_Strategy = p_GBInitSimple_ByGVolumCtl
                call this%Load_GB_Simple_Distribution_ByGVolumCtl(hBoxFile,*100)
            case default
                write(*,*) "MCPSCUERROR: The Illegal flag: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "At box file Line: ",LINE
                pause
                stop
        end select

    END DO

    return

    100 return 1
  end subroutine Load_GB_Simple

  !*********************************************
  subroutine Load_GB_Simple_Distribution_ByGSeedCtl(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRTEMP(10)
    !---Body---
    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&MINCUTOFF")
                call EXTRACT_NUMB(STR,1,N,STRTEMP)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: Too few parameters for the minial (cut-off) distance between seeds."
                    write(*,*) "At line: ",LINE
                    write(*,*) "You should special: &MINCUTOFF   The minial (cut-off) distance between seeds ="
                    pause
                    stop
                end if
                this%m_GrainBoundary%Cutoff(1) = DRSTR(STRTEMP(1))*this%LatticeLength
            case("&MAXCUTOFF")
                call EXTRACT_NUMB(STR,1,N,STRTEMP)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: Too few parameters for the max (cut-off) distance between seeds."
                    write(*,*) "At line: ",LINE
                    write(*,*) "You should special: &MAXCUTOFF   The max (cut-off) distance between seeds ="
                    pause
                    stop
                end if
                this%m_GrainBoundary%Cutoff(2) = DRSTR(STRTEMP(1))*this%LatticeLength
            case("&DISTANCE_GAUSS")
                call EXTRACT_NUMB(STR,2,N,STRTEMP)
                if(N .LT. 2) then
                    write(*,*) "MCPSCUERROR: Too few parameters for the distance distribution between seeds."
                    write(*,*) "At line: ",LINE
                    write(*,*) "You should special: &DISTANCE_GAUSS THE GAUSS DISTRIBUTION CENTRAL =, THE HALF WIDTH ="
                    pause
                    stop
                end if
                this%m_GrainBoundary%SeedsDistINI = DRSTR(STRTEMP(1))*this%LatticeLength
                this%m_GrainBoundary%SeedsDistSD = DRSTR(STRTEMP(2))*this%LatticeLength
            case default
                write(*,*) "MCPSCUERROR: The Illegal flag: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "At box file Line: ",LINE
                pause
                stop
        end select

    END DO

    if(this%m_GrainBoundary%Cutoff(2) .LT. this%m_GrainBoundary%Cutoff(1)) then
        write(*,*) "MCPSCUERROR: the Cut-off distance setting error."
        write(*,*) "Min cut-off: ",this%m_GrainBoundary%Cutoff(1)
        write(*,*) "Max cut-off: ",this%m_GrainBoundary%Cutoff(2)
        pause
        stop
    end if

    return
    100 return 1
  end subroutine Load_GB_Simple_Distribution_ByGSeedCtl

  !*********************************************
  subroutine Load_GB_Simple_Distribution_ByGVolumCtl(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*32::STRTEMP(10)
    !---Body---
    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&MINCUTOFF")
                call EXTRACT_NUMB(STR,1,N,STRTEMP)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: Too few parameters for the minial (cut-off) volum between grains."
                    write(*,*) "At line: ",LINE
                    write(*,*) "You should special: &MINCUTOFF The minial(cut-off) volum for grain ="
                    pause
                    stop
                end if
                this%m_GrainBoundary%Cutoff(1) = DRSTR(STRTEMP(1))*(this%LatticeLength**3)
            case("&MAXCUTOFF")
                call EXTRACT_NUMB(STR,1,N,STRTEMP)
                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: Too few parameters for the max (cut-off) volum between grains."
                    write(*,*) "At line: ",LINE
                    write(*,*) "You should special: &MAXCUTOFF The max(cut-off) volum for grain ="
                    pause
                    stop
                end if
                this%m_GrainBoundary%Cutoff(2) = DRSTR(STRTEMP(2))*(this%LatticeLength**3)
            case("&VOLUM_GAUSS")
                call EXTRACT_NUMB(STR,2,N,STRTEMP)
                if(N .LT. 2) then
                    write(*,*) "MCPSCUERROR: Too few parameters for the volum distribution between grains."
                    write(*,*) "At line: ",LINE
                    write(*,*) "You should special: &VOLUM_GAUSS THE GAUSS DISTRIBUTION CENTRAL =, THE HALF WIDTH ="
                    pause
                    stop
                end if
                this%m_GrainBoundary%GVolumINI = DRSTR(STRTEMP(1))*(this%LatticeLength**3)
                this%m_GrainBoundary%GVolumSD  = DRSTR(STRTEMP(2))*(this%LatticeLength**3)
            case default
                write(*,*) "MCPSCUERROR: The Illegal flag: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "At box file Line: ",LINE
                pause
                stop
        end select

    END DO

    if(this%m_GrainBoundary%Cutoff(2) .LT. this%m_GrainBoundary%Cutoff(1)) then
        write(*,*) "MCPSCUERROR: the Cut-off distance setting error."
        write(*,*) "Min cut-off: ",this%m_GrainBoundary%Cutoff(1)
        write(*,*) "Max cut-off: ",this%m_GrainBoundary%Cutoff(2)
        pause
        stop
    end if

    return
    100 return 1
  end subroutine Load_GB_Simple_Distribution_ByGVolumCtl

  !*********************************************
  subroutine Load_GB_SpecialDistFromFile(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    integer::LINE
    integer::N
    character*256::STR
    character*32::KEYWORD
    character*10::STRTMP(10)
    !---Body---

    DO While(.true.)
        call GETINPUTSTRLINE(hBoxFile,STR,LINE,"!",*100)
        call RemoveComments(STR,"!")
        STR = adjustl(STR)
        call GETKEYWORD("&",STR,KEYWORD)
        call UPCASE(KEYWORD)

        select case(KEYWORD(1:LENTRIM(KEYWORD)))
            case("&ENDSUBCTL")
                exit
            case("&FGBDIST")
                call EXTRACT_SUBSTR(STR,1,N,STRTMP)

                if(N .LT. 1) then
                    write(*,*) "MCPSCUERROR: You must special the grain boundary configuration file path !"
                    write(*,*) "At line : ",LINE
                    pause
                    stop
                end if

                if(LENTRIM(STRTMP(1)) .LE. 0) then
                    write(*,*) "MCPSCUERROR: The grain boundary configuration file path is null !"
                    pause
                    stop
                end if

                this%m_GrainBoundary%GBCfgFileName = adjustl((trim(STRTMP(1))))
            case default
                write(*,*) "MCPSCUERROR: The Illegal flag: ",KEYWORD(1:LENTRIM(KEYWORD))
                write(*,*) "At box file Line: ",LINE
                pause
                stop
        end select

    END DO

    return
    100 return 1
  end subroutine Load_GB_SpecialDistFromFile

  !*********************************************
  subroutine Load_GB_SpecialDistFromExteFunc(this,hBoxFile,*)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    integer,intent(in)::hBoxFile
    !---Local Vars---
    return
    100 return 1
  end subroutine Load_GB_SpecialDistFromExteFunc

  !**********************************************
  subroutine CleanSimulationBoxes(this)
    implicit none
    !---Dummy Vars---
    CLASS(SimulationBoxes)::this
    !---Body---


    return
  end subroutine CleanSimulationBoxes

  !**********************************************
  subroutine DestorySimulationBoxes(this)
    implicit none
    !---Dummy Vars---
    type(SimulationBoxes)::this
    !---Body---

    call this%Clean()

    call this%m_DiffusorTypesMap%Clean()

    call this%m_ReactionsMap%Clean()

    call this%m_GrainBoundary%Clean_Grainboundary()

    this%LatticeLength = 3.14D-8
    this%BOXBOUNDARY = 0.D0
    this%BOXSIZE = 0.D0
    this%HBOXSIZE = 0.D0
    this%BOXVOLUM = 0.D0

    call this%MatrixAtom%CleanAtom()

    call this%Atoms_list%CleanAtomsList()

    this%IniConfig = ""

    this%ImpFile = ""

    call this%ReadDiffusorProp_List%Clean_ReadDiffusorPropList()

    call this%ReadReactionProp_List%Clean_ReadReactionPropList()

    return
  end subroutine DestorySimulationBoxes

end module MCMF_TYPEDEF_SIMULATIONBOXARRAY

