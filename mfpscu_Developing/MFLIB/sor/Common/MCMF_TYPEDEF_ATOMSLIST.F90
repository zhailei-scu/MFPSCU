!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MCMF_TYPEDEF_ATOMSLIST
    use MCMF_CONSTANTS
    use MCMF_UTILITIES_FORMER
    implicit none

    type,public::ATOM
        integer::m_ID = 0                      ! the inner index used in program
        character(len=10)::m_Symbol = ""
        integer::m_ElementIndex = 0            ! the index in peridoic table of elements

        real::m_AtomMass = 0.D0

        real(KIND=KMCDF)::m_Volum = 0.D0      ! the volum of the atom (only need by matrix atoms)

        contains
        procedure,non_overridable,private,pass::CopyAtomFromOther
        procedure,non_overridable,public,pass::CleanAtom
        Generic::Assignment(=)=>CopyAtomFromOther
        Final::Clean_Atom
    end type

    type,public::AtomsList
        type(ATOM)::m_Atom
        integer::m_AtomsNumber = 0
        integer,private::m_ListCount = 0
        type(AtomsList),pointer::next=>null()

        contains
        procedure,non_overridable,public,pass::Get_ListCount
        procedure,non_overridable,public,pass::AppendOne
        procedure,non_overridable,public,pass::CleanAtomsList
        procedure,non_overridable,private,pass::CopyAtomsListFromOther
        procedure,non_overridable,public,pass::FindIndexBySymbol
        procedure,non_overridable,public,pass::GetSymbolByIndex
        Generic::Assignment(=)=>CopyAtomsListFromOther
        Final::Clean_AtomsList
    end type

    private::Get_ListCount
    private::CopyAtomFromOther
    private::AppendOne
    private::CleanAtom
    private::Clean_Atom
    private::CopyAtomsListFromOther
    private::FindIndexBySymbol
    private::GetSymbolByIndex
    private::CleanAtomsList
    private::Clean_AtomsList
    contains

    !*****************************************
    subroutine CopyAtomFromOther(this,others)
        implicit none
        !---Dummy Vars---
        CLASS(ATOM),intent(out)::this
        type(ATOM),intent(in)::others
        !---Body---

        this%m_ID = others%m_ID
        this%m_Symbol = others%m_Symbol
        this%m_ElementIndex = others%m_ElementIndex
        this%m_AtomMass = others%m_AtomMass
        this%m_Volum = others%m_Volum

        return
    end subroutine

    !****************************************
    subroutine CleanAtom(this)
        implicit none
        CLASS(ATOM)::this

        this%m_ID = 0
        this%m_Symbol = ""
        this%m_ElementIndex = 0
        this%m_AtomMass = 0.D0
        this%m_Volum = 0.D0

        return
    end subroutine CleanAtom

    !*******************************************
    subroutine Clean_Atom(this)
        implicit none
        type(ATOM)::this
        !---Body---
        call this%CleanAtom()

        return
    end subroutine Clean_Atom

    !*******************************************
    function Get_ListCount(this) result(TheCount)
        implicit none
        !---Dummy Vars---
        CLASS(AtomsList),target::this
        integer,intent(out)::TheCount
        !---Local Vars---
        type(AtomsList),pointer::cursorP=>null()
        !---Body---
        cursorP=>this
        if(.not. associated(cursorP)) then
            write(*,*) "MFPSCUERROR: You should allocate the AtomsList first!"
            pause
            stop
        end if

        TheCount = this%m_ListCount

        return
    end function Get_ListCount

    !*****************************************
    subroutine AppendOne(this,newOne,atomsNumb)
        implicit none
        !----Dummy Vars---
        CLASS(AtomsList),target::this
        type(ATOM)::newOne
        integer,intent(in)::atomsNumb
        !---Local Vars---
        type(AtomsList),pointer::cursor=>null()
        type(AtomsList),pointer::cursorP=>null()
        !---Body---
        cursorP=>this
        if(.not. associated(cursorP)) then
            write(*,*) "MFPSCUERROR: You should allocate the AtomsList first!"
            pause
            stop
        end if

        if(this%m_ListCount .LE. 0) then
            ! the (=) had been overrided
            newOne%m_ID = 1
            this%m_Atom = newOne
            this%m_AtomsNumber = atomsNumb
        else
            cursorP=>this
            cursor=>this%next

            if(IsStrEqual(cursorP%m_Atom%m_Symbol,newOne%m_Symbol)) then
                write(*,*) "MFPSCUERROR: The element symbol is redifined: ",newOne%m_Symbol
                pause
                stop
            end if

            DO While(associated(cursor))
                if(IsStrEqual(cursor%m_Atom%m_Symbol,newOne%m_Symbol)) then
                    write(*,*) "MFPSCUERROR: The element symbol is redifined: ",newOne%m_Symbol
                    pause
                    stop
                end if

                cursor=>cursor%next
                cursorP=>cursorP%next
            END DO

            allocate(cursor)
            cursor%next=>null()

            ! the (=) had been overrided
            newOne%m_ID = this%m_ListCount + 1
            cursor%m_Atom = newOne
            cursorP%next=>cursor
            cursor%m_AtomsNumber = atomsNumb
        end if

        this%m_ListCount = this%m_ListCount + 1

        if(this%m_ListCount .GT. p_ATOMS_GROUPS_NUMBER) then
            write(*,*) "MFPSCUERROR: The defined elements group number is greater than defined max atoms kinds: ",p_ATOMS_GROUPS_NUMBER
            pause
            stop
        end if

        Nullify(cursorP)
        cursorP=>null()
        NullifY(cursor)
        cursor=>null()
        return
    end subroutine AppendOne

    !*****************************************
    subroutine CleanAtomsList(this)
        implicit none
        !---Dummy Vars---
        CLASS(AtomsList),target::this
        !---Local Vars---
        type(AtomsList),pointer::cursor=>null()
        type(AtomsList),pointer::next=>null()
        !---Body---
        cursor=>this
        if(.not. associated(cursor)) then
            return
        end if

        if(cursor%Get_ListCount() .LE. 0) then
            return
        end if

        cursor=>this%next

        call CleanAtom(cursor%m_Atom)
        cursor%m_ListCount = 0
        cursor%m_AtomsNumber = 0

        DO While(associated(cursor))
            next=>cursor%next
            call CleanAtom(cursor%m_Atom)
            cursor%m_ListCount = 0
            cursor%m_AtomsNumber = 0

            deallocate(cursor)
            Nullify(cursor)
            cursor=>next
        END DO

        Nullify(cursor)
        cursor=>null()
        Nullify(next)
        next=>null()

        return
    end subroutine CleanAtomsList

    !*********************************
    function FindIndexBySymbol(this,Symbol) result(TheIndex)
        use MiniUtilities
        implicit none
        !---Dummy Vars---
        CLASS(AtomsList),intent(in),target::this
        character*(*),intent(in)::Symbol
        integer,intent(out)::TheIndex
        !---Local Vars---
        type(AtomsList),pointer::cursor=>null()
        integer::tempIndex
        character*256::tempSymbol
        !---Body---
        TheIndex = 0

        tempIndex = 1

        cursor=>this

        tempSymbol = ""

        tempSymbol(1:LENTRIM(Symbol)) = Symbol(1:LENTRIM(Symbol))

        call UPCASE(tempSymbol)

        DO While(associated(cursor))

            if(IsStrEqual(tempSymbol,cursor%m_Atom%m_Symbol)) then
                TheIndex = tempIndex
                if(TheIndex .ne. cursor%m_Atom%m_ID) then
                    write(*,*) "MFPSCUERROR: The elements id is not stored correct: ",Symbol
                    write(*,*) "The definded id is ",cursor%m_Atom%m_ID
                    write(*,*) "However, it is located in the atoms defined list for postion: ",TheIndex
                    pause
                    stop
                end if

                exit
            end if

            tempIndex = tempIndex + 1

            cursor=>cursor%next

        END DO

        nullify(cursor)

        if(TheIndex .LE. 0) then
            write(*,*) "MFPSCUERROR: The element is not defined: ",Symbol
            pause
            stop
        end if

        return
    end function FindIndexBySymbol

    !*********************************
    subroutine GetSymbolByIndex(this,TheIndex,TheResult)
        implicit none
        !---Dummy Vars---
        CLASS(AtomsList),target,intent(in)::this
        integer::TheIndex
        character*(*)::TheResult
        !---Local Vars---
        TYPE(AtomsList),pointer::cursor=>null()
        integer::tempCount
        !---Body---
        TheResult = ""

        cursor=>this
        tempCount = 0

        if(this%Get_ListCount() .GT. 0) then
            DO while(associated(cursor))
                tempCount = tempCount + 1

                if(cursor%m_Atom%m_ID .eq. TheIndex) then
                    TheResult = cursor%m_Atom%m_Symbol
                    exit
                end if
                cursor=>cursor%next
            END DO
        end if

        if(tempCount .eq. 0) then
            write(*,*) "The elements index is not existed in the defined elements."
            pause
            stop
        end if

        return
    end subroutine

    !*********************************
    subroutine CopyAtomsListFromOther(this,others)
        implicit none
        !---Dummy Vars---
        CLASS(AtomsList),intent(out),target::this
        CLASS(AtomsList),intent(in),target::others
        !---Local Vars---
        type(AtomsList),pointer::othersCursor=>null()
        type(AtomsList),pointer::cursor=>null()
        type(AtomsList),pointer::pCursor=>null()
        !---Body---
        pCursor=>this

        if(.not. associated(pCursor)) then
            write(*,*) "MFPSCUERROR: You need to allocate the AtomsList first!"
            pause
            stop
        end if

        othersCursor=>others
        if(.not. associated(othersCursor)) then
            Nullify(pCursor)
            return
        end if

        call CleanAtomsList(this)

        pCursor=>this
        othersCursor=>others%next
        cursor=this%next

        !---the (=) had been overrided
        this%m_Atom = others%m_Atom
        this%m_ListCount = others%m_ListCount
        this%m_AtomsNumber = others%m_AtomsNumber

        DO While(associated(othersCursor))
            allocate(cursor)
            cursor%m_Atom = othersCursor%m_Atom
            cursor%m_ListCount = othersCursor%m_ListCount
            cursor%m_AtomsNumber = othersCursor%m_AtomsNumber
            pCursor%next=>cursor
            cursor=>cursor%next
            pCursor=>cursor%next
            othersCursor=>othersCursor%next
        END DO

        Nullify(cursor)
        Nullify(pCursor)
        Nullify(othersCursor)

        return
    end subroutine CopyAtomsListFromOther

    !*****************************************
    subroutine Clean_AtomsList(this)
        implicit none
        !---Dummy Vars---
        type(AtomsList)::this
        !---Body---
        call this%CleanAtomsList()

        return
    end subroutine Clean_AtomsList

end module MCMF_TYPEDEF_ATOMSLIST

