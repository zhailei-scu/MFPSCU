#include "../../../Macro"

module MCMF_TYPEDEF_REACTIONSVALUE
    use MCMF_CONSTANTS
    use iso_c_binding
    use MCMF_TYPEDEF_ACLUSTER
    use MCMF_UTILITIES
    implicit none

    type,public::ReadReactionPair
        character(kind=c_char,len=20)::SubjectSymbol = ""
        character(kind=c_char,len=20)::ObjectSymbol = ""

        integer(c_int)::ReactionCoefficientType = p_ReactionCoefficient_ByValue

        real(c_double)::ReactionCoefficient_Value = -1.D0 ! < 0 means not occur, >=1 means must occur
        real(c_double)::PreFactor = 0.D0
        real(c_double)::ActEnergy = 0.D0

        integer(c_int)::ECRValueType = p_ECR_ByValue
        real(c_double)::ECR = 0.D0

        contains
        procedure,private,non_overridable,pass::CopyReadedReactionPairFromOther
        procedure,public,non_overridable,pass::Convert2ReactionValue
        GENERIC::Assignment(=)=>CopyReadedReactionPairFromOther
        Final::CleanReadedReactionPair

    end type ReadReactionPair


    !**********Based by our test, if we want to combine C, the same data structure need to be defined in fortran and C
    !**********More important, we can not use the inherit in fotran types define, or it would meeting mismatch while running the code
    !**********So we must define the type ReadReactionPair in fortran to match C type to reslove the reaction from C.
    !**********Howeve, the member "symbol" is not required after the reactions is used for calculation, so we re-define type ReactionValue
    !**********Which contains all same members except the "symbol".
    type,public::ReactionValue
        integer::ReactionCoefficientType PREASSIGN p_ReactionCoefficient_ByValue

        real(kind=KMCDF)::ReactionCoefficient_Value PREASSIGN -1.D0 ! < 0 means not occur, >=1 means must occur
        real(kind=KMCDF)::PreFactor PREASSIGN 0.D0
        real(kind=KMCDF)::ActEnergy PREASSIGN 0.D0

        integer::ECRValueType PREASSIGN p_ECR_ByValue
        real(kind=KMCDF)::ECR PREASSIGN 0.D0

        contains
        procedure,private,non_overridable,pass::CopyReactionValueFromOther
        GENERIC::Assignment(=)=>CopyReactionValueFromOther
        !*********Important note: The PGI CUDA Fortran not support the Final symbol, the final symbol would cause the
        !*********Compiler error************
        !Final::CleanDiffusorValue
    end type ReactionValue

    type,public::ReactionEntity
        integer(kind=KMCLINT)::SubjectCode PREASSIGN 0
        integer(kind=KMCLINT)::ObjectCode PREASSIGN 0

        type(ReactionValue)::TheValue

        integer::NextIndex PREASSIGN 0

        contains
        procedure,public,non_overridable,pass::CopyReactionEntityFromOther
        Generic::Assignment(=)=>CopyReactionEntityFromOther
        !---Similarly, the final procedure cannot be used here---
    end type ReactionEntity

    type,public::ReactionsMap

        integer::MaxDivideGroups_SingleElement = ISHFT(1,9) - 1

        integer::MapBitLength = 18
        integer::MapLength = ISHFT(1,18)

        integer,dimension(:,:),allocatable::SingleAtomsDivideArrays

        type(ReactionEntity),dimension(:),allocatable::RecordsEntities

        contains
        procedure,public,non_overridable,pass::put=>putToReactionsMap
        procedure,public,non_overridable,pass::get=>getValueFromReactionsMap
        procedure,public,non_overridable,pass::constructor=>ReactionsMapConstructor

        procedure,private,non_overridable,pass::getCode
        procedure,private,non_overridable,pass::hash
        procedure,private,non_overridable,pass::GetIndexFor
        procedure,private,non_overridable,pass::CopyReactionsMapFromOther
        procedure,public,non_overridable,pass::Clean=>Clean_ReactionsMap
        Generic::Assignment(=)=>CopyReactionsMapFromOther
        Final::CleanReactionsMap
    end type ReactionsMap


    private::CopyReadedReactionPairFromOther
    private::Convert2ReactionValue
    !private::CleanReadedReactionPair
    private::CopyReactionValueFromOther
    private::CopyReactionEntityFromOther
    private::putToReactionsMap
    private::getValueFromReactionsMap
    private::ReactionsMapConstructor
    private::getCode
    private::hash
    private::GetIndexFor
    private::CopyReactionsMapFromOther
    private::Clean_ReactionsMap
    private::CleanReactionsMap

    contains

    !*************************************
    subroutine CopyReadedReactionPairFromOther(this,Others)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPair),intent(out)::this
        type(ReadReactionPair),intent(in)::Others
        !---Body---

        this%SubjectSymbol = Others%SubjectSymbol
        this%ObjectSymbol = Others%ObjectSymbol

        this%ReactionCoefficientType = Others%ReactionCoefficientType

        this%ReactionCoefficient_Value = Others%ReactionCoefficient_Value
        this%PreFactor = Others%PreFactor
        this%ActEnergy = Others%ActEnergy

        this%ECRValueType = Others%ECRValueType
        this%ECR = Others%ECR

        return
    end subroutine CopyReadedReactionPairFromOther

    !*******************************************
    function Convert2ReactionValue(this) result(TheReactionValue)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPair),intent(in)::this
        type(ReactionValue),intent(out)::TheReactionValue
        !---Body---

        TheReactionValue%ReactionCoefficientType = this%ReactionCoefficientType

        TheReactionValue%ReactionCoefficient_Value = this%ReactionCoefficient_Value

        TheReactionValue%PreFactor = this%PreFactor

        TheReactionValue%ActEnergy = this%ActEnergy

        TheReactionValue%ECRValueType = this%ECRValueType

        TheReactionValue%ECR = this%ECR

        return
    end function Convert2ReactionValue

    !*******************************************
    subroutine CleanReadedReactionPair(this)
        implicit none
        !---Dummy Vars---
        type(ReadReactionPair)::this
        !---Body---
        this%SubjectSymbol = ""
        this%ObjectSymbol = ""

        this%ReactionCoefficientType = p_ReactionCoefficient_ByValue

        this%ReactionCoefficient_Value = -1.D0 ! < 0 means not occur, >=1 means must occur
        this%PreFactor = 0.D0
        this%ActEnergy = 0.D0

        this%ECRValueType = p_ECR_ByValue
        this%ECR = 0.D0
        return
    end subroutine

    !*********************************************
    subroutine CopyReactionValueFromOther(this,Others)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionValue),intent(out)::this
        type(ReactionValue),intent(in)::Others
        !---Body---

        this%ReactionCoefficientType = Others%ReactionCoefficientType

        this%ReactionCoefficient_Value = Others%ReactionCoefficient_Value
        this%PreFactor = Others%PreFactor
        this%ActEnergy = Others%ActEnergy

        this%ECRValueType = Others%ECRValueType
        this%ECR = Others%ECR

        return
    end subroutine CopyReactionValueFromOther


    !*******************************************
    subroutine CleanReactionValue(this)
        implicit none
        !---Dummy Vars---
        type(ReactionValue)::this
        !---Body---

        this%ReactionCoefficientType = p_ReactionCoefficient_ByValue

        this%ReactionCoefficient_Value = -1.D0 ! < 0 means not occur, >=1 means must occur
        this%PreFactor = 0.D0
        this%ActEnergy = 0.D0

        this%ECRValueType = p_ECR_ByValue
        this%ECR = 0.D0
        return
    end subroutine CleanReactionValue

    !*******************************************
    subroutine CopyReactionEntityFromOther(this,others)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionEntity),intent(out)::this
        type(ReactionEntity),intent(in)::others
        !---Body---

        this%SubjectCode = others%SubjectCode
        this%ObjectCode = others%ObjectCode

        !---The (=) hand been overrided
        this%TheValue = others%TheValue

        this%NextIndex = others%NextIndex

        return
    end subroutine CopyReactionEntityFromOther

    !********************************************
    subroutine ReactionsMapConstructor(this,SingleAtomsDivideArrays,MapLength)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionsMap)::this
        integer,dimension(:,:),allocatable::SingleAtomsDivideArrays
        integer::MapLength
        !---Local Vars---
        integer::tempLength
        integer::MapBitLength
        !---Body---

        MapBitLength = 0

        tempLength = MapLength

        if(size(SingleAtomsDivideArrays,DIM=1) .ne. p_ATOMS_GROUPS_NUMBER) then
            write(*,*) "MFPSCUERROR: In current version, the max elements group is ",p_ATOMS_GROUPS_NUMBER
            write(*,*) "However, the SingleAtomsDivideArrays own elements group: ",size(SingleAtomsDivideArrays,DIM=1)
            pause
            stop
        end if

        this%MaxDivideGroups_SingleElement = size(SingleAtomsDivideArrays,DIM=2)


        this%MapLength = 0


        DO While(tempLength .GT. 0)
            this%MapLength = ISHFT(this%MapLength,1)
            this%MapLength = this%MapLength + 1
            tempLength = ISHFT(tempLength,-1)
            MapBitLength = MapBitLength + 1

        END DO

        this%MapBitLength = MapBitLength

        call AllocateArray_Host(this%SingleAtomsDivideArrays,p_ATOMS_GROUPS_NUMBER,this%MaxDivideGroups_SingleElement,"SingleAtomsDivideArrays")

        this%SingleAtomsDivideArrays = SingleAtomsDivideArrays

        allocate(this%RecordsEntities(this%MapLength))

        return
    end subroutine ReactionsMapConstructor

    !********************************************
    subroutine CopyReactionsMapFromOther(this,Others)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionsMap),intent(out)::this
        type(ReactionsMap),intent(in)::Others
        !---Body---
        this%MaxDivideGroups_SingleElement = Others%MaxDivideGroups_SingleElement
        this%MapLength = Others%MapLength
        this%MapBitLength = Others%MapBitLength

        this%SingleAtomsDivideArrays = Others%SingleAtomsDivideArrays
        this%RecordsEntities = Others%RecordsEntities

        return
    end subroutine CopyReactionsMapFromOther

    !**********************************************
    subroutine Clean_ReactionsMap(this)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionsMap)::this
        !---Body---
        this%MaxDivideGroups_SingleElement = 0

        this%MapBitLength = 0

        this%MapLength = 0

        call DeAllocateArray_Host(this%SingleAtomsDivideArrays,"SingleAtomsDivideArrays")

        if(allocated(this%RecordsEntities)) then
            deallocate(this%RecordsEntities)
        end if

        return
    end subroutine Clean_ReactionsMap

    !**********************************************
    subroutine CleanReactionsMap(this)
        implicit none
        !---Dummy Vars---
        type(ReactionsMap)::this
        !---Body---

        call this%Clean()

        return
    end subroutine

    !**********************************
    subroutine putToReactionsMap(this,KeySubject,KeyObject,TheValue)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionsMap)::this
        type(ACluster),intent(in)::KeySubject
        type(ACluster),intent(in)::KeyObject
        CLASS(ReactionValue),intent(in)::TheValue
        !---Local Vars---
        integer(kind=KMCLINT)::SubjectCode
        integer(kind=KMCLINT)::ObjectCode
        integer(kind=KMCLINT)::reSparedCode
        integer(kind=KMCLINT)::IndexFor
        integer(kind=KMCLINT)::NextIndex
        integer::ICount
        !---Body---
        SubjectCode = this%getCode(KeySubject%m_Atoms)
        ObjectCode = this%getCode(KeyObject%m_Atoms)

        reSparedCode = this%hash(SubjectCode,ObjectCode)

        IndexFor = this%GetIndexFor(reSparedCode)

        ! Handle the conflictions
        if((this%RecordsEntities(IndexFor)%SubjectCode .NE. 0 .or. this%RecordsEntities(IndexFor)%ObjectCode .NE. 0))then
            if(SubjectCode .NE. this%RecordsEntities(IndexFor)%SubjectCode .or. ObjectCode .NE. this%RecordsEntities(IndexFor)%ObjectCode) then

                NextIndex = this%RecordsEntities(IndexFor)%NextIndex

                ICount = 1
                DO While(NextIndex .GT. 0)

                    if(SubjectCode .eq. this%RecordsEntities(IndexFor)%SubjectCode .AND. ObjectCode .eq. this%RecordsEntities(IndexFor)%ObjectCode) then
                        write(*,*) "MFPSCUERROR: The reaction is redefined !"
                        write(*,*) "Subject : ",KeySubject%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA
                        write(*,*) "Subject : ",KeySubject%m_Atoms(1:p_ATOMS_GROUPS_NUMBER)%m_NA
                        pause
                        stop
                    end if

                    ICount = ICount + 1

                    if(ICount .GT. this%MapLength) then
                        write(*,*) "MFPSCUERROR: The reactions map is not sufficient to store all kinds of reaction pairs."
                        pause
                        stop
                    end if

                    IndexFor = NextIndex

                    NextIndex = this%RecordsEntities(IndexFor)%NextIndex

                END DO

                ICount = 1
                DO While(.true.)
                    ICount = ICount + 1
                    if(ICount .GT. this%MapLength) then
                        write(*,*) "MFPSCUERROR: The reactions map is not sufficient to store all kinds of reaction pairs."
                        pause
                        stop
                    end if

                    reSparedCode = reSparedCode + 1

                    NextIndex = this%GetIndexFor(reSparedCode)

                    if(NextIndex .LE. 0 .or. NextIndex .GT. this%MapLength) then
                        reSparedCode = 1
                    else if(this%RecordsEntities(NextIndex)%SubjectCode .eq. 0 .AND. this%RecordsEntities(NextIndex)%ObjectCode .eq. 0) then
                        exit
                    end if

                END DO

                this%RecordsEntities(IndexFor)%NextIndex = NextIndex

                IndexFor = NextIndex
            end if
        end if

        this%RecordsEntities(IndexFor)%SubjectCode = SubjectCode
        this%RecordsEntities(IndexFor)%ObjectCode = ObjectCode
        this%RecordsEntities(IndexFor)%TheValue = TheValue

        return
    end subroutine putToReactionsMap

    !**********************************
    function getValueFromReactionsMap(this,KeySubject,KeyObject) result(TheValue)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionsMap)::this
        type(ACluster),intent(in)::KeySubject
        type(ACluster),intent(in)::KeyObject
        type(ReactionValue)::TheValue
        !---Local Vars---
        integer(kind=KMCLINT)::SubjectCode
        integer(kind=KMCLINT)::ObjectCode
        integer(kind=KMCLINT)::reSparedCode
        integer(kind=KMCLINT)::IndexFor
        integer(kind=KMCLINT)::NextIndex
        !---Body---
        SubjectCode = this%getCode(KeySubject%m_Atoms)
        ObjectCode = this%getCode(KeyObject%m_Atoms)

        reSparedCode = this%hash(SubjectCode,ObjectCode)

        IndexFor = this%GetIndexFor(reSparedCode)

        DO While(IndexFor .GT. 0)

            if(this%RecordsEntities(IndexFor)%SubjectCode .eq. SubjectCode .AND. this%RecordsEntities(IndexFor)%ObjectCode .eq. ObjectCode) then
                TheValue = this%RecordsEntities(IndexFor)%TheValue
                exit
            end if

            !---here, we consider that the reaction is symmetrical which means if A can react with B , vice versa
            if(this%RecordsEntities(IndexFor)%SubjectCode .eq. ObjectCode .AND. this%RecordsEntities(IndexFor)%ObjectCode .eq. SubjectCode) then
                TheValue = this%RecordsEntities(IndexFor)%TheValue
                exit
            end if

            IndexFor = this%RecordsEntities(IndexFor)%NextIndex
        END DO


        return
    end function getValueFromReactionsMap

    !**********************************
    function getCode(this,Atoms) result(Code)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionsMap)::this
        type(Single_AtomsSet)::Atoms(:)
        integer(kind=KMCLINT)::Code
        !---Local Vars---
        integer::I
        integer::J
        !---Body---
        Code = 0


        DO I = 1,p_ATOMS_GROUPS_NUMBER

            J = BinarySearch_GE(Atoms(I)%m_NA,this%SingleAtomsDivideArrays(I,:),1,this%MaxDivideGroups_SingleElement)

            Code = ISHFT(Code,this%MapBitLength) + J

        END DO

        return
    end function getCode

    !**********************************
    function hash(this,SubjectCode,ObjectCode) result(reSparedCode)
        implicit none
        ! Purpose: to spare the code to be more uniform
        !---Dummy Vars---
        CLASS(ReactionsMap)::this
        integer(kind=KMCLINT)::SubjectCode
        integer(kind=KMCLINT)::ObjectCode
        integer(kind=KMCLINT)::reSparedCode
        !---Local Vars---
        integer(kind=KMCLINT)::TempCode
        !---Body---
        reSparedCode = IOR(SubjectCode,ObjectCode)

        TempCode = reSparedCode
        TempCode = ISHFT(TempCode,-this%MapBitLength)
        DO While(TempCode .GT. 0)

            reSparedCode = IOR(reSparedCode,IBITS(TempCode,0,this%MapBitLength-1))

            TempCode = ISHFT(TempCode,-this%MapBitLength)

        END DO

        return
    end function hash

    !************************************
    function GetIndexFor(this,Code) result(IndexFor)
        implicit none
        !---Dummy Vars---
        CLASS(ReactionsMap)::this
        integer(kind=KMCLINT)::Code
        integer(kind=KMCLINT)::IndexFor
        !---Body---
        IndexFor = IAND(Code,this%MapLength)

        return
    end function GetIndexFor



end module MCMF_TYPEDEF_REACTIONSVALUE
