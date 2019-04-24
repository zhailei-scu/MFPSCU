module MCMF_TYPEDEF_ReactionPropList
    use MCMF_TYPEDEF_REACTIONSVALUE
    use MCMF_TYPEDEF_ATOMSLIST
    use MCMF_TYPEDEF_ACLUSTER
    use MCMF_UTILITIES

    implicit none

    type,public::ReadReactionPropList
        type(ReadReactionPair)::Reaction
        type(ReadReactionPropList),pointer::next=>null()

        integer,private::ListCount = 0

        contains
        procedure,non_overridable,pass,public::AppendOne_ReadReactionPropList
        procedure,non_overridable,pass,public::AppendArray_ReadReactionPropList
        procedure,non_overridable,pass,public::GetReadReactionByListIndex
        procedure,non_overridable,pass,public::ConvertToReactionsMap
        procedure,non_overridable,pass,public::GetList_Count=>GetReadReactionPropList_Count
        procedure,non_overridable,pass,private::CopyReadReactionPropListFromOther
        procedure,non_overridable,pass,public::Clean_ReadReactionPropList
        procedure,non_overridable,pass,public::PrintOutCheckingResult
        procedure,non_overridable,pass,public::WhetherFreeDiffusion
        Generic::Assignment(=)=>CopyReadReactionPropListFromOther
        Final::CleanReadReactionPropList

    end type ReadReactionPropList

    interface InterpCScript_ReactionsDef
        function InterpCScript_ReactionsDef(scriptStr) result(ArraySize) bind(c,name="InterpCScript_ReactionsDef")
            use iso_c_binding
            implicit none
            !---Dummy Vars---
            character(kind=c_char,len=10000)::scriptStr
            integer(kind=c_int)::ArraySize
        end function
    end interface

    interface GetInterpedReactionsArray
        subroutine GetInterpedReactionsArray(FReactionsDefArray) bind(c,name="GetInterpedReactionsArray")
            use iso_c_binding
            use MCMF_TYPEDEF_REACTIONSVALUE
            implicit none
            !---Dummy Vars---
            type(ReadReactionPair)::FReactionsDefArray(*)
        end subroutine

    end interface


    private::AppendOne_ReadReactionPropList
    private::AppendArray_ReadReactionPropList
    private::GetReadReactionByListIndex
    private::ConvertToReactionsMap
    private::GetReadReactionPropList_Count
    private::CopyReadReactionPropListFromOther
    private::Clean_ReadReactionPropList
    private::PrintOutCheckingResult
    private::CleanReadReactionPropList

    contains

    !*******************************************************
    subroutine ConvertToReactionsMap(this,BasicAtomsList,TheReactionsMap)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),intent(in),target::this
        type(AtomsList),intent(in)::BasicAtomsList
        type(ReactionsMap)::TheReactionsMap
        !---Local Vars---
        type(AtomsSetRange),dimension(:),allocatable::TheAtomsSetsRangesArray
        type(ReadReactionPropList),pointer::cursor=>null()
        integer,dimension(:,:),allocatable::RangesArray
        integer,dimension(:,:),allocatable::tempSingleAtomsDivideArrays
        integer,dimension(:,:),allocatable::SingleAtomsDivideArrays
        integer,dimension(:,:),allocatable::AtomsSetsRangesMarkArray
        type(AClusterList),dimension(:),pointer::ConstructClusterListsArray=>null()
        type(AClusterList),pointer::SubjectClusterListCursor=>null()
        type(AClusterList),pointer::ObjectClusterListCursor=>null()
        integer::UsedAtomsType(p_ATOMS_GROUPS_NUMBER)
        type(ReadReactionPair)::Reaction
        integer::temp_MaxDivideGroups
        integer::tempMax
        integer::tempMin
        integer::tempMinLoc
        integer::DivideGroups_SingleElement(p_ATOMS_GROUPS_NUMBER)
        integer::tempDivideGroups_SingleElement(p_ATOMS_GROUPS_NUMBER)
        integer::tempIndex
        integer::ListCount
        integer::IElement
        integer::I
        integer::J
        integer::TheID
        integer::coverageCount
        logical::coverageOnElement
        integer::Maplength
        !---Body---
        cursor=>this

        if(.not. associated(cursor)) then
            write(*,*) "MCPSCUERROR: You should allocate the ReadReactionPropList first!"
            pause
            stop
        end if

        ListCount = this%GetList_Count()

        if(ListCount .LE. 0) then
            return
        else
            allocate(TheAtomsSetsRangesArray(2*ListCount))
        end if

        tempIndex = 0

        DO While(associated(cursor))

            tempIndex = tempIndex + 1

            call TheAtomsSetsRangesArray(tempIndex)%ReleaseSetsRange()

            TheAtomsSetsRangesArray(tempIndex) = ResolveSymbol2AtomsSetRange(cursor%Reaction%SubjectSymbol,BasicAtomsList)

            tempIndex = tempIndex + 1

            call TheAtomsSetsRangesArray(tempIndex)%ReleaseSetsRange()

            TheAtomsSetsRangesArray(tempIndex) = ResolveSymbol2AtomsSetRange(cursor%Reaction%ObjectSymbol,BasicAtomsList)

            cursor=>cursor%next

        END DO

        if(2*ListCount .ne. tempIndex) then
            write(*,*) "MCPSCUERROR: The reaction define List count is error."
            write(*,*) "The record list count is ",2*ListCount
            write(*,*) "In fact, the actual list count is: ",2*tempIndex
            pause
            stop
        end if

        call AllocateArray_Host(RangesArray,p_ATOMS_GROUPS_NUMBER,2*2*ListCount,"RangesArray")
        RangesArray = 0

        call AllocateArray_Host(tempSingleAtomsDivideArrays,p_ATOMS_GROUPS_NUMBER,2*2*ListCount,"tempSingleAtomsDivideArrays")
        tempSingleAtomsDivideArrays = 0

        UsedAtomsType = 0

        DO I = 1,2*ListCount
            DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                TheID = TheAtomsSetsRangesArray(I)%m_SetsRange(IElement)%m_ID
                if(TheID .GT. 0) then
                    RangesArray(TheID,(I -1)*2 + 1) = TheAtomsSetsRangesArray(I)%m_SetsRange(TheID)%m_NA_From
                    RangesArray(TheID,(I -1)*2 + 2) = TheAtomsSetsRangesArray(I)%m_SetsRange(TheID)%m_NA_To
                    UsedAtomsType(TheID) = 1
                end if
            END DO
        END DO

        !---Check the coverage---
        DO I = 1,ListCount
            DO J = I+1,ListCount
                coverageCount = 0

                DO IElement = 1,p_ATOMS_GROUPS_NUMBER

                    if(UsedAtomsType(IElement) .GT. 0) then

                        coverageOnElement = .false.

                        if(IsRangeCoverage(RangesArray(IElement,(I -1)*4 + 1),RangesArray(IElement,(I -1)*4 + 2),       &
                                           RangesArray(IElement,(J -1)*4 + 1),RangesArray(IElement,(J -1)*4 + 2)) .AND. &
                           IsRangeCoverage(RangesArray(IElement,(I -1)*4 + 3),RangesArray(IElement,(I -1)*4 + 4),       &
                                           RangesArray(IElement,(J -1)*4 + 3),RangesArray(IElement,(J -1)*4 + 4))) then

                            coverageOnElement = .true.
                        end if

                        if(IsRangeCoverage(RangesArray(IElement,(I -1)*4 + 1),RangesArray(IElement,(I -1)*4 + 2),       &
                                           RangesArray(IElement,(J -1)*4 + 3),RangesArray(IElement,(J -1)*4 + 4)) .AND. &
                           IsRangeCoverage(RangesArray(IElement,(I -1)*4 + 3),RangesArray(IElement,(I -1)*4 + 4),       &
                                           RangesArray(IElement,(J -1)*4 + 1),RangesArray(IElement,(J -1)*4 + 2))) then

                            coverageOnElement = .true.
                        end if

                        if(coverageOnElement .eq. .true.) then
                            coverageCount = coverageCount + 1
                        end if
                    else
                        coverageCount = coverageCount + 1
                    end if
                END DO

                if(coverageCount .GE. p_ATOMS_GROUPS_NUMBER) then
                    Reaction = this%GetReadReactionByListIndex(I)

                    write(*,*) "MCPSCUERROR: The reaction define is overlapping between reactions pair: ",Reaction%SubjectSymbol,Reaction%ObjectSymbol
                    Reaction = this%GetReadReactionByListIndex(J)
                    write(*,*) "and reaction pair: ",Reaction%SubjectSymbol,Reaction%ObjectSymbol
                    pause
                    stop
                end if
            END DO
        END DO

        !---Construct the range array---
        DO IElement = 1,p_ATOMS_GROUPS_NUMBER

            temp_MaxDivideGroups = 0

            tempMax = maxval(RangesArray(IElement,:))

            if(tempMax .GT. 0) then

                DO I = 1,4*ListCount

                    tempMin = minval(RangesArray(IElement,:),MASK=(RangesArray(IElement,:) .GT. 0))

                    tempMinLoc = minloc(RangesArray(IElement,:),DIM=1,mask=(RangesArray(IElement,:) .GT. 0))

                    RangesArray(IElement,tempMinLoc) = tempMax

                    if(temp_MaxDivideGroups .eq. 0) then
                            temp_MaxDivideGroups = temp_MaxDivideGroups + 1
                            tempSingleAtomsDivideArrays(IElement,temp_MaxDivideGroups) = tempMin
                    else if(tempMin .GT. tempSingleAtomsDivideArrays(IElement,temp_MaxDivideGroups)) then
                            temp_MaxDivideGroups = temp_MaxDivideGroups + 1
                            tempSingleAtomsDivideArrays(IElement,temp_MaxDivideGroups) = tempMin
                    end if
                END DO

            end if

        END DO

        DivideGroups_SingleElement = 0

        DO IElement = 1,p_ATOMS_GROUPS_NUMBER
            DivideGroups_SingleElement(IElement) = count(tempSingleAtomsDivideArrays(IElement,:) .GT. 0)
        END DO

        !----Start to set barrier for the void range ---------------------
        if(maxval(DivideGroups_SingleElement) .GT. 1) then
            call AllocateArray_Host(AtomsSetsRangesMarkArray,p_ATOMS_GROUPS_NUMBER,maxval(DivideGroups_SingleElement)-1,"AtomsSetsRangesMarkArray")

            AtomsSetsRangesMarkArray = 0

            DO tempIndex = 1,2*ListCount

                DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                    DO J = 1,DivideGroups_SingleElement(IElement)-1
                        if((tempSingleAtomsDivideArrays(IElement,J) + 1) .LT. tempSingleAtomsDivideArrays(IElement,J+1)) then
                            if(TheAtomsSetsRangesArray(tempIndex)%m_SetsRange(IElement)%m_NA_From .LT. (tempSingleAtomsDivideArrays(IElement,J) + 1) .AND. &
                                TheAtomsSetsRangesArray(tempIndex)%m_SetsRange(IElement)%m_NA_To .GT. (tempSingleAtomsDivideArrays(IElement,J) + 1)) then
                                AtomsSetsRangesMarkArray(IElement,J) = 1
                            end if
                        else
                            AtomsSetsRangesMarkArray(IElement,J) = 1
                        end if
                    END DO

                END DO

            END DO

            tempDivideGroups_SingleElement = DivideGroups_SingleElement

            DO IElement = 1,p_ATOMS_GROUPS_NUMBER
                DO J = 1,DivideGroups_SingleElement(IElement)-1
                    if(AtomsSetsRangesMarkArray(IElement,J) .LE. 0) then
                        if((tempSingleAtomsDivideArrays(IElement,J) + 2) .LT. tempSingleAtomsDivideArrays(IElement,J+1)) then
                            tempDivideGroups_SingleElement(IElement) = tempDivideGroups_SingleElement(IElement) + 2
                        else
                            tempDivideGroups_SingleElement(IElement) = tempDivideGroups_SingleElement(IElement) + 1
                        end if
                    end if
                END DO

            END DO

            !---The left start of the range array should be considered----
            tempDivideGroups_SingleElement = tempDivideGroups_SingleElement + 1

            if(maxval(DivideGroups_SingleElement) .GT. 0) then

                call AllocateArray_Host(SingleAtomsDivideArrays,p_ATOMS_GROUPS_NUMBER,maxval(tempDivideGroups_SingleElement) ,"SingleAtomsDivideArrays")

                DO IElement = 1,p_ATOMS_GROUPS_NUMBER

                    tempIndex = 1
                    !---The left start of the range array should be considered----
                    SingleAtomsDivideArrays(IElement,tempIndex) = max(tempSingleAtomsDivideArrays(IElement,1)-1,0)

                    DO I = 1,DivideGroups_SingleElement(IElement)-1

                        if(tempSingleAtomsDivideArrays(IElement,I) .LE. 0) then
                            exit
                        end if
                        tempIndex = tempIndex + 1
                        SingleAtomsDivideArrays(IElement,tempIndex) = tempSingleAtomsDivideArrays(IElement,I)

                        if(AtomsSetsRangesMarkArray(IElement,I) .LE. 0) then
                            if((tempSingleAtomsDivideArrays(IElement,I) + 2) .LT. tempSingleAtomsDivideArrays(IElement,I+1)) then
                                SingleAtomsDivideArrays(IElement,tempIndex+1) = tempSingleAtomsDivideArrays(IElement,I) + 1
                                SingleAtomsDivideArrays(IElement,tempIndex+2) = tempSingleAtomsDivideArrays(IElement,I+1) - 1
                                tempIndex = tempIndex + 2
                            else
                                SingleAtomsDivideArrays(IElement,tempIndex+1) = tempSingleAtomsDivideArrays(IElement,I) + 1
                                tempIndex = tempIndex + 1
                            end if
                        end if
                    END DO

                    SingleAtomsDivideArrays(IElement,tempIndex+1:maxval(tempDivideGroups_SingleElement)) = maxval(tempSingleAtomsDivideArrays(IElement,:))
                END DO
            end if

        else if(maxval(DivideGroups_SingleElement) .eq. 1) then
            call AllocateArray_Host(SingleAtomsDivideArrays,p_ATOMS_GROUPS_NUMBER,3,"SingleAtomsDivideArrays")

            DO IElement = 1,p_ATOMS_GROUPS_NUMBER

                if(UsedAtomsType(IElement) .GT. 0) then
                    SingleAtomsDivideArrays(IElement,1) = max(tempSingleAtomsDivideArrays(IElement,1)-1,0)
                    SingleAtomsDivideArrays(IElement,2) = tempSingleAtomsDivideArrays(IElement,1)
                    SingleAtomsDivideArrays(IElement,3) = tempSingleAtomsDivideArrays(IElement,1) + 1
                end if
            END DO
        else
            call AllocateArray_Host(SingleAtomsDivideArrays,p_ATOMS_GROUPS_NUMBER,1,"SingleAtomsDivideArrays")
            SingleAtomsDivideArrays = 0
        end if

        !---Start to construct the reactions map---
        allocate(ConstructClusterListsArray(2*ListCount))

        Maplength = 0

        DO tempIndex = 1,2*ListCount

            call ConstructClusterListsArray(tempIndex)%Clean_ClusterList()

            ConstructClusterListsArray(tempIndex) = TheAtomsSetsRangesArray(tempIndex)%AtomsSetRange2ClusterList(SingleAtomsDivideArrays)

            Maplength = Maplength + ConstructClusterListsArray(tempIndex)%GetList_Count()
        END DO

        call TheReactionsMap%constructor(SingleAtomsDivideArrays,Maplength)

        !*********Note: The array TheAtomsSetsRangesArray need to response to this(ReadReactionPropList) two by one, so, please ensure*********
        !*********The ConstructClusterListsArray is not modified after resolved from ReadReactionPropList***********
        cursor=>this
        tempIndex = 0
        DO While(associated(cursor))

            if(ConstructClusterListsArray(tempIndex*2+1)%GetList_Count() .GT. 0 .AND. ConstructClusterListsArray(tempIndex*2+2)%GetList_Count() .GT. 0) then
                SubjectClusterListCursor=>ConstructClusterListsArray(tempIndex*2+1)
                DO While(associated(SubjectClusterListCursor))

                    ObjectClusterListCursor=>ConstructClusterListsArray(tempIndex*2+2)

                    DO While(associated(ObjectClusterListCursor))

                        call TheReactionsMap%put(SubjectClusterListCursor%TheCluster,ObjectClusterListCursor%TheCluster,cursor%Reaction%Convert2ReactionValue())

                        ObjectClusterListCursor=>ObjectClusterListCursor%next
                    END DO

                    SubjectClusterListCursor=>SubjectClusterListCursor%next
                END DO
            end if

            tempIndex = tempIndex + 1

            cursor=>cursor%next
        END DO

        call DeAllocateArray_Host(RangesArray,"RangesArray")

        call DeAllocateArray_Host(tempSingleAtomsDivideArrays,"tempSingleAtomsDivideArrays")

        call DeAllocateArray_Host(SingleAtomsDivideArrays,"SingleAtomsDivideArrays")

        call DeAllocateArray_Host(AtomsSetsRangesMarkArray,"AtomsSetsRangesMarkArray")

        if(allocated(TheAtomsSetsRangesArray)) then
            deallocate(TheAtomsSetsRangesArray)
        end if

        if(allocated(ConstructClusterListsArray)) then
            deallocate(ConstructClusterListsArray)
        end if
        Nullify(ConstructClusterListsArray)
        ConstructClusterListsArray=>null()

        Nullify(cursor)
        cursor=>null()

        Nullify(SubjectClusterListCursor)
        SubjectClusterListCursor=>null()
        Nullify(ObjectClusterListCursor)
        ObjectClusterListCursor=>null()

        return
    end subroutine ConvertToReactionsMap

    !***************************************
    subroutine CopyReadReactionPropListFromOther(this,otherOne)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),intent(out),target::this
        type(ReadReactionPropList),target,intent(in)::otherOne
        !---Local Vars---
        type(ReadReactionPropList),pointer::cursorOfOthers=>null()
        type(ReadReactionPropList),pointer::cursorOfSelf=>null()
        type(ReadReactionPropList),pointer::cursorOfSelfP=>null()
        !---Body---
        cursorOfSelf=>this
        if(.not. associated(cursorOfSelf)) then
            write(*,*) "MCPSCUERROR: You need to allocate the ReadReactionPropList first!"
            pause
            stop
        end if

        cursorOfOthers=>otherOne
        if(.not. associated(cursorOfOthers)) then
            Nullify(cursorOfSelf)
            return
        end if

        call this%Clean_ReadReactionPropList()

        this%Reaction = otherOne%Reaction

        cursorOfOthers=>otherOne%next
        cursorOfSelfP=>this
        cursorOfSelf=>this%next
        DO While(associated(cursorOfOthers))
            allocate(cursorOfSelf)
            cursorOfSelf%Reaction = cursorOfOthers%Reaction
            cursorOfSelfP%next=>cursorOfSelf

            cursorOfOthers=>cursorOfOthers%next
            cursorOfSelfP=>cursorOfSelfP%next
            cursorOfSelf=>cursorOfSelf%next
        END DO
        this%ListCount = otherOne%GetList_Count()

        Nullify(cursorOfSelfP)
        Nullify(cursorOfSelf)
        Nullify(cursorOfOthers)
        return
    end subroutine CopyReadReactionPropListFromOther

    !***************************************
    subroutine AppendOne_ReadReactionPropList(this,newOne)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),target::this
        type(ReadReactionPair)::newOne
        !---Local Vars---
        type(ReadReactionPropList),pointer::cursor=>null(),cursorP=>null()
        !---Body---
        cursorP=>this
        if(.not. associated(cursorP)) then
            write(*,*) "MCPSCUERROR: You need to allocate the ReadReactionPropList first!"
            pause
            stop
        end if

        if(this%GetList_Count() .LE. 0) then
            this%ListCount = 1
            this%Reaction = newOne
        else
            cursor=>this%next
            cursorP=>this
            if(IsStrEqual(cursorP%Reaction%SubjectSymbol,newOne%SubjectSymbol) .AND. IsStrEqual(cursorP%Reaction%ObjectSymbol,newOne%ObjectSymbol)) then
                write(*,*) "MCPSCUERROR: The Reaction is dumplicated: Subject: ",newOne%SubjectSymbol," object: ",newOne%ObjectSymbol
                pause
                stop
            end if

            if(IsStrEqual(cursorP%Reaction%SubjectSymbol,newOne%ObjectSymbol) .AND. IsStrEqual(cursorP%Reaction%ObjectSymbol,newOne%SubjectSymbol)) then
                write(*,*) "MCPSCUERROR: The Reaction is dumplicated for pairs Subject : ",cursorP%Reaction%SubjectSymbol," object: ",cursorP%Reaction%ObjectSymbol
                write(*,*) "and : Subject: ",newOne%SubjectSymbol," object: ",newOne%ObjectSymbol
                pause
                stop
            end if

            DO while(associated(cursor))
                cursor=>cursor%next
                cursorP=>cursorP%next

                if(IsStrEqual(cursorP%Reaction%SubjectSymbol,newOne%SubjectSymbol) .AND. IsStrEqual(cursorP%Reaction%ObjectSymbol,newOne%ObjectSymbol)) then
                    write(*,*) "MCPSCUERROR: The Reaction is dumplicated: Subject: ",newOne%SubjectSymbol," object: ",newOne%ObjectSymbol
                    pause
                    stop
                end if

                if(IsStrEqual(cursorP%Reaction%SubjectSymbol,newOne%ObjectSymbol) .AND. IsStrEqual(cursorP%Reaction%ObjectSymbol,newOne%SubjectSymbol)) then
                    write(*,*) "MCPSCUERROR: The Reaction is dumplicated for pairs Subject : ",cursorP%Reaction%SubjectSymbol," object: ",cursorP%Reaction%ObjectSymbol
                    write(*,*) "and : Subject: ",newOne%SubjectSymbol," object: ",newOne%ObjectSymbol
                    pause
                    stop
                end if

            END DO

            this%ListCount = this%ListCount + 1

            allocate(cursor)
            NUllify(cursor%next)
            ! The assignment(=) had been overrided
            cursor%Reaction = newOne
            cursorP%next=>cursor
        end if

        Nullify(cursorP)
        cursorP=>null()
        Nullify(cursor)
        cursor=>null()
        return
    end subroutine AppendOne_ReadReactionPropList


    !***************************************
    subroutine AppendArray_ReadReactionPropList(this,ReactionsArray,ArraySize)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),target::this
        type(ReadReactionPair),allocatable::ReactionsArray(:)
        integer,intent(in)::ArraySize
        !---Local Vars---
        integer::I
        type(ReadReactionPropList),pointer::cursor=>null(),cursorP=>null()
        !---Body---
        cursorP=>this
        if(.not. associated(cursorP)) then
            write(*,*) "MCPSCUERROR: You need to allocate the ReadReactionPropList first!"
            pause
            stop
        end if

        if(ArraySize  .LE. 0 .or. size(ReactionsArray) .LE. 0) then
            write(*,*) "MCPSCUWARNING: No reaction pair would be appended to ReadReactionPropList"
            return
        end if

        if(ArraySize .GT. size(ReactionsArray)) then
            write(*,*) "MCPSCUERROR: The aimmed size to appended to the ReadReactionPropList is greater than the Array size",ArraySize,size(ReactionsArray)
            pause
            stop
        end if


        DO I=1,ArraySize
            call this%AppendOne_ReadReactionPropList(ReactionsArray(I))
        END DO

        return
    end subroutine AppendArray_ReadReactionPropList

    !**************************************
    function GetReadReactionByListIndex(this,ListIndex) result(Reaction)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),target::this
        integer,intent(in)::ListIndex
        type(ReadReactionPair),intent(out)::Reaction
        !---Local Vars---
        integer::tempIndex
        type(ReadReactionPropList),pointer::cursor=>null()
        !---Body---
        cursor=>this
        if(.not. associated(cursor)) then
            write(*,*) "MCPSCUERROR: You need to allocate the ReadReactionPropList first!"
            pause
            stop
        end if

        tempIndex = 0

        DO while(associated(cursor))

            tempIndex = tempIndex + 1

            if(tempIndex .eq. ListIndex) then
                Reaction = cursor%Reaction
                exit
            end if

            cursor=>cursor%next

        END DO

        if(ListIndex .ne. tempIndex) then
            write(*,*) "MCPSCUERROR: Cannot get the reaction pair form reaction list by index: ",ListIndex
            pause
            stop
        end if

        Nullify(cursor)
        cursor=>null()
        return
    end function GetReadReactionByListIndex


    !**************************************
    integer function GetReadReactionPropList_Count(this)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),target::this
        !---Local Vars---
        type(ReadReactionPropList),pointer::cursor=>null()
        !---Body---
        cursor=>this
        if(.not. associated(cursor)) then
            write(*,*) "MCPSCUERROR: You need to allocate the ReadReactionPropList first!"
            pause
            stop
        end if


        GetReadReactionPropList_Count = this%ListCount

        return
    end function GetReadReactionPropList_Count

    !**************************************
    subroutine PrintOutCheckingResult(this,hFile,BasicAtomsList,TheReactionsMap)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),intent(in),target::this
        integer,intent(in)::hFile
        type(AtomsList),intent(in)::BasicAtomsList
        type(ReactionsMap),intent(in)::TheReactionsMap
        !---Local Vars---
        type(AtomsSetRange),dimension(:),allocatable::TheAtomsSetsRangesArray
        type(ReadReactionPropList),pointer::cursor=>null()
        type(AClusterList),target::SubjectConstructClusterList
        type(AClusterList),pointer::SubjectClusterListCursor=>null()
        type(AClusterList),target::ObjectConstructClusterList
        type(AClusterList),pointer::ObjectClusterListCursor=>null()
        type(ReactionValue)::TheValue
        integer::ListCount
        integer::tempIndex
        character*32::SubjectSymbol
        character*32::ObjectSymbol
        character*32::CNum
        character*32::CElements
        integer::IAtomsGroup
        integer::I
        !---Body---
        cursor=>this

        if(.not. associated(cursor)) then
            return
        end if

        ListCount = this%GetList_Count()

        if(ListCount .LE. 0) then
            return
        else
            allocate(TheAtomsSetsRangesArray(2*ListCount))
        end if

        tempIndex = 0

        DO While(associated(cursor))

            tempIndex = tempIndex + 1

            call TheAtomsSetsRangesArray(tempIndex)%ReleaseSetsRange()

            TheAtomsSetsRangesArray(tempIndex) = ResolveSymbol2AtomsSetRange(cursor%Reaction%SubjectSymbol,BasicAtomsList)

            tempIndex = tempIndex + 1

            call TheAtomsSetsRangesArray(tempIndex)%ReleaseSetsRange()

            TheAtomsSetsRangesArray(tempIndex) = ResolveSymbol2AtomsSetRange(cursor%Reaction%ObjectSymbol,BasicAtomsList)

            cursor=>cursor%next
        END DO

        cursor=>this
        tempIndex = 0
        DO While(associated(cursor))

            call SubjectConstructClusterList%Clean_ClusterList()

            SubjectConstructClusterList = TheAtomsSetsRangesArray(tempIndex*2+1)%AtomsSetRange2ClusterList(TheReactionsMap%SingleAtomsDivideArrays)

            call ObjectConstructClusterList%Clean_ClusterList()

            ObjectConstructClusterList = TheAtomsSetsRangesArray(tempIndex*2+2)%AtomsSetRange2ClusterList(TheReactionsMap%SingleAtomsDivideArrays)

            write(*,*) "######################################################################"

            write(hFile,fmt="('!','The reaction subject symbol =',A20,2x,  &
                              '!','The reaction object symbol =',A20,2x,  &
                              '!','CoefficentsGenerate way =',I1,2x, &
                              '!','Reaction Coefficents value =',1F10.4,2x, &
                              '!','PreFactor = ',1PE10.4,2x, &
                              '!','ActEnergy = ',1PE10.4,2x, &
                              '!','ECR Generate way =',I1,2x, &
                              '!','ECR Value =',1PE10.4)")                     cursor%Reaction%SubjectSymbol, &
                                                                               cursor%Reaction%ObjectSymbol, &
                                                                               cursor%Reaction%ReactionCoefficientType, &
                                                                               cursor%Reaction%ReactionCoefficient_Value,  &
                                                                               cursor%Reaction%PreFactor, &
                                                                               cursor%Reaction%ActEnergy, &
                                                                               cursor%Reaction%ECRValueType, &
                                                                               cursor%Reaction%ECR

            if(SubjectConstructClusterList%GetList_Count() .GT. 0 .AND. ObjectConstructClusterList%GetList_Count()) then
                SubjectClusterListCursor=>SubjectConstructClusterList
                DO While(associated(SubjectClusterListCursor))

                    ObjectClusterListCursor=>ObjectConstructClusterList

                    DO While(associated(ObjectClusterListCursor))

                        TheValue = TheReactionsMap%get(SubjectClusterListCursor%TheCluster,ObjectClusterListCursor%TheCluster)

                        SubjectSymbol = ""
                        ObjectSymbol = ""
                        DO IAtomsGroup = 1,p_ATOMS_GROUPS_NUMBER
                            CNum = ""
                            call BasicAtomsList%GetSymbolByIndex(IAtomsGroup,CElements)
                            CElements = adjustl(CElements)
                            write(CNum,*) SubjectClusterListCursor%TheCluster%m_Atoms(IAtomsGroup)%m_NA
                            CNum = adjustl(CNum)
                            SubjectSymbol = SubjectSymbol(1:LENTRIM(SubjectSymbol))//CElements(1:LENTRIM(CElements))//"#"//CNum(1:LENTRIM(CNum))

                            write(CNum,*) ObjectClusterListCursor%TheCluster%m_Atoms(IAtomsGroup)%m_NA
                            CNum = adjustl(CNum)
                            ObjectSymbol = ObjectSymbol(1:LENTRIM(ObjectSymbol))//CElements(1:LENTRIM(CElements))//"#"//CNum(1:LENTRIM(CNum))

                            if(IAtomsGroup .LT. p_ATOMS_GROUPS_NUMBER) then
                                SubjectSymbol = SubjectSymbol(1:LENTRIM(SubjectSymbol))//"@"
                                ObjectSymbol = ObjectSymbol(1:LENTRIM(ObjectSymbol))//"@"
                            end if
                        END DO

                        SubjectSymbol = adjustl(SubjectSymbol)
                        ObjectSymbol = adjustl(ObjectSymbol)

                        write(hFile,fmt="('! |--','The subject symbol =',A20,2x,  &
                                          '!','The object symbol =',A20,2x,  &
                                          '!','CoefficentsGenerate way =',I1,2x, &
                                          '!','Reaction Coefficents value =',1F10.4,2x, &
                                          '!','PreFactor = ',1PE10.4,2x, &
                                          '!','ActEnergy = ',1PE10.4,2x, &
                                          '!','ECR Generate way =',I1,2x, &
                                          '!','ECR Value =',1PE10.4)")         SubjectSymbol,                       &
                                                                               ObjectSymbol,                        &
                                                                               TheValue%ReactionCoefficientType,    &
                                                                               TheValue%ReactionCoefficient_Value,  &
                                                                               TheValue%PreFactor,                  &
                                                                               TheValue%ActEnergy,                  &
                                                                               TheValue%ECRValueType,               &
                                                                               TheValue%ECR

                        ObjectClusterListCursor=>ObjectClusterListCursor%next
                    END DO

                    SubjectClusterListCursor=>SubjectClusterListCursor%next
                END DO
            end if

            tempIndex = tempIndex + 1

            cursor=>cursor%next
        END DO

        if(allocated(TheAtomsSetsRangesArray)) then
            deallocate(TheAtomsSetsRangesArray)
        end if

        call SubjectConstructClusterList%Clean_ClusterList()
        Call ObjectConstructClusterList%Clean_ClusterList()

        Nullify(cursor)
        Nullify(SubjectClusterListCursor)
        Nullify(ObjectClusterListCursor)
        return
    end subroutine PrintOutCheckingResult

    !**************************************
    subroutine Clean_ReadReactionPropList(this)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),target::this
        !---Local Vars---
        type(ReadReactionPropList),pointer::cursor=>null()
        type(ReadReactionPropList),pointer::next=>null()
        !---Body---
        cursor=>this

        if(.not. associated(cursor)) then
            return
        end if

        if(cursor%GetList_Count() .LE. 0) then
            return
        end if

        cursor=>this%next

        call CleanReadedReactionPair(this%Reaction)

        DO While(associated(cursor))
            next=>cursor%next
            call CleanReadedReactionPair(cursor%Reaction)
            deallocate(cursor)
            Nullify(cursor)
            cursor=>next
        END DO

        this%next=>null()

        this%ListCount = 0

        Nullify(cursor)
        Nullify(next)
        cursor=>null()
        next=>null()

        return
    end subroutine Clean_ReadReactionPropList

    !************************************
    subroutine CleanReadReactionPropList(this)
        implicit none
        !---Dummy Vars---
        type(ReadReactionPropList)::this
        !---Body---

        call this%Clean_ReadReactionPropList()

        return
    end subroutine CleanReadReactionPropList

    !**************************************************
    subroutine ResloveReactionsFromCScript(scriptStr,ReactionProp_List)
        implicit none
        !---Dummy Vars---
        character(kind=c_char,len=10000)::scriptStr
        type(ReadReactionPropList)::ReactionProp_List
        !---Local Vars---
        type(ReadReactionPair),allocatable::FReactionsDefArray(:)
        integer(kind=c_int)::ArraySize
        !---Body---

        ArraySize = InterpCScript_ReactionsDef(scriptStr)

        if(allocated(FReactionsDefArray)) then
            deallocate(FReactionsDefArray)
        end if
        allocate(FReactionsDefArray(ArraySize))

        call GetInterpedReactionsArray(FReactionsDefArray)


        call AppendArray_ReadReactionPropList(ReactionProp_List,FReactionsDefArray,ArraySize)


        return
    end subroutine ResloveReactionsFromCScript

    !*********************************************************
    function WhetherFreeDiffusion(this,TKB) result(TheResult)
        implicit none
        !---Dummy Vars---
        CLASS(ReadReactionPropList),target::this
        real(kind=KMCDF),intent(in)::TKB
        logical::TheResult
        !---Local Vars---
        type(ReadReactionPropList),pointer::cursor=>null()
        real(kind=KMCDF)::TheValue
        !---Body---

        TheResult = .false.

        cursor=>this

        DO While(associated(cursor))

            select case(cursor%Reaction%ReactionCoefficientType)
                case(p_ReactionCoefficient_ByValue)
                    TheValue = cursor%Reaction%ReactionCoefficient_Value
                case(p_ReactionCoefficient_ByArrhenius)
                    TheValue = cursor%Reaction%PreFactor*exp(-C_EV2ERG*cursor%Reaction%ActEnergy/TKB)
                case default
                    write(*,*) "MCPSCUERROR: undefined type of reaction value.",cursor%Reaction%ReactionCoefficientType
                    pause
                    stop
            end select

            if(TheValue .LE. 0.D0) then
                TheResult = .true.
            else
                TheResult = .false.
            end if

            cursor=>cursor%next

        END DO


        cursor=>null()

        return
    end function



end module MCMF_TYPEDEF_ReactionPropList
