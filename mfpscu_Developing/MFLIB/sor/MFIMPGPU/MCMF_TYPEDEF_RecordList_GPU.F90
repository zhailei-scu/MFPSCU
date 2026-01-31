!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MCMF_TYPEDEF_RECORDLIST_GPU
!    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
!    use MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
!    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY_GPU
!    implicit none
!
!    abstract interface
!        subroutine Procedure_OneStep(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TSTEP,SimRecord)
!            use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
!            use MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
!            use MFLIB_TYPEDEF_SIMULATIONBOXARRAY_GPU
!            implicit none
!            type(SimulationBoxes)::Host_Boxes
!            type(SimulationCtrlParam)::Host_SimuCtrlParam
!            type(SimulationBoxes_GPU)::Dev_Boxes
!            real(kind=KMCDF)::TSTEP
!            CLASS(SimulationRecord)::SimRecord
!        end subroutine
!    end interface
!
!    type,public::OneStepProcedureList
!        procedure(Procedure_OneStep),pointer::m_Procedure=>null()
!
!        type(OneStepProcedureList),pointer::next=>null()
!
!        integer::m_ListCount = 0
!        contains
!        procedure,non_overridable,public,pass::DoProcedure=>DoProcedure_OneStepProcedureList
!        procedure,non_overridable,public,pass::CopyFromOther=>CopyOneStepProcedureListFromOther
!        procedure,non_overridable,public,pass::AppendOne=>AppendOne_OneStepProcedureList
!        procedure,non_overridable,public,pass::GetList_Count=>GetOneStepProcedureList_Count
!        procedure,non_overridable,public,pass::Clean=>Clean_OneStepProcedureList
!
!        Final::CleanOneStepProcedureList
!
!    end type OneStepProcedureList
!
!    contains
!
!    !***************************************
!    subroutine DoProcedure_OneStepProcedureList(this,Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TSTEP,SimRecord)
!        implicit none
!        !---Dummy Vars---
!        CLASS(OneStepProcedureList),intent(out),target::this
!        type(SimulationBoxes)::Host_Boxes
!        type(SimulationCtrlParam)::Host_SimuCtrlParam
!        type(SimulationBoxes_GPU)::Dev_Boxes
!        real(kind=KMCDF)::TSTEP
!        CLASS(SimulationRecord)::SimRecord
!        !---Local Vars---
!        type(OneStepProcedureList),pointer::cursor=>null()
!        !---Body---
!
!        cursor=>this
!
!        Do while(associated(cursor))
!            if(associated(cursor%m_Procedure)) then
!                call cursor%m_Procedure(Host_Boxes,Host_SimuCtrlParam,Dev_Boxes,TSTEP,SimRecord)
!            end if
!            cursor=>cursor%next
!        End Do
!
!        Nullify(cursor)
!
!        return
!    end subroutine
!
!    !***************************************
!    subroutine CopyOneStepProcedureListFromOther(this,otherOne)
!        implicit none
!        !---Dummy Vars---
!        CLASS(OneStepProcedureList),intent(out),target::this
!        type(OneStepProcedureList),intent(in)::otherOne
!        !---Local Vars---
!        type(OneStepProcedureList),pointer::cursorOfOthers=>null()
!        type(OneStepProcedureList),pointer::cursorOfSelf=>null()
!        type(OneStepProcedureList),pointer::cursorOfSelfP=>null()
!        !---Body---
!
!        this%m_Procedure=>null()
!
!        this%m_Procedure = otherOne%m_Procedure
!
!        cursorOfOthers=>otherOne%next
!        cursorOfSelfP=>this
!        cursorOfSelf=>this%next
!        DO While(associated(cursorOfOthers))
!            allocate(cursorOfSelf)
!            cursorOfSelf%m_Procedure = cursorOfOthers%m_Procedure
!            cursorOfSelfP%next=>cursorOfSelf
!
!            cursorOfOthers=>cursorOfOthers%next
!            cursorOfSelfP=>cursorOfSelfP%next
!            cursorOfSelf=>cursorOfSelf%next
!        END DO
!        this%m_ListCount = otherOne%GetList_Count()
!
!        Nullify(cursorOfSelfP)
!        Nullify(cursorOfSelf)
!        Nullify(cursorOfOthers)
!        return
!    end subroutine CopyOneStepProcedureListFromOther
!
!    !***************************************
!    subroutine AppendOne_OneStepProcedureList(this,newOne)
!        implicit none
!        !---Dummy Vars---
!        CLASS(OneStepProcedureList),target::this
!        procedure(Procedure_OneStep)::newOne
!        !---Local Vars---
!        type(OneStepProcedureList),pointer::cursor=>null(),cursorP=>null()
!        !---Body---
!        if(this%GetList_Count() .LE. 0) then
!            this%m_ListCount = 1
!            this%m_Procedure=>newOne
!        else
!            cursor=>this%next
!            cursorP=>this
!
!            DO while(associated(cursor))
!                cursor=>cursor%next
!                cursorP=>cursorP%next
!            END DO
!
!            this%m_ListCount = this%m_ListCount + 1
!
!            allocate(cursor)
!            NUllify(cursor%next)
!            cursor%m_Procedure=>newOne
!            cursorP%next=>cursor
!        end if
!
!        Nullify(cursorP)
!        cursorP=>null()
!        Nullify(cursor)
!        cursor=>null()
!        return
!    end subroutine AppendOne_OneStepProcedureList
!
!
!    !**************************************
!    integer function GetOneStepProcedureList_Count(this)
!        implicit none
!        !---Dummy Vars---
!        CLASS(OneStepProcedureList)::this
!        !---Body---
!        GetOneStepProcedureList_Count = this%m_ListCount
!
!        return
!    end function
!
!    !**************************************
!    subroutine Clean_OneStepProcedureList(this)
!        implicit none
!        !---Dummy Vars---
!        CLASS(OneStepProcedureList),target::this
!        !---Local Vars---
!        type(OneStepProcedureList),pointer::cursor=>null()
!        type(OneStepProcedureList),pointer::next=>null()
!        !---Body---
!
!        Nullify(this%m_Procedure)
!
!        cursor=>this%next
!
!        DO While(associated(cursor))
!            next=>cursor%next
!            Nullify(cursor%m_Procedure)
!            deallocate(cursor)
!            Nullify(cursor)
!            cursor=>next
!        END DO
!
!        this%next=>null()
!
!        this%m_ListCount = 0
!
!        Nullify(cursor)
!        Nullify(next)
!        cursor=>null()
!        next=>null()
!
!        return
!    end subroutine Clean_OneStepProcedureList
!
!    !************************************
!    subroutine CleanOneStepProcedureList(this)
!        implicit none
!        !---Dummy Vars---
!        type(OneStepProcedureList)::this
!        !---Body---
!
!        call this%Clean()
!
!        return
!    end subroutine CleanOneStepProcedureList
!
!
!
end module MCMF_TYPEDEF_RECORDLIST_GPU
