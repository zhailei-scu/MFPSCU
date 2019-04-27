module MF_MethodClass_Factory_CPU
!    use MFLIB_TYPEDEF_SIMULATIONBOXARRAY
!    use MFLIB_TYPEDEF_SIMULATIONCTRLPARAM
!    implicit none
!
!    abstract interface
!        subroutine For_One_Test(Host_SimBoxes,Host_SimuCtrlParam,Dev_Boxes,JobIndex)
!            use MCLIB_TYPEDEF_SIMULATIONBOXARRAY
!            use MCLIB_TYPEDEF_SIMULATIONCTRLPARAM
!            use MCLIB_TYPEDEF_SIMULATIONBOXARRAY_GPU
!            use MCLIB_TYPEDEF_USUAL
!            implicit none
!            type(SimulationBoxes)::Host_SimBoxes
!            type(SimulationCtrlParam),target::Host_SimuCtrlParam
!            type(SimulationBoxes_GPU)::Dev_Boxes
!            integer,intent(in)::JobIndex
!        end subroutine
!    end interface
!
!
!    type,public::MFMethodClassCPU
!        character(len=20)::name
!        type(SimulationBoxes),pointer::pSimulationBoxes=>null()
!        type(SimulationCtrlParam),pointer::pSimulationCtrlParam=>null()
!        procedure(For_One_Test),pointer,nopass::ForOneTest=>null()
!
!        contains
!
!        procedure,public,pass,non_overridable::Register_Method_Class
!        procedure,public,pass,non_overridable::Clean_MethodClass
!        Final::CleanMethodClass
!    end type
!
!    private::Register_Method_Class
!    private::Clean_MethodClass
!    private::CleanMethodClass
!
!    contains
!    !*********************************************
!    subroutine Register_Method_Class(this,className,SimBoxes,SimCtrlParams)
!        use MC_Method_MIGCOALE_CLUSTER_GPU, only:For_One_Test_MIGCOALE_CLUSTER_GPU => For_One_Test
!        implicit none
!        !---Dummy Vars---
!        CLASS(MFMethodClassCPU)::this
!        character*(*)::className
!        type(SimulationBoxes),target::SimBoxes
!        type(SimulationCtrlParam),target::SimCtrlParams
!        !---Local Vars---
!        !---Body---
!        this%pSimulationBoxes=>SimBoxes
!        this%pSimulationCtrlParam=>SimCtrlParams
!
!        select case(className(1:LENTRIM(className)))
!            case("MIGCOALE_CLUSTER_GPU")
!                this%name = "MIGCOALE_CLUSTER_GPU"
!                this%ForOneTest=>For_One_Test_MIGCOALE_CLUSTER_GPU
!            case default
!                write(*,*) "MCPSCUERROR: The unknown method name: ",className
!                pause
!                stop
!        end select
!
!        return
!    end subroutine Register_Method_Class
!
!    !**************************************************
!    subroutine Clean_MethodClass(this)
!        implicit none
!        !---Dummy Vars---
!        class(MFMethodClassCPU)::this
!        !---Body---
!        Nullify(this%ForOneTest)
!        Nullify(this%pSimulationBoxes)
!        Nullify(this%pSimulationCtrlParam)
!        this%ForOneTest=>null()
!        this%pSimulationBoxes=>null()
!        this%pSimulationCtrlParam=>null()
!        this%name = ""
!
!        return
!    end subroutine
!
!    !**************************************************
!    subroutine CleanMethodClass(this)
!        implicit none
!        !---Dummy Vars---
!        type(MFMethodClassCPU)::this
!        !---Body---
!        call this%Clean_MethodClass()
!        return
!    end subroutine
!
!
end module MF_MethodClass_Factory_CPU
