!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MIGCOALE_TYPEDEF_STATISTICINFO

    use MCMF_CONSTANTS
    use MCMF_UTILITIES

    implicit none

    type,public::MigCoaleStatisticOneBox

        integer::ICMAX(p_NUMBER_OF_STATU) = 0                               ! the ID of the largest cluster
        real(kind=KMCDF)::RMAX(p_NUMBER_OF_STATU) = 0.D0                    ! the radius of the largest cluster at current time intervale
        real(kind=KMCDF)::RMIN(p_NUMBER_OF_STATU) = 0.D0                    ! the radius of the smallest cluster at current time intervale, to be used to determine the time step
        real(kind=KMCDF)::RAVA(p_NUMBER_OF_STATU) = 0.D0                    ! the average radius of the clusters at current time intervale

        real(kind=KMCDF)::DiffusorValueMax(p_NUMBER_OF_STATU) = 0.D0

        contains
        procedure,non_overridable,public,pass::InitStatisticInfo
        procedure,non_overridable,public,pass::Clean_StatisticInfo
        procedure,non_overridable,private,pass::CopyMigCoaleStatisticOneBoxFromOther
        Generic::Assignment(=)=>CopyMigCoaleStatisticOneBoxFromOther
        Final::CleanStatisticInfo

    end type MigCoaleStatisticOneBox

    type,public::MigCoaleStatisticInfo
        type(MigCoaleStatisticOneBox),dimension(:),allocatable::statistic_SingleBoxes
        type(MigCoaleStatisticOneBox)::statistic_IntegralBox

        contains
        procedure,non_overridable,public,pass::Init=>InitMigCoaleStatisticInfo
        procedure,non_overridable,public,pass::Clean=>Clean_MigCoaleStatisticInfo
        procedure,non_overridable,private,pass::CopyMigCoaleStatisticInfoFromOther
        Generic::Assignment(=)=>CopyMigCoaleStatisticInfoFromOther
        Final::CleanMigCoaleStatisticInfo
    end type MigCoaleStatisticInfo

    type,public,extends(MigCoaleStatisticInfo)::MigCoaleStatisticInfo_Used

    end type

    type,public,extends(MigCoaleStatisticInfo)::MigCoaleStatisticInfo_Expd
        contains
        procedure,non_overridable,public,pass::ConverFromUsed=>ConvertUsedToExpd
    end type

    type,public,extends(MigCoaleStatisticInfo)::MigCoaleStatisticInfo_Virtual
        contains
        procedure,non_overridable,public,pass::ConverFromUsed=>ConvertUsedToVirtual
    end type


    type,public::MigCoaleStatInfoWrap
        type(MigCoaleStatisticInfo_Used)::m_MigCoaleStatisticInfo_Used
        type(MigCoaleStatisticInfo_Expd)::m_MigCoaleStatisticInfo_Expd
        type(MigCoaleStatisticInfo_Virtual)::m_MigCoaleStatisticInfo_Virtual

        contains
        procedure,non_overridable,public,pass::Init=>InitMigCoaleStatInfoWrap
        procedure,non_overridable,public,pass::Clean=>Clean_MigCoaleStatInfoWrap
        Final::CleanMigCoaleStatInfoWrap
    end type


    private::InitStatisticInfo
    private::CopyMigCoaleStatisticOneBoxFromOther
    private::Clean_StatisticInfo
    private::CleanStatisticInfo
    private::InitMigCoaleStatisticInfo
    private::CopyMigCoaleStatisticInfoFromOther
    private::Clean_MigCoaleStatisticInfo
    private::CleanMigCoaleStatisticInfo
    private::ConvertUsedToExpd
    private::ConvertUsedToVirtual
    private::InitMigCoaleStatInfoWrap
    private::Clean_MigCoaleStatInfoWrap
    private::CleanMigCoaleStatInfoWrap

    contains

    !****************Type MigCoaleStatisticOneBox*****************************
    subroutine InitStatisticInfo(this)
        implicit none
        !---Dummy Vars---
        Class(MigCoaleStatisticOneBox)::this
        !---Body---

        this%ICMAX = 0
        this%RMAX = 0.D0
        this%RMIN = 0.D0
        this%RAVA = 0.D0

        this%DiffusorValueMax = 0.D0

        return
    end subroutine InitStatisticInfo

    !***********************************************
    subroutine CopyMigCoaleStatisticOneBoxFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoaleStatisticOneBox),intent(out)::this
        TYPE(MigCoaleStatisticOneBox),intent(in)::other
        !---Body---
        this%RMAX = other%RMAX
        this%RMIN = other%RMIN
        this%RAVA = other%RAVA

        this%ICMAX = other%ICMAX

        this%DiffusorValueMax = other%DiffusorValueMax

        return
    end subroutine

    !**********************************************
    subroutine Clean_StatisticInfo(this)
        implicit none
        !---Dummy Vars---
        Class(MigCoaleStatisticOneBox)::this
        !---Body---

        this%ICMAX = 0
        this%RMAX = 0.D0
        this%RMIN = 0.D0
        this%RAVA = 0.D0

        this%DiffusorValueMax = 0.D0

        return
    end subroutine Clean_StatisticInfo

    subroutine CleanStatisticInfo(this)
        implicit none
        !---Dummy Vars---
        type(MigCoaleStatisticOneBox)::this
        !---Body---

        call this%Clean_StatisticInfo()

        return
    end subroutine CleanStatisticInfo

    !****************Type MigCoaleStatisticOneBox*****************************
    subroutine InitMigCoaleStatisticInfo(this,MultiBox)
        implicit none
        !---Dummy Vars---
        Class(MigCoaleStatisticInfo)::this
        integer,intent(in)::MultiBox
        !---Local Vars---
        integer::IBox
        !---Body---
        if(allocated(this%statistic_SingleBoxes)) then
            deallocate(this%statistic_SingleBoxes)
        end if

        allocate(this%statistic_SingleBoxes(MultiBox))

        DO IBox = 1,MultiBox
            call this%statistic_SingleBoxes(IBox)%InitStatisticInfo()
        END DO

        call this%statistic_IntegralBox%InitStatisticInfo()

        return
    end subroutine

    !***********************************************
    subroutine CopyMigCoaleStatisticInfoFromOther(this,other)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoaleStatisticInfo),intent(out)::this
        TYPE(MigCoaleStatisticInfo),intent(in)::other
        !---Local Vars---
        integer::I
        !---Body---
        this%statistic_IntegralBox = other%statistic_IntegralBox

        if(allocated(this%statistic_SingleBoxes)) then
            deallocate(this%statistic_SingleBoxes)
        end if

        if(size(other%statistic_SingleBoxes) .GT. 0) then
             allocate(this%statistic_SingleBoxes(size(other%statistic_SingleBoxes)))

             DO I = 1,size(other%statistic_SingleBoxes)
                this%statistic_SingleBoxes(I) = other%statistic_SingleBoxes(I)
             END DO
        end if

        return
    end subroutine

    subroutine Clean_MigCoaleStatisticInfo(this)
        implicit none
        !---Dummy Vars---
        Class(MigCoaleStatisticInfo)::this
        !---Body---
        if(allocated(this%statistic_SingleBoxes)) then
            deallocate(this%statistic_SingleBoxes)
        end if

        call this%statistic_IntegralBox%Clean_StatisticInfo()

        return
    end subroutine Clean_MigCoaleStatisticInfo

    subroutine CleanMigCoaleStatisticInfo(this)
        implicit none
        !---Dummy Vars---
        type(MigCoaleStatisticInfo)::this
        !---Body---
        call this%Clean()

        return
    end subroutine CleanMigCoaleStatisticInfo

    !****************For type MigCoaleStatisticInfo_Expd*************************
    subroutine ConvertUsedToExpd(this,Used)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoaleStatisticInfo_Expd)::this
        TYPE(MigCoaleStatisticInfo_Used)::Used
        !---Local Vars---
        integer::I
        !---Body---
        this%statistic_IntegralBox = Used%statistic_IntegralBox

        if(allocated(this%statistic_SingleBoxes)) then
            deallocate(this%statistic_SingleBoxes)
        end if

        if(size(Used%statistic_SingleBoxes) .GT. 0) then
             allocate(this%statistic_SingleBoxes(size(Used%statistic_SingleBoxes)))

             DO I = 1,size(Used%statistic_SingleBoxes)
                this%statistic_SingleBoxes(I) = Used%statistic_SingleBoxes(I)
             END DO
        end if

        return
    end subroutine

    !****************For type MigCoaleStatisticInfo_Virtual*************************
    subroutine ConvertUsedToVirtual(this,Used)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoaleStatisticInfo_Virtual)::this
        TYPE(MigCoaleStatisticInfo_Used)::Used
        !---Local Vars---
        integer::I
        !---Body---
        this%statistic_IntegralBox = Used%statistic_IntegralBox

        if(allocated(this%statistic_SingleBoxes)) then
            deallocate(this%statistic_SingleBoxes)
        end if

        if(size(Used%statistic_SingleBoxes) .GT. 0) then
             allocate(this%statistic_SingleBoxes(size(Used%statistic_SingleBoxes)))

             DO I = 1,size(Used%statistic_SingleBoxes)
                this%statistic_SingleBoxes(I) = Used%statistic_SingleBoxes(I)
             END DO
        end if

        return
    end subroutine

    !****************For type MigCoaleStatInfoWrap************************
    subroutine InitMigCoaleStatInfoWrap(this,MultiBox)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoaleStatInfoWrap)::this
        integer,intent(in)::MultiBox
        !---Body---
        call this%m_MigCoaleStatisticInfo_Expd%Init(MultiBox)

        call this%m_MigCoaleStatisticInfo_Used%Init(MultiBox)

        call this%m_MigCoaleStatisticInfo_Virtual%Init(MultiBox)

        return
    end subroutine

    subroutine Clean_MigCoaleStatInfoWrap(this)
        implicit none
        !---Dummy Vars---
        CLASS(MigCoaleStatInfoWrap)::this
        !---Body---
        call this%m_MigCoaleStatisticInfo_Expd%Clean()

        call this%m_MigCoaleStatisticInfo_Used%Clean()

        call this%m_MigCoaleStatisticInfo_Virtual%Clean()

        return
    end subroutine

    subroutine CleanMigCoaleStatInfoWrap(this)
        implicit none
        !---Dummy Vars---
        type(MigCoaleStatInfoWrap)::this
        !---Body---
        call this%Clean()

        return
    end subroutine

end module MIGCOALE_TYPEDEF_STATISTICINFO
