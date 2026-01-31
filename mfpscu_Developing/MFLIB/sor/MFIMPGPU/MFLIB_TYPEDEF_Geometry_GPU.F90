!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MFLIB_TYPEDEF_GEOMETRY_GPU
    use MFLIB_TYPEDEF_GEOMETRY
    use MCMF_Utilities_GPU
    implicit none

    type,public::GrainBoundary_Dev
        type(GrainSeed),device,dimension(:),allocatable::dm_GrainSeeds

        contains

        procedure,non_overridable,public,pass::InitGrainBoundary_Dev
        procedure,non_overridable,public,pass::CopyInGrainBoundaryFromHost
        procedure,non_overridable,public,pass::CopyOutGrainBoundaryToHost
        procedure,non_overridable,public,pass::Clean_GrainBoundary_Dev
        Final::CleanGB_Dev

    end type GrainBoundary_Dev

    private::InitGrainBoundary_Dev
    private::CopyInGrainBoundaryFromHost
    private::CopyOutGrainBoundaryToHost
    private::Clean_GrainBoundary_Dev
    private::CleanGB_Dev

    contains

    !**************************************
    subroutine InitGrainBoundary_Dev(this,SeedNum)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary_Dev)::this
        integer::SeedNum
        !---Body---
        if(SeedNum .GT. 0) then
            allocate(this%dm_GrainSeeds(SeedNum))
        else
            ! the device array should be at least for 1 size while used in kernel
            allocate(this%dm_GrainSeeds(1))
        end if
        return
    end subroutine InitGrainBoundary_Dev

    !**************************************
    subroutine CopyInGrainBoundaryFromHost(this,Host_GB)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary_Dev)::this
        type(GrainBoundary)::Host_GB
        !---Local Vars---
        integer::SeedNum
        !---Body---

        SeedNum = Host_GB%GrainNum

        if(SeedNum .GT. 0) then
            call copyInGrainSeedsSync(Host_GB%GrainSeeds,this%dm_GrainSeeds,SeedNum)
        end if

        return
    end subroutine

    !**************************************
    subroutine CopyOutGrainBoundaryToHost(this,Host_GB)
        implicit none
        !---Dummy Vars---
        CLASS(GrainBoundary_Dev)::this
        type(GrainBoundary)::Host_GB
        !---Local Vars---
        integer::SeedNum
        !---Body---

        SeedNum = Host_GB%GrainNum

        if(SeedNum .GT. 0) then
            call copyOutGrainSeedsSync(Host_GB%GrainSeeds,this%dm_GrainSeeds,SeedNum)
        end if

        return
    end subroutine CopyOutGrainBoundaryToHost

    !***************************************
    subroutine Clean_GrainBoundary_Dev(this)
        !---Dummy Vars---
        CLASS(GrainBoundary_Dev)::this
        !---Body---

        if(allocated(this%dm_GrainSeeds)) then
            deallocate(this%dm_GrainSeeds)
        end if

        return
    end subroutine Clean_GrainBoundary_Dev

    !***************************************
    subroutine CleanGB_Dev(this)
        !---Dummy Vars---
        type(GrainBoundary_Dev)::this
        !---Body---

        call this%Clean_GrainBoundary_Dev()

        return
    end subroutine CleanGB_Dev

    !********************************************************
    attributes(device) function GrainBelongsTo_Dev(NSeeds,Dev_GrainSeeds,POS) result(SeedID)
        !---Dummy Vars---
        integer,value::NSeeds
        type(GrainSeed),device::Dev_GrainSeeds(*) ! When the nollvm compiler option is used, the attributes(device) dummy vars array should write as (*) for one dimension, cannot be (:)
        real(kind=KMCDF)::POS(3)
        integer::SeedID
        !---Local Vars---
        integer::ISeed
        real(kind=KMCSF)::SEP(3)
        real(kind=KMCDF)::DIST2
        real(kind=KMCDF)::MinDist2
        !---Body--

        MinDist2 = 1.D32

        SeedID = 0

        DO ISeed = 1,NSeeds
            SEP = POS - Dev_GrainSeeds(ISeed)%m_POS

            if(ABS(SEP(1)) .GT. dm_HBOXSIZE(1) .AND. dm_PERIOD(1) .GT. 0) then
                SEP(1) = SEP(1) - SIGN(dm_BOXSIZE(1),SEP(1))
            end if

            if(ABS(SEP(2)) .GT. dm_HBOXSIZE(2) .AND. dm_PERIOD(2) .GT. 0) then
                SEP(2) = SEP(2) - SIGN(dm_BOXSIZE(2),SEP(2))
            end if

            if(ABS(SEP(3)) .GT. dm_HBOXSIZE(3) .AND. dm_PERIOD(3) .GT. 0) then
                SEP(3) = SEP(3) - SIGN(dm_BOXSIZE(3),SEP(3))
            end if

            DIST2 = SEP(1)*SEP(1) + SEP(2)*SEP(2) + SEP(3)*SEP(3)

            if(DIST2 .LT. MinDist2) then
                MinDist2 = DIST2
                SeedID = ISeed
            end if

        END DO

        return
    end function GrainBelongsTo_Dev

    !********************************************************
    attributes(device) subroutine CalCrossPointInGB_Dev(Dev_GrainSeeds,Seed1,Seed2,POS1,POS2,CrossPoint)
        !---Dummy Vars---
        type(GrainSeed),device::Dev_GrainSeeds(*) ! When the nollvm compiler option is used, the attributes(device) dummy vars array should write as (*) for one dimension,cannot be (:)
        integer,intent(in)::Seed1
        integer,intent(in)::Seed2
        real(kind=KMCDF),intent(in)::POS1(3)
        real(kind=KMCDF),intent(in)::POS2(3)
        real(kind=KMCDF),intent(out)::CrossPoint(3) !Note: If the nollvm compiler option is used, we cannot apply the array as a function result, so this array result
                                                    !      must be an dummy array of the subroutine
        !---Local Vars---
        real(kind=KMCSF)::GrainSeed1(3)
        real(kind=KMCSF)::GrainSeed2(3)
        real(kind=KMCSF)::SEP(3)
        real(kind=KMCDF)::halfPoint(3)
        real(kind=KMCDF)::GBNormVector(3)
        real(kind=KMCDF)::K1
        real(kind=KMCDF)::K2
        !---Body--
        GrainSeed1 = Dev_GrainSeeds(Seed1)%m_POS

        SEP = GrainSeed1 - POS1
        if(ABS(SEP(1)) .GT. dm_HBOXSIZE(1) .AND. dm_PERIOD(1) .GT. 0) then
            GrainSeed1(1) = GrainSeed1(1) - SIGN(dm_BOXSIZE(1),SEP(1))
        end if
        if(ABS(SEP(2)) .GT. dm_HBOXSIZE(2) .AND. dm_PERIOD(2) .GT. 0) then
            GrainSeed1(2) = GrainSeed1(2) - SIGN(dm_BOXSIZE(2),SEP(2))
        end if
        if(ABS(SEP(3)) .GT. dm_HBOXSIZE(3) .AND. dm_PERIOD(3) .GT. 0) then
            GrainSeed1(3) = GrainSeed1(3) - SIGN(dm_BOXSIZE(3),SEP(3))
        end if

        GrainSeed2 = Dev_GrainSeeds(Seed2)%m_POS
        SEP = GrainSeed2 - POS1
        if(ABS(SEP(1)) .GT. dm_HBOXSIZE(1) .AND. dm_PERIOD(1) .GT. 0) then
            GrainSeed2(1) = GrainSeed2(1) - SIGN(dm_BOXSIZE(1),SEP(1))
        end if
        if(ABS(SEP(2)) .GT. dm_HBOXSIZE(2) .AND. dm_PERIOD(2) .GT. 0) then
            GrainSeed2(2) = GrainSeed2(2) - SIGN(dm_BOXSIZE(2),SEP(2))
        end if
        if(ABS(SEP(3)) .GT. dm_HBOXSIZE(3) .AND. dm_PERIOD(3) .GT. 0) then
            GrainSeed2(3) = GrainSeed2(3) - SIGN(dm_BOXSIZE(3),SEP(3))
        end if

        halfPoint = (GrainSeed1 + GrainSeed2)/2

        GBNormVector = GrainSeed1 - GrainSeed2

        K1 = DABS(GBNormVector(1)*(halfPoint(1) - POS1(1)) + GBNormVector(2)*(halfPoint(2) - POS1(2)) + GBNormVector(3)*(halfPoint(3) - POS1(3)))

        K2 = K1 + DABS(GBNormVector(1)*(halfPoint(1) - POS2(1)) + GBNormVector(2)*(halfPoint(2) - POS2(2)) + GBNormVector(3)*(halfPoint(3) - POS2(3)))


        CrossPoint = POS1 + (POS2 - POS1)*(K1/K2)

        return
    end subroutine CalCrossPointInGB_Dev



end module MFLIB_TYPEDEF_GEOMETRY_GPU
