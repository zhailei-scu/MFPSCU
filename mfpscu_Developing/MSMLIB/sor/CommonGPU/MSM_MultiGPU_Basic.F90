!*********************************************************************************!
!--- Description:
!--- Author : Qing Hou, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : qhou@scu.edu.cn
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
  module MSM_MultiGPU_Basic
  !***  DESCRIPTION: this module is to define a number of functions for the use
  !                  of multiple GPU. This version is eaxtracted from old version
  !                  of MD_Globle_Variables_GPU.F90. 
  !                  ______________________________________________________
  !
  ! **** HOSTORY:
  !       * Sep.    2018(HOU Qing): Seperated from MD_Globle_Variables_GPU.F90, 
  !                                 Keep some generic functions 
  !
  use CUDAFOR
  use MSM_CONSTANTS
  implicit none

     !--- The list of devices to be used
          integer,parameter::m_MXDEVICE=6                             !The number limitation of devices to be used.
                                                                      !without the major device counted

          integer::m_NDEVICE=0                                        !The actual number of devices to be used
          integer::m_DEVICES(0:m_MXDEVICE)=(/0,0,1,2,3,4,5/)          !The index of devices to be used

      !--- The GPU random number generator
           !--- uniform generator
           integer(kind=int_ptr_kind())::dm_URanGen(m_MXDEVICE) = 0
           !--- gaussian generator
           integer(kind=int_ptr_kind())::dm_GRanGen(m_MXDEVICE) = 0



  !--- interface of routines in this module -------------------
  !---------------------------------------------------------
  !       to check if device have been initialized
          public:: Check_DEVICES

  !---------------------------------------------------------
          public :: Clear_DeviceSwapData

  !---------------------------------------------------------
  !       to exit devices to release all resource on device
          public:: End_DEVICES

  !---------------------------------------------------------
  !       to set the number and index of devices to be used
          public:: Initialize_DEVICES

  !---------------------------------------------------------
  !       to synchronize all devices
          public:: SynchronizeDevices

  !---------------------------------------------------------
  !       some parameter and generic array operations
  !       
          integer, parameter, private::mp_ArrayOp_2Power   = 6
          integer, parameter, private::mp_ArrayOp_Blocksize = 2**mp_ArrayOp_2Power
          integer, parameter, private::mp_ArrayOp_Gridsize  = 1024
          !--- the data type used for storing temp. data ---------------
          private DeviceSwapData
          type::DeviceSwapData
                integer                                          ::IDEV
                real(KINDDF), device, dimension(:),  allocatable::Swap2D
          end   type DeviceSwapData

          type(DeviceSwapData), dimension(:), allocatable, private::dm_Array2DSwap
          private:: Clear_DeviceSwapData
          private:: &
                    Dot_product2d_KERNEL0, &
                    Dot_product2d_KERNEL1
          private:: &
                    Dot_product2d_template0,           &
                    Dot_product2d_template1

          public::  Dot_product_template
          interface Dot_product_template
                    module procedure Dot_product2d_template0
                    module procedure Dot_product2d_template1
          end interface Dot_product_template

          private:: Normal2d_template0
          public::  Norm_template
          interface Norm_template
                    module procedure Normal2d_template0
          end interface Norm_template
          
          private:: &
                    Sep_product2d_KERNEL0, &
                    Sep_product2d_KERNEL1
          private:: &
                    Sep_product2d_template0,   &
                    Sep_product2d_template0_a, &
                    Sep_product2d_template1,   &
                    Sep_product2d_template1_a

          public::  Sep_product_template
          interface Sep_product_template
                    module procedure Sep_product2d_template0
                    module procedure Sep_product2d_template0_a
                    module procedure Sep_product2d_template1
                    module procedure Sep_product2d_template1_a
          end interface Sep_product_template
  contains


  !****************************************************************************
  subroutine Initialize_DEVICES( FIRSTDEV, NDEV)
  !***  PURPOSE:  to set the number and index of devices to be used
  !     INPUT     FIRSTDEV,  the first device ID used by this host process
  !               NDEV,      the number of devices involved in the host process
  !
      implicit none
      !--- dummy variables
      integer::FIRSTDEV,NDEV
      !--- Local vairables
      integer::ERR, K

             if(NDEV .GT. m_MXDEVICE) then
                write(*,*) "MDPSCU Error: the number of devices larger than permitted value ", m_MXDEVICE
                stop
              end if

              m_DEVICES(0) = FIRSTDEV
              m_NDEVICE    = NDEV
              m_DEVICES(1) = m_DEVICES(0)
              do K=2, m_NDEVICE
                 m_DEVICES(K) = m_DEVICES(K-1)+1
              end do

              !--- allocate memery for temp arrays
              if(allocated(dm_Array2DSwap)) then
                 do K=1, size(dm_Array2DSwap)
                    call Clear_DeviceSwapData(dm_Array2DSwap(K))
                 end do
                 deallocate(dm_Array2DSwap)
              end if

              allocate(dm_Array2DSwap(m_NDEVICE))
              do K=1, m_NDEVICE
                 dm_Array2DSwap(K)%IDEV = m_DEVICES(K)
                 ERR = cudaSetDevice(m_DEVICES(K))
                 allocate(dm_Array2DSwap(K)%Swap2D(mp_ArrayOp_Gridsize))
              end do
              ERR = cudaSetDevice(m_DEVICES(0) )

            return
  end subroutine Initialize_DEVICES
  !****************************************************************************

  !****************************************************************************
  subroutine Check_DEVICES( )
  !***  PURPOSE:  to check if device have been initialized
  !     INPUT
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::ERR, K


               if(m_NDEVICE .le. 0) then
                  write(*,fmt="(A)")  ' MDPSCU Error: GPU version of MDPSCU is to be used. But devices are not initialized.'
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               end if
            return
  end subroutine Check_DEVICES
  !****************************************************************************

  !****************************************************************************
  subroutine End_DEVICES( )
  !***  PURPOSE:  to exit devices to release all resource on device
  !     INPUT
  !
      use CudaRandomC2F_M
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::ERR, K

              !--- to destroy a sequence random number generator on devices
              do K=1, m_NDEVICE
                 ERR = cudaSetDevice(m_DEVICES(K) )
                 if(dm_URanGen(K) .gt. 0) ERR = curandDestroyGenerator(dm_URanGen(K))
                 if(dm_GRanGen(K) .gt. 0) ERR = curandDestroyGenerator(dm_GRanGen(K))
              end do
              dm_URanGen = 0
              dm_GRanGen = 0

              do K=0, m_NDEVICE
                 ERR = cudaSetDevice(m_DEVICES(K) )
                 ERR = cudaDeviceReset()
              end do

            return
  end subroutine End_DEVICES
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_Rand_DEVICES()
  !***  PURPOSE:  to intialize the random numbers on devices
  !     INPUT     FIRSTDEV,  the first device ID used by this host process
  !               NDEV,      the number of devices involved in the host process
  !
      use RAND32SEEDLIB_MODULE
      use RAND32_MODULE
      use CudaRandomC2F_M
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::ERR, CURDEV
      integer(kind=8)::SEED
      integer(kind=4)::ISEED(2)
      equivalence(SEED, ISEED(1))
      integer::I
      !--- for rand number test
      !integer::N = 32*10000
      !real(kind=8), device, dimension(:),allocatable :: dRandomArray1
      !real(kind=8), device, dimension(:),allocatable :: dRandomArray2
      !real(kind=8), dimension(:),allocatable :: hRandomArray1
      !real(kind=8), dimension(:),allocatable :: hRandomArray2


              ERR = cudaGetDevice(CURDEV )
              !--- to create a sequence random number generator on devices
              do I=1, m_NDEVICE
                 ERR = cudaSetDevice(m_DEVICES(I) )
                 if(dm_URanGen(I) .gt. 0) ERR = curandDestroyGenerator(dm_URanGen(I))
                 ERR = curandCreateGenerator(dm_URanGen(I),CURAND_RNG_PSEUDO_XORWOW)
                 ERR = DRand32()*RAND32SEEDLIB_SIZE+1
                 call GetSeed_RAND32SEEDLIB(ERR, ISEED(1), ISEED(2))
                 ERR = curandSetPseudoRandomGeneratorSeed(dm_URanGen(I),SEED)

                 !--- to crerate gaussian generator
                 if(dm_GRanGen(I) .gt. 0) ERR = curandDestroyGenerator(dm_GRanGen(I))
                 ERR = curandCreateGenerator(dm_GRanGen(I),CURAND_RNG_QUASI_SCRAMBLED_SOBOL64)
                 ERR = DRand32()*RAND32SEEDLIB_SIZE+1
                 call GetSeed_RAND32SEEDLIB(ERR, ISEED(1), ISEED(2))
                 ERR = curandSetPseudoRandomGeneratorSeed(dm_GRanGen(I),SEED)
              end do

              ERR = cudaSetDevice(CURDEV)
            return
  end subroutine Initialize_Rand_DEVICES
  !****************************************************************************

  !****************************************************************************
  subroutine SynchronizeDevices()
  !***  PURPOSE:   to synchronize all devices
  !
      implicit none
      !--- dummy variables
           integer::CURDEV,I,ERR

           ERR = cudaGetDevice(CURDEV)
           DO I=0, m_NDEVICE
              err = cudaSetDevice(m_DEVICES(I))
              err = cudaThreadSynchronize();
           END DO
           ERR = cudaSetDevice(CURDEV)

  end subroutine SynchronizeDevices
  !****************************************************************************
  
  !*****************************************************************************
  !*********************************************************************
  subroutine Clear_DeviceSwapData(SwapData)
        implicit none
        type(DeviceSwapData)::SwapData
        !--- local varaiables
        integer::CURDEV, ERR

            if( SwapData%IDEV .lt. 0) return

               ERR = cudaGetDevice(CURDEV)
               ERR = cudaSetDevice(SwapData%IDEV)
               SwapData%IDEV  = -1 
               if(allocated(SwapData%Swap2D))  deallocate( SwapData%Swap2D)
               ERR = cudaSetDevice(CURDEV)
       return
  end subroutine Clear_DeviceSwapData
  !*********************************************************************

  !****************************************************************************
  attributes(global) subroutine Dot_product2d_KERNEL0(N1, N2, V1, V2, NPB, BRES)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array: sum(V1(N1:N2,1:3)*V2(N1:N2,1:3))
  !
  !
  implicit none
  !----   DUMMY Variables
          integer,     value                 :: NPB, N1, N2
          real(KINDDF),device, dimension(:,:):: V1, V2
          real(KINDDF),device, dimension(:)  :: BRES
  
  !----   Local variables
          integer        :: IT, IB, I, J, IM, OFFSET    
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = 0.D0
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + N1 -1
                 if(J.le.N2) then
                    S(IT) = S(IT) + V1(J,1)*V2(J,1) + V1(J,2)*V2(J,2) + V1(J,3)*V2(J,3)
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = S(IT) + S(IT+IM)
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
        return
  end subroutine Dot_product2d_KERNEL0
  !****************************************************************************

  !****************************************************************************
  attributes(global) subroutine Dot_product2d_KERNEL1(N1, N2, V1, V2, GID, NPB, BRES)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(N1:N2,:)*V2(GID(N1:N2),:)
  !
  !
  implicit none
  !----   DUMMY Variables
          integer,     value                 :: NPB, N1, N2
          integer,     device, dimension(:)  :: GID
          real(KINDDF),device, dimension(:,:):: V1, V2
          real(KINDDF),device, dimension(:)  :: BRES
  
  !----   Local variables
          integer             :: IT, IB, I, J, IM, OFFSET    
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = 0.D0
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + N1- 1
                 if(J.le.N2) then
                    S(IT) = S(IT) + V1(J,1)*V2(GID(J),1) + V1(J,2)*V2(GID(J),2) + V1(J,3)*V2(GID(J),3)
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = S(IT) + S(IT+IM)
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
        return
  end subroutine Dot_product2d_KERNEL1
  !**************************************************************************

  !**************************************************************************
  subroutine Dot_product2d_template0(IDEV, N1, N2, V1, V2, PRDCT)
  !***  PURPOSE:              to calculate the dot between two array
  !
  !
  implicit none
      !----   DUMMY Variables
       integer,      intent(in)            :: IDEV, N1, N2
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: PRDCT
       !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I
         real(KINDDF)::hRES(mp_ArrayOp_Gridsize)

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (N2 - N1 + 1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, size(dm_Array2DSwap)
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call Dot_product2d_KERNEL0<<<blocks, threads>>>(N1, N2, V1, V2, NPB, dm_Array2DSwap(I)%Swap2D)
                  hRES = dm_Array2DSwap(I)%Swap2D
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             PRDCT = sum(hRES)

             return
   end subroutine Dot_product2d_template0
  !***************************************************************************

  !**************************************************************************
  subroutine Dot_product2d_template1(IDEV, N1, N2, V1, V2, GID, PRDCT)
  !***  PURPOSE:             to calculate the dot between two array
  !
  !
  implicit none
      !----   DUMMY Variables
       integer,      intent(in)            :: IDEV, N1, N2
       real(KINDDF), device, dimension(:,:):: V1, V2
       integer,      device, dimension(:)  :: GID
       real(KINDDF)                        :: PRDCT
       !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I 
         real(KINDDF)::hRES(mp_ArrayOp_Gridsize)

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (N2 - N1 + 1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, size(dm_Array2DSwap)
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call Dot_product2d_KERNEL1<<<blocks, threads>>>(N1, N2, V1, V2, GID, NPB, dm_Array2DSwap(I)%Swap2D)
                  hRES = dm_Array2DSwap(I)%Swap2D
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             PRDCT = sum(hRES)

             return
   end subroutine Dot_product2d_template1
  !***************************************************************************

  !**************************************************************************
  subroutine Normal2d_template0(IDEV, N1, N2, V, NORM)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array
  !
  !
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV, N1, N2
       real(KINDDF), device, dimension(:,:):: V
       real(KINDDF)                        :: NORM
 
                 call Dot_product2d_template0(IDEV, N1, N2, V, V, NORM)
                 return
  end subroutine Normal2d_template0
  !**************************************************************************

  !****************************************************************************
  attributes(global) subroutine Sep_product2d_KERNEL0(N1, N2, V1, V2 , BSX, BSY, BSZ, NPB, BRES)
  !***  PURPOSE:   KERNEL    to calculate the dot between  difference of two array
  !
  !
  implicit none
  !----   DUMMY Variables
          integer,      value                 :: NPB, N1, N2
          real(KINDDF), value                 :: BSX, BSY, BSZ
          real(KINDDF), device, dimension(:,:):: V1, V2
          real(KINDDF), device, dimension(:)  :: BRES
  
  !----   Local variables
          integer             :: IT, IB, I, J, J1, J2, IM, OFFSET    
          real(KINDDF)        :: SEPX, SEPY, SEPZ, HBSX, HBSY, HBSZ
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB

              HBSX = 0.5D0*BSX
              HBSY = 0.5D0*BSY
              HBSZ = 0.5D0*BSZ
              
              S(IT) = 0.D0
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + N1 -1
                 if(J.le.N2) then
                    SEPX  = V1(J,1) - V2(J,1)
                    SEPY  = V1(J,1) - V2(J,1)
                    SEPZ  = V1(J,1) - V2(J,1)
                    if(dabs(SEPX) .gt. HBSX) SEPX = SEPX - dsign(BSX, SEPX)
                    if(dabs(SEPY) .gt. HBSY) SEPY = SEPY - dsign(BSY, SEPY)
                    if(dabs(SEPZ) .gt. HBSZ) SEPZ = SEPZ - dsign(BSZ, SEPY)
                    S(IT) = S(IT) + SEPX*SEPX + SEPY*SEPZ + SEPY*SEPZ
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = S(IT) + S(IT+IM)
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
        return
  end subroutine Sep_product2d_KERNEL0
  !****************************************************************************

  !****************************************************************************
  attributes(global) subroutine Sep_product2d_KERNEL1(N1, N2, V1, V2, GID, BSX, BSY, BSZ, NPB, BRES)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array
  !
  !
  !
  implicit none
  !----   DUMMY Variables
          integer,      value                 :: NPB, N1, N2
          real(KINDDF), value                 :: BSX, BSY, BSZ
          real(KINDDF), device, dimension(:,:):: V1, V2
          integer,      device, dimension(:)  :: GID
          real(KINDDF), device, dimension(:)  :: BRES

  
  !----   Local variables
          integer             :: IT, IB, I, J, IM, OFFSET    
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = 0.D0
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + N1 - 1
                 if(J.le.N2) then
                    S(IT) = S(IT) + (V1(J,1)-V2(GID(J),1))*(V1(J,1)-V2(GID(J),1)) &
                                  + (V1(J,2)-V2(GID(J),2))*(V1(J,2)-V2(GID(J),2)) &
                                  + (V1(J,3)-V2(GID(J),3))*(V1(J,3)-V2(GID(J),3))
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = S(IT) + S(IT+IM)
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
        return
  end subroutine Sep_product2d_KERNEL1
  !**************************************************************************

  !**************************************************************************
  subroutine Sep_product2d_template0(IDEV, N1, N2, V1, V2, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:   KERNEL    to calculate the dot between subset of two array:
  !                          sum(V1(N1:N2)*V2(N1,N2)  
  !
  !     INPUT:   IDEV,       the ID of the device  
  !              N1, N2,    the bound of the subset
  !              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !     OUTPUT   PRDCT,     the results
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV, N1, N2
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       real(KINDDF)                        :: PRDCT
       !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: N, NPB, CURDEV, ERR, I 
         real(KINDDF)::hRES(mp_ArrayOp_Gridsize)

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (N2-N1+1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, size(dm_Array2DSwap)
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call Sep_product2d_KERNEL0<<<blocks, threads>>>(N1, N2, V1, V2, BSX, BSY, BSZ, NPB, dm_Array2DSwap(I)%Swap2D)
                  hRES = dm_Array2DSwap(I)%Swap2D
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             PRDCT = sum(hRES)

             return
   end subroutine Sep_product2d_template0
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_product2d_template1(IDEV, N1, N2, V1, V2, GID, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array
  !
  !
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV, N1, N2
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       integer,      device, dimension(:)  :: GID
       real(KINDDF)                        :: PRDCT
       !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I
         real(KINDDF)::hRES(mp_ArrayOp_Gridsize)

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (N2 - N1 +1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, size(dm_Array2DSwap)
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call Sep_product2d_KERNEL1<<<blocks, threads>>>(N1, N2, V1, V2, GID, BSX, BSY, BSZ, NPB, dm_Array2DSwap(I)%Swap2D)
                  hRES = dm_Array2DSwap(I)%Swap2D
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             PRDCT = sum(hRES)

             return
   end subroutine Sep_product2d_template1
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_product2d_template0_a(IDEV, V1, V2, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array
  !
  !
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       real(KINDDF)                        :: PRDCT
       !--- Device variables and variables to be used in GPU

             call Sep_product2d_template0(IDEV, 1, size(V1,dim=1), V1, V2, BSX, BSY, BSZ, PRDCT)
             return
   end subroutine Sep_product2d_template0_a
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_product2d_template1_a(IDEV, V1, V2, GID, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array
  !
  !
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       integer,      device, dimension(:)  :: GID
       real(KINDDF)                        :: PRDCT
       !--- Device variables and variables to be used in GPU

             call Sep_product2d_template1(IDEV, 1, size(V1,dim=1), V1, V2, GID, BSX, BSY, BSZ, PRDCT)
             return
   end subroutine Sep_product2d_template1_a
  !***************************************************************************
  !**************************************************************************
  end module MSM_MultiGPU_Basic

