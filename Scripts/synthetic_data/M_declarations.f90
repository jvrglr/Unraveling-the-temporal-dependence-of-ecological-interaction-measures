module declarations_module
    !Module for public variables to be used in the rest of the modules
    implicit none
    public
    
    integer*4 :: ios
    integer*8 :: n_R,n_C
    double precision :: cr
    double precision :: t
    double precision, dimension(:),allocatable :: R,C
    double precision, dimension(:),allocatable :: e,d,K,nu,Lims,flow,D_R,D_C
    double precision, dimension(:,:),allocatable :: a,e_cross
    double precision, parameter :: pi = &
     3.1415926535897932384626433832795028841972 !4.0d0*atan(1.0d0)
  
  
  contains
  
    subroutine assingments()
      implicit none
      n_R = 2 !Number of resources
      n_C = 2 !Number of consumners
      allocate(R(n_R),nu(n_R),Lims(n_R),flow(n_R))
      allocate(e(n_C),d(n_C),C(n_C),K(n_R),D_R(n_R),D_C(n_C)) 
      allocate(a(n_C,n_R),e_cross(n_C,n_R))
      R = 0.0d0 !Consumers
      C = 0.0d0 !Resources
      cr = 3.0d0 !Intrinsic renewal rate
      !---------attack-consumption rate-------
      a(1,1) = 1.0d0 
      a(1,2) = 0.5d0 
      a(2,1) = 0.0d0 
      a(2,2) = 1.0d0
      !---------------------------------------
      e_cross = 0.0d0 !Efficiency of cross-feeding
      e = 0.08d0 
      nu = 1.0d0
      K(1) = 20.0d0 !Carrying capacity of resource
      K(2) = 10.0d0 !Carrying capacity of resource
      d = 0.1d0 !death rate
      flow = 100.0d0 !Chemostat flow rate

      t = 0.0d0 !time
    end subroutine assingments
  
  end module declarations_module