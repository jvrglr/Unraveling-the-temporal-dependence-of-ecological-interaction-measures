module functions
    !Module with public functions.
    !Depends on declarations_module.
    use declarations_module
    implicit none

  contains
  
    function f(Rs,Cs,mode_renewal) result(output)
      !Renewal function (production of resources)
      implicit none
      double precision, dimension (n_R) :: output
      double precision, dimension (n_R) :: Rs
      double precision, dimension (n_C) :: Cs,consumption
      integer*8, intent(in) :: mode_renewal 
      integer*8 :: ii,ll
      double precision :: dummy
      
      select case (mode_renewal)
        case (1)! Chemostat: constant influx of resources 
          output = flow
        case (2)! Batch: closed system with no renewal (resource decay)
          output = Rs*0.0d0
        case (3)! Biotic resources
          output = cr*Rs*(1.0d0 - Rs/K)
        case (4)! Cross-feeding
          
          do ii = 1, n_c, 1 !Think on how to do this vectorial, see how I access a
            dummy = 0.0d0
              do ll = 1, n_R, 1
                dummy = dummy+g(Cs,Rs,ii,ll,mode_consumption=int8(2))
              end do
            consumption(ii) = dummy
          end do

          do ll = 1, n_R, 1 !From what was consumed, how much new resources are produced
            dummy = 0.0d0
            do ii = 1, n_C, 1
              dummy = dummy + e_cross(ii,ll)*consumption(ii)
            enddo
            output(ll) = dummy
          end do
        case default
          print *, "Error: Invalid mode_renewal in f function"
          stop
      end select
    end function f

    function g(Cs,Rs,ii,ll,mode_consumption) result(output)
    !Consumption function of resource ll from species ii
        implicit none
        integer*8, intent(in) :: mode_consumption 
        double precision :: output
        double precision, dimension (n_R) :: Rs
        double precision, dimension (n_C) :: Cs
        integer*8 :: ii,ll

        select case (mode_consumption)
          case (1)! Linear consumption
              output = a(ii,ll)*Cs(ii)*Rs(ll)
          case (2)! Monod consumption
              output = a(ii,ll)*Cs(ii)*Rs(ll)/(Lims(ll)+Rs(ll))
          case default
              print *, "Error: Invalid mode_consumption in g function"
              stop
        end select
    end function g

    function gen(Cs,Rs,ii,ll,mode_generation) result(output)
      !Consumer productions
      implicit none
      integer*8, intent(in), optional :: mode_generation 
      double precision :: output
      double precision, dimension (n_R) :: Rs
      double precision, dimension (n_C) :: Cs
      integer*8 :: ii,ll

      select case (mode_generation)
        case (1)! Efficiency from consumer
          output = g(Cs,Rs,ii,ll,mode_generation)*e(ii)
        case (2)! Efficiency from resource (Picciani Mori et al. 2020)
          output = g(Cs,Rs,ii,ll,mode_generation)*nu(ll)
        case default
          print *, "Error: Invalid mode_generation in gen function"
          stop
      end select

    end function gen

    function h(Cs) result(output)
      implicit none
      !Death function
      double precision, dimension (n_C) :: output
      double precision, dimension (n_C) :: Cs
        output = d*Cs
    end function h

    function der_R(Cs,Rs,mode_consumption,mode_renewal) result(output)
        implicit none
        !Time derivative of Resources
        double precision, dimension (n_R) :: output,consumption
        double precision, dimension (n_R) :: Rs
        double precision, dimension (n_C) :: Cs
        integer*8, intent(in) :: mode_consumption, mode_renewal  ! Optional mode parameter
        double precision :: sum_over_Cs
        integer*8 :: jj,ll,mode

        do ll = 1, n_R, 1 !Think on how to do this vectorial, see how I access a
          sum_over_Cs = 0.0d0
          do jj = 1, n_C, 1
            sum_over_Cs = sum_over_Cs+g(Cs,Rs,jj,ll,mode_consumption)
          end do
          consumption(ll)=sum_over_Cs
        end do

        output = f(Rs,Cs,mode_renewal)-consumption
    end function der_R

    function der_C(Cs,Rs,mode) result(output)
      implicit none
      !Time derivative of Consumers
      double precision, dimension (n_C) :: output,consumption
      double precision, dimension (n_R) :: Rs
      double precision, dimension (n_C) :: Cs
      integer*8, intent(in) :: mode  ! Optional mode parameter
      double precision :: sum_over_Rs
      integer*8 :: ii,ll

      do ii = 1, n_c, 1 !Think on how to do this vectorial, see how I access a
        sum_over_Rs = 0.0d0
          do ll = 1, n_R, 1
            sum_over_Rs = sum_over_Rs+gen(Cs,Rs,ii,ll,mode)
          end do
        consumption(ii) = sum_over_Rs
      end do
      output = consumption - h(Cs)
    end function der_C

    character(len=30) function str(k) !fUNCTION
      implicit none
  !   "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k !write to a string
      str = adjustl(str)
    end function str

    character(len=30) function str8(k) !fUNCTION
      implicit none
  !   "Convert an integer to string."
      integer*8, intent(in) :: k
      write (str8, *) k !write to a string
      str8 = adjustl(str8)
    end function str8
  
    function dble_mean(x) result(a)
      !compute average of double precision vector
      implicit none
      double precision :: a
      double precision, dimension (:), intent(in) :: x
      a=sum(x)/size(x)
    end function dble_mean
  
    function dble_var(x) result(a)
      !compute variance of double precision vector
      implicit none
      double precision :: a
      double precision, dimension (:), intent(in) :: x
      double precision, dimension (:), allocatable :: dummy
      integer*4 :: i, len
  
      len=size(x)
      allocate(dummy(len))
  
      dummy=x-dble_mean(x)
      a=sum(dummy*dummy)/len
      deallocate(dummy)
    end function dble_var
  
    function dble_err(x) result(a)
      !compute standard error of double precision vector
      implicit none
      double precision :: a
      double precision, dimension (:), intent(in) :: x
      a=sqrt(dble_var(x)/size(x))
    end function dble_err
    
    function Normal_dist(x,mu,sig) result(a)
      !Gaussian function
      implicit none
      double precision :: a
      double precision, intent (in):: x,mu,sig
      a=(1.0/(sqrt(2.0*pi)*sig))*exp(-(x - mu)*(x - mu) / (2.0d0 *sig*sig))
    end function Normal_dist
  
  end module functions