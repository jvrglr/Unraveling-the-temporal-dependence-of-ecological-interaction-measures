module Subroutines

    !Module with public subroutines.
    !Depends on declarations_module and functions.
  
    use declarations_module
    use functions
  
  contains
    subroutine Compute_interactions_over_time_1R1C(v_ref,v_per,dt,steps,M)
      implicit none
      integer*8, intent(in) :: steps
      double precision, intent(in) :: dt
      double precision, intent(in), dimension (0:steps,*) :: v_ref,v_per
      double precision, intent(inout), dimension (steps) :: M
      
      M = (v_per(1:,1)-v_per(:steps-1,1))/v_per(:steps-1,1)/dt - (v_ref(1:,1)-v_ref(:steps-1,1))/v_ref(:steps-1,1)/dt
    end subroutine Compute_interactions_over_time_1R1C

    subroutine Compute_interactions_over_time_general_CR(v_ref,v_per,dt,steps,M,i_target)
      !Compute interactions over time for a general CR model
      !v_ref and v_per are the reference and perturbed trajectories, respectively. They are matrices with dimensions (0:steps,n_C+n_R) where the first n_C columns correspond to consumers and the last n_R columns correspond to resources. 
      !i_target is the index of the target species (consumer or resource) for which interactions are computed.
      implicit none
      integer*8, intent(in) :: steps,i_target
      double precision, intent(in) :: dt
      double precision, intent(in), dimension (0:steps,*) :: v_ref,v_per
      double precision, intent(inout), dimension (steps) :: M
      
      M = (v_per(1:,i_target)-v_per(:steps-1,i_target))/v_per(:steps-1,i_target)/dt - (v_ref(1:,i_target)-v_ref(:steps-1,i_target))/v_ref(:steps-1,i_target)/dt
    end subroutine Compute_interactions_over_time_general_CR

    subroutine Compute_g_empirical(C0,R0,tf,dt_traj,g_C,g_R,dt_estimation_g,mode_consumption,mode_renewal,mode_generation)
      !Estimate growth rate per capita at time tf in 1C 1R model
      !dt for estimation of g is assumed to be equal to the dt for the trajectory integration unless explicirly stated.
      implicit none
      double precision, intent(in) :: C0,R0,dt_traj
      double precision, intent(inout) :: tf
      double precision, optional, intent(in) :: dt_estimation_g
      double precision, intent(out) :: g_C,g_R
      integer*8, intent(in) :: mode_consumption,mode_renewal,mode_generation
      double precision :: dt_g,C_t,R_t,C_tpdt,R_tpdt,dum_t

      
      if ( present(dt_estimation_g) ) then
          dt_g = dt_estimation_g
      else
          dt_g = dt_traj
      end if
      dum_t = tf+dt_g

      R(1) = R0 
      C(1) = C0
      t = 0.0d0
      call Propagate_IC(tf ,mode_consumption,mode_renewal,mode_generation,temp_dt=dt_traj)

      C_t = C(1)
      R_t = R(1)
      call Propagate_IC_one_step(dt_g,mode_consumption,mode_renewal,mode_generation)
      C_tpdt = C(1)
      R_tpdt = R(1)
      g_C = (C_tpdt-C_t)/C_t/dt_g
      g_R = (R_tpdt-R_t)/R_t/dt_g
    ! print *, "g_C=",g_C,e*a(1,1)*R_t-d
    ! print *, "g_R=",g_R,-a(1,1)*C_t+ cr*(1.0-R_t/K)

    end subroutine Compute_g_empirical

      subroutine Compute_g_empirical_general_consumer(tf,dt_traj,g_C,jj,dt_estimation_g, mode_consumption,mode_renewal,mode_generation)
      !Estimate growth rate per capita at time tf in CR model over
      !dt for estimation of g is assumed to be equal to the dt for the trajectory integration unless explicirly stated.
      implicit none
      double precision, intent(in) :: dt_traj
      double precision, intent(inout) :: tf
      double precision, optional, intent(in) :: dt_estimation_g
      integer*8, intent(in) :: mode_consumption,mode_renewal,mode_generation
      double precision, intent(out) :: g_C
      integer*8, intent(in) :: jj
      double precision :: dt_g,C_t,C_tpdt,dumt
      
      if ( present(dt_estimation_g) ) then
          dt_g = dt_estimation_g
      else
          dt_g = dt_traj
      end if
      
      t = 0.0d0
      call Propagate_IC(tf ,mode_consumption,mode_renewal,mode_generation,temp_dt=dt_traj)

      C_t = C(jj)
      if ( dt_g < 2*dt_traj ) then !Integration = one iteration of RK4 method
        call Propagate_IC_one_step(dt_g,mode_consumption,mode_renewal,mode_generation)
      else 
        dumt = tf + dt_g
        call Propagate_IC(dumt,mode_consumption,mode_renewal,mode_generation,temp_dt=dt_traj)
      end if
      
      C_tpdt = C(jj)

      g_C = (C_tpdt-C_t)/C_t/dt_g


    end subroutine Compute_g_empirical_general_consumer


    subroutine Propagate_IC(tf,mode_consumption,mode_renewal,mode_generation,temp_dt)
        !Evolve initial condition (C, R) from t to tf.
        !C,R are overwritten with state of the system at time tf.
        implicit none
        double precision,intent(inout) :: tf
        integer*8, intent(in) :: mode_consumption,mode_renewal,mode_generation
        double precision,optional,intent(in) :: temp_dt
        double precision, dimension (n_C) :: temp_C,k1_C,k2_C,k3_C,k4_C
        double precision, dimension (n_R) :: temp_R,k1_R,k2_R,k3_R,k4_R
        integer*8 :: steps,ii
        double precision :: dt

        if ( present(temp_dt) ) then
            dt = temp_dt
        else
            dt = 0.001d0
        end if
        steps = int8((tf-t)/dt)
        if ( steps<1 ) then
            print *, "tf<dt so tf is redefined as t+dt"
            print *
            steps = 1
            tf = t + dt
        end if
        do ii = 1, steps
            t=t+dt
            ! Euler method-------------------
            ! R = R + dt*der_R(C,R,mode_consumption,mode_renewal)
            ! C = C + dt*der_C(C,R,mode_generation)
            ! -------------------------------
            ! RK4 method---------------------
            k1_R = der_R(C,R,mode_consumption,mode_renewal)
            k1_C = der_C(C,R,mode_generation)
            k2_R = der_R(C+dt*k1_C/2.0d0,R+dt*k1_R/2.0d0,mode_consumption,mode_renewal)
            k2_C = der_C(C+dt*k1_C/2.0d0,R+dt*k1_R/2.0d0,mode_generation)
            k3_R = der_R(C+dt*k2_C/2.0d0,R+dt*k2_R/2.0d0,mode_consumption,mode_renewal)
            k3_C = der_C(C+dt*k2_C/2.0d0,R+dt*k2_R/2.0d0,mode_generation)
            k4_R = der_R(C+dt*k3_C,R+dt*k3_R,mode_consumption,mode_renewal)
            k4_C = der_C(C+dt*k3_C,R+dt*k3_R,mode_generation)
            R = R+dt*(k1_R+2.0d0*k2_R+2.0d0*k3_R+k4_R)/6.0d0
            C = C+dt*(k1_C+2.0d0*k2_C+2.0d0*k3_C+k4_C)/6.0d0
        end do
        
    end subroutine Propagate_IC

    subroutine Propagate_IC_stochastic(tf,mode_consumption,mode_renewal,mode_generation,temp_dt)
        !Evolve initial condition (C, R) from t to tf.
        !C,R are overwritten with state of the system at time tf.
        implicit none
        double precision,intent(inout) :: tf
        integer*8, intent(in) :: mode_consumption,mode_renewal,mode_generation
        double precision,optional,intent(in) :: temp_dt
        double precision, dimension (n_C) :: temp_C,k1_C,k2_C,k3_C,k4_C
        double precision, dimension (n_R) :: temp_R,k1_R,k2_R,k3_R,k4_R
        integer*8 :: steps,ii
        double precision :: dt,sqdt,dran_g,u

        if ( present(temp_dt) ) then
            dt = temp_dt
        else
            dt = 0.001d0
        end if
        sqdt = sqrt(dt)
        steps = int8((tf-t)/dt)
        if ( steps<1 ) then
            print *, "tf<dt so tf is redefined as t+dt"
            print *
            steps = 1
            tf = t + dt
        end if
        do ii = 1, steps
            t=t+dt
            ! Euler method-------------------
            u = dran_g()
            R = R + dt*der_R(C,R,mode_consumption,mode_renewal) + sqdt*D_R*u
            C = C + dt*der_C(C,R,mode_generation) + sqdt*D_C*dran_g()
        end do
        
    end subroutine Propagate_IC_stochastic

    
    subroutine Propagate_IC_one_step(temp_dt,mode_consumption,mode_renewal,mode_generation)
        !Evolve initial condition (C, R) from t to t+dt.
        !C,R are overwritten with state of the system at time tf.
        implicit none
        double precision,optional,intent(in) :: temp_dt
        integer*8, intent(in) :: mode_consumption,mode_renewal,mode_generation
        double precision, dimension (n_C) :: temp_C,k1_C,k2_C,k3_C,k4_C
        double precision, dimension (n_R) :: temp_R,k1_R,k2_R,k3_R,k4_R
        double precision :: dt

        if ( present(temp_dt) ) then
            dt = temp_dt
        else
            dt = 0.001d0
        end if

        t=t+dt
        ! Euler method-------------------
        ! temp_R = R + dt*der_R(C,R)
        ! temp_C = C + dt*der_C(C,R)
        ! -------------------------------
        ! RK4 method---------------------
        k1_R = der_R(C,R,mode_consumption,mode_renewal)
        k1_C = der_C(C,R,mode_generation)
        k2_R = der_R(C+dt*k1_C/2.0d0,R+dt*k1_R/2.0d0,mode_consumption,mode_renewal)
        k2_C = der_C(C+dt*k1_C/2.0d0,R+dt*k1_R/2.0d0,mode_generation)
        k3_R = der_R(C+dt*k2_C/2.0d0,R+dt*k2_R/2.0d0,mode_consumption,mode_renewal)
        k3_C = der_C(C+dt*k2_C/2.0d0,R+dt*k2_R/2.0d0,mode_generation)
        k4_R = der_R(C+dt*k3_C,R+dt*k3_R,mode_consumption,mode_renewal)
        k4_C = der_C(C+dt*k3_C,R+dt*k3_R,mode_generation)
        R = R+dt*(k1_R+2.0d0*k2_R+2.0d0*k3_R+k4_R)/6.0d0
        C = C+dt*(k1_C+2.0d0*k2_C+2.0d0*k3_C+k4_C)/6.0d0

    end subroutine Propagate_IC_one_step

    subroutine Draw_trajectory(from_t,until_t,Delta_t,filename,temp_dt,mode_consumption,mode_renewal,mode_generation)
        !Concatenate runs of "Propagate_IC" to generate a trajectory, then saves the trajectory into "filename"
        !temp_dt is an optional variable that sets discretization for "Propagate_IC"
        !temp_addname is an optional variable that is used to modi
        implicit none
        double precision, intent(in) :: from_t,until_t,Delta_t
        double precision,optional, intent(in) :: temp_dt
        character(len=*),intent(in) :: filename
        integer*8, intent(in) :: mode_consumption,mode_renewal,mode_generation
        double precision :: tf
        integer*8 :: steps,ii
        double precision, dimension (:,:), allocatable :: temp_Cs,temp_Rs
        double precision, dimension (:), allocatable :: ts
        
        steps = int8((until_t-from_t)/Delta_t)
        if ( steps<1 ) then
            print *, "until_t-from_t<Delta_t so no evolution is performed"
            print *
        end if
        
        allocate(temp_Cs(0:steps,n_C),temp_Rs(0:steps,n_R),ts(0:steps))

        t= from_t
        temp_Cs(0,:) = C
        temp_Rs(0,:) = R
        ts(0) = t
        do ii = 1, steps
            tf = from_t+ii*(until_t-from_t)/steps
            if ( present(temp_dt) ) then
                call Propagate_IC(tf,mode_consumption,mode_renewal,mode_generation,temp_dt=temp_dt)
            else 
                call Propagate_IC(tf,mode_consumption,mode_renewal,mode_generation)
            end if
            t = tf
            temp_Cs(ii,:) = C
            temp_Rs(ii,:) = R
            ts(ii) = t
        end do

        open(unit=1001, file=filename, iostat=ios, status="unknown", action="write")
        if ( ios /= 0 ) stop "Error opening file "//filename
        do ii = 0, steps, 1
            write(1001, *) ts(ii),temp_Cs(ii,:),temp_Rs(ii,:)
        end do

        close(1001)

        deallocate(ts,temp_Cs,temp_Rs)
    end subroutine Draw_trajectory

    subroutine Draw_trajectory_in_array(temp_Rs,temp_Cs,ts,steps,from_t,until_t,Delta_t,temp_dt,mode_consumption,mode_renewal,mode_generation)
        !Concatenate runs of "Propagate_IC" to generate a trajectory in CR model, then saves the trajectory into array
        !temp_dt is an optional variable that sets discretization for "Propagate_IC"
        implicit none
        double precision, intent(in) :: from_t,until_t,Delta_t
        double precision,optional, intent(in) :: temp_dt
        integer*8, intent(in) :: steps,mode_consumption,mode_renewal,mode_generation
        double precision :: tf
        integer*8 :: ii
        double precision, intent(out), dimension (:,:), allocatable :: temp_Cs,temp_Rs
        double precision, intent(out), dimension (:), allocatable :: ts

        allocate(temp_Cs(0:steps,n_C),temp_Rs(0:steps,n_R),ts(0:steps))

        t= from_t
        temp_Cs(0,:) = C
        temp_Rs(0,:) = R
        ts(0) = t
        do ii = 1, steps
            tf = from_t+ii*(until_t-from_t)/steps
            if ( present(temp_dt) ) then
                call Propagate_IC(tf,mode_consumption,mode_renewal,mode_generation,temp_dt=temp_dt)
            else 
                call Propagate_IC(tf,mode_consumption,mode_renewal,mode_generation)
            end if
            t = tf
            temp_Cs(ii,:) = C
            temp_Rs(ii,:) = R
            ts(ii) = t
        end do

    end subroutine Draw_trajectory_in_array

    subroutine Draw_trajectory_in_array_stochastic(temp_Rs,temp_Cs,ts,steps,from_t,until_t,Delta_t,temp_dt,mode_consumption,mode_renewal,mode_generation)
        !Concatenate runs of "Propagate_IC" to generate a trajectory in CR model, then saves the trajectory into array
        !temp_dt is an optional variable that sets discretization for "Propagate_IC"
        implicit none
        double precision, intent(in) :: from_t,until_t,Delta_t
        double precision,optional, intent(in) :: temp_dt
        integer*8, intent(in) :: steps,mode_consumption,mode_renewal,mode_generation
        double precision :: tf
        integer*8 :: ii
        double precision, intent(out), dimension (:,:), allocatable :: temp_Cs,temp_Rs
        double precision, intent(out), dimension (:), allocatable :: ts

        allocate(temp_Cs(0:steps,n_C),temp_Rs(0:steps,n_R),ts(0:steps))

        t= from_t
        temp_Cs(0,:) = C
        temp_Rs(0,:) = R
        ts(0) = t
        do ii = 1, steps
            tf = from_t+ii*(until_t-from_t)/steps
            if ( present(temp_dt) ) then
                call Propagate_IC_stochastic(tf,mode_consumption,mode_renewal,mode_generation,temp_dt=temp_dt)
            else 
                call Propagate_IC_stochastic(tf,mode_consumption,mode_renewal,mode_generation)
            end if
            t = tf
            temp_Cs(ii,:) = C
            temp_Rs(ii,:) = R
            ts(ii) = t
        end do
        
    end subroutine Draw_trajectory_in_array_stochastic
    
    subroutine read_xF_tf()
      !Example of subroutine, read data file
      implicit none
      integer*4 :: ios,i,dum,datapoints
      double precision :: dummy
  
      open(unit=1001, file="data/F_Target_test3.dat", iostat=ios, status="old", action="read")
      if ( ios /= 0 ) stop "Error opening file "
      do i = 1, datapoints, 1
        read(1001,*) dummy,dummy,dummy,dum
      end do
  
  
    end subroutine read_xF_tf
    
    subroutine search_list_binary_algorithm(list,position,p)
      !Rafle event:
      !Given a list of probabilities called "list" such that sum(list)=1 and a probability p.
      !Look for "position" such that C(position)>=p and C(j)<p for all j in [1,position[.
      !Where C is the cumulative of list: C(i)=list(1)+list(2)+...+list(i)-
      !REFERENCE:Brainerd, W. S. (2015). Guide to Fortran 2008 programming (p. 141). Berlin: Springer.
  
      implicit none
      double precision, dimension(:), intent (in) :: list
      double precision, intent (in) :: p
      integer*4, intent(out) :: position
      double precision, dimension(size(list)) :: C
      integer*4 ii,N,first,last,half
  
      N=size(list) !It would be cool to define this as a parameter (constant), I don't know how...
      C(1)=list(1)
      do ii = 2, N, 1 !Compute cumulative of list
        c(ii)=C(ii-1)+list(ii)
      end do
  
      first=1;last=N
      do while ( first.ne.last )
        half=(first+last)/2
        if ( p>C(half) ) then
          first=half+1
        else
          last=half
        end if
      end do
  
      position=first
  
    end subroutine search_list_binary_algorithm
  
    subroutine float_to_string(a, n, result)
      !Convert float number a to string with format int(a)//"d"//dec(a).
      !E.g. 0.00543--> "0d00543"
      !For dec(a) select first n decimals
      !If number of digits in dec(a)>n then add zeros to the left
      double precision, intent(in) :: a
      integer, intent(in) :: n
      character(len=*), intent(out) :: result
      character(len=n) :: dumc
      integer :: digit,zeros
  
    
      integer :: int_part
      double precision :: dec_part
    
      ! Get the integer and decimal parts
      int_part = int(a)
      dec_part = a - int_part
  
      ! Convert the decimal part to a string with n digits
      zeros=0
      digit=int(dec_part*10**(zeros+1))
      do while ((digit==0).and.(zeros<n))
        zeros=zeros+1
        digit=int(dec_part*10**(zeros+1))
      enddo
      dumc=repeat("0",zeros)//trim(str(int(dec_part*10**(n)))) 
  
  
      result=trim(str(int_part))//"d"//trim(dumc)
  
    
    end subroutine float_to_string
  
  end module subroutines