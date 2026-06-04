program main_program
    !Main program, depends on the rest of modules (declarations_module, functions and subroutines))
    use declarations_module
    use functions
    use subroutines
    implicit none
    double precision :: t0,tf,dt_samples,dt_integration,delta
    integer*8 :: steps,ii,points_e,jj,kk,points_r,count,i_target,j_perturbed,seed, mode,mode_consumption,mode_renewal,mode_generation
    double precision :: d_min,d_max,nu_min,nu_max,a_min,a_max,lims_min,lims_max,C0_min,C0_max,R0_min,R0_max,overall_e_cross,ZZ,cr_max,cr_min,flow_min,flow_max
    double precision :: dran_u
    double precision,dimension(:,:),allocatable :: Cs,Rs,Cs_dr,Rs_dr,Cs_dc,Rs_dc
    double precision,dimension(:),allocatable :: M_CC,M_CR,M_RR,M_RC
    double precision,dimension(:),allocatable :: C_I,R_I,ts
    double precision :: dt
    character (len=100) :: filename,path,name,addname

    seed = 1994 !time()
    call dran_ini(seed)
    call assingments()
    allocate(C_I(n_C),R_I(n_R))
    delta = 1.0D-2
    
    i_target = 1
    j_perturbed = 2
    mode_consumption = int8(1) !linear
    mode_generation = int8(1) !efficiency from consumer

    t0 = 0.0d0
    tf = 50.0d0
    dt_samples = 1.0D-3
    dt_integration = 1.0D-4
    do ii = 1, N_C, 1
      do jj = 1, N_C, 1
        a(ii,jj) = dran_u()*1
      end do
    end do
    cr = 1.0d0 !1.0 !Intrinsic renewal rate
    e_cross = 0.0d0 ! Efficiency of cross-feeding
    e = 0.1d0 !1.0d0 !0.12 !0.05
    nu = 1.0d0
    K = 5.0d0  !20
    d = 0.1d0
    flow = 1d0
    do ii = 1, N_C, 1
      C_I(ii) =  dran_u()*10d-3
    end do
    do ii = 1, N_R, 1
      R_I(ii) =  dran_u()*7.0d0 !1.0d0
    end do

    print*, "dt_samples=",dt_samples,"dt_integration=",dt_integration
    steps = int8((tf-t0)/dt_samples)
    allocate(Cs(0:steps,n_C),Cs_dc(0:steps,n_C),ts(0:steps)) 
    allocate(M_CC(steps))
    if ( steps<1 ) then
        print *, "tf-t0<dt_samples so no evolution is performed"
        print *
    end if

    do mode_renewal = 1, 3, 1
      R = R_I 
      C = C_I
      call Draw_trajectory_in_array(Rs,Cs,ts,steps = steps,from_t = t0,until_t = tf ,Delta_t = dt_samples,temp_dt = dt_integration,mode_consumption = mode_consumption ,mode_renewal = mode_renewal,mode_generation = mode_generation)
      R = R_I
      C = C_I
      C(j_perturbed) = C(j_perturbed)+delta
      call Draw_trajectory_in_array(Rs_dc,Cs_dc,ts,steps = steps,from_t = t0,until_t = tf ,Delta_t = dt_samples,temp_dt = dt_integration,mode_consumption = mode_consumption ,mode_renewal = mode_renewal,mode_generation = mode_generation)

      call Compute_interactions_over_time_general_CR(v_ref = Cs,v_per = Cs_dc ,dt = dt_samples,steps = steps,M = M_CC,i_target = i_target) 
      
      select case (mode_renewal)
        case (1)
            addname = trim("_chemostat")
        case (2)
            addname = trim("_batch")
        case (3)
            addname = trim("_biotic")
        case default
            addname = trim("_biotic")
      end select

      path = trim("data/interactions_general_CR/")
      filename = trim(trim(path)//"MCC"//trim(addname)//".dat")
      print *, "  Save MCCs in: ",filename
      open(unit=1001, file=trim(filename), iostat=ios, status="unknown", action="write")
      if ( ios /= 0 ) stop "Error opening file "//filename
        do ii = 1, steps, 1
          write(1001,*) ts(ii),M_CC(ii),-e(1)*sum(a(1,:)*a(2,:)*Rs(ii,:))*ts(ii)
        end do
      close(1001)
      print *,"sum=", -e(1)*sum(a(1,:)*a(2,:)*Rs(1,:))

      filename = trim(trim(path)//"traj"//trim(addname)//".dat")
      print *, "  Save trajs in: ",filename
      open(unit=1001, file=trim(filename), iostat=ios, status="unknown", action="write")
      if ( ios /= 0 ) stop "Error opening file "//filename
        do ii = 1, steps, 1
          write(1001,*) ts(ii),Cs(ii,1),Cs(ii,2)
        end do
      close(1001)

      filename = trim(trim(path)//"traj_R"//trim(addname)//".dat")
      print *, "  Save trajs in: ",filename
      open(unit=1001, file=trim(filename), iostat=ios, status="unknown", action="write")
      if ( ios /= 0 ) stop "Error opening file "//filename
        do ii = 1, steps, 1
          write(1001,*) ts(ii),Rs(ii,1),Rs(ii,2)
        end do
      close(1001)
    end do
    

    filename = trim(trim(path)//"params"//trim(addname)//".dat")
    print *, "  Save params in: ",filename
    open(unit=2002, file=trim(filename), iostat=ios, status="unknown", action="write")
    if ( ios /= 0 ) stop "Error opening file "//filename
      write(2002,*) C_I
      write(2002,*) R_I
      do ii = 1, n_C, 1
        write(2002,*) a(ii,:)
      end do
      write(2002,*) cr
      write(2002,*) K
      write(2002,*) e
      write(2002,*) d
      write(2002,*) delta
      write(2002,*) seed
    close(2002)

    deallocate(C_I,R_I,Cs,Rs,Cs_dc,ts)
    deallocate(M_CC)
    
101 end program main_program