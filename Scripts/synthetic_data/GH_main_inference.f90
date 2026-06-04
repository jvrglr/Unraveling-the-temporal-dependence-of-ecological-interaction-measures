program main_program
    !Main program, depends on the rest of modules (declarations_module, functions and subroutines))
    use declarations_module
    use functions
    use subroutines
    implicit none
    double precision :: tf,dc,g_C_dC,g_C_base,alpha,dt_integration,dt_samples
    integer*8 :: ii,points,jj
    character (len=100) :: filename,path,name,addname
    character(8)  :: date
    double precision :: dt
    double precision, dimension(:),allocatable :: C_I
    double precision, dimension(:),allocatable :: R_I
    double precision, dimension(:,:),allocatable :: estimation_alpha
    double precision :: dran_u

    call assingments()
    call dran_ini(1994)

    allocate(C_I(n_C),R_I(n_R),estimation_alpha(n_C,n_R))
    dc = 0.001d0
    tf = 0.001d0
    dt_integration = 0.001d0
    dt_samples = dt_integration
    C_I = 1.0d0
    R_I = 1.0d0
    do ii = 1, n_C, 1
      do jj = 1, n_R, 1
        a(ii,jj) = 2*dran_u()
      end do
    end do

    do ii = 1, n_C, 1
      do jj = 1, n_R, 1
      
        t = 0.0d0
        C = C_I
        R = R_I
        call Compute_g_empirical_general_consumer(tf,dt_integration,g_C_base,ii)
        C = C_I
        R = R_I
        R(jj) = R(jj) + dc
        
        call Compute_g_empirical_general_consumer(tf,dt_integration,g_C_dc,ii)
        alpha = (g_C_dc-g_C_base)/dc
        estimation_alpha(ii,jj) = alpha
        
      enddo
    enddo

     ! call date_and_time(DATE=date)
    date = "20250212"
    addname = "" !"_transient"
    path = trim("data/"//trim(date)//"/")
    name = trim("tf_0d001.dat")

    filename = trim(path)//"as_"//trim(name)
    print *, "  alphas saved in: ",filename
    open(unit=2002, file=filename, iostat=ios, status="unknown", action="write")
    if ( ios /= 0 ) stop "Error opening file "//filename
      do ii = 1, n_C, 1
        do jj = 1, n_R, 1
          write(2002,*) ii,jj,estimation_alpha(ii,jj),e(ii)*a(ii,jj)
        enddo
      enddo
    close(2002)

    deallocate(C_I,R_I,estimation_alpha)
    
  101 end program main_program
  