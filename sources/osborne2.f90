 Program osborne2
    use sort

    implicit none 
    
    integer :: allocerr,samples,noutliers,q,iterations,n_eval
    real(kind=8) :: fxk,fxtrial,ti,sigma
    real(kind=8), allocatable :: xtrial(:),faux(:),indices(:),nu_l(:),nu_u(:),opt_cond(:),&
                                 xinit(:),y(:),data(:,:),t(:),xbest(:)
    integer, allocatable :: Idelta(:),outliers(:),outliers_best(:)
    real(kind=8) :: fovo,delta,sigmin,gamma,start,finish
    
    ! LOCAL SCALARS
    logical :: checkder
    integer :: hnnzmax,inform,jcnnzmax,m,n,nvparam
    real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt,f,nlpsupn,snorm

    ! LOCAL ARRAYS
    character(len=80) :: specfnm,outputfnm,vparam(10) 
    logical :: coded(11)
    real(kind=8),   pointer :: l(:),u(:),x(:),xk(:),grad(:,:)

    logical,        pointer :: equatn(:),linear(:)
    real(kind=8),   pointer :: lambda(:)

    integer :: i
    real(kind=8), dimension(3,3) :: solutions

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    ! Reading data and storing it in the variables t and y
    Open(Unit = 100, File = trim(pwd)//"/../data/osborne2.txt", ACCESS = "SEQUENTIAL")

    ! Set parameters
    read(100,*) samples

    n = 12

    allocate(t(samples),y(samples),x(n),xk(n-1),xbest(n-1),xtrial(n-1),l(n),u(n),xinit(n-1),data(2,samples),faux(samples),&
    indices(samples),Idelta(samples),nu_l(n-1),nu_u(n-1),opt_cond(n-1),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    do i = 1, samples
        read(100,*) data(:,i)
    enddo

    close(100)

    ! Coded subroutines
    coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
    coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp

    ! Upper bounds on the number of sparse-matrices non-null elements
    jcnnzmax = 10000
    hnnzmax  = 10000

    ! Checking derivatives?
    checkder = .false.

    ! Parameters setting
    epsfeas   = 1.0d-08
    epsopt    = 1.0d-08
  
    efstain   = sqrt( epsfeas )
    eostain   = epsopt ** 1.5d0
  
    efacc     = sqrt( epsfeas )
    eoacc     = sqrt( epsopt )

    outputfnm = ''
    specfnm   = ''

    nvparam   = 0
    vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 0' 

    l(1:n) = -1.0d+20
    u(1:n-1) = 1.0d+20; u(n) = 0.0d0

    ! Number of days
    t(:) = data(1,:)
    y(:) = data(2,:)

    Open(Unit= 100, file = 'param.txt')
    read(100,*) delta,sigmin,gamma,noutliers
    close(100)

    q = samples - noutliers

    allocate(outliers(noutliers),outliers_best(noutliers),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    Open(Unit = 100, File = trim(pwd)//"/../output/sol_ls_osborne2.txt", ACCESS = "SEQUENTIAL")

    do i = 1, n-1
        read(100,*) xinit(i)
    enddo

    close(100)

    outliers(:) = 0
    xk(:) = xinit(:)
    call cpu_time(start)

    call ovo_algorithm(q,noutliers,t,y,indices,Idelta,samples,m,n,xtrial,&
        delta,sigmin,gamma,outliers,.false.,fovo,iterations,n_eval)

    call cpu_time(finish)
    
    write(*,100) "esta", noutliers,"&",fovo,"&",iterations,"&",n_eval,"&",finish-start,"\\"
    100 format (A5,1X,I2,1X,A1,1X,ES10.3,1X,A1,1X,I3,1X,A1,1X,I4,1X,A1,1X,ES10.3,1X,A2)

    Open(Unit = 98, File = trim(pwd)//"/../output/solution_osborne2.txt", ACCESS = "SEQUENTIAL")
    write(98,"(11F7.3)") xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7),xk(8),xk(9),xk(10),xk(11)

    Open(Unit = 99, File = trim(pwd)//"/../output/outliers_osborne2.txt", ACCESS = "SEQUENTIAL")
    write(99,"(I2)") noutliers

    do i = 1, noutliers
        write(99,"(I2)") outliers(i)
    enddo
    
    close(98)
    close(99)

    CONTAINS

    !==============================================================================
    ! MAIN ALGORITHM
    !==============================================================================
    subroutine ovo_algorithm(q,noutliers,t,y,indices,Idelta,samples,m,n,xtrial, &
                             delta,sigmin,gamma,outliers,print_iter,fovo,iterations,n_eval)
        implicit none

        logical,        intent(in) :: print_iter
        integer,        intent(in) :: q,noutliers,samples,n
        real(kind=8),   intent(in) :: t(samples),y(samples),delta,sigmin,gamma
        integer,        intent(inout) :: Idelta(samples),m
        real(kind=8),   intent(inout) :: indices(samples),xtrial(n-1),fovo
        integer,        intent(inout) :: outliers(noutliers),iterations,n_eval

        integer, parameter  :: max_iter = 10000, max_iter_sub = 100, kflag = 2
        integer             :: iter,iter_sub,i,j
        real(kind=8)        :: gaux,terminate,alpha,epsilon

        alpha   = 1.0d-8
        epsilon = 1.0d-3
        iter    = 0 
        
        indices(:) = (/(i, i = 1, samples)/)
    
        ! Scenarios
        do i = 1, samples
            call fi(xk,i,n,t,y,samples,faux(i))
        end do
    
        ! Sorting
        call DSORT(faux,indices,samples,kflag)

        ! q-Order-Value function 
        fxk = faux(q)
        n_eval = 1

        call mount_Idelta(faux,delta,q,indices,samples,Idelta,m)

        if (print_iter) then

            print*,"-----------------------------------------------------------------------------"
            write(*,10) "Iterations","Inter. Iter.","Objective func.","Optimality cond.","Idelta","Sum LM"
            10 format (2X,A11,2X,A12,2X,A15,2X,A16,2X,A6,2X,A6)
            print*,"-----------------------------------------------------------------------------"

            write(*,20)  0,"-",fxk,"-",m,"-"
            20 format (7X,I1,13X,A1,6X,ES14.6,12X,A1,11X,I2,6X,A1)

        endif

        do
            iter = iter + 1
    
            allocate(equatn(m),linear(m),lambda(m),grad(m,n-1),stat=allocerr)
    
            if ( allocerr .ne. 0 ) then
                write(*,*) 'Allocation error in main program'
                stop
            end if
    
            equatn(:) = .false.
            linear(:) = .false.
            lambda(:) = 0.0d0

            do i = 1, m
                ti = t(Idelta(i))

                call model(xk,Idelta(i),n,t,samples,gaux)

                gaux = gaux - y(Idelta(i))
    
                grad(i,1) = exp(-ti * xk(5))
                grad(i,2) = exp(-xk(6) * (ti - xk(9))**2)
                grad(i,3) = exp(-xk(7) * (ti - xk(10))**2)
                grad(i,4) = exp(-xk(8) * (ti - xk(11))**2)
                grad(i,5) = -ti * x(1) * exp(-ti * x(5))
                grad(i,6) = (-x(2) * (ti - x(9))**2) * exp(-x(6) * (ti - x(9))**2)
                grad(i,7) = (-x(3) * (ti - x(10))**2) * exp(-x(7) * (ti - x(10))**2)
                grad(i,8) = (-x(4) * (ti - x(11))**2) * exp(-x(8) * (ti - x(11))**2)
                grad(i,9) = 2.d0 * x(2) * x(6) * (ti - x(9)) * exp(-x(6) * (ti - x(9))**2)
                grad(i,10) = 2.d0 * x(3) * x(7) * (ti - x(10)) * exp(-x(7) * (ti - x(10))**2)
                grad(i,11) = 2.d0 * x(4) * x(8) * (ti - x(11)) * exp(-x(8) * (ti - x(11))**2)
    
                grad(i,:) = gaux * grad(i,:)
            end do
            
            sigma = sigmin
            ! if (iter .eq. 1) then
            !     sigma = sigmin
            ! else
            !     sigma = sigma / gamma
            ! endif
            
            iter_sub = 1
            x(:) = (/xk(:),0.0d0/)
    
            ! Minimizing using ALGENCAN
            do 
                call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,   &
                    myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax,    &
                    hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm,   &
                    specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,    &
                    checkder,f,cnorm,snorm,nlpsupn,inform)

                xtrial(1:n-1) = x(1:n-1)
                indices(:) = (/(i, i = 1, samples)/)
    
                ! Scenarios
                do i = 1, samples
                    call fi(xtrial,i,n,t,y,samples,faux(i))
                end do
    
                ! Sorting
                call DSORT(faux,indices,samples,kflag)
        
                fxtrial = faux(q)
                n_eval = n_eval + 1
        
                ! Test the sufficient descent condition
                if (fxtrial .le. (fxk - alpha * norm2(xtrial(1:n-1) - xk(1:n-1))**2)) exit
                if (iter_sub .ge. max_iter_sub) exit
    
                sigma = gamma * sigma
                iter_sub = iter_sub + 1
            end do ! End of internal iterations
    
            opt_cond(:) = 0.0d0
            nu_l(:) = 0.0d0
            nu_u(:) = 0.0d0
    
            do j = 1, n-1
                if (xtrial(j) .le. l(j)) then 
                    do i = 1, m
                        nu_l(j) = nu_l(j) + lambda(i) * grad(i,j)
                    end do
                else if (xtrial(j) .ge. u(j)) then
                    do i = 1, m
                        nu_u(j) = nu_u(j) - lambda(i) * grad(i,j)
                    end do
                end if
            end do
    
            do i = 1, m
                opt_cond(:) = opt_cond(:) + lambda(i) * grad(i,:)
            enddo
    
            opt_cond(:) = opt_cond(:) + nu_u(:) - nu_l(:)
            terminate = norm2(opt_cond)
            ! terminate = norm2(xk-xtrial)

            if (print_iter) then

                write(*,30)  iter,iter_sub,fxtrial,terminate,m,sum(lambda(:))
                30 format (2X,I6,10X,I4,6X,ES14.6,4X,ES14.6,6X,I2,5X,F3.1)

            endif

            deallocate(lambda,equatn,linear,grad)
            fxk = fxtrial
            xk(1:n-1) = xtrial(1:n-1)

            if (terminate .lt. epsilon) exit
            if (iter .ge. max_iter) exit
    
            call mount_Idelta(faux,delta,q,indices,samples,Idelta,m)
            
        end do ! End of Main Algorithm
        print*,"-----------------------------------------------------------------------------"

        outliers(:) = int(indices(samples - noutliers + 1:))
        fovo = fxk
        iterations = iter
        
    end subroutine ovo_algorithm

    function drand(ix)

        implicit none
      
        ! This is the random number generator of Schrage:
        !
        ! L. Schrage, A more portable Fortran random number generator, ACM
        ! Transactions on Mathematical Software 5 (1979), 132-138.
      
        ! FUNCTION TYPE
        real(kind=8) :: drand
      
        ! SCALAR ARGUMENT
        real(kind=8), intent(inout) :: ix
      
        ! LOCAL ARRAYS
        real(kind=8) :: a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      
        data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/
      
        xhi= ix/b16
        xhi= xhi - dmod(xhi,1.d0)
        xalo= (ix-xhi*b16)*a
        leftlo= xalo/b16
        leftlo= leftlo - dmod(leftlo,1.d0)
        fhi= xhi*a + leftlo
        k= fhi/b15
        k= k - dmod(k,1.d0)
        ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
        if (ix.lt.0) ix= ix + p
        drand= ix*4.656612875d-10
      
        return
      
      end function drand

    !==============================================================================
    ! MOUNT THE SET OF INDICES I(x,delta)
    !==============================================================================
    subroutine mount_Idelta(f,delta,q,indices,samples,Idelta,m)
        implicit none

        integer,        intent(in) :: samples,q
        real(kind=8),   intent(in) :: delta,f(samples),indices(samples)
        integer,        intent(out) :: Idelta(samples),m
        integer :: i
        real(kind=8) :: fq

        Idelta(:) = 0
        fq = f(q)
        m = 0

        do i = 1, samples
            if (abs(fq - f(i)) .le. delta) then
                m = m + 1
                Idelta(m) = int(indices(i))
            end if
        end do

    end subroutine

    !==============================================================================
    ! QUADRATIC ERROR OF EACH SCENARIO
    !==============================================================================
    subroutine fi(x,i,n,t,y,samples,res)
        implicit none

        integer,        intent(in) :: n,i,samples
        real(kind=8),   intent(in) :: x(n-1),t(samples),y(samples)
        real(kind=8),   intent(out) :: res
        
        call model(x,i,n,t,samples,res)
        res = res - y(i)
        res = 0.5d0 * (res**2)

    end subroutine fi

    !==============================================================================
    ! MODEL TO BE FITTED TO THE DATA
    !==============================================================================
    subroutine model(x,i,n,t,samples,res)
        implicit none 

        integer,        intent(in) :: n,i,samples
        real(kind=8),   intent(in) :: x(n-1),t(samples)
        real(kind=8),   intent(out) :: res

        res = x(1) * exp(-t(i) * x(5)) + x(2) * exp(-x(6) * (t(i) - x(9))**2) + &
        x(3) * exp(-x(7) * (t(i) - x(10))**2) + x(4) * exp(-x(8) * (t(i) - x(11))**2)

    end subroutine model

    !==============================================================================
    ! SUBROUTINES FOR ALGENCAN
    !==============================================================================

    !******************************************************************************
    ! OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalf(n,x,f,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: f

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)

        ! Compute objective function

        flag = 0

        f = x(n)

    end subroutine myevalf

    !******************************************************************************
    ! GRADIENT OF THE OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalg(n,x,g,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        integer, intent(out) :: flag

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n)

        ! Compute gradient of the objective function

        flag = 0

        g(1:n-1) = 0.0d0
        g(n)     = 1.0d0

    end subroutine myevalg

    !******************************************************************************
    ! HESSIAN FOR THE OBJECTIVE FUNCTION
    !******************************************************************************
    subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)
        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,n
        integer, intent(out) :: flag,hnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hcol(lim),hrow(lim)
        real(kind=8), intent(in)  :: x(n)
        real(kind=8), intent(out) :: hval(lim)

        ! Compute (lower triangle of the) Hessian of the objective function
        flag = 0
        lmem = .false.
        hnnz = 0
    end subroutine myevalh

    !******************************************************************************
    ! CONSTRAINTS
    !******************************************************************************
    subroutine myevalc(n,x,ind,c,flag)
        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: ind,n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: c

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)

        ! Compute ind-th constraint
        flag = 0

        c = dot_product(x(1:n-1) - xk(1:n-1),grad(ind,1:n-1)) + &
        (sigma * 0.5d0) * dot_product(x(1:n-1) - xk(1:n-1),x(1:n-1) - xk(1:n-1)) - x(n)

    end subroutine myevalc

    !******************************************************************************
    ! JACOBIAN OF THE CONSTRAINTS
    !******************************************************************************
    subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: ind,lim,n
        integer, intent(out) :: flag,jcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: jcval(lim)

        integer :: i

        flag = 0
        lmem = .false.

        jcnnz = n

        if ( jcnnz .gt. lim ) then
            lmem = .true.
            return
        end if

        jcvar(1:n) = (/(i, i = 1, n)/)
        jcval(1:n) = (/grad(ind,1:n-1) + sigma * (x(1:n-1) - xk(1:n-1)),-1.0d0/)
        ! jcval(1:n) = (/(grad(ind,i) + sigma * (x(i) - xk(i)), i = 1, n-1), -1.0d0/)

    end subroutine myevaljac

    !******************************************************************************
    ! HESSIAN OF THE CONSTRAINTS
    !******************************************************************************
    subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: ind,lim,n
        integer, intent(out) :: flag,hcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hccol(lim),hcrow(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: hcval(lim)

        flag = 0
        lmem = .false.
    
        hcnnz = n - 1
    
        if ( hcnnz .gt. lim ) then
            lmem = .true.
            return
        end if
    
        hcrow(1:n-1) = (/(i, i = 1, n-1)/)
        hccol(1:n-1) = (/(i, i = 1, n-1)/)
        hcval(1:n-1) = sigma

    end subroutine myevalhc

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalfc(n,x,f,m,c,flag)

        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: f

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: c(m)

        flag = - 1

    end subroutine myevalfc

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,m,n
        integer, intent(out) :: flag,jcnnz

        ! ARRAY ARGUMENTS
        integer, intent(out) :: jcfun(lim),jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n),jcval(lim)

        flag = - 1

    end subroutine myevalgjac

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(inout) :: gotj
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        character, intent(in) :: work

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(inout) :: p(m),q(n)
        real(kind=8), intent(out) :: g(n)

        flag = - 1

    end subroutine myevalgjacp

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,m,n
        integer, intent(out) :: flag,hlnnz
        real(kind=8), intent(in) :: sf

        ! ARRAY ARGUMENTS
        integer, intent(out) :: hlcol(lim),hlrow(lim)
        real(kind=8), intent(in) :: lambda(m),sc(m),x(n)
        real(kind=8), intent(out) :: hlval(lim)

        flag = - 1

    end subroutine myevalhl

    ! ******************************************************************
    ! ******************************************************************

    subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

        implicit none

        ! SCALAR ARGUMENTS
        logical, intent(inout) :: goth
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: sf

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: lambda(m),p(n),sc(m),x(n)
        real(kind=8), intent(out) :: hp(n)

        flag = - 1

    end subroutine myevalhlp
end Program osborne2
