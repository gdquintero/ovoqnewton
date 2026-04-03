 Program main
    use sort

    implicit none 
    
    integer :: allocerr,samples,inf,sup
    real(kind=8) :: fxk,fxtrial,ti,sigma,delta,sigmin,gamma
    real(kind=8), allocatable :: xtrial(:),faux(:),indices(:),nu_l(:),nu_u(:),opt_cond(:),&
                                 xstar(:),y(:),data(:,:),t(:),xinit(:,:)
    integer, allocatable :: Idelta(:),outliers(:)
    
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
    Open(Unit = 100, File = trim(pwd)//"/../data/seropositives.txt", ACCESS = "SEQUENTIAL")

    ! Set parameters
    read(100,*) samples

    n = 4

    allocate(t(samples),y(samples),x(n),xk(n-1),xtrial(n-1),l(n),u(n),xstar(n-1),data(5,samples),faux(samples),&
    indices(samples),Idelta(samples),nu_l(n-1),nu_u(n-1),opt_cond(n-1),xinit(3,3),stat=allocerr)

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

    l(1:n-1) = 0.0d0; l(n) = -1.0d+20
    u(1:n-1) = 1.0d+20; u(n) = 0.0d0

    ! Number of days
    t(:) = data(1,:) ! Initial point
    ! t(:) = data(5,:) ! Midpoint

    Open(Unit= 100, file = 'param.txt')
    read(100,*) delta,sigmin,gamma
    close(100)

    inf = 0
    sup = 10

    allocate(outliers(3*samples*(sup-inf+1)),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    outliers(:) = 0

    call mixed_test(inf,sup,delta,sigmin,gamma,outliers,t,y,indices,Idelta,samples,m,n,xtrial)

    call export(xtrial,outliers,sup)

    CONTAINS

    subroutine mixed_test(out_inf,out_sup,delta,sigmin,gamma,outliers,t,y,indices,Idelta,samples,m,n,xtrial)
        implicit none

        integer,        intent(in) :: samples,n,out_inf,out_sup
        real(kind=8),   intent(in) :: t(samples),delta,sigmin,gamma
        integer,        intent(inout) :: Idelta(samples),outliers(3*samples*(sup-inf+1)),m
        real(kind=8),   intent(inout) :: indices(samples),xtrial(n-1),y(samples)

        integer :: noutliers,q,ind,iterations,n_eval
        real(kind=8) :: fovo,start,finish

        print*
        Print*, "OVO Algorithm for Measles"

        y(:) = data(2,:)

        Open(Unit = 100, file =trim(pwd)//"/../output/sol_ls_farrington.txt")
        read(100,*) xinit(1,:)
        read(100,*) xinit(2,:)
        read(100,*) xinit(3,:)
        close(100)


        do noutliers = out_inf,out_sup
            call cpu_time(start)
            q = samples - noutliers

            print*
            write(*,1100) "Number of outliers: ",noutliers
            ! xk(:) = (/0.379029d0,0.500859d0,0.016986d0/)
            xk(:) = xinit(1,:)

            ind = 1
            
            call ovo_algorithm(q,noutliers,t,y,indices,Idelta,samples,m,n,xtrial, &
            delta,sigmin,gamma,outliers(ind:ind+noutliers-1),fovo,iterations,n_eval)

            call cpu_time(finish)

            Open(Unit = 100, File =trim(pwd)//"/../output/solutions_mixed_measles.txt", ACCESS = "SEQUENTIAL")
            Open(Unit = 600, File =trim(pwd)//"/../output/outliers.txt", ACCESS = "SEQUENTIAL")

            write(100,1000) xtrial(1), xtrial(2), xtrial(3)

            write(*,1300) "esta", noutliers,"&",fovo,"&",iterations,"&",n_eval,"&",finish-start,"\\"
            
        enddo

        print*
        Print*, "OVO Algorithm for Mumps"

        y(:) = data(3,:)

        do noutliers = out_inf,out_sup
            call cpu_time(start)
            q = samples - noutliers
            print*
            write(*,1100) "Number of outliers: ",noutliers
            ! xk(:) = (/0.285745d0,0.424520d0,0.005894d0/)
            xk(:) = xinit(2,:)

            ind = ind + noutliers

            call ovo_algorithm(q,noutliers,t,y,indices,Idelta,samples,m,n,xtrial, &
            delta,sigmin,gamma,outliers(ind:ind+noutliers-1),fovo,iterations,n_eval)

            call cpu_time(finish)

            Open(Unit = 110, File =trim(pwd)//"/../output/solutions_mixed_mumps.txt", ACCESS = "SEQUENTIAL")

            write(110,1000) xtrial(1), xtrial(2), xtrial(3)

            write(*,1300) "esta", noutliers,"&",fovo,"&",iterations,"&",n_eval,"&",finish-start,"\\"
        enddo

        print*
        Print*, "OVO Algorithm for Rubella"

        y(:) = data(4,:)

        do noutliers = out_inf,out_sup
            call cpu_time(start)
            q = samples - noutliers
            print*
            write(*,1100) "Number of outliers: ",noutliers
            ! xk(:) = (/0.117309d0,0.341322d0,0.026605d0/)
            xk(:) = xinit(3,:)

            ind = ind + noutliers

            call ovo_algorithm(q,noutliers,t,y,indices,Idelta,samples,m,n,xtrial, &
            delta,sigmin,gamma,outliers(ind:ind+noutliers-1),fovo,iterations,n_eval)

            call cpu_time(finish)

            Open(Unit = 120, File =trim(pwd)//"/../output/solutions_mixed_rubella.txt", ACCESS = "SEQUENTIAL")

            write(120,1000) xtrial(1), xtrial(2), xtrial(3)

            write(*,1300) "esta", noutliers,"&",fovo,"&",iterations,"&",n_eval,"&",finish-start,"\\"
    
        enddo
        
        Open(Unit = 500, File =trim(pwd)//"/../output/num_mixed_test.txt", ACCESS = "SEQUENTIAL")
        write(500,1200) out_inf
        write(500,1200) out_sup
        
        1000 format (ES12.6,1X,ES12.6,1X,ES12.6)
        1100 format (1X,A20,I2)
        1200 format (I2)
        1300 format (A16,1X,I2,1X,A1,1X,ES10.3,1X,A1,1X,I3,1X,A1,1X,I3,1X,A1,1X,ES10.3,1X,A2)

        
        
    end subroutine mixed_test

    !==============================================================================
    ! MAIN ALGORITHM
    !==============================================================================
    subroutine ovo_algorithm(q,noutliers,t,y,indices,Idelta,samples,m,n,xtrial, &
                             delta,sigmin,gamma,outliers,fovo,iterations,n_eval)
        implicit none

        integer,        intent(in) :: q,noutliers,samples,n
        real(kind=8),   intent(in) :: t(samples),y(samples),delta,sigmin,gamma
        integer,        intent(inout) :: Idelta(samples),m
        real(kind=8),   intent(inout) :: indices(samples),xtrial(n-1),fovo
        integer,        intent(inout) :: outliers(noutliers),iterations,n_eval

        integer, parameter  :: max_iter = 1000, max_iter_sub = 100, kflag = 2
        integer             :: iter,iter_sub,i,j
        real(kind=8)        :: gaux1,gaux2,a,b,c,ebt,terminate,alpha,epsilon

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

        print*,"-----------------------------------------------------------------------------"
        write(*,10) "Iterations","Inter. Iter.","Objective func.","Optimality cond.","Idelta","Sum LM"
        10 format (2X,A11,2X,A12,2X,A15,2X,A16,2X,A6,2X,A6)
        print*,"-----------------------------------------------------------------------------"

        write(*,20)  0,"-",fxk,"-",m,"-"
        20 format (7X,I1,13X,A1,6X,ES14.6,12X,A1,11X,I2,6X,A1)

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
    
            a = xk(1)
            b = xk(2)
            c = xk(3)

            do i = 1, m
                ti = t(Idelta(i))

                ebt = exp(-b * ti)

                call model(xk,Idelta(i),n,t,samples,gaux1)

                gaux1 = y(Idelta(i)) - gaux1

                gaux2 = (a / b) * ti * ebt + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) - c * ti

                gaux2 = exp(gaux2)
    
                grad(i,1) = (1.0d0 / b**2) * (ebt * (ti * b + 1.0d0) - 1.0d0)
    
                grad(i,2) = ebt * ((-2.0d0 * a * ti / b**2) - ((a * ti**2) / b) - (2.0d0 * a / b**3) + &
                            (c / b**2) + (c * ti / b)) + (2.0d0 * a / b**3) - (c / b**2)
    
                grad(i,3) = (1.0d0 / b) * (1.0d0 - ebt) - ti
    
                grad(i,:) = gaux1 * gaux2 * grad(i,:)
            end do
    
            sigma = sigmin
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

            write(*,30)  iter,iter_sub,fxtrial,terminate,m,sum(lambda(:))
            30 format (2X,I6,10X,I4,6X,ES14.6,4X,ES14.6,6X,I2,5X,F3.1)

            deallocate(lambda,equatn,linear,grad)
            fxk = fxtrial
            xk(1:n-1) = xtrial(1:n-1)

            if (terminate .le. epsilon) exit
            if (iter .ge. max_iter) exit
    
            call mount_Idelta(faux,delta,q,indices,samples,Idelta,m)
            
        end do ! End of Main Algorithm
        print*,"-----------------------------------------------------------------------------"

        outliers(:) = int(indices(samples - noutliers + 1:))
        fovo = fxk
        iterations = iter
        
    end subroutine ovo_algorithm

    !==============================================================================
    ! EXPORT RESULT TO PLOT
    !==============================================================================
    subroutine export(xsol,outliers,noutliers)
        implicit none

        integer,        intent(in) :: noutliers,outliers(3*samples)
        real(kind=8),   intent(in) :: xsol(n-1,n-1)

        integer :: i

        Open(Unit = 100, File =trim(pwd)//"/../output/solutions_ovo.txt", ACCESS = "SEQUENTIAL")

        write(100,110) xsol(1,1), xsol(1,2), xsol(1,3)
        write(100,110) xsol(2,1), xsol(2,2), xsol(2,3)
        write(100,110) xsol(3,1), xsol(3,2), xsol(3,3)

        110 format (ES12.6,1X,ES12.6,1X,ES12.6)
    
        close(100)

        Open(Unit = 200, File =trim(pwd)//"/../output/outliers.txt", ACCESS = "SEQUENTIAL")

        write(200,210) noutliers

        do i = 1, 3*noutliers
            write(200,210) outliers(i)
        enddo

        210 format (I2)

        close(200)

    end subroutine export

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
        real(kind=8) :: a,b,c,ti,ebt

        a = x(1)
        b = x(2)
        c = x(3)
        ti = t(i)
        ebt = exp(-b * ti)

        res = (a / b) * ti * ebt
        res = res + (1.0d0 / b) * ((a / b) - c) * (ebt - 1.0d0) 
        res = 1.0d0 - exp(res - c * ti)

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
            (sigma * 0.5d0) * (norm2(x(1:n-1) - xk(1:n-1))**2) - x(n)

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
        jcval(1:n) = (/(grad(ind,i) + sigma * (x(i) - xk(i)), i = 1, n-1), -1.0d0/)

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
end Program main
