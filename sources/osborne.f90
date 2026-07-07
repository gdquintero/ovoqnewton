 Program osborne
    use sort

    implicit none

    ! Structure created for the construction of the Bkj matrices
    type :: pdata_type
        real(kind=8), allocatable :: Bkj(:,:),identity(:,:),eig_hess(:),M
        character(len=1) :: JOBZ,UPLO ! lapack variables
        integer :: LDA,LWORK,INFO,NRHS,LDB ! lapack variables
        real(kind=8), allocatable :: WORK(:),IPIV(:) ! lapack variables
    end type pdata_type

    integer :: allocerr,samples,noutliers,q,iterations,n_eval,ntrials,itrial,iter_best,neval_best,cauchy_flag
    integer :: ndata,obar,o,omin,omax,ostep,nvals,k,kstar
    real(kind=8) :: fxk,fxtrial,ti,sigma,seed,dseed,fovo_best,time_best,tiempo
    real(kind=8) :: E1,E2,ytrue,p1,p2,ratio,ratio_best,rperturb
    real(kind=8), allocatable :: xtrial(:),faux(:),indices(:),nu_l(:),nu_u(:),opt_cond(:),&
                                 xinit(:),y(:),t(:),xbest(:),xtrue(:)
    integer, allocatable :: Idelta(:),outliers(:),outliers_best(:)
    real(kind=8) :: fovo,delta,sigmin,gamma,start,finish

    ! Arrays indexed by the swept value of o (drop detection)
    integer :: o_arr(4000),it_arr(4000),ne_arr(4000)
    real(kind=8) :: fo_arr(4000),tm_arr(4000)

    ! LOCAL SCALARS
    logical :: checkder
    integer :: hnnzmax,inform,jcnnzmax,m,n,nvparam
    real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt,f,nlpsupn,snorm

    ! LOCAL ARRAYS
    character(len=80) :: specfnm,outputfnm,vparam(10)
    character(len=8)  :: tag
    logical :: coded(11)
    real(kind=8),   pointer :: l(:),u(:),x(:),xk(:),grad(:,:),hess(:,:,:)

    logical,        pointer :: equatn(:),linear(:)
    real(kind=8),   pointer :: lambda(:)

    type(pdata_type), target :: pdata
    integer :: i
    logical, pointer :: cauchy

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    ! ------------------------------------------------------------------
    ! Read parameters: delta, sigmin, gamma, cauchy_flag, ndata
    ! ------------------------------------------------------------------
    Open(Unit = 100, file = 'param.txt')
    read(100,*) delta,sigmin,gamma,cauchy_flag,ndata
    close(100)

    n = 6                 ! 5 model parameters + 1 epigraph variable w
    samples = ndata

    allocate(cauchy)
    cauchy = (cauchy_flag .eq. 1)
    if (cauchy) then
        tag = 'fo'
    else
        tag = 'qn'
    end if

    ! lapack variables
    pdata%JOBZ = 'N'
    pdata%UPLO = 'U'
    pdata%LDA = n-1
    pdata%LDB = n-1
    pdata%LWORK = 3*(n-1) - 1
    pdata%NRHS = 1

    allocate(t(samples),y(samples),x(n),xk(n-1),xbest(n-1),xtrial(n-1),l(n),u(n),xinit(n-1),xtrue(n-1),faux(samples),&
    indices(samples),Idelta(samples),nu_l(n-1),nu_u(n-1),opt_cond(n-1),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    allocate(pdata%WORK(pdata%LWORK),pdata%IPIV(n-1),pdata%Bkj(n-1,n-1),pdata%eig_hess(n-1),pdata%identity(n-1,n-1),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    ! Coded subroutines
    coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
    coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp

    ! Upper bounds on the number of sparse-matrices non-null elements
    jcnnzmax = samples*n + 100
    hnnzmax  = samples*n*n + 100

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

    nvparam   = 1
    vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 0'

    ! Box constraints. Osborne 1: keep x4,x5 > 0 so the exponentials decay; x3<0
    ! at the true solution, so its box straddles zero.
    l(1) = -10.0d0;  u(1) = 10.0d0
    l(2) =   0.1d0;  u(2) = 10.0d0
    l(3) = -10.0d0;  u(3) = 10.0d0
    l(4) = 0.001d0;  u(4) =  1.0d0
    l(5) = 0.001d0;  u(5) =  1.0d0
    l(n) = -1.0d+20; u(n) =  0.0d0

    ! Ground-truth parameters (standard Osborne 1 solution, MGH #17)
    xtrue(:) = (/0.3754d0, 1.9358d0, -1.4647d0, 0.01287d0, 0.02212d0/)

    ! Standard MGH starting point for Osborne 1
    xinit(:) = (/0.5d0, 1.5d0, -1.0d0, 0.01d0, 0.02d0/)

    ! ------------------------------------------------------------------
    ! Generate synthetic data: t(i) in [0,320], y(i) = model(t,xtrue) + noise,
    ! each observation an outlier with probability 0.1 (above with prob. 0.8,
    ! below with prob. 0.2), on the scale of the model (y in ~[0.4,0.95]).
    ! ------------------------------------------------------------------
    dseed = 777.0d0
    obar  = 0
    do i = 1, samples
        t(i)  = dble(i-1) / dble(samples-1) * 320.0d0
        ytrue = xtrue(1) + xtrue(2)*exp(-xtrue(4)*t(i)) + xtrue(3)*exp(-xtrue(5)*t(i))
        p1 = drand(dseed)
        if (p1 .lt. 0.1d0) then
            obar = obar + 1
            p2 = drand(dseed)
            if (p2 .lt. 0.8d0) then
                y(i) = ytrue + drand(dseed) * (3.0d0 - ytrue)      ! outlier above: U[ytrue,3]
            else
                y(i) = -1.0d0 + drand(dseed) * (ytrue + 1.0d0)     ! outlier below: U[-1,ytrue]
            end if
        else
            y(i) = ytrue + (2.0d0*drand(dseed) - 1.0d0) * 0.05d0   ! clean: y + U[-0.05,0.05]
        end if
    end do

    ! ------------------------------------------------------------------
    ! Range of o to sweep, following the original scalability experiment
    ! (expected number of outliers is 0.1*ndata).
    ! ------------------------------------------------------------------
    if (ndata .eq. 100) then
        omin = 5;   omax = 15;   ostep = 1
    else if (ndata .eq. 1000) then
        omin = 50;  omax = 150;  ostep = 1
    else if (ndata .eq. 10000) then
        omin = 500; omax = 1500; ostep = 10
    else
        omin = ndata/100 - ndata/20; omax = ndata/100 + ndata/20; ostep = max(1,ndata/1000)
        if (omin .lt. 1) omin = 1
    end if

    ! Write the (o, f(x*)) curve to a file for plotting (Figure-8 style)
    Open(Unit = 97, File = trim(pwd)//"/../output/osborne_"//trim(tag)//"_curve.txt", ACCESS = "SEQUENTIAL")

    ntrials = 100
    nvals   = 0

    ! ==================================================================
    ! Sweep over o
    ! ==================================================================
    do o = omin, omax, ostep
        nvals = nvals + 1
        noutliers = o
        q = samples - noutliers

        allocate(outliers(noutliers),outliers_best(noutliers),stat=allocerr)
        if ( allocerr .ne. 0 ) then
            write(*,*) 'Allocation error (outliers)'
            stop
        end if
        outliers(:) = 0

        ! Multistart: same 100 starting points for every o and both methods,
        ! obtained by a relative perturbation of xinit projected onto the box.
        seed = 123456.0d0
        fovo_best = huge(1.0d0)

        do itrial = 1, ntrials
            xk(:) = xinit(:)
            do i = 1, n-1
                rperturb = (drand(seed) - 0.5d0)          ! r in [-0.5,0.5]
                xk(i) = xk(i) + rperturb * abs(xinit(i))
                xk(i) = max(l(i), min(u(i), xk(i)))       ! project onto the box
            end do

            call cpu_time(start)
            call ovo_algorithm(pdata,q,noutliers,t,y,indices,Idelta,samples,m,n,xtrial,&
            delta,sigmin,gamma,outliers,.false.,fovo,iterations,n_eval)
            call cpu_time(finish)
            tiempo = finish - start

            if (fovo .lt. fovo_best) then
                time_best = tiempo
                fovo_best = fovo
                xbest(:) = xk(:)
                outliers_best(:) = outliers(:)
                iter_best = iterations
                neval_best = n_eval
            endif
        enddo

        o_arr(nvals)  = o
        fo_arr(nvals) = fovo_best
        it_arr(nvals) = iter_best
        ne_arr(nvals) = neval_best
        tm_arr(nvals) = time_best

        write(97,"(I8,1X,ES14.6)") o, fovo_best

        deallocate(outliers,outliers_best)
    end do

    close(97)

    ! ------------------------------------------------------------------
    ! Detect the drop: the o at which f(x*) falls most sharply (largest
    ! ratio between consecutive best values). Report it as the number of
    ! outliers identified.
    ! ------------------------------------------------------------------
    kstar = 1
    ratio_best = 0.0d0
    do k = 2, nvals
        ratio = fo_arr(k-1) / max(fo_arr(k), 1.0d-30)
        if (ratio .gt. ratio_best) then
            ratio_best = ratio
            kstar = k
        end if
    end do

    ! ------------------------------------------------------------------
    ! Table-5-style row for the detected o:
    ! method  ndata  obar  o  f(x*)  #it  #fcnt  Time  Time/#fcnt
    ! ------------------------------------------------------------------
    write(*,200) trim(tag), samples, obar, o_arr(kstar), fo_arr(kstar), &
                 it_arr(kstar), ne_arr(kstar), tm_arr(kstar), tm_arr(kstar)/dble(max(1,ne_arr(kstar)))
    200 format (A3,1X,"&",1X,I8,1X,"&",1X,I7,1X,"&",1X,I7,1X,"&",1X,ES11.3,1X,"&",1X,&
                I4,1X,"&",1X,I5,1X,"&",1X,ES11.3,1X,"&",1X,ES11.3,1X,"\\")

    CONTAINS

    !==============================================================================
    ! MAIN ALGORITHM
    !==============================================================================
    subroutine ovo_algorithm(pdata,q,noutliers,t,y,indices,Idelta,samples,m,n,xtrial, &
                             delta,sigmin,gamma,outliers,print_iter,fovo,iterations,n_eval)
        implicit none

        logical,        intent(in) :: print_iter
        integer,        intent(in) :: q,noutliers,samples,n
        real(kind=8),   intent(in) :: t(samples),y(samples),delta,sigmin,gamma
        integer,        intent(inout) :: Idelta(samples),m
        real(kind=8),   intent(inout) :: indices(samples),xtrial(n-1),fovo
        integer,        intent(inout) :: outliers(noutliers),iterations,n_eval

        type(pdata_type), intent(inout) :: pdata

        integer, parameter  :: max_iter = 10000, max_iter_sub = 100, kflag = 2
        integer             :: iter,iter_sub,i,j
        real(kind=8)        :: terminate,alpha,epsilon,lambda_min,aux_iden,lambda_max,ri,E1,E2
        real(kind=8)        :: gr1,gr2,gr3,gr4,gr5

        alpha   = 1.0d-8
        epsilon = 1.0d-4
        iter    = 0
        pdata%M = 1

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

        do
            iter = iter + 1

            allocate(equatn(m),linear(m),lambda(m),grad(m,n-1),hess(m,n-1,n-1),stat=allocerr)

            if ( allocerr .ne. 0 ) then
                write(*,*) 'Allocation error in main program'
                stop
            end if

            equatn(:) = .false.
            linear(:) = .false.
            lambda(:) = 0.0d0

            do i = 1, m
                ti = t(Idelta(i))
                E1 = exp(-xk(4) * ti)
                E2 = exp(-xk(5) * ti)
                ri = xk(1) + xk(2)*E1 + xk(3)*E2 - y(Idelta(i))   ! residual r = model - y

                ! Gradient of the residual, grad_r = (1, E1, E2, -x2*t*E1, -x3*t*E2)
                gr1 = 1.0d0
                gr2 = E1
                gr3 = E2
                gr4 = -xk(2) * ti * E1
                gr5 = -xk(3) * ti * E2

                ! Gradient of f_i = 0.5*r^2 is r*grad_r
                grad(i,1) = ri * gr1
                grad(i,2) = ri * gr2
                grad(i,3) = ri * gr3
                grad(i,4) = ri * gr4
                grad(i,5) = ri * gr5

                ! Hessian of f_i = grad_r*grad_r^T + r*Hess_r
                hess(i,1,1) = gr1*gr1
                hess(i,1,2) = gr1*gr2
                hess(i,1,3) = gr1*gr3
                hess(i,1,4) = gr1*gr4
                hess(i,1,5) = gr1*gr5
                hess(i,2,2) = gr2*gr2
                hess(i,2,3) = gr2*gr3
                hess(i,2,4) = gr2*gr4 + ri * (-ti * E1)
                hess(i,2,5) = gr2*gr5
                hess(i,3,3) = gr3*gr3
                hess(i,3,4) = gr3*gr4
                hess(i,3,5) = gr3*gr5 + ri * (-ti * E2)
                hess(i,4,4) = gr4*gr4 + ri * (xk(2) * ti*ti * E1)
                hess(i,4,5) = gr4*gr5
                hess(i,5,5) = gr5*gr5 + ri * (xk(3) * ti*ti * E2)

                ! Symmetrize
                hess(i,2,1) = hess(i,1,2)
                hess(i,3,1) = hess(i,1,3)
                hess(i,4,1) = hess(i,1,4)
                hess(i,5,1) = hess(i,1,5)
                hess(i,3,2) = hess(i,2,3)
                hess(i,4,2) = hess(i,2,4)
                hess(i,5,2) = hess(i,2,5)
                hess(i,4,3) = hess(i,3,4)
                hess(i,5,3) = hess(i,3,5)
                hess(i,5,4) = hess(i,4,5)

                ! Copy the Hessian, since DSYEV overwrites the input matrix
                pdata%Bkj(:,:) = hess(i,:,:)

                ! Eigenvalues of the Hessian via DSYEV
                call dsyev(pdata%JOBZ,pdata%UPLO,n-1,pdata%Bkj,pdata%LDA,&
                pdata%eig_hess,pdata%WORK,pdata%LWORK,pdata%INFO)

                lambda_min = minval(pdata%eig_hess)

                ! Shift to make the curvature matrix positive definite
                aux_iden =  max(0.d0,-lambda_min + 1.d-8)
                call dlaset('A',n-1,n-1,0.0d0,aux_iden,pdata%identity,n-1)
                hess(i,:,:) = hess(i,:,:) + pdata%identity(:,:)

                lambda_max = maxval(pdata%eig_hess) + aux_iden

                ! Enforce ||B_{k,j}|| <= M
                if (lambda_max > pdata%M) then
                    hess(i,:,:) = (pdata%M / lambda_max) * hess(i,:,:)
                end if

                if (cauchy .eqv. .true.) then
                    hess(:,:,:) = 0.d0
                end if
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

                ! Sufficient descent condition
                if (fxtrial .le. (fxk - alpha * dot_product(xtrial(1:n-1) - xk(1:n-1),xtrial(1:n-1) - xk(1:n-1)))) exit
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

            deallocate(lambda,equatn,linear,grad,hess)
            fxk = fxtrial
            xk(1:n-1) = xtrial(1:n-1)

            if (terminate .lt. epsilon) exit
            if (iter .ge. max_iter) exit

            call mount_Idelta(faux,delta,q,indices,samples,Idelta,m)

        end do ! End of Main Algorithm

        outliers(:) = int(indices(samples - noutliers + 1:))
        fovo = fxk
        iterations = iter

    end subroutine ovo_algorithm

    function drand(ix)

        implicit none

        ! Schrage's portable random number generator.
        real(kind=8) :: drand
        real(kind=8), intent(inout) :: ix
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
        real(kind=8),   intent(out) :: res
        real(kind=8),   intent(in) :: x(n-1),t(samples)

        ! Osborne 1 model: x1 + x2*exp(-x4*t) + x3*exp(-x5*t)
        res = x(1) + x(2)*exp(-x(4)*t(i)) + x(3)*exp(-x(5)*t(i))

    end subroutine model

    !==============================================================================
    ! SUBROUTINES FOR ALGENCAN
    !==============================================================================

    subroutine myevalf(n,x,f,flag)
        implicit none
        integer, intent(in) :: n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: f
        real(kind=8), intent(in) :: x(n)

        flag = 0
        f = x(n)
    end subroutine myevalf

    subroutine myevalg(n,x,g,flag)
        implicit none
        integer, intent(in) :: n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n)

        flag = 0
        g(1:n-1) = 0.0d0
        g(n)     = 1.0d0
    end subroutine myevalg

    subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)
        implicit none
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,n
        integer, intent(out) :: flag,hnnz
        integer, intent(out) :: hcol(lim),hrow(lim)
        real(kind=8), intent(in)  :: x(n)
        real(kind=8), intent(out) :: hval(lim)

        flag = 0
        lmem = .false.
        hnnz = 0
    end subroutine myevalh

    subroutine myevalc(n,x,ind,c,flag)
        implicit none
        integer, intent(in) :: ind,n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: c
        real(kind=8), intent(in) :: x(n)

        flag = 0

        c = dot_product(x(1:n-1) - xk(1:n-1),grad(ind,1:n-1)) + 0.5d0 * &
            (dot_product(x(1:n-1) - xk(1:n-1),matmul(hess(ind,1:n-1,1:n-1),x(1:n-1) - xk(1:n-1))) + &
            sigma * dot_product(x(1:n-1) - xk(1:n-1),x(1:n-1) - xk(1:n-1))) - x(n)
    end subroutine myevalc

    subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)
        implicit none
        logical, intent(out) :: lmem
        integer, intent(in) :: ind,lim,n
        integer, intent(out) :: flag,jcnnz
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
        jcval(1:n) = (/grad(ind,1:n-1) + &
                    matmul(hess(ind,1:n-1,1:n-1),x(1:n-1) - xk(1:n-1)) + &
                    sigma * (x(1:n-1) - xk(1:n-1)), -1.0d0/)
    end subroutine myevaljac

    subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)
        implicit none
        logical, intent(out) :: lmem
        integer, intent(in) :: ind,lim,n
        integer, intent(out) :: flag,hcnnz
        integer, intent(out) :: hccol(lim),hcrow(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: hcval(lim)
        integer :: i,j

        flag = 0
        lmem = .false.

        hcnnz = 0
        do j = 1, n-1
            do i = j, n-1
                hcnnz = hcnnz + 1

                if ( hcnnz .gt. lim ) then
                    lmem = .true.
                    return
                end if

                hcrow(hcnnz) = i
                hccol(hcnnz) = j

                if ( i .eq. j ) then
                    hcval(hcnnz) = hess(ind,i,j) + sigma
                else
                    hcval(hcnnz) = hess(ind,i,j)
                end if
            end do
        end do
    end subroutine myevalhc

    subroutine myevalfc(n,x,f,m,c,flag)
        implicit none
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        real(kind=8), intent(out) :: f
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: c(m)

        flag = - 1
    end subroutine myevalfc

    subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)
        implicit none
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,m,n
        integer, intent(out) :: flag,jcnnz
        integer, intent(out) :: jcfun(lim),jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n),jcval(lim)

        flag = - 1
    end subroutine myevalgjac

    subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)
        implicit none
        logical, intent(inout) :: gotj
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        character, intent(in) :: work
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(inout) :: p(m),q(n)
        real(kind=8), intent(out) :: g(n)

        flag = - 1
    end subroutine myevalgjacp

    subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)
        implicit none
        logical, intent(out) :: lmem
        integer, intent(in) :: lim,m,n
        integer, intent(out) :: flag,hlnnz
        real(kind=8), intent(in) :: sf
        integer, intent(out) :: hlcol(lim),hlrow(lim)
        real(kind=8), intent(in) :: lambda(m),sc(m),x(n)
        real(kind=8), intent(out) :: hlval(lim)

        flag = - 1
    end subroutine myevalhl

    subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)
        implicit none
        logical, intent(inout) :: goth
        integer, intent(in) :: m,n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: sf
        real(kind=8), intent(in) :: lambda(m),p(n),sc(m),x(n)
        real(kind=8), intent(out) :: hp(n)

        flag = - 1
    end subroutine myevalhlp
end Program osborne
