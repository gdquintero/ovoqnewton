!==============================================================================
! Algorithm 2.1 (regularized Newton-type method for box-constrained OVO) on a
! test instance, with a multistart over perturbed starting points.
!
! Usage (key=value arguments, all optional except inst):
!   ovo_bench inst=FILE o=K B=0|H|GN|QN M=1 mmode=scale|clip tau=1e-8
!             delta=0.1 sigmin=0.1 gamma=5 eps=1e-4 alpha=1e-8
!             ntrials=100 seed=123456 start=1|2 pert=0.5 pscale=abs|rel
!             project=1 maxit=10000 maxtime=1800 first=1 out=FILE.csv
!
! Curvature families (Step 2 of Algorithm 2.1):
!   B=0   first-order method (B_{k,j,i} = 0)
!   B=H   exact component Hessian, shifted to positive definiteness
!   B=GN  Gauss-Newton matrix grad r_i grad r_i^T
!   B=QN  a single BFGS matrix shared by all components, as in Andreani,
!         Martinez, Salvatierra and Yano (2006)
! Norm control ||B_{k,j,i}|| <= M (M <= 0 disables it):
!   mmode=scale  rescale the whole matrix by M/lambda_max
!   mmode=clip   replace eigenvalues larger than M by M
! Adaptive bound (cad=c > 0, overrides M): ||B_{k,j,i}|| <= c*sigma_{k,j},
! the matrix being rescaled for each value of sigma in the inner loop.
!
! Each trial appends one line to the CSV file `out`; the last two fields are
! the mean and the maximum size of the active set over the iterations.
! Status: 0 stopped at
! Step 5, 1 iteration limit, 2 CPU-time limit (maxtime seconds per run).
!==============================================================================
program ovo_bench

    use sort
    use hyperdual
    use models

    implicit none

    ! Instance data
    character(len=256) :: inst, outfile, name
    integer :: model_id, np, samples, nout_true
    real(kind=8), allocatable :: lo(:), up(:), start1(:), start2(:), xref(:), t(:), y(:)
    integer, allocatable :: true_out(:)

    ! Parameters
    character(len=8) :: bmode, mmode, pscale
    integer :: noutliers, ntrials, start_id, maxit, project, first
    real(kind=8) :: maxtime, tstart
    real(kind=8) :: cad
    real(kind=8) :: bigM, tau, delta, sigmin, gamma, epsilon, alpha, pert, seed

    ! Algencan
    logical :: checkder
    integer :: hnnzmax, inform, jcnnzmax, m, n, nvparam
    real(kind=8) :: cnorm, efacc, efstain, eoacc, eostain, epsfeas, epsopt, f, nlpsupn, snorm
    character(len=80) :: specfnm, outputfnm, vparam(10)
    logical :: coded(11)
    real(kind=8), pointer :: l(:), u(:), x(:), lambda(:)
    logical,      pointer :: equatn(:), linear(:)

    ! Algorithm state shared with the Algencan callbacks
    real(kind=8), allocatable :: xk(:), xtrial(:), faux(:), indices(:), xinit(:), x0(:)
    real(kind=8), allocatable :: grad(:,:), hess(:,:,:), bqn(:,:)
    ! Adaptive bound: unscaled matrices and their largest eigenvalues
    real(kind=8), allocatable :: bfull(:,:,:), lmaxv(:)
    integer, allocatable :: Idelta(:)
    real(kind=8) :: sigma, fxk, fxtrial

    ! Lapack workspace
    real(kind=8), allocatable :: work(:), eigval(:), eigvec(:,:)
    integer :: lwork, info

    ! Results
    integer :: itrial, iterations, n_eval, status, subfail, ncorrect, q, i, iu
    ! Size of the active set I_delta(x^k): sum and maximum over the iterations
    integer :: msum, mmax
    real(kind=8) :: fovo, start, finish

    call read_arguments()
    call read_instance()

    np = size(xref)
    n = np + 1
    q = samples - noutliers

    allocate(l(n), u(n), x(n), xk(np), xtrial(np), xinit(np), x0(np), faux(samples), &
             indices(samples), Idelta(samples), bqn(np,np), eigval(np), eigvec(np,np))
    lwork = max(1, 3*np)
    allocate(work(lwork))

    l(1:np) = lo; u(1:np) = up
    l(n) = -1.0d+20; u(n) = 0.0d0

    if (start_id == 2) then
        xinit = start2
    else
        xinit = start1
    end if

    ! Algencan settings (as in the original implementation)
    coded(1:6)  = .true.
    coded(7:11) = .false.
    jcnnzmax = 10000 + 10*n*samples
    hnnzmax  = 10000 + 10*n*n*samples
    checkder = .false.
    epsfeas  = 1.0d-08
    epsopt   = 1.0d-08
    efstain  = sqrt(epsfeas)
    eostain  = epsopt ** 1.5d0
    efacc    = sqrt(epsfeas)
    eoacc    = sqrt(epsopt)
    outputfnm = ''
    specfnm   = ''
    nvparam   = 1
    vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 0'

    open(newunit=iu, file=trim(outfile), position='append', action='write')

    do itrial = 1, ntrials
        do i = 1, np
            if (pscale == 'rel') then
                x0(i) = xinit(i) + (2.0d0*drand(seed) - 1.0d0) * pert * abs(xinit(i))
            else
                x0(i) = xinit(i) + (2.0d0*drand(seed) - 1.0d0) * pert * max(1.0d0, abs(xinit(i)))
            end if
        end do
        if (project == 1) x0 = max(lo, min(up, x0))
        ! Trials before `first` only advance the random sequence
        if (itrial < first) cycle

        call cpu_time(start)
        tstart = start
        call ovo_algorithm(x0, fovo, iterations, n_eval, status, subfail)
        call cpu_time(finish)

        ncorrect = count_correct()

        write(iu,'(A,",",I0,",",A,",",ES10.3,",",A,",",ES10.3,",",ES10.3,",",ES10.3,",",I0,",",I0,",",I0,",",I0,",",I0,",",ES12.5,",",ES16.9,",",I0,",",I0,",",F10.2,",",I0)') &
            trim(name), noutliers, trim(bmode), bigM, trim(mmode), tau, delta, sigmin, start_id, itrial, &
            status, subfail, iterations, finish - start, fovo, n_eval, ncorrect, &
            dble(msum) / dble(max(1, iterations)), mmax
        flush(iu)
    end do

    close(iu)

contains

    !==========================================================================
    ! Algorithm 2.1 from the starting point x0
    !==========================================================================
    subroutine ovo_algorithm(x0, fovo, iter, n_eval, status, subfail)
        real(kind=8), intent(in)  :: x0(np)
        real(kind=8), intent(out) :: fovo
        integer,      intent(out) :: iter, n_eval, status, subfail

        integer, parameter :: max_iter_sub = 100
        integer :: iter_sub, i, ip
        real(kind=8) :: theta, s(np), yv(np), ys, bs(np)
        type(hd) :: r

        xk = x0
        iter = 0
        msum = 0
        mmax = 0
        status = 1
        subfail = 0
        if (bmode == 'QN') then
            bqn = 0.0d0
            do i = 1, np
                bqn(i,i) = 1.0d0
            end do
        end if

        call eval_all(xk, faux, indices)
        fxk = faux(q)
        n_eval = 1
        call mount_Idelta()

        do
            iter = iter + 1
            msum = msum + m
            mmax = max(mmax, m)

            allocate(equatn(m), linear(m), lambda(m), grad(m,np), hess(m,np,np), &
                     bfull(m,np,np), lmaxv(m))
            equatn = .false.
            linear = .false.
            lambda = 0.0d0

            ! Gradients and curvature matrices of the active components (Step 2)
            hd_order = 2
            do i = 1, m
                r = residual(model_id, t(Idelta(i)), y(Idelta(i)), xk, np)
                grad(i,:) = r%v * r%g(1:np)
                call curvature(r, hess(i,:,:), lmaxv(i))
            end do
            if (cad > 0.0d0) bfull = hess

            sigma = sigmin
            iter_sub = 1
            x(1:np) = xk
            x(n) = 0.0d0

            do
                ! Adaptive bound ||B_{k,j,i}|| <= cad * sigma_{k,j}: rescale the
                ! unscaled matrices for the current regularization parameter
                if (cad > 0.0d0) then
                    do i = 1, m
                        hess(i,:,:) = min(1.0d0, cad * sigma / max(lmaxv(i), tiny(1.0d0))) &
                                      * bfull(i,:,:)
                    end do
                end if
                call algencan(myevalf, myevalg, myevalh, myevalc, myevaljac, myevalhc,  &
                    myevalfc, myevalgjac, myevalgjacp, myevalhl, myevalhlp, jcnnzmax,    &
                    hnnzmax, epsfeas, epsopt, efstain, eostain, efacc, eoacc, outputfnm, &
                    specfnm, nvparam, vparam, n, x, l, u, m, lambda, equatn, linear,     &
                    coded, checkder, f, cnorm, snorm, nlpsupn, inform)

                xtrial = x(1:np)
                call eval_all(xtrial, faux, indices)
                fxtrial = faux(q)
                n_eval = n_eval + 1

                ! Sufficient decrease (Step 3)
                if (fxtrial <= fxk - alpha * sum((xtrial - xk)**2)) exit
                if (elapsed() > maxtime) exit
                if (iter_sub >= max_iter_sub) then
                    subfail = subfail + 1
                    exit
                end if
                sigma = gamma * sigma
                iter_sub = iter_sub + 1
            end do

            ! Stopping quantity (Step 5)
            theta = 0.0d0
            do i = 1, m
                theta = max(theta, norm2(matmul(hess(i,:,:), xtrial - xk) + sigma * (xtrial - xk)))
            end do

            ! Shared BFGS update, with y built from the component that attains
            ! the order value at x^k (as in Andreani et al., 2006)
            if (bmode == 'QN') then
                hd_order = 2
                s = xtrial - xk
                ip = i_p_at_xk()
                r = residual(model_id, t(ip), y(ip), xtrial, np)
                yv = r%v * r%g(1:np)
                r = residual(model_id, t(ip), y(ip), xk, np)
                yv = yv - r%v * r%g(1:np)
                ys = dot_product(yv, s)
                if (ys > 0.0d0) then
                    bs = matmul(bqn, s)
                    do i = 1, np
                        bqn(:,i) = bqn(:,i) + yv * yv(i) / ys - bs * bs(i) / dot_product(s, bs)
                    end do
                end if
            end if

            deallocate(lambda, equatn, linear, grad, hess, bfull, lmaxv)
            fxk = fxtrial
            xk = xtrial

            if (theta <= epsilon) then
                status = 0
                exit
            end if
            if (iter >= maxit) exit
            if (elapsed() > maxtime) then
                status = 2
                exit
            end if

            call mount_Idelta()
        end do

        fovo = fxk

    end subroutine ovo_algorithm

    real(kind=8) function elapsed()
        real(kind=8) :: now
        call cpu_time(now)
        elapsed = now - tstart
    end function elapsed

    !==========================================================================
    ! Curvature matrix of an active component, with the norm control
    !==========================================================================
    subroutine curvature(r, B, lmax)
        type(hd),     intent(in)  :: r
        real(kind=8), intent(out) :: B(np,np), lmax
        real(kind=8) :: shift
        integer :: j
        character(len=1) :: jobz

        select case (trim(bmode))
        case ('0')
            B = 0.0d0
            lmax = 0.0d0
            return
        case ('GN')
            do j = 1, np
                B(:,j) = r%g(1:np) * r%g(j)
            end do
        case ('H')
            do j = 1, np
                B(:,j) = r%g(1:np) * r%g(j) + r%v * r%h(1:np,j)
            end do
        case ('QN')
            B = bqn
        case default
            write(*,*) 'Unknown curvature family ', trim(bmode)
            stop
        end select

        ! Eigenvalues (and eigenvectors if they are needed for clipping)
        jobz = 'N'
        if (mmode == 'clip') jobz = 'V'
        eigvec = B
        call dsyev(jobz, 'U', np, eigvec, np, eigval, work, lwork, info)

        ! Shift to positive definiteness, as in (23) (Hessians only)
        shift = 0.0d0
        if (bmode == 'H') shift = max(0.0d0, -minval(eigval) + tau)
        do j = 1, np
            B(j,j) = B(j,j) + shift
        end do
        eigval = eigval + shift
        lmax = maxval(eigval)

        ! With the adaptive bound the matrix is rescaled later, for each sigma
        if (cad > 0.0d0) return

        ! Norm control ||B|| <= M
        if (bigM > 0.0d0 .and. lmax > bigM) then
            if (mmode == 'clip') then
                eigval = min(eigval, bigM)
                B = 0.0d0
                do j = 1, np
                    B = B + eigval(j) * spread(eigvec(:,j), 2, np) * spread(eigvec(:,j), 1, np)
                end do
            else
                B = (bigM / lmax) * B
            end if
        end if

    end subroutine curvature

    !==========================================================================
    ! Values f_i(x) of all components, sorted increasingly (one evaluation)
    !==========================================================================
    subroutine eval_all(xv, fv, idx)
        real(kind=8), intent(in)  :: xv(np)
        real(kind=8), intent(out) :: fv(samples), idx(samples)
        type(hd) :: r
        integer :: i

        hd_order = 0
        do i = 1, samples
            r = residual(model_id, t(i), y(i), xv, np)
            fv(i) = 0.5d0 * r%v**2
            idx(i) = dble(i)
        end do
        hd_order = 2
        call DSORT(fv, idx, samples, 2)
    end subroutine eval_all

    integer function i_p_at_xk()
        ! Index attaining the order value at x^k (indices is sorted at x^k
        ! before the trial points are evaluated, so recompute it here)
        real(kind=8) :: fv(samples), idx(samples)
        call eval_all(xk, fv, idx)
        i_p_at_xk = int(idx(q))
    end function i_p_at_xk

    !==========================================================================
    ! Active set I_delta(x^k) from the sorted values in faux/indices
    !==========================================================================
    subroutine mount_Idelta()
        integer :: i
        real(kind=8) :: fq

        fq = faux(q)
        m = 0
        do i = 1, samples
            if (abs(fq - faux(i)) <= delta) then
                m = m + 1
                Idelta(m) = int(indices(i))
            end if
        end do
    end subroutine mount_Idelta

    !==========================================================================
    ! Number of true outliers among the noutliers discarded components at x^k
    !==========================================================================
    integer function count_correct()
        integer :: i, j
        count_correct = 0
        call eval_all(xk, faux, indices)
        do i = samples - noutliers + 1, samples
            do j = 1, nout_true
                if (int(indices(i)) == true_out(j)) count_correct = count_correct + 1
            end do
        end do
    end function count_correct

    !==========================================================================
    ! Input
    !==========================================================================
    subroutine read_arguments()
        character(len=256) :: arg, key, val
        integer :: k, p

        inst = ''; outfile = 'results.csv'
        noutliers = -1; bmode = '0'; mmode = 'scale'; pscale = 'abs'
        bigM = 1.0d0; tau = 1.0d-8; delta = 1.0d-1; sigmin = 1.0d-1; gamma = 5.0d0
        epsilon = 1.0d-4; alpha = 1.0d-8; pert = 0.5d0; seed = 123456.0d0
        ntrials = 100; start_id = 1; maxit = 10000; project = 1; maxtime = 1800.0d0; first = 1; cad = 0.0d0

        do k = 1, command_argument_count()
            call get_command_argument(k, arg)
            p = index(arg, '=')
            if (p == 0) cycle
            key = arg(1:p-1)
            val = arg(p+1:)
            select case (trim(key))
            case ('inst');    inst = val
            case ('out');     outfile = val
            case ('o');       read(val,*) noutliers
            case ('B');       bmode = val
            case ('M');       read(val,*) bigM
            case ('cad');     read(val,*) cad
            case ('mmode');   mmode = val
            case ('tau');     read(val,*) tau
            case ('delta');   read(val,*) delta
            case ('sigmin');  read(val,*) sigmin
            case ('gamma');   read(val,*) gamma
            case ('eps');     read(val,*) epsilon
            case ('alpha');   read(val,*) alpha
            case ('ntrials'); read(val,*) ntrials
            case ('seed');    read(val,*) seed
            case ('start');   read(val,*) start_id
            case ('pert');    read(val,*) pert
            case ('pscale');  pscale = val
            case ('project'); read(val,*) project
            case ('maxit');   read(val,*) maxit
            case ('maxtime'); read(val,*) maxtime
            case ('first');   read(val,*) first
            case default
                write(*,*) 'Unknown argument ', trim(key)
                stop
            end select
        end do
        if (len_trim(inst) == 0) then
            write(*,*) 'Missing inst=FILE'
            stop
        end if
    end subroutine read_arguments

    subroutine read_instance()
        integer :: iu, i, nn

        open(newunit=iu, file=trim(inst), status='old', action='read')
        read(iu,'(A)') name
        read(iu,*) model_id, nn, samples, nout_true
        allocate(lo(nn), up(nn), start1(nn), start2(nn), xref(nn), t(samples), y(samples), &
                 true_out(nout_true))
        read(iu,*) lo
        read(iu,*) up
        read(iu,*) start1
        read(iu,*) start2
        read(iu,*) xref
        do i = 1, samples
            read(iu,*) t(i), y(i)
        end do
        do i = 1, nout_true
            read(iu,*) true_out(i)
        end do
        close(iu)
        if (noutliers < 0) noutliers = nout_true
    end subroutine read_instance

    !==========================================================================
    ! Schrage's portable random number generator (as in the original code)
    !==========================================================================
    function drand(ix)
        real(kind=8) :: drand
        real(kind=8), intent(inout) :: ix
        real(kind=8) :: a, p, b15, b16, xhi, xalo, leftlo, fhi, k
        data a/16807.d0/, b15/32768.d0/, b16/65536.d0/, p/2147483647.d0/

        xhi = ix/b16
        xhi = xhi - dmod(xhi,1.d0)
        xalo = (ix-xhi*b16)*a
        leftlo = xalo/b16
        leftlo = leftlo - dmod(leftlo,1.d0)
        fhi = xhi*a + leftlo
        k = fhi/b15
        k = k - dmod(k,1.d0)
        ix = (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
        if (ix < 0) ix = ix + p
        drand = ix*4.656612875d-10
    end function drand

    !==========================================================================
    ! Algencan callbacks for the subproblem
    !   minimize w  s.t.  grad_i^T d + (d^T B_i d + sigma ||d||^2)/2 <= w,
    !   i in I_delta(x^k), with d = x - x^k and x in the box
    !==========================================================================
    subroutine myevalf(n, x, f, flag)
        integer, intent(in) :: n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: f
        flag = 0
        f = x(n)
    end subroutine myevalf

    subroutine myevalg(n, x, g, flag)
        integer, intent(in) :: n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n)
        flag = 0
        g(1:n-1) = 0.0d0
        g(n) = 1.0d0
    end subroutine myevalg

    subroutine myevalh(n, x, hrow, hcol, hval, hnnz, lim, lmem, flag)
        logical, intent(out) :: lmem
        integer, intent(in) :: lim, n
        integer, intent(out) :: flag, hnnz
        integer, intent(out) :: hcol(lim), hrow(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: hval(lim)
        flag = 0
        lmem = .false.
        hnnz = 0
    end subroutine myevalh

    subroutine myevalc(n, x, ind, c, flag)
        integer, intent(in) :: ind, n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: c
        real(kind=8) :: d(n-1)
        flag = 0
        d = x(1:n-1) - xk
        c = dot_product(d, grad(ind,:)) + 0.5d0 * (dot_product(d, matmul(hess(ind,:,:), d)) &
            + sigma * dot_product(d, d)) - x(n)
    end subroutine myevalc

    subroutine myevaljac(n, x, ind, jcvar, jcval, jcnnz, lim, lmem, flag)
        logical, intent(out) :: lmem
        integer, intent(in) :: ind, lim, n
        integer, intent(out) :: flag, jcnnz
        integer, intent(out) :: jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: jcval(lim)
        real(kind=8) :: d(n-1)
        integer :: i
        flag = 0
        lmem = .false.
        jcnnz = n
        if (jcnnz > lim) then
            lmem = .true.
            return
        end if
        d = x(1:n-1) - xk
        jcvar(1:n) = (/(i, i = 1, n)/)
        jcval(1:n-1) = grad(ind,:) + matmul(hess(ind,:,:), d) + sigma * d
        jcval(n) = -1.0d0
    end subroutine myevaljac

    subroutine myevalhc(n, x, ind, hcrow, hccol, hcval, hcnnz, lim, lmem, flag)
        logical, intent(out) :: lmem
        integer, intent(in) :: ind, lim, n
        integer, intent(out) :: flag, hcnnz
        integer, intent(out) :: hccol(lim), hcrow(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: hcval(lim)
        integer :: i, j
        flag = 0
        lmem = .false.
        hcnnz = 0
        do j = 1, n-1
            do i = j, n-1
                hcnnz = hcnnz + 1
                if (hcnnz > lim) then
                    lmem = .true.
                    return
                end if
                hcrow(hcnnz) = i
                hccol(hcnnz) = j
                hcval(hcnnz) = hess(ind,i,j)
                if (i == j) hcval(hcnnz) = hcval(hcnnz) + sigma
            end do
        end do
    end subroutine myevalhc

    subroutine myevalfc(n, x, f, m, c, flag)
        integer, intent(in) :: m, n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: f, c(m)
        flag = -1
    end subroutine myevalfc

    subroutine myevalgjac(n, x, g, m, jcfun, jcvar, jcval, jcnnz, lim, lmem, flag)
        logical, intent(out) :: lmem
        integer, intent(in) :: lim, m, n
        integer, intent(out) :: flag, jcnnz
        integer, intent(out) :: jcfun(lim), jcvar(lim)
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(out) :: g(n), jcval(lim)
        flag = -1
    end subroutine myevalgjac

    subroutine myevalgjacp(n, x, g, m, p, q, work, gotj, flag)
        logical, intent(inout) :: gotj
        integer, intent(in) :: m, n
        integer, intent(out) :: flag
        character, intent(in) :: work
        real(kind=8), intent(in) :: x(n)
        real(kind=8), intent(inout) :: p(m), q(n)
        real(kind=8), intent(out) :: g(n)
        flag = -1
    end subroutine myevalgjacp

    subroutine myevalhl(n, x, m, lambda, sf, sc, hlrow, hlcol, hlval, hlnnz, lim, lmem, flag)
        logical, intent(out) :: lmem
        integer, intent(in) :: lim, m, n
        integer, intent(out) :: flag, hlnnz
        real(kind=8), intent(in) :: sf
        integer, intent(out) :: hlcol(lim), hlrow(lim)
        real(kind=8), intent(in) :: lambda(m), sc(m), x(n)
        real(kind=8), intent(out) :: hlval(lim)
        flag = -1
    end subroutine myevalhl

    subroutine myevalhlp(n, x, m, lambda, sf, sc, p, hp, goth, flag)
        logical, intent(inout) :: goth
        integer, intent(in) :: m, n
        integer, intent(out) :: flag
        real(kind=8), intent(in) :: sf
        real(kind=8), intent(in) :: lambda(m), p(n), sc(m), x(n)
        real(kind=8), intent(out) :: hp(n)
        flag = -1
    end subroutine myevalhlp

end program ovo_bench
