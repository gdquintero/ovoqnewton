! Checks the hyper-dual derivatives of every model against central finite
! differences, at the reference parameters of each instance.
program test_derivs
    use hyperdual
    use models
    implicit none
    character(len=64) :: files(12) = [character(len=64) :: 'Poly_ref','Bard_ref','Misra1a_05', &
        'Chwirut2_05','Rat43_05','MGH17_05','Lanczos3_05','Thurber_05','Kirby2_05','Bennett5_05', &
        'ENSO_05','Gauss1_05']
    integer :: k, iu, mid, n, m, nout, i, j, jt
    real(kind=8), allocatable :: lo(:), up(:), s1(:), s2(:), b(:), t(:), y(:), bp(:), bm(:)
    real(kind=8) :: h, eg, eh, errg, errh, fp, fm
    type(hd) :: r, rp, rm
    character(len=64) :: name
    do k = 1, 12
        open(newunit=iu, file='../instances/'//trim(files(k))//'.txt', status='old')
        read(iu,'(A)') name
        read(iu,*) mid, n, m, nout
        allocate(lo(n), up(n), s1(n), s2(n), b(n), t(m), y(m), bp(n), bm(n))
        read(iu,*) lo; read(iu,*) up; read(iu,*) s1; read(iu,*) s2; read(iu,*) b
        do i = 1, m
            read(iu,*) t(i), y(i)
        end do
        close(iu)
        errg = 0; errh = 0
        do jt = 1, m, max(1, m/7)
            hd_order = 2
            r = residual(mid, t(jt), y(jt), b, n)
            do j = 1, n
                h = 1.0d-6 * max(abs(b(j)), 1.0d-12)
                bp = b; bm = b; bp(j) = b(j) + h; bm(j) = b(j) - h
                rp = residual(mid, t(jt), y(jt), bp, n)
                rm = residual(mid, t(jt), y(jt), bm, n)
                eg = abs((rp%v - rm%v)/(2*h) - r%g(j)) / max(1.0d0, abs(r%g(j)))
                errg = max(errg, eg)
                do i = 1, n
                    eh = abs((rp%g(i) - rm%g(i))/(2*h) - r%h(i,j)) / max(1.0d0, abs(r%h(i,j)))
                    errh = max(errh, eh)
                end do
            end do
        end do
        write(*,'(A12," model ",I2,"  max rel. error: gradient ",ES9.2,"  Hessian ",ES9.2)') trim(name), mid, errg, errh
        deallocate(lo, up, s1, s2, b, t, y, bp, bm)
    end do
end program test_derivs
