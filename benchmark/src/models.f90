!==============================================================================
! Parametric models y(t,x) used in the experiments. Each model is written once
! in hyper-dual arithmetic, so the residual r_i(x) = y(t_i,x) - y_i comes with
! its exact gradient and Hessian.
!==============================================================================
module models

    use hyperdual

    implicit none

    real(kind=8), parameter :: pi = 3.14159265358979323846d0

contains

    !--------------------------------------------------------------------------
    ! Model value y(t,b) for model identifier id
    !--------------------------------------------------------------------------
    function model_eval(id, t, b) result(y)
        integer,      intent(in) :: id
        real(kind=8), intent(in) :: t
        type(hd),     intent(in) :: b(:)
        type(hd) :: y
        real(kind=8) :: u, v, w

        select case (id)

        case (1)   ! Cubic polynomial (Andreani, Dunder and Martinez)
            y = b(1) + b(2)*t + b(3)*t**2 + b(4)*t**3

        case (2)   ! Bard (More, Garbow and Hillstrom)
            u = t
            v = 16.0d0 - t
            w = min(u, v)
            y = b(1) + u / (v*b(2) + w*b(3))

        case (3)   ! Misra1a (NIST)
            y = b(1) * (1.0d0 - exp(-b(2)*t))

        case (4)   ! Chwirut2 (NIST)
            y = exp(-b(1)*t) / (b(2) + b(3)*t)

        case (5)   ! Rat43 (NIST)
            y = b(1) / ((1.0d0 + exp(b(2) - b(3)*t)) ** (1.0d0 / b(4)))

        case (6)   ! MGH17 = Osborne 1 (NIST)
            y = b(1) + b(2)*exp(-t*b(4)) + b(3)*exp(-t*b(5))

        case (7)   ! Lanczos3 (NIST)
            y = b(1)*exp(-b(2)*t) + b(3)*exp(-b(4)*t) + b(5)*exp(-b(6)*t)

        case (8)   ! Thurber (NIST)
            y = (b(1) + b(2)*t + b(3)*t**2 + b(4)*t**3) / &
                (1.0d0 + b(5)*t + b(6)*t**2 + b(7)*t**3)

        case (9)   ! Kirby2 (NIST)
            y = (b(1) + b(2)*t + b(3)*t**2) / (1.0d0 + b(4)*t + b(5)*t**2)

        case (10)  ! Bennett5 (NIST)
            y = b(1) * (b(2) + t) ** (-1.0d0 / b(3))

        case (11)  ! ENSO (NIST)
            y = b(1) + b(2)*cos(2.0d0*pi*t/12.0d0) + b(3)*sin(2.0d0*pi*t/12.0d0) &
                + b(5)*cos(2.0d0*pi*t/b(4)) + b(6)*sin(2.0d0*pi*t/b(4))        &
                + b(8)*cos(2.0d0*pi*t/b(7)) + b(9)*sin(2.0d0*pi*t/b(7))

        case (12)  ! Gauss1 (NIST)
            y = b(1)*exp(-b(2)*t) + b(3)*exp(-(t - b(4))**2 / b(5)**2) &
                + b(6)*exp(-(t - b(7))**2 / b(8)**2)

        case default
            write(*,*) 'Unknown model identifier ', id
            stop
        end select

    end function model_eval

    !--------------------------------------------------------------------------
    ! Residual r = y(t,x) - yobs, with gradient and Hessian if hd_order = 2
    !--------------------------------------------------------------------------
    function residual(id, t, yobs, x, n) result(r)
        integer,      intent(in) :: id, n
        real(kind=8), intent(in) :: t, yobs, x(n)
        type(hd) :: r
        type(hd) :: b(n)
        integer :: i

        hd_n = n
        do i = 1, n
            b(i) = hd_var(x(i), i)
        end do
        r = model_eval(id, t, b) - yobs

    end function residual

end module models
