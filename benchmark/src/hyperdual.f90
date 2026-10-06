!==============================================================================
! Forward-mode automatic differentiation of second order ("hyper-dual"
! numbers). A variable of type hd carries a value, its gradient and its
! Hessian with respect to the first hd_n parameters, so that each model is
! written once and its exact first and second derivatives follow from the
! overloaded arithmetic below.
!
! hd_order = 0 propagates values only (cheap function evaluations);
! hd_order = 2 propagates values, gradients and Hessians.
!==============================================================================
module hyperdual

    implicit none

    integer, parameter :: hd_nmax = 12

    integer :: hd_n = 0
    integer :: hd_order = 2

    type :: hd
        real(kind=8) :: v = 0.0d0
        real(kind=8) :: g(hd_nmax) = 0.0d0
        real(kind=8) :: h(hd_nmax,hd_nmax) = 0.0d0
    end type hd

    interface operator(+)
        module procedure add_hh, add_hr, add_rh
    end interface

    interface operator(-)
        module procedure sub_hh, sub_hr, sub_rh, neg_h
    end interface

    interface operator(*)
        module procedure mul_hh, mul_hr, mul_rh
    end interface

    interface operator(/)
        module procedure div_hh, div_hr, div_rh
    end interface

    interface operator(**)
        module procedure pow_hr, pow_hi, pow_hh
    end interface

    interface exp
        module procedure exp_h
    end interface

    interface log
        module procedure log_h
    end interface

    interface sin
        module procedure sin_h
    end interface

    interface cos
        module procedure cos_h
    end interface

    private :: add_hh, add_hr, add_rh, sub_hh, sub_hr, sub_rh, neg_h,   &
               mul_hh, mul_hr, mul_rh, div_hh, div_hr, div_rh,          &
               pow_hr, pow_hi, pow_hh, exp_h, log_h, sin_h, cos_h, chain

contains

    !--------------------------------------------------------------------------
    ! Constructors
    !--------------------------------------------------------------------------
    function hd_var(x, i) result(c)
        real(kind=8), intent(in) :: x
        integer,      intent(in) :: i
        type(hd) :: c

        c%v = x
        if (hd_order > 0) c%g(i) = 1.0d0
    end function hd_var

    function hd_const(x) result(c)
        real(kind=8), intent(in) :: x
        type(hd) :: c

        c%v = x
    end function hd_const

    !--------------------------------------------------------------------------
    ! Univariate chain rule: c = phi(a), with phi'(a%v) = d1, phi''(a%v) = d2
    !--------------------------------------------------------------------------
    function chain(a, v, d1, d2) result(c)
        type(hd),     intent(in) :: a
        real(kind=8), intent(in) :: v, d1, d2
        type(hd) :: c
        integer :: j, n

        c%v = v
        if (hd_order == 0) return
        n = hd_n
        c%g(1:n) = d1 * a%g(1:n)
        do j = 1, n
            c%h(1:n,j) = d1 * a%h(1:n,j) + d2 * a%g(1:n) * a%g(j)
        end do
    end function chain

    !--------------------------------------------------------------------------
    ! Addition and subtraction
    !--------------------------------------------------------------------------
    function add_hh(a, b) result(c)
        type(hd), intent(in) :: a, b
        type(hd) :: c
        integer :: n

        c%v = a%v + b%v
        if (hd_order == 0) return
        n = hd_n
        c%g(1:n) = a%g(1:n) + b%g(1:n)
        c%h(1:n,1:n) = a%h(1:n,1:n) + b%h(1:n,1:n)
    end function add_hh

    function add_hr(a, r) result(c)
        type(hd),     intent(in) :: a
        real(kind=8), intent(in) :: r
        type(hd) :: c

        c = a
        c%v = a%v + r
    end function add_hr

    function add_rh(r, a) result(c)
        real(kind=8), intent(in) :: r
        type(hd),     intent(in) :: a
        type(hd) :: c

        c = a
        c%v = a%v + r
    end function add_rh

    function sub_hh(a, b) result(c)
        type(hd), intent(in) :: a, b
        type(hd) :: c
        integer :: n

        c%v = a%v - b%v
        if (hd_order == 0) return
        n = hd_n
        c%g(1:n) = a%g(1:n) - b%g(1:n)
        c%h(1:n,1:n) = a%h(1:n,1:n) - b%h(1:n,1:n)
    end function sub_hh

    function sub_hr(a, r) result(c)
        type(hd),     intent(in) :: a
        real(kind=8), intent(in) :: r
        type(hd) :: c

        c = a
        c%v = a%v - r
    end function sub_hr

    function sub_rh(r, a) result(c)
        real(kind=8), intent(in) :: r
        type(hd),     intent(in) :: a
        type(hd) :: c

        c = neg_h(a)
        c%v = r - a%v
    end function sub_rh

    function neg_h(a) result(c)
        type(hd), intent(in) :: a
        type(hd) :: c
        integer :: n

        c%v = -a%v
        if (hd_order == 0) return
        n = hd_n
        c%g(1:n) = -a%g(1:n)
        c%h(1:n,1:n) = -a%h(1:n,1:n)
    end function neg_h

    !--------------------------------------------------------------------------
    ! Multiplication and division
    !--------------------------------------------------------------------------
    function mul_hh(a, b) result(c)
        type(hd), intent(in) :: a, b
        type(hd) :: c
        integer :: j, n

        c%v = a%v * b%v
        if (hd_order == 0) return
        n = hd_n
        c%g(1:n) = a%v * b%g(1:n) + b%v * a%g(1:n)
        do j = 1, n
            c%h(1:n,j) = a%v * b%h(1:n,j) + b%v * a%h(1:n,j) &
                       + a%g(1:n) * b%g(j) + b%g(1:n) * a%g(j)
        end do
    end function mul_hh

    function mul_hr(a, r) result(c)
        type(hd),     intent(in) :: a
        real(kind=8), intent(in) :: r
        type(hd) :: c
        integer :: n

        c%v = a%v * r
        if (hd_order == 0) return
        n = hd_n
        c%g(1:n) = r * a%g(1:n)
        c%h(1:n,1:n) = r * a%h(1:n,1:n)
    end function mul_hr

    function mul_rh(r, a) result(c)
        real(kind=8), intent(in) :: r
        type(hd),     intent(in) :: a
        type(hd) :: c

        c = mul_hr(a, r)
    end function mul_rh

    function div_hh(a, b) result(c)
        type(hd), intent(in) :: a, b
        type(hd) :: c

        c = mul_hh(a, chain(b, 1.0d0/b%v, -1.0d0/b%v**2, 2.0d0/b%v**3))
    end function div_hh

    function div_hr(a, r) result(c)
        type(hd),     intent(in) :: a
        real(kind=8), intent(in) :: r
        type(hd) :: c

        c = mul_hr(a, 1.0d0/r)
    end function div_hr

    function div_rh(r, a) result(c)
        real(kind=8), intent(in) :: r
        type(hd),     intent(in) :: a
        type(hd) :: c

        c = mul_hr(chain(a, 1.0d0/a%v, -1.0d0/a%v**2, 2.0d0/a%v**3), r)
    end function div_rh

    !--------------------------------------------------------------------------
    ! Powers
    !--------------------------------------------------------------------------
    function pow_hr(a, p) result(c)
        type(hd),     intent(in) :: a
        real(kind=8), intent(in) :: p
        type(hd) :: c

        c = chain(a, a%v**p, p * a%v**(p-1.0d0), p * (p-1.0d0) * a%v**(p-2.0d0))
    end function pow_hr

    function pow_hi(a, k) result(c)
        type(hd), intent(in) :: a
        integer,  intent(in) :: k
        type(hd) :: c
        real(kind=8) :: p

        p = dble(k)
        if (k >= 2) then
            c = chain(a, a%v**k, p * a%v**(k-1), p * (p-1.0d0) * a%v**(k-2))
        else
            c = pow_hr(a, p)
        end if
    end function pow_hi

    function pow_hh(a, b) result(c)
        type(hd), intent(in) :: a, b
        type(hd) :: c

        c = exp_h(mul_hh(b, log_h(a)))
    end function pow_hh

    !--------------------------------------------------------------------------
    ! Elementary functions
    !--------------------------------------------------------------------------
    function exp_h(a) result(c)
        type(hd), intent(in) :: a
        type(hd) :: c
        real(kind=8) :: e

        e = exp(a%v)
        c = chain(a, e, e, e)
    end function exp_h

    function log_h(a) result(c)
        type(hd), intent(in) :: a
        type(hd) :: c

        c = chain(a, log(a%v), 1.0d0/a%v, -1.0d0/a%v**2)
    end function log_h

    function sin_h(a) result(c)
        type(hd), intent(in) :: a
        type(hd) :: c

        c = chain(a, sin(a%v), cos(a%v), -sin(a%v))
    end function sin_h

    function cos_h(a) result(c)
        type(hd), intent(in) :: a
        type(hd) :: c

        c = chain(a, cos(a%v), -sin(a%v), -cos(a%v))
    end function cos_h

end module hyperdual
