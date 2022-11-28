module root_finding_mod

    use stdlib_linalg, only: outer_product

    use shared_mod, only: wp, Func00, Func11, Func12
    use utils_mod, only: solve

    implicit none

contains

    subroutine bound_search(f, lb_n, ub_n, lb, ub, &
                            max_iter, stat, &
                            mon_incr)

        procedure(Func00) :: f

        real(wp), intent(in) :: lb_n, ub_n
        real(wp), intent(out) :: lb, ub

        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        logical, optional, intent(in) :: mon_incr

        real(wp) :: lb_scale

        logical :: mask
        logical :: comp

        integer :: i

        if (lb_n >= ub_n) &
            error stop "bound_search: lower bound must be less than upper &
                       &bound."

        ! The error will only occur if one is possitive and the other
        ! is negative or if either is zero.
        if ((lb_n*ub_n)<=0) &
            error stop "bound_search: domain is inconsistent."

        lb = lb_n
        ub = ub_n

        ! create mask of monotonically-increasing ranges
        if (present(mon_incr)) then
            mask = mon_incr
        else
            mask = f(lb) < f(ub)
        end if

        if (lb_n > 0) then
            lb_scale = 2.0
        else
            lb_scale = 0.5
        end if

        stat = 0

        do i=1,max_iter
            ! flip bit of monotonically-decreasing ranges
            comp = (f(lb)<0).eqv.mask
            if (comp) exit
            if (.not.comp) lb = lb/lb_scale
        end do

        if (i > max_iter) stat = 1

        do i=1,max_iter
            ! flip bit of monotonically-decreasing ranges
            comp = (f(ub)>0).eqv.mask
            if (comp) exit
            if (.not.comp) ub = ub*lb_scale
        end do

        if (i > max_iter) stat = 1

    end subroutine bound_search

    subroutine batch_bound_search(f, lb_n, ub_n, lb, ub, &
                                  max_iter, stat, &
                                  mon_incr)

        procedure(Func11) :: f

        real(wp), intent(in) :: lb_n, ub_n
        real(wp), dimension(:), intent(out) :: lb, ub

        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        logical, optional, intent(in) :: mon_incr

        real(wp) :: lb_scale

        logical, dimension(size(lb)) :: mask
        logical, dimension(size(lb)) :: comp

        integer :: i

        if (lb_n >= ub_n) &
            error stop "bound_search: lower bound must be less than upper &
                       &bound."

        if ((size(lb)/=size(ub))) &
            error stop "bound_search: lower and upper bound must have &
                       &the same shape."

        ! The error will only occur if one is possitive and the other
        ! is negative or if either is zero.
        if ((lb_n*ub_n)<=0) &
            error stop "bound_search: domain is inconsistent."

        lb(:) = lb_n
        ub(:) = ub_n

        ! create mask of monotonically-increasing ranges
        if (present(mon_incr)) then
            mask(:) = mon_incr
        else
            mask = f(lb) < f(ub)
        end if

        if (lb_n > 0) then
            lb_scale = 2.0
        else
            lb_scale = 0.5
        end if

        stat = 0

        do i=1,max_iter
            ! flip bit of monotonically-decreasing ranges
            comp = (f(lb)<0).eqv.mask
            if (all(comp)) exit
            where (.not.comp) lb = lb/lb_scale
        end do

        if (i > max_iter) stat = 1

        do i=1,max_iter
            ! flip bit of monotonically-decreasing ranges
            comp = (f(ub)>0).eqv.mask
            if (all(comp)) exit
            where (.not.comp) ub = ub*lb_scale
        end do

        if (i > max_iter) stat = 1

    end subroutine batch_bound_search

    subroutine bisection_search(f, x, lb, ub, &
                                tol, max_iter, stat)

        procedure(Func00) :: f

        real(wp), intent(inout) :: x
        real(wp), intent(inout) :: lb, ub

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        logical :: mask
        real(wp) :: fx

        integer :: i

        if (lb > ub) &
            error stop "bisection_search: all lower bounds must be lower &
                       &than upper bounds."

        stat = 0

        ! create mask of monotonically-increasing ranges
        mask = f(lb) < f(ub)

        do i=1,max_iter
            ! This won't overflow
            x = lb + (ub-lb)/2
            fx = f(x)
            if (abs(fx) < tol) exit

            ! flip bit of monotonically-decreasing ranges
            if ((fx>0).eqv.mask) then
                ub = x
            else
                lb = x
            end if
        end do

        if (i > max_iter) stat = 1

    end subroutine bisection_search

    subroutine batch_bisection_search(f, x, lb, ub, &
                                      tol, max_iter, stat)

        procedure(Func11) :: f

        real(wp), dimension(:), intent(inout) :: x
        real(wp), dimension(size(x)), intent(inout) :: lb, ub

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        logical, dimension(size(x)) :: mask
        real(wp), dimension(size(x)) :: fx

        integer :: i

        if (any(lb > ub)) &
            error stop "bisection_search: all lower bounds must be lower &
                       &than upper bounds."

        stat = 0

        ! create mask of monotonically-increasing ranges
        mask = f(lb) < f(ub)

        do i=1,max_iter
            ! This won't overflow
            x = lb + (ub-lb)/2
            fx = f(x)
            if (maxval(abs(fx)) < tol) exit

            ! flip bit of monotonically-decreasing ranges
            where ((fx>0).eqv.mask)
                ub = x
            elsewhere
                lb = x
            end where
        end do

        if (i > max_iter) stat = 1

    end subroutine batch_bisection_search

    subroutine newton_method(f, df, x, tol, max_iter, stat)

        procedure(Func00) :: f
        procedure(Func00) :: df
        real(wp), intent(inout) :: x

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp) :: fx
        integer :: i

        stat = 0

        do i=1,max_iter
            fx = f(x)
            if (abs(fx) < tol) exit
            x = x - (fx / df(x))
        end do

        if (i > max_iter) stat = 1

    end subroutine newton_method

    subroutine batch_newton_method(f, df, x, tol, max_iter, stat)

        procedure(Func11) :: f
        procedure(Func11) :: df
        real(wp), dimension(:), intent(inout) :: x

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp), dimension(size(x)) :: fx
        integer :: i

        stat = 0

        do i=1,max_iter
            fx = f(x)
            if (maxval(abs(fx)) < tol) exit
            x = x - (fx / df(x))
        end do

        if (i > max_iter) stat = 1

    end subroutine batch_newton_method

    subroutine newton_method_mult(f, J, x, tol, max_iter, stat)

        procedure(Func11) :: f
        procedure(Func12) :: J
        real(wp), dimension(:), intent(inout) :: x

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp), dimension(size(x)) :: sb
        integer :: n, i

        stat = 0

        n = size(x)

        do i=1,max_iter
            sb = -f(x) ! starts as b
            if (maxval(abs(sb)) < tol) exit
            call solve(J(x), sb, n) ! becomes the x
            x = x + sb
        end do

        if (i > max_iter) stat = 1

    end subroutine newton_method_mult

    subroutine secant_method(f, x1, x2, tol, max_iter, stat)

        procedure(Func00) :: f
        real(wp), intent(inout) :: x1
        real(wp), intent(inout) :: x2

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp) :: x_temp, fx1, fx2
        integer :: i

        stat = 0

        fx2 = f(x2)
        if (abs(fx2) < tol) then
            x1 = x2
            return
        end if

        do i=1,max_iter
            fx1 = f(x1)
            if (abs(fx1) < tol) exit

            x_temp = x1
            x1 = x1 - fx1*((x1-x2) / (fx1-fx2))
            x2 = x_temp
            fx2 = fx1
        end do

        if (i > max_iter) stat = 1

    end subroutine secant_method

    subroutine batch_secant_method(f, x1, x2, tol, max_iter, stat)

        procedure(Func11) :: f
        real(wp), dimension(:), intent(inout) :: x1
        real(wp), dimension(size(x1)), intent(inout) :: x2

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp), dimension(size(x1)) :: x_temp, fx1, fx2
        integer :: i

        stat = 0

        fx2 = f(x2)
        if (maxval(abs(fx2)) < tol) then
            x1 = x2
            return
        end if

        do i=1,max_iter
            fx1 = f(x1)
            if (maxval(abs(fx1)) < tol) exit

            x_temp = x1
            x1 = x1 - fx1*((x1-x2) / (fx1-fx2))
            x2 = x_temp
            fx2 = fx1
        end do

        if (i > max_iter) stat = 1

    end subroutine batch_secant_method

    subroutine broyden_method(f, x, B, tol, max_iter, stat)

        procedure(Func11) :: f
        real(wp), dimension(:), intent(inout) :: x
        real(wp), dimension(size(x), size(x)), intent(inout) :: B

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp), dimension(size(x)) :: fx1, fx2, sb

        integer :: n, i

        n = size(x)

        fx1 = f(x)
        do i=1,max_iter
            sb = -fx1 ! starts as b
            if (maxval(abs(sb)) < tol) exit
            call solve(B, sb, n) ! becomes the x
            x = x + sb
            fx2 = fx1
            fx1 = f(x)
            B = B + (outer_product(fx1 - fx2 - matmul(B,sb), sb)/sum(sb**2))
        end do

        if (i > max_iter) stat = 1

    end subroutine broyden_method

end module
