module static_opt_mod

    use shared_mod, only: wp, golden_ratio, Func00, Func11

    implicit none

contains

    subroutine golden_section_search(f, lb, ub, &
                                     tol, max_iter, stat)

        procedure(Func00) :: f

        real(wp), intent(inout) :: lb, ub

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp) :: x1, x2, fx1, fx2

        integer :: i

        stat = 0

        associate( tau => golden_ratio )

        x1 = lb + (1-tau)*(ub-lb)
        x2 = lb + tau*(ub-lb)

        fx1 = f(x1)
        fx2 = f(x2)

        do i=1,max_iter
            if ((ub-lb) < tol) exit

            if (fx1 > fx2) then
                lb = x1
                x1 = x2
                fx1 = fx2
                x2 = lb + tau*(ub-lb)
                fx2 = f(x2)
            else
                ub = x2
                x2 = x1
                fx2 = fx1
                x1 = lb + (1-tau)*(ub-lb)
                fx1 = f(x1)
            end if
        end do

        if (fx1 > fx2) then
            lb = x2
        else
            lb = x1
        end if
        ub = lb

        if (i > max_iter) stat = 1

        end associate

    end subroutine golden_section_search

    subroutine steepest_descent(df, x, alpha, tol, max_iter, stat)

        procedure(Func00) :: df
        real(wp), intent(inout) :: x
        real(wp), intent(in) :: alpha

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp) :: dfx
        integer :: i

        stat = 0

        do i=1,max_iter
            dfx = df(x)
            if (abs(dfx) < tol) exit
            x = x - alpha*dfx
        end do

        if (i > max_iter) stat = 1

    end subroutine steepest_descent

    subroutine steepest_descent_mult(df, x, alpha, tol, max_iter, stat)

        procedure(Func11) :: df
        real(wp), dimension(:), intent(inout) :: x
        real(wp), intent(in) :: alpha

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp), dimension(size(x)) :: dfx
        integer :: i

        stat = 0

        do i=1,max_iter
            dfx = df(x)
            if (maxval(abs(dfx)) < tol) exit
            x = x - alpha*dfx
        end do

        if (i > max_iter) stat = 1

    end subroutine steepest_descent_mult

end module static_opt_mod
