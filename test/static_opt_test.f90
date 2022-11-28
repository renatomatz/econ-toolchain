program root_finding_test

    use shared_mod, only: wp, Func00
    use static_opt_mod

    use unit_testing_mod

    implicit none

    real(wp) :: x
    real(wp) :: lb, ub
    real(wp) :: alpha

    real(wp) :: tol

    integer :: max_iter
    integer :: stat

    type(UnitTester) :: ut

    call ut%set_verbose(.true.)

    call ut%init_test("Test 1: y = x^2 - 2x - 4")

    tol = 0.001
    max_iter = 1000

    lb = 0
    ub = 2

    call golden_section_search(quadratic, lb, ub, &
                               tol, max_iter, stat)

    call ut%int_scalar("1.1.1: stat 0", 0, stat)
    call ut%real_scalar("1.1.2: golden section search", 1.0_wp, lb, tol)

    x = 0
    alpha = 0.01

    call steepest_descent(quadratic_df, x, alpha, tol, max_iter, stat)

    call ut%int_scalar("2.1.1: stat 0", 0, stat)
    call ut%real_scalar("2.1.2: golden section search", 1.0_wp, x, tol)

    call ut%end_test()

contains


    function quadratic(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y

        y = x**2 - 2*x - 4
    end function quadratic

    function quadratic_df(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y

        y = 2*x - 2
    end function quadratic_df

end program root_finding_test
