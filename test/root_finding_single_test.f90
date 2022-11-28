program root_finding_test

    use shared_mod, only: wp, Func00
    use utils_mod, only: fixed_point
    use root_finding_mod

    use unit_testing_mod

    implicit none

    real(wp) :: lb_n, ub_n
    real(wp) :: x, x2
    real(wp) :: lb, ub

    real(wp) :: tol

    integer :: max_iter
    integer :: stat

    type(UnitTester) :: ut

    call ut%set_verbose(.true.)

    call ut%init_test("Test 1: Bound Search")

    tol = 0.01
    max_iter = 100

    ! Normal linear
    lb_n = 6
    ub_n = 14

    call bound_search(linear, lb_n, ub_n, lb, ub, &
                      max_iter, stat)

    call ut%int_scalar("1.1.1: stat 0", 0, stat)
    call ut%real_scalar("1.1.2: lb no change", lb_n, lb)
    call ut%real_scalar("1.1.3: ub no change", ub_n, ub)

    ! Lower bound search linear
    lb_n = 12
    ub_n = 14

    call bound_search(linear, lb_n, ub_n, lb, ub, &
                      max_iter, stat)

    call ut%int_scalar("1.2.1: stat 0", 0, stat)
    call ut%real_scalar("1.2.2: lb correction", lb_n/2, lb)

    ! Upper bound search linear
    lb_n = 6
    ub_n = 7

    call bound_search(linear, lb_n, ub_n, lb, ub, &
                      max_iter, stat)

    call ut%int_scalar("1.3.1: stat 0", 0, stat)
    call ut%real_scalar("1.3.2: ub correction", ub_n*2, ub)

    call ut%end_test()

    tol = 0.001
    max_iter = 1000

    call ut%init_test("Test 2: y = x - 8")

    lb = 4
    ub = 14

    call bisection_search(linear, x, lb, ub, &
                          tol, max_iter, stat)

    call ut%int_scalar("2.1: stat 0", 0, stat)
    call ut%real_scalar("2.2: bisection search", 8.0_wp, x, tol)

    call ut%end_test()

    call ut%init_test("Test 3: y = x^2 - x - 2")

    lb = 0
    ub = 4

    call bisection_search(quadratic, x, lb, ub, &
                          tol, max_iter, stat)

    call ut%int_scalar("3.1.1: stat 0", 0, stat)
    call ut%real_scalar("3.1.2: bisection search", 2.0_wp, x, tol)

    x = 1

    call fixed_point(quadratic_fp, x, tol, max_iter, stat)

    call ut%int_scalar("3.2.1: stat 0", 0, stat)
    call ut%real_scalar("3.2.2: fixed point search", 2.0_wp, x, tol)

    x = 1

    call newton_method(quadratic, quadratic_df, x, tol, max_iter, stat)

    call ut%int_scalar("3.3.1: stat 0", 0, stat)
    call ut%real_scalar("3.3.2: newton method search", 2.0_wp, x, tol)

    x = 1
    x2 = 3

    call secant_method(quadratic, x, x2, tol, max_iter, stat)

    call ut%int_scalar("3.3.1: stat 0", 0, stat)
    call ut%real_scalar("3.3.2: newton method search", 2.0_wp, x, tol)

    call ut%end_test()

contains

    function linear(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y

        y = x - 8
    end function linear

    function quadratic(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y

        y = x**2 - x - 2
    end function quadratic

    function quadratic_fp(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y

        y = sqrt(x+2)
    end function quadratic_fp

    function quadratic_df(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y

        y = 2*x - 1
    end function quadratic_df

end program root_finding_test
