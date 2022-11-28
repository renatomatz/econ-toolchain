program root_finding_test

    use shared_mod, only: wp
    use root_finding_mod

    use unit_testing_mod

    implicit none

    real(wp), dimension(:), allocatable :: x
    real(wp), dimension(:,:), allocatable :: B

    real(wp) :: tol

    integer :: max_iter
    integer :: stat

    real(wp), dimension(:), allocatable :: wanted

    type(UnitTester) :: ut

    call ut%set_verbose(.true.)

    tol = 0.001
    max_iter = 1000

    call ut%init_test("Test 1")

    allocate(x(2), B(2,2), wanted(2))

    wanted = [0.0_wp, 1.0_wp]

    x = [1.0_wp, 2.0_wp]

    call newton_method_mult(test1_f, test1_J, x, tol, max_iter, stat)

    call ut%int_scalar("1.1.1: stat 0", 0, stat)
    call ut%real_vector("1.1.2: newton method search", wanted, x, tol)

    x = [1.0_wp, 2.0_wp]
    B = test1_J(x)

    call broyden_method(test1_f, x, B, tol, max_iter, stat)

    call ut%int_scalar("1.2.1: stat 0", 0, stat)
    call ut%real_vector("1.2.2: broyden method", wanted, x, tol)

    call ut%end_test()

contains

    function test1_f(x) result(y)
        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(size(x)) :: y

        y(1) = x(1) + 2*x(2) - 2
        y(2) = x(1)**2 + 4*x(2)**2 - 4
    end function test1_f

    function test1_J(x) result(y)
        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(size(x),size(x)) :: y

        y(1,1) = 1
        y(1,2) = 2
        y(2,1) = 2*x(1)
        y(2,2) = 8*x(2)
    end function test1_J

end program root_finding_test
