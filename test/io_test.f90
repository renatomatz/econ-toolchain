program io_test

    use shared_mod, only: wp
    use utils_mod, only: write_res, read_res

    use unit_testing_mod

    implicit none

    real(wp), allocatable, dimension(:) :: x_in
    real(wp), allocatable, dimension(:) :: x_out
    real(wp) :: y_in
    real(wp) :: y_out
    integer :: n

    integer :: i

    type(UnitTester) :: ut

    call ut%init_test("Test 1")

    n = 5
    allocate(x_in(n), x_out(n))

    x_in = [(1.0_wp/i, i=1, n)]
    y_in = 2.0_wp
    call write_res(x_in, y_in, "test_file")
    call read_res(x_out, y_out, "test_file")

    call ut%real_vector("1.1: x", x_in, x_out)
    call ut%real_scalar("1.2: y", y_in, y_out)

    call ut%end_test()

end program io_test
