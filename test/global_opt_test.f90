program global_opt_test

    use shared_mod, only: wp
    use global_opt_mod, only: rand_global_search, sa_wrapper

    implicit none

    integer, parameter :: n = 3

    real(wp), dimension(n) :: x
    real(wp) :: y

    real(wp), dimension(n) :: lb
    real(wp), dimension(n) :: ub

    integer :: n_evals

    integer :: state_size
    integer, allocatable, dimension(:) :: state

    n_evals = 50

    lb = -5.0_wp
    ub = 10.0_wp

    call random_seed(size=state_size)
    allocate(state(state_size))
    ! for tests only
    state = 42 + this_image()
    ! call random_seed(get=state)
    call random_seed(put=state)
    deallocate(state)

    call rand_global_search(rosenbrock, sa_wrapper, x, n, y, lb, ub, n_evals)

    sync all

    if (this_image() == 1) then
        print *, "Best Input"
        print *, x
        print *, "Best Output"
        print *, y
    end if

contains

    function rosenbrock(x, n) result(y)
        real(wp), dimension(n), intent(in) :: x
        integer, intent(in) :: n

        real(wp) :: y

        real(kind=wp), dimension(size(x)-1) :: x1, x2

        real(wp) :: burn

        burn = 1
        do while (burn > 1e-0)
            call random_number(burn)
        end do

        x1 = x(:n-1)
        x2 = x(2:)

        y = sum(100*(x2-x1**2)**2+(1-x1)**2) + 1
    end function rosenbrock

end program global_opt_test
