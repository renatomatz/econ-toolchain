program tauchen_test

    use shared_mod, only: wp
    use dist_mod, only: tauchen

    use unit_testing_mod

    implicit none

    integer :: n
    real(wp) :: mu
    real(wp) :: rho
    real(wp) :: sigma
    real(wp) :: m

    real(wp), allocatable, dimension(:) :: z
    real(wp), allocatable, dimension(:,:) :: zprob

    type(UnitTester) :: ut
    integer :: i

    call ut%init_test("Test 1: Simple dist")

    n = 3

    allocate(z(n), zprob(n,n))

    mu = log(0.15)
    rho = 0.8
    sigma = 2
    m = 3

    call tauchen(n, mu, rho, sigma, m, z, zprob)

    call ut%real_vector("1: sum to one", &
                        [(1.0_wp, i=1, n)], &
                        sum(zprob, 2))

    call ut%end_test()

    deallocate(z, zprob)

end program tauchen_test
