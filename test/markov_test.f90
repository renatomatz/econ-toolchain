program tauchen_test

    use shared_mod, only: wp
    use dist_mod, only: markov_steady_state

    use unit_testing_mod

    implicit none

    real(wp), allocatable, dimension(:,:) :: adj_mat
    real(wp), allocatable, dimension(:) :: s_dist

    real(wp) :: tol
    integer :: max_iter
    integer :: stat

    integer :: n_stt

    type(UnitTester) :: ut

    call ut%init_test("Test 1: Simple chain")

    n_stt = 2

    allocate(adj_mat(n_stt, n_stt), s_dist(n_stt))

    adj_mat(1,:) = [0.2, 0.8]
    adj_mat(2,:) = [0.8, 0.2]

    tol = 0.001
    max_iter = 100

    call markov_steady_state(adj_mat, s_dist, tol, max_iter, stat)

    call ut%int_scalar("1.1: stat 0", 0, stat)
    call ut%real_scalar("1.2: sum to one", 1.0_wp, sum(s_dist), tol)

    call ut%end_test()

    deallocate(adj_mat, s_dist)

end program tauchen_test
