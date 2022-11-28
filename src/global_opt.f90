module global_opt_mod

    use shared_mod, only: wp, dp, inf, Func1m0
    use coarray_utils_mod, only: co_min_with_vec, partition_params

    use simulated_annealing_module, only: simulated_annealing_type

    implicit none

    abstract interface

        subroutine LocSolver(f, x, n, y, lb, ub)
            import wp, Func1m0

            procedure(Func1m0) :: f
            real(wp), dimension(n), intent(inout) :: x
            real(wp), intent(out) :: y
            integer, intent(in) :: n
            real(wp), dimension(n), intent(in) :: lb
            real(wp), dimension(n), intent(in) :: ub
        end subroutine LocSolver

    end interface

contains

    subroutine rand_global_search(f, loc_solver, x, n, y, &
                                  lb, ub, n_evals, best_img)
        procedure(Func1m0) :: f
        procedure(LocSolver) :: loc_solver
        real(wp), dimension(n), intent(inout) :: x
        integer, intent(in) :: n
        real(wp), intent(inout) :: y
        real(wp), dimension(n), intent(in) :: lb
        real(wp), dimension(n), intent(in) :: ub
        integer, intent(in) :: n_evals
        integer, optional, intent(out) :: best_img

        real(wp), dimension(n) :: cur_x
        real(wp) :: cur_y

        integer :: i

        y = inf
        do i=1,n_evals

            call random_number(cur_x)
            cur_x = lb + (ub - lb)*cur_x
            call loc_solver(f, cur_x, n, cur_y, lb, ub)

            if (cur_y < y) then
                x = cur_x
                y = cur_y
            end if

        end do

        call co_min_with_vec(x, y, best_img)

    end subroutine rand_global_search

    subroutine cust_global_search(f, loc_solver, x, n, y, &
                                  cust_x, lb, ub, best_img)
        procedure(Func1m0) :: f
        procedure(LocSolver) :: loc_solver
        real(wp), dimension(n), intent(inout) :: x
        integer, intent(in) :: n
        real(wp), intent(inout) :: y
        real(wp), dimension(:,:), intent(in) :: cust_x
        real(wp), dimension(n), intent(in) :: lb
        real(wp), dimension(n), intent(in) :: ub
        integer, optional, intent(out) :: best_img

        real(wp), dimension(n) :: cur_x
        real(wp) :: cur_y

        integer :: i, max_loc_len, loc_len, loc_base

        if (size(cust_x,2) /= n) &
            error stop "custom x parameters should be along dimension 2"

        call partition_params(size(cust_x,1), max_loc_len, loc_len, loc_base)

        y = inf
        do i=loc_base,loc_base+loc_len

            cur_x = cust_x(i,:)
            call loc_solver(f, cur_x, n, cur_y, lb, ub)

            if (cur_y < y) then
                x = cur_x
                y = cur_y
            end if

        end do

        call co_min_with_vec(x, y, best_img)

    end subroutine cust_global_search

    subroutine sa_wrapper(f, x, n, y, lb, ub)
        procedure(Func1m0) :: f
        real(wp), dimension(n), intent(inout) :: x
        integer, intent(in) :: n
        real(wp), intent(out) :: y
        real(wp), dimension(n), intent(in) :: lb
        real(wp), dimension(n), intent(in) :: ub

        real(dp) :: c(n), vm(n), lb_dp(n), ub_dp(n)
        real(dp) :: eps, vms
        integer :: maxevl, ns, nt, neps, iseed1, iseed2, iprint, &
                   step_mode, iunit
        logical :: mx

        real(dp) :: xin(n), xopt(n)
        real(dp) :: t, rt, fopt
        integer :: nacc, nfcnev, ier

        type(simulated_annealing_type) :: sa

        c = 2.0_dp

        ! set input parameters.
        lb_dp = real(lb, dp)
        ub_dp = real(ub, dp)

        ! set input parameters.
        mx = .false.
        eps = 1.0e-6_dp
        iseed1 = 1 + this_image()
        iseed2 = 2 + this_image()
        ns = 20
        nt = 10
        neps = 4
        maxevl = 20
        iprint = 0

        ! others:
        step_mode = 1
        vms = 0.1_dp

        call sa%initialize(f_wrap, n, lb_dp, ub_dp, c, &
                           mx,eps, ns, nt, neps, maxevl, &
                           iprint, iseed1, iseed2, step_mode, vms, iunit)

        xin = real(x, dp)
        rt = 0.5_dp
        t = 5.0_dp
        vm = 1.0_dp

        call sa%optimize(xin, rt, t, vm, xopt, fopt, nacc, nfcnev, ier)

        x = real(xopt, wp)
        y = real(fopt, wp)

    contains

        subroutine f_wrap(self, x, y, istat)
            class(simulated_annealing_type),intent(inout) :: self
            real(dp),dimension(:),intent(in) :: x
            real(dp), intent(out) :: y
            integer, intent(out) :: istat

            y = f(x, n)
            istat = 0
        end subroutine f_wrap

    end subroutine sa_wrapper

end module global_opt_mod
