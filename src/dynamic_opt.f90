module dynamic_opt_mod

    use shared_mod, only: wp, inf
    use utils_mod, only: check_rtv_

    abstract interface

        subroutine BellmanOpInterface(r, t_exg, t_end, v_in, v_out, beta, &
                                      n_exg, n_end, n_act)
            import wp

            real(wp), dimension(n_exg,n_end,n_act), intent(in) :: r
        real(wp), dimension(n_exg,n_exg), intent(in) :: t_exg
        integer, dimension(n_end,n_act), intent(in) :: t_end

            real(wp), dimension(n_exg,n_end), intent(in) :: v_in
            real(wp), dimension(n_exg,n_end), intent(out) :: v_out

            real(wp), intent(in) :: beta

            integer, intent(in) :: n_exg
            integer, intent(in) :: n_end
            integer, intent(in) :: n_act
        end subroutine BellmanOpInterface

    end interface

contains

    subroutine bellman_op_loop(r, t_exg, t_end, v_in, v_out, beta, &
                               n_exg, n_end, n_act)

        real(wp), dimension(n_exg,n_end,n_act), intent(in) :: r
        real(wp), dimension(n_exg,n_exg), intent(in) :: t_exg
        integer, dimension(n_end,n_act), intent(in) :: t_end

        real(wp), dimension(n_exg,n_end), intent(in) :: v_in
        real(wp), dimension(n_exg,n_end), intent(out) :: v_out

        real(wp), intent(in) :: beta

        integer, intent(in) :: n_exg
        integer, intent(in) :: n_end
        integer, intent(in) :: n_act

        integer :: xs, ns, a, xsp, nsp
        real(wp) :: max_a, a_eu

        do xs=1,n_exg
            do ns=1,n_end
                ! this removes non-initialized variable warning
                max_a = -inf
                do a=1,n_act
                    nsp = t_end(ns,a)
                    a_eu = r(xs,ns,a)
                    do xsp=1,n_exg
                        a_eu = a_eu + beta*t_exg(xs,xsp)*v_in(xsp,nsp)
                    end do

                    ! update best action
                    max_a = max(max_a, a_eu)
                end do
                v_out(xs,ns) = max_a

            end do
        end do

    end subroutine bellman_op_loop

    function get_eu(r, t_exg, t_end, v, beta, n_exg, n_end, n_act) result(a_eu)

        real(wp), dimension(n_exg,n_end,n_act), intent(in) :: r
        real(wp), dimension(n_exg,n_exg), intent(in) :: t_exg
        integer, dimension(n_end,n_act), intent(in) :: t_end

        real(wp), dimension(n_exg,n_end), intent(in) :: v

        real(wp), intent(in) :: beta

        integer, intent(in) :: n_exg
        integer, intent(in) :: n_end
        integer, intent(in) :: n_act

        real(wp), dimension(n_exg,n_end,n_act) :: a_eu

        real(wp), dimension(n_exg,n_end) :: s_eu

        integer i, j

        s_eu = beta*matmul(t_exg, v)

        a_eu = 0
        do i=1,n_act
            do j=1,n_end
                a_eu(:,j,i) = r(:,j,i) + s_eu(:,t_end(j,i))
            end do
        end do

    end function get_eu

    subroutine bellman_op_vectorized(r, t_exg, t_end, v_in, v_out, beta, &
                                     n_exg, n_end, n_act)

        real(wp), dimension(n_exg,n_end,n_act), intent(in) :: r
        real(wp), dimension(n_exg,n_exg), intent(in) :: t_exg
        integer, dimension(n_end,n_act), intent(in) :: t_end

        real(wp), dimension(n_exg,n_end), intent(in) :: v_in
        real(wp), dimension(n_exg,n_end), intent(out) :: v_out

        real(wp), intent(in) :: beta

        integer, intent(in) :: n_exg
        integer, intent(in) :: n_end
        integer, intent(in) :: n_act

        v_out = maxval(get_eu(r, t_exg, t_end, v_in, &
                              beta, n_exg, n_end, n_act), 3)

    end subroutine bellman_op_vectorized

    subroutine test_vfi(r, t_exg, t_end, v, beta, &
                        n_exg, n_end, n_act, &
                        bellman_op, tol, max_iter, stat)

        real(wp), dimension(n_exg,n_end,n_act), intent(in) :: r
        real(wp), dimension(n_exg,n_exg), intent(in) :: t_exg
        integer, dimension(n_end,n_act), intent(in) :: t_end
        real(wp), dimension(n_exg,n_end), intent(inout) :: v

        real(wp), intent(in) :: beta

        integer, intent(in) :: n_exg, n_end, n_act

        procedure(BellmanOpInterface) :: bellman_op

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter

        integer, intent(out) :: stat

        integer :: i
        real(wp), dimension(n_exg, n_end) :: last

        stat = 0

        do i=1,max_iter
            last = v
            call bellman_op(r, t_exg, t_end, last, v, beta, &
                            n_exg, n_end, n_act)
            if (maxval(abs(v - last)) < tol) exit
        end do

        if (i > max_iter) stat = 1

    end subroutine test_vfi

    subroutine vfi(r, t_exg, t_end, v, beta, &
                          n_exg, n_end, n_act, &
                          tol, max_iter, stat)

        real(wp), dimension(n_exg,n_end,n_act), intent(in) :: r
        real(wp), dimension(n_exg,n_exg), intent(in) :: t_exg
        integer, dimension(n_end,n_act), intent(in) :: t_end
        real(wp), dimension(n_exg,n_end), intent(inout) :: v

        real(wp), intent(in) :: beta

        integer, intent(in) :: n_exg, n_end, n_act

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter

        integer, intent(out) :: stat

        real(wp), dimension(n_exg, n_end) :: last
        real(wp), dimension(n_exg,n_end,n_act) :: a_eu

        integer i

        stat = 0

        do i=1,max_iter
            last = v
            a_eu = get_eu(r, t_exg, t_end, v, &
                          beta, n_exg, n_end, n_act)
            v = maxval(a_eu, 3)

            if (maxval(abs(v - last)) < tol) exit
        end do

        if (i > max_iter) stat = 1

    end subroutine vfi

    function get_policy(r, t_exg, t_end, v, beta) result (pi)

        real(wp), dimension(:,:,:), intent(in) :: r
        real(wp), dimension(:,:), intent(in) :: t_exg
        integer, dimension(:,:), intent(in) :: t_end
        real(wp), dimension(:,:), intent(out) :: v

        real(wp), intent(in) :: beta
        integer, dimension(size(v,1),size(v,2)) :: pi

        integer :: n_exg, n_end, n_act

        n_exg = size(r, 1)
        n_end = size(r, 2)
        n_act = size(r, 3)

        call check_rtv_(r, t_exg, t_end, v, n_exg, n_end, n_act)

        pi = maxloc(get_eu(r, t_exg, t_end, v, beta, n_exg, n_end, n_act), 3)

    end function get_policy

end module dynamic_opt_mod
