module dist_mod

    use shared_mod, only: wp

    implicit none

contains

    subroutine tauchen(n, mu, rho, sigma, m, z, zprob)
        !! Finds a Markov chain whose sample paths approximate those of the
        !! AR(1) process:
        !! z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
        !! where eps are normal with stddev sigma.

        integer, intent(in) :: n
            !! number of nodes for Z
        real(wp), intent(in) :: mu
            !! unconditional mean of process
        real(wp), intent(in) :: rho
            !! rho parameter in AR process
        real(wp), intent(in) :: sigma
            !! standard deviation of epsilons
        real(wp), intent(in) :: m
            !! max +- standard deviation

        real(wp), dimension(n), intent(out) :: z
            !! nodes for Z
        real(wp), dimension(n,n), intent(out) :: zprob
            !! transition probabilities

        real(wp) :: a, step
        integer :: i

        z = 0
        zprob = 0
        a = (1-rho)*mu

        z(n) = m * sqrt(sigma**2/(1-rho**2))
        z(1) = -z(n)
        step = (z(n) - z(1)) / (n-1)

        z(2:(n-1)) = z(1) + [(i, i=1, (n-2))]*step
        z = z + mu

        ! these assignments cannot be vectorized because of rho*z, as z is
        ! a vector.
        zprob(:,1) =     normal_cdf((z(1) - a - rho*z + step/2)  / sigma)
        do concurrent(i=2:(n-1))
            zprob(:,i) = normal_cdf((z(i) - a - rho*z + step/2)  / sigma) - &
                         normal_cdf((z(i) - a - rho*z - step/2)  / sigma)
        end do
        zprob(:,n) = 1 - normal_cdf((z(n) - a - rho*z - step/2)  / sigma)
    end subroutine tauchen

    pure elemental function normal_cdf(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y

        y = 0.5_wp * erfc(-x/sqrt(2.0_wp))
    end function normal_cdf

    subroutine markov_steady_state(adj_mat, s_dist, tol, max_iter, stat)

        real(wp), dimension(:,:), intent(in) :: adj_mat
        real(wp), dimension(:), intent(inout) :: s_dist

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp), dimension(size(s_dist)) :: last
        integer :: i, n

        if (size(adj_mat,1) /= size(adj_mat,2)) &
            error stop "markov_steady_state: adj_mat must be square"

        if (size(adj_mat,1) /= size(s_dist,1)) &
            error stop "markov_steady_state: s_dist must have the same number &
                       &of states as adj_mat"

        stat = 0

        n = size(s_dist, 1)
        last = 1.0_wp/n

        do i=1,max_iter
            s_dist = matmul(adj_mat, last)
            if (sqrt(sum((s_dist-last)**2)/n) < tol) exit
            last = s_dist
        end do

        if (i > max_iter) stat = 1
    end subroutine markov_steady_state

    function make_disc_unif_dist(n) result(dist)

        integer, intent(in) :: n

        real(wp), dimension(n) :: dist

        integer :: i

        dist = [(1.0_wp, i=1, n)] / n

    end function make_disc_unif_dist

end module dist_mod
