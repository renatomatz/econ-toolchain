module simulate_mod

    use shared_mod, only: wp
    use dist_mod, only: make_disc_unif_dist

    implicit none

contains

    subroutine simulate_states(t_exg, tn_pi, sim_idx, &
                               n_exg, n_end, n_sim, n_t, n_burn, &
                               init_exg_dist, init_end_dist)

        real(wp), dimension(n_exg,n_exg), intent(in) :: t_exg
        integer, dimension(n_exg,n_end), intent(in) :: tn_pi

        integer, dimension(n_sim,n_t,2), intent(out) :: sim_idx

        integer, intent(in) :: n_exg, n_end
        integer, intent(in) :: n_sim, n_t, n_burn

        real(wp), dimension(n_exg), intent(inout) :: init_exg_dist
        real(wp), dimension(n_end), intent(inout) :: init_end_dist

        real(wp), dimension(n_sim) :: harvest

        real(wp), dimension(n_exg) :: init_exg_cdf
        real(wp), dimension(n_end) :: init_end_cdf
        real(wp), dimension(n_exg,n_exg) :: tx_cdf

        integer, dimension(n_sim, 2) :: burn_idx

        integer :: i, j

        ! cdf initializations
        init_exg_cdf = init_exg_dist
        init_end_cdf = init_end_dist
        tx_cdf = t_exg

        ! cdf cummulative sums
        do i=2,n_exg
            init_exg_cdf(i) = init_exg_cdf(i) + init_exg_cdf(i-1)
            tx_cdf(:,i) = tx_cdf(:,i-1) + tx_cdf(:,i)
        end do
        do i=2,n_end
            init_end_cdf(i) = init_end_cdf(i) + init_end_cdf(i-1)
        end do

        call random_number(harvest)
        burn_idx(:,1) = rnd_idx_arr(harvest, init_exg_cdf)
        call random_number(harvest)
        burn_idx(:,2) = rnd_idx_arr(harvest, init_end_cdf)

        do i=1,n_burn
            call random_number(harvest)
            do j=1,n_sim
                burn_idx(j,2) = tn_pi(burn_idx(j,1), burn_idx(j,2))
                burn_idx(j,1) = rnd_idx(harvest(j), &
                                        tx_cdf(burn_idx(j,1),:))
            end do
        end do

        call random_number(harvest)
        sim_idx(:,1,1) = burn_idx(:,1)
        sim_idx(:,1,2) = burn_idx(:,2)
        do i=2,n_t
            call random_number(harvest)
            do j=1,n_sim
                sim_idx(j,i,2) = tn_pi(sim_idx(j,i-1,1), sim_idx(j,i-1,2))
                sim_idx(j,i,1) = rnd_idx(harvest(j), &
                                         tx_cdf(sim_idx(j,i-1,1),:))
            end do
        end do

    end subroutine simulate_states

    function rnd_idx_arr(rnd, cdf) result(idx)
        real(wp), dimension(:), intent(in) :: rnd
        real(wp), dimension(:), intent(in) :: cdf

        integer, dimension(size(rnd)) :: idx

        logical, dimension(size(cdf),size(rnd)) :: mask

        integer :: i

        do i=1,size(rnd)
            mask(:,i) = cdf <= rnd(i)
        end do

        idx = count(mask, 1) + 1
    end function rnd_idx_arr

    function rnd_idx(rnd, cdf) result(idx)
        real(wp), intent(in) :: rnd
        real(wp), dimension(:), intent(in) :: cdf
        integer :: idx

        idx = count(cdf <= rnd) + 1
    end function rnd_idx

    function get_tn_pi(t_end, pi, n_exg, n_end, n_act) result(tn_pi)

        integer, dimension(n_end,n_act), intent(in) :: t_end
        integer, dimension(n_end,n_exg), intent(in) :: pi

        integer, intent(in) :: n_exg, n_end, n_act

        integer, dimension(n_end,n_exg) :: tn_pi

        integer :: i, j

        do j=1,size(t_end,2) ! endogenous states
            do i=1,size(t_end,1) ! exogenous states
                tn_pi(i,j) = t_end(j,pi(i,j))
            end do
        end do

    end function get_tn_pi

    function make_sim_mat(sim_idx, mat) result(sim_mat)
        integer,  dimension(:,:,:), intent(in) :: sim_idx
        real(wp), dimension(:,:), intent(in) :: mat

        real(wp), dimension(size(sim_idx,1),size(sim_idx,2)) :: sim_mat

        integer :: i, j

        do j=1,size(sim_idx,2) ! time axis
            do i=1,size(sim_idx,1) ! sim axis
                sim_mat(i,j) = mat(sim_idx(i,j,1),sim_idx(i,j,2))
            end do
        end do
    end function make_sim_mat

    function make_exg_sim_mat(sim_idx, vec) result(sim_mat)
        integer,  dimension(:,:,:), intent(in) :: sim_idx
        real(wp), dimension(:), intent(in) :: vec

        real(wp), dimension(size(sim_idx,1),size(sim_idx,2)) :: sim_mat

        integer :: i, j

        do j=1,size(sim_idx,2) ! time axis
            do i=1,size(sim_idx,1) ! sim axis
                sim_mat(i,j) = vec(sim_idx(i,j,1))
            end do
        end do
    end function make_exg_sim_mat

    function make_end_sim_mat(sim_idx, vec) result(sim_mat)
        integer,  dimension(:,:,:), intent(in) :: sim_idx
        real(wp), dimension(:), intent(in) :: vec

        real(wp), dimension(size(sim_idx,1),size(sim_idx,2)) :: sim_mat

        integer :: i, j

        do j=1,size(sim_idx,2) ! time axis
            do i=1,size(sim_idx,1) ! sim axis
                sim_mat(i,j) = vec(sim_idx(i,j,2))
            end do
        end do
    end function make_end_sim_mat

    function make_end_pol_sim_mat(sim_idx, vec, tn_pi) result(sim_mat)
        integer,  dimension(:,:,:), intent(in) :: sim_idx
        real(wp), dimension(:), intent(in) :: vec
        integer, dimension(:,:), intent(in) :: tn_pi

        real(wp), dimension(size(sim_idx,1),size(sim_idx,2)) :: sim_mat

        integer :: i, j

        do j=1,size(sim_idx,2) ! time axis
            do i=1,size(sim_idx,1) ! sim axis
                sim_mat(i,j) = vec(tn_pi(sim_idx(i,j,1), sim_idx(i,j,2)))
            end do
        end do
    end function make_end_pol_sim_mat

end module simulate_mod
