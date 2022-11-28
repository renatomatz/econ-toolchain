module extra_stats_mod

    use stdlib_stats, only: mean, var, corr, cov

    use shared_mod, only: wp,  eps
    use utils_mod, only: col_norm_squeeze, col_concat

    implicit none

    interface

        subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
            import :: wp

            character :: trans
            integer :: info
            integer :: lda
            integer :: ldb
            integer :: lwork
            integer :: m
            integer :: n
            integer :: nrhs
            real(wp) :: a(lda,*)
            real(wp) :: b(ldb,*)
            real(wp) :: work(*)
        end subroutine dgels

    end interface

contains

    function xy_cov(x, y) result(cov_res)
        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y

        real(wp) :: cov_res
        real(wp), dimension(2,2) :: cov_mat

        cov_mat = cov(col_concat([x, y], 2), 1)
        cov_res = cov_mat(1,2)
    end function xy_cov

    function xy_corr(x, y) result(corr_res)
        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y

        real(wp) :: corr_res
        real(wp), dimension(2,2) :: corr_mat

        corr_mat = corr(col_concat([x, y], 2), 1)
        corr_res = corr_mat(1,2)
    end function xy_corr

    subroutine simple_regress(x, y, coef, intr, r, stat)
        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(size(x)), intent(in) :: y
        real(wp), intent(out) :: coef
        real(wp), intent(out) :: intr
        real(wp), dimension(size(x)), intent(out) :: r
        integer, intent(out) :: stat

        real(wp) :: tmp_coef

        stat = 0

        tmp_coef =  var(x)
        if (abs(tmp_coef) < eps) stat = 1
        tmp_coef = xy_cov(x, y) / tmp_coef
        coef = tmp_coef
        intr = mean(y) - tmp_coef*mean(x)
        r = y - tmp_coef*x - intr
    end subroutine simple_regress

    subroutine multi_regress(x, y, coef, intr, r, stat)
        real(wp), dimension(:,:), intent(in) :: x
        real(wp), dimension(size(x,1)), intent(in) :: y
        real(wp), dimension(size(x,2)), intent(out) :: coef
        real(wp), intent(out) :: intr
        real(wp), dimension(size(x,1)), intent(out) :: r
        integer, intent(out) :: stat

        real(wp), dimension(size(x,1),size(x,2)+1) :: A
        real(wp), dimension(size(x,1)) :: b
        real(wp), dimension(min(size(x,1),size(x,2)+1) &
                            + max(1,size(x,1),size(x,2)+1)) :: work
        integer :: m, n, lwork, info

        real(wp), dimension(size(x,2)+1) :: tmp_coef

        stat = 0

        m = size(x,1)
        ! we add one to account for the constant column
        n = size(x,2)+1
        lwork = size(work)

        A(:,1) = 1
        A(:,2:) = x
        b = y

        call dgels('N', m, n, 1, A, m, b, m, work, lwork, info)
        tmp_coef = b(:n)
        coef = tmp_coef(2:)
        intr = tmp_coef(1)
        r = y - matmul(A,tmp_coef)

        if (info /= 0) stat = 1
    end subroutine multi_regress

    subroutine panel_ar1_regress(x, coef, intr, r, stat)
        real(wp), dimension(:,:), intent(in) :: x
        real(wp), intent(out) :: coef
        real(wp), intent(out) :: intr
        real(wp), dimension(size(x)), intent(out) :: r
        integer, intent(out) :: stat

        real(wp), dimension(size(x)-size(x,1)) :: x_lag, x_fut

        stat = 0

        x_lag = col_norm_squeeze(x(:,:size(x,2)-1))
        x_fut = col_norm_squeeze(x(:,2:))

        call simple_regress(x_lag, x_fut, coef, intr, r, stat)
    end subroutine panel_ar1_regress

end module extra_stats_mod
