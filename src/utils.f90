module utils_mod

    use stdlib_stats, only: mean

    use shared_mod, only: wp, real_fmt0, Func00, Func11, is_err

    implicit none

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import wp

            integer :: n
            integer nrhs
            real(wp), dimension(lda,*) :: a
            integer :: lda
            integer, dimension(*) :: ipiv
            real(wp), dimension(ldb,*) :: b
            integer :: ldb
            integer :: info
        end subroutine dgesv
    end interface

contains

    subroutine check_rtv_(r, t_exg, t_end, v, n_exg, n_end, n_act)

        real(wp), dimension(:,:,:), intent(in) :: r
        real(wp), dimension(:,:), intent(in) :: t_exg
        integer, dimension(:,:), intent(in) :: t_end
        real(wp), dimension(:,:), intent(out) :: v

        integer, intent(in) :: n_exg
        integer, intent(in) :: n_end
        integer, intent(in) :: n_act

        if ((size(r, 1) /= n_exg).or.&
            (size(r, 2) /= n_end).or.&
            (size(r, 3) /= n_act)) &
            error stop "check_rtv: incorrect r (rewards) shape"

        if ((size(t_exg, 1) /= n_exg).or.&
            (size(t_exg, 2) /= n_exg)) &
            error stop "check_rtv: incorrect t_exg (exogenous transition) &
                       &shape"

        if ((size(t_end, 1) /= n_end).or.&
            (size(t_end, 2) /= n_act)) &
            error stop "check_rtv: incorrect t_end (endogenous transition) &
                       &shape"

        if ((size(v, 1) /= n_exg).or.&
            (size(v, 2) /= n_end)) &
            error stop "check_rtv: incorrect v (values) shape"

    end subroutine check_rtv_

    function row_concat(x, n) result(x_concat)
        real(wp), dimension(:), intent(in) :: x
        integer, intent(in) :: n

        real(wp), dimension(n,size(x)/n) :: x_concat

        x_concat = reshape(x, [n,size(x)/n])
    end function row_concat

    function col_concat(x, n) result(x_concat)
        real(wp), dimension(:), intent(in) :: x
        integer, intent(in) :: n

        real(wp), dimension(size(x)/n,n) :: x_concat

        x_concat = reshape(x, [size(x)/n,n])
    end function col_concat

    function norm_cols(x) result(x_norm)
        real(wp), dimension(:,:), intent(in) :: x

        real(wp), dimension(size(x,1),size(x,2)) :: x_norm

        real(wp), dimension(size(x,1)) :: x_mean

        integer :: i

        x_mean = mean(x, 2)

        x_norm = x
        do i=1,size(x,1)
            x_norm(i,:) = x_norm(i,:) - x_mean(i)
        end do
    end function norm_cols

    function squeeze(x) result(x_squeeze)
        real(wp), dimension(:,:), intent(in) :: x

        real(wp), dimension(size(x)) :: x_squeeze

        x_squeeze = reshape(x, [size(x)])
    end function squeeze

    function col_norm_squeeze(x) result(x_normed_vec)
        real(wp), dimension(:,:), intent(in) :: x

        real(wp), dimension(size(x)) :: x_normed_vec

        x_normed_vec = squeeze(norm_cols(x))
    end function col_norm_squeeze

    function proportion(x) result(prop)
        logical, dimension(:), intent(in) :: x

        real(wp) :: prop

        prop = real(count(x), wp) / size(x)
    end function proportion

    subroutine fixed_point(f, x, tol, max_iter, stat)

        procedure(Func00) :: f
        real(wp), intent(inout) :: x

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp) :: x_old
        integer :: i

        stat = 0

        x_old = x
        do i=1,max_iter
            x = f(x_old)
            if (abs(x_old - x) < tol) exit
            x_old = x
        end do

        if (i > max_iter) stat = 1

    end subroutine fixed_point

    subroutine fixed_point_mult(f, x, tol, max_iter, stat)

        procedure(Func11) :: f
        real(wp), dimension(:), intent(inout) :: x

        real(wp), intent(in) :: tol
        integer, intent(in) :: max_iter
        integer, intent(out) :: stat

        real(wp), dimension(size(x)) :: x_old
        integer :: i

        stat = 0

        x_old = x
        do i=1,max_iter
            x = f(x_old)
            if (maxval(abs(x_old - x)) < tol) exit
            x_old = x
        end do

        if (i > max_iter) stat = 1

    end subroutine fixed_point_mult

    subroutine solve(A, xb, n)

        real(wp), dimension(n,n), intent(in) :: A
        real(wp), dimension(n), intent(inout) :: xb
        integer, intent(in) :: n

        integer, dimension(n) :: ipiv
        integer :: info

        call dgesv(n, 1, A, n, ipiv, xb, n, info)

    end subroutine solve

    subroutine write_vec(x, file_name)
        real(wp), dimension(:), intent(in) :: x
        character(len=*), intent(in) :: file_name

        character(80) :: err_msg
        integer :: unt, err_stat

        open(newunit=unt, file=file_name, &
             action="write", status="unknown", iostat=err_stat, iomsg=err_msg)

        if (err_stat /= 0) then
            write (*,'(A,I6)') "Error opening file, stat: ", err_stat
            write (*,'(A)') trim(err_msg)
            error stop
        end if

        write(unt, *) x
        close(unt)
    end subroutine write_vec

    subroutine read_vec(x, file_name)
        real(wp), dimension(:), intent(out) :: x
        character(len=*), intent(in) :: file_name

        character(80) :: err_msg
        integer :: unt, err_stat

        open(newunit=unt, file=file_name, &
             action="read", status="unknown", iostat=err_stat, iomsg=err_msg)

        if (err_stat /= 0) then
            write (*,'(A,I6)') "Error opening file, stat: ", err_stat
            write (*,'(A)') trim(err_msg)
            error stop
        end if

        read(unt, *) x
        close(unt)
    end subroutine read_vec

    subroutine write_scalar(x, file_name)
        real(wp), intent(in) :: x
        character(len=*), intent(in) :: file_name

        character(80) :: err_msg
        integer :: unt, err_stat

        open(newunit=unt, file=file_name, &
             action="write", status="unknown", iostat=err_stat, iomsg=err_msg)

        if (err_stat /= 0) then
            write (*,'(A,I6)') "Error opening file, stat: ", err_stat
            write (*,'(A)') trim(err_msg)
            error stop
        end if

        write(unt, *) x
        close(unt)
    end subroutine write_scalar

    subroutine read_scalar(x, file_name)
        real(wp), intent(out) :: x
        character(len=*), intent(in) :: file_name

        character(80) :: err_msg
        integer :: unt, err_stat

        open(newunit=unt, file=file_name, &
             action="read", status="unknown", iostat=err_stat, iomsg=err_msg)

        if (err_stat /= 0) then
            write (*,'(A,I6)') "Error opening file, stat: ", err_stat
            write (*,'(A)') trim(err_msg)
            error stop
        end if

        read(unt, *) x
        close(unt)
    end subroutine read_scalar

    subroutine write_res(x, y, base_name)
        real(wp), dimension(:), intent(in) :: x
        real(wp), intent(in) :: y
        character(len=*), intent(in) :: base_name

        logical :: file_exists

        inquire(file=base_name // "_x", exist=file_exists)
        if (file_exists) then
            call execute_command_line("rm " // base_name // "_*")
        end if

        call write_vec(x, base_name // "_x")
        call write_scalar(y, base_name // "_y")
    end subroutine write_res

    subroutine read_res(x, y, base_name)
        real(wp), dimension(:), intent(out) :: x
        real(wp), intent(out) :: y
        character(len=*), intent(in) :: base_name

        call read_vec(x, base_name // "_x")
        call read_scalar(y, base_name // "_y")
    end subroutine read_res

    function time_vec_to_sec(x) result(y)
        real(wp), dimension(3), intent(in) :: x

        real(wp) :: y

        y = sum(x * [3600, 60, 1])
    end function time_vec_to_sec

    function time_sec_to_vec(x) result(y)
        real(wp), intent(in) :: x
        real(wp), dimension(3) :: y

        y(1) = floor(x / 3600)
        y(2) = floor((x-y(1)) / 60)
        y(3) = x-y(1)-y(2)
    end function time_sec_to_vec

    function lclip(x, c) result(y)
        real(wp), dimension(:), intent(in) :: x
        real(wp), intent(in) :: c

        real(wp), dimension(size(x)) :: y

        real(wp), dimension(size(x),2) :: comp_vec

        comp_vec(:,1) = x
        comp_vec(:,2) = c
        y = maxval(comp_vec, 2)
    end function lclip

    function uclip(x, c) result(y)
        real(wp), dimension(:), intent(in) :: x
        real(wp), intent(in) :: c

        real(wp), dimension(size(x)) :: y

        real(wp), dimension(size(x),2) :: comp_vec

        comp_vec(:,1) = x
        comp_vec(:,2) = c
        y = minval(comp_vec, 2)
    end function uclip

    function not_err_elems(x, n) result(y)
        real(wp), dimension(:), intent(in) :: x
        integer, optional, intent(out) :: n

        logical, dimension(size(x)) :: y

        integer :: i

        do i=1,size(x)
            y(i) = .not.is_err(x(i))
        end do
        if (present(n)) n = count(.not.y)
    end function not_err_elems

end module utils_mod
