module coarray_utils_mod

    use shared_mod, only: wp

contains

    subroutine co_min_with_vec(x, y, best_img)
        real(wp), dimension(:), intent(inout) :: x
        real(wp), intent(inout) :: y
        integer, optional, intent(out) :: best_img

        real(wp), allocatable, dimension(:), codimension[:] :: best_x
        real(wp), allocatable, codimension[:] :: best_y

        integer :: L, cur_best_img

        integer:: p, me

        p = num_images()
        me = this_image()

        allocate(best_x(size(x))[*], best_y[*])

        best_x = x
        best_y = y
        cur_best_img = me

        ! Reduce Routine
        L = 1
        do while (L < p)
            if ((me+L<=p).and.(mod(me-1,2*L)==0)) then
                if (best_y[me+L] < best_y) then
                    best_x = best_x(:)[me+L]
                    best_y = best_y[me+L]
                    cur_best_img = me+L
                end if
            end if
            L = 2*L
            sync all
        end do

        call co_broadcast(best_x, 1)
        call co_broadcast(best_y, 1)

        x = best_x
        y = best_y
        if (present(best_img)) best_img = cur_best_img

        deallocate(best_x, best_y)
    end subroutine co_min_with_vec

    subroutine partition_params(tot_len, max_loc_len, loc_len, loc_base)

        integer, intent(in) :: tot_len
        integer, intent(out) :: max_loc_len
        integer, intent(out) :: loc_len
        integer, intent(out) :: loc_base

        integer :: r

        integer:: p, me

        p = num_images()
        me = this_image()

        r = mod(tot_len, p)
        max_loc_len = (tot_len-1)/p + 1

        loc_len = max_loc_len
        if ((r>0).and.(me>r)) loc_len = loc_len-1

        loc_base = (me-1)*loc_len
        if ((r>0).and.(me>r)) loc_base = loc_base+r

    end subroutine partition_params

end module coarray_utils_mod
