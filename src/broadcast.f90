module broadcast_mod

    use shared_mod, only: wp

    implicit none

contains

    function bc_add(m, v) result(mv)
        real(wp), dimension(:,:), intent(in) :: m
        real(wp), dimension(size(m,2)), intent(in) :: v

        real(wp), dimension(size(m,1),size(m,2)) :: mv

        integer :: i

        do i=1,size(m,2)
            mv(:,i) = m(:,i) + v(i)
        end do
    end function bc_add

    function bc_sub(m, v) result(mv)
        real(wp), dimension(:,:), intent(in) :: m
        real(wp), dimension(size(m,2)), intent(in) :: v

        real(wp), dimension(size(m,1),size(m,2)) :: mv

        integer :: i

        do i=1,size(m,2)
            mv(:,i) = m(:,i) - v(i)
        end do
    end function bc_sub

    function bc_mul(m, v) result(mv)
        real(wp), dimension(:,:), intent(in) :: m
        real(wp), dimension(size(m,2)), intent(in) :: v

        real(wp), dimension(size(m,1),size(m,2)) :: mv

        integer :: i

        do i=1,size(m,2)
            mv(:,i) = m(:,i) * v(i)
        end do
    end function bc_mul

    function bc_div(m, v) result(mv)
        real(wp), dimension(:,:), intent(in) :: m
        real(wp), dimension(size(m,2)), intent(in) :: v

        real(wp), dimension(size(m,1),size(m,2)) :: mv

        integer :: i

        do i=1,size(m,2)
            mv(:,i) = m(:,i) / v(i)
        end do
    end function bc_div

    function bc_pow(m, v) result(mv)
        real(wp), dimension(:,:), intent(in) :: m
        real(wp), dimension(size(m,2)), intent(in) :: v

        real(wp), dimension(size(m,1),size(m,2)) :: mv

        integer :: i

        do i=1,size(m,2)
            mv(:,i) = m(:,i) ** v(i)
        end do
    end function bc_pow

end module broadcast_mod
