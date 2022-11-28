module shared_mod

    use iso_fortran_env, only: real64, real32

    implicit none

    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
    integer, parameter :: wp = dp

    real(wp), parameter :: zero = 0.0_wp
    real(wp), parameter :: one = 1.0_wp
    real(wp), parameter :: inf = HUGE(zero)

    real(wp), parameter :: eps = 1.0e-5_wp
    real(wp), parameter :: golden_ratio = (sqrt(5.0_wp) - one) / 2.0_wp

    character(len=*), parameter :: real_fmt0 = "(F24.12)"

    abstract interface

        function Func0() result(y)
            import wp

            real(wp) :: y
        end function Func0

        function Func00(x) result(y)
            import wp

            real(wp), intent(in) :: x
            real(wp) :: y
        end function Func00

        function Func11(x) result(y)
            import wp

            real(wp), dimension(:), intent(in) :: x
            real(wp), dimension(size(x)) :: y
        end function Func11

        function Func12(x) result(y)
            import wp

            real(wp), dimension(:), intent(in) :: x
            real(wp), dimension(size(x), size(x)) :: y
        end function Func12

        function Func1m0(x, n) result(y)
            import wp

            real(wp), dimension(n), intent(in) :: x
            integer, intent(in) :: n

            real(wp) :: y
        end function Func1m0

        function Func1m1m(x, n, m) result(y)
            import wp

            real(wp), dimension(n), intent(in) :: x
            integer, intent(in) :: n
            integer, intent(in) :: m

            real(wp), dimension(m) :: y
        end function Func1m1m

        function Func1m2mn(x, m, n) result(y)
            import wp

            real(wp), dimension(m), intent(in) :: x
            integer, intent(in) :: n
            integer, intent(in) :: m

            real(wp), dimension(m, n) :: y
        end function Func1m2mn

    end interface

contains

    function is_err(x) result(y)
        real(wp), intent(in) :: x

        logical :: y

        y = .not.(abs(x-x)<1e-10_wp).or.(x<-inf).or.(x>inf)
    end function is_err

end module shared_mod
