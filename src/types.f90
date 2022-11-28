module types_mod

    use shared_mod, only: wp

    implicit none

    type, abstract, public :: Problem

        integer :: n_exg
        integer :: n_end
        integer :: n_act

        integer :: n_mom

    contains

        procedure(GetRInterface), deferred, pass :: get_r
        procedure(GetTExgInterface), deferred, pass :: get_t_exg
        procedure(GetTEndInterface), deferred, pass :: get_t_end
        procedure(GetMmtInterface), deferred, pass :: get_mmt

    end type Problem

    abstract interface

        function GetRInterface(self) result(r)
            import Problem, wp

            class(Problem), intent(in) :: self
            real(wp), dimension(self%n_exg,self%n_end,self%n_act) :: r
        end function GetRInterface

        function GetTExgInterface(self) result(t_exg)
            import Problem, wp

            class(Problem), intent(in) :: self
            real(wp), dimension(self%n_exg,self%n_exg) :: t_exg
        end function GetTExgInterface

        function GetTEndInterface(self) result(t_end)
            import Problem, wp

            class(Problem), intent(in) :: self
            integer, dimension(self%n_end,self%n_act) :: t_end
        end function GetTEndInterface

        function GetMmtInterface(self, t_exg, tn_pi, &
                                 n_sim, n_t, n_burn, &
                                 init_exg_dist, init_end_dist) result(mmt)
            import Problem, wp

            class(Problem), intent(inout) :: self
            real(wp), dimension(self%n_exg,self%n_exg), intent(in) :: t_exg
            real(wp), dimension(self%n_exg), intent(inout) :: init_exg_dist
            real(wp), dimension(self%n_end), intent(inout) :: init_end_dist
            integer, dimension(self%n_exg,self%n_end), intent(in) :: tn_pi

            integer, intent(in) :: n_sim
            integer, intent(in) :: n_t
            integer, intent(in) :: n_burn

            real(wp), dimension(self%n_mom) :: mmt
        end function GetMmtInterface

    end interface

    type, abstract, public :: Model

        class(Problem), pointer :: problem

    contains

        procedure(ModelOnlyInterface), deferred, public, pass :: run

    end type Model

    abstract interface

        subroutine ModelOnlyInterface(self)
            import Model, wp

            class(Model), intent(inout) :: self
        end subroutine ModelOnlyInterface

    end interface

end module types_mod
