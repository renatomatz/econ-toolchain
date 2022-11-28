module unit_testing_mod
    use iso_fortran_env, only: error_unit
    use shared_mod, only: wp

    implicit none

    character(len=*), parameter :: passed_msg = "TEST PASSED"
    character(len=*), parameter :: final_msg = "TEST RESULTS:"
    character(len=*), parameter :: wanted_msg = "WANTED:"
    character(len=*), parameter :: got_msg = "GOT:"

    character(len=*), parameter :: sep1 = "---------------------------------&
                                          &---------------------------------"
    character(len=*), parameter :: sep2 = "=================================&
                                          &================================="
    character(len=*), parameter :: sep3 = "+++++++++++++++++++++++++++++++++&
                                          &+++++++++++++++++++++++++++++++++"

    character(len=*), parameter :: passed_symb = "(:"
    character(len=*), parameter :: failed_symb = "XX"

    real(wp), parameter :: default_eps = 1e-10

    type, public :: UnitTester

        integer :: n_tests
        integer :: n_passed

        logical :: verbose = .true.

    contains

        procedure, public, pass ::  init_test
        procedure, public, pass ::  end_test

        procedure, public, pass ::  set_verbose

        procedure, private, pass ::  sub_t_helper1

        procedure, public, pass ::  int_scalar

        procedure, public, pass ::  real_scalar
        procedure, public, pass ::  real_vector

    end type UnitTester

contains

    subroutine init_test(self, test_name)

        class(UnitTester), intent(inout) :: self
        character(len=*), intent(in) :: test_name

        write(error_unit,*) sep2
        write(error_unit,*) test_name
        write(error_unit,*) sep2

        self%n_tests = 0
        self%n_passed = 0

    end subroutine init_test

    subroutine set_verbose(self, verbose)

        class(UnitTester), intent(inout) :: self
        logical, intent(in) :: verbose

        self%verbose = verbose

    end subroutine set_verbose

    subroutine end_test(self)

        class(UnitTester), intent(in) :: self

        write(error_unit,*)
        write(error_unit,*) final_msg
        write(error_unit,*) self%n_passed, '/', self%n_tests
        write(error_unit,*) sep2

    end subroutine end_test

    function sub_t_helper1(self, res, sub_t_name) result(not_res)

        class(UnitTester), intent(inout) :: self
        logical, intent(in) :: res
        character(len=*), intent(in) :: sub_t_name

        logical :: not_res

        character(len=64) :: msg
        character(len=3) :: symb

        self%n_tests = self%n_tests + 1

        msg = adjustl(sub_t_name)

        if (res) then
            symb = passed_symb
        else
            symb = failed_symb
        end if

        write(error_unit,'(X, A64, A3)') msg, symb

        if (res) self%n_passed = self%n_passed + 1

        ! return whether the test didn't pass.
        ! this is purely to make calling function code cleaner when errors happen
        ! see examples of why we use this below.
        not_res = (.not.res).and.(self%verbose)

    end function sub_t_helper1

    subroutine int_scalar(self, sub_t_name, wanted, got)

        class(UnitTester), intent(inout) :: self

        character(len=*), intent(in) :: sub_t_name
        integer, intent(in) :: wanted
        integer, intent(in) :: got

        if (self%sub_t_helper1(wanted==got, sub_t_name)) then
            write(error_unit,*) sep1
            write(error_unit,*) wanted_msg
            write(error_unit,*) wanted
            write(error_unit,*) got_msg
            write(error_unit,*) got
            write(error_unit,*) sep1
        end if
    end subroutine int_scalar

    subroutine real_scalar(self, sub_t_name, wanted, got, opt_eps)

        class(UnitTester), intent(inout) :: self

        character(len=*), intent(in) :: sub_t_name

        real(wp), intent(in) :: wanted
        real(wp), intent(in) :: got

        real(wp), optional, intent(in) :: opt_eps

        real(wp) :: eps

        if (present(opt_eps)) then
            eps = opt_eps
        else
            eps = default_eps
        end if


        if (self%sub_t_helper1(abs(wanted-got)<eps, sub_t_name)) then
            write(error_unit,*) sep1
            write(error_unit,*) wanted_msg
            write(error_unit,*) wanted
            write(error_unit,*) got_msg
            write(error_unit,*) got
            write(error_unit,*) sep1
        end if
    end subroutine real_scalar

    subroutine real_vector(self, sub_t_name, wanted, got, opt_eps)

        class(UnitTester), intent(inout) :: self

        character(len=*), intent(in) :: sub_t_name

        real(wp), dimension(:), intent(in) :: wanted
        real(wp), dimension(:), intent(in) :: got

        real(wp), optional, intent(in) :: opt_eps

        real(wp) :: eps

        if (present(opt_eps)) then
            eps = opt_eps
        else
            eps = default_eps
        end if

        if (self%sub_t_helper1(maxval(abs(wanted-got))<eps, sub_t_name)) then
            write(error_unit,*) sep1
            write(error_unit,*) wanted_msg
            write(error_unit,*) wanted
            write(error_unit,*) got_msg
            write(error_unit,*) got
            write(error_unit,*) sep1
        end if
    end subroutine real_vector

end module unit_testing_mod
