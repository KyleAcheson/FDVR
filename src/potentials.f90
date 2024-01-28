module potentials
	use types
	use grid_gen
	public :: get_potential
	private :: harmonic

	contains

	function get_potential(grids, nelems, grid_sizes, masses, ndims, potential_type) result(potential)
		use types
		use grid_gen

		implicit none
		abstract interface
			function func_type(grid_values, masses, ndim_total)
				use types
				real(dp), intent(in) :: grid_values(:)
				real(dp), intent(in) :: masses(:)
				integer, intent(in) :: ndim_total
				real(dp) :: func_type
			end function func_type
		end interface

		integer, intent(in) :: nelems
		integer, intent(in) :: ndims
		integer, intent(in) :: potential_type
		real(dp), intent(in) :: grids(nelems)
		real(dp), intent(in) :: masses(ndims)
		integer, intent(in) :: grid_sizes(ndims)
		procedure(func_type), pointer :: function_pointer
		real(dp), allocatable :: potential(:)
		integer, allocatable :: strides(:)
		integer, allocatable :: indicies(:)
		integer :: i

		select case(potential_type)
			case (0)
				function_pointer => harmonic
			case default
				print *, 'Error: Undefined potential type'
				stop
		end select

		allocate(indicies(ndims))
		allocate(strides(ndims))

		strides(1) = 0
		do i = 1, ndims-1
			strides(i+1) = strides(i) + grid_sizes(i)
		end do

		allocate(potential(nelems))

		call meshgrid_recursive(potential, function_pointer, grids, masses, indicies, strides, grid_sizes, ndims, ndims)

		deallocate(indicies)
		deallocate(strides)

	end function get_potential


	function harmonic(grid_values, masses, ndims) result(func_value)
		use types
		implicit none
		real(dp), intent(in) :: grid_values(:)
		real(dp), intent(in) :: masses(:)
		integer, intent(in) :: ndims
		real(dp) :: func_value
		integer :: i

		func_value = 0

		do i = 1, ndims
			func_value = func_value + ((grid_values(i)**2) / masses(i))
		end do

		func_value = func_value * 0.5_dp

	end function harmonic


end module potentials