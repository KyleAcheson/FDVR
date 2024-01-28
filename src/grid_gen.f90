module grid_gen
	use types
	implicit none

	public :: uniform_grids, meshgrid_recursive

	contains

		function uniform_grids(params) result(grids)
			use types
			implicit none
			type(Parameters), intent(in) :: params
			real(dp), allocatable :: grids(:)
			integer :: grid_sum
			integer :: i, j, n, ngrid
			real(dp) :: dx
			real(dp) :: min_val, max_val

			grid_sum = sum(params%grid_sizes)
			allocate(grids(grid_sum))

			n = 1
			do i = 1, params%ndims
				ngrid = params%grid_sizes(i)
				min_val = params%min_grid_vals(i)
				max_val = params%max_grid_vals(i)
				dx = (max_val - min_val) / (ngrid - 1)
				do j = 1, ngrid
					grids(n) = (j-1) * dx + min_val
					n = n + 1
				end do
			end do

		end function uniform_grids

		recursive subroutine meshgrid_recursive(f, potential_func, grids, masses, indices, strides, sizes, n, ndim_total)
			use types
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

			real(dp), intent(inout) :: f(:)
			integer, intent(inout) :: indices(:)
			real(dp), intent(in) :: masses(:)
			real(dp), intent(in) :: grids(:)
			integer, intent(in) :: sizes(:)
			integer, intent(in) :: strides(:)
			integer, intent(in) :: n
			integer, intent(in) :: ndim_total
			procedure(func_type), pointer, intent(in) :: potential_func
			integer :: linear_index, i
			real(dp) :: grid_values(ndim_total)

			if (n > 0) then
				do i = 1, sizes(n)
					indices(n) = i
					call meshgrid_recursive(f, potential_func, grids, masses, indices, strides, sizes, n - 1, ndim_total)
				end do
			else
				! Base case: calculate the array value
				linear_index = 1
				do i = size(sizes), 1, -1
					linear_index = linear_index + (indices(i) - 1) * product(sizes(i+1:))
				end do
				grid_values = grids(indices + strides)
				f(linear_index) = potential_func(grid_values, masses, ndim_total)
			end if
			!print *, 'strided indicies:', (indices + strides), 'grid_vals:', grids(indices + strides)
		end subroutine meshgrid_recursive

end module grid_gen