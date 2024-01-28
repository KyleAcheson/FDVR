module types
	implicit none
	integer, parameter :: sp = selected_real_kind(6, 37)
	integer, parameter :: dp = selected_real_kind(15, 307)
	real(sp) :: r_sp = 1.0
	real(dp) :: r_dp = 1.0_dp

	type :: Parameters
		integer :: ndims
		integer :: neig
		integer :: potential
		integer :: algorithm
		integer, allocatable :: grid_sizes(:)
		real(dp), allocatable :: min_grid_vals(:)
		real(dp), allocatable :: max_grid_vals(:)
		real(dp), allocatable :: masses(:)
	end type Parameters

	type :: list
		real(dp), allocatable :: elements(:, :)
	end type list

end module types