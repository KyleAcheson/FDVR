program DVR
	use types
	use io
	use grid_gen
	use potentials
	use dvr_solvers
	!use lapack

	implicit none

	type(Parameters) :: input_params
	character(len=256) :: input_fname
	real(dp), allocatable :: grids(:), potential(:)
	integer :: nelems

	if (COMMAND_ARGUMENT_COUNT() /= 1) then
		print *, 'USAGE: ./dvr input.txt'
		stop
	end if

	call GET_COMMAND_ARGUMENT(1, input_fname)
	input_fname = trim(input_fname)
	call read_inputs(input_fname, input_params)
	grids = uniform_grids(input_params)
	nelems = product(input_params%grid_sizes)
	potential = get_potential(grids, nelems, input_params%grid_sizes, input_params%masses, input_params%ndims, input_params%potential)
	call dvr_base(potential, grids, input_params%grid_sizes, nelems, input_params%masses, &
			      input_params%ndims, input_params%neig, input_params%algorithm)

	print *, 'CALCULATION FINISHED'

end program
