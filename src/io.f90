module io
	use types
	implicit none
	public :: read_inputs
	private :: split_line, process_keyword

	contains

		subroutine read_inputs(file_name, params)
			use types
			implicit none
			character(len=*), intent(in) :: file_name
			type(Parameters), intent(out) :: params
			character(100) :: parameter_name
			character(100) :: emsg
			integer :: iostatus

			open(unit=1, file=file_name, status='old', iostat=iostatus, iomsg=emsg)

			if (iostatus /= 0) then
				print *, 'Error opening file: ', trim(file_name)
				print *, trim(emsg)
				stop
			end if

			! Read parameters from the file
			do
				read(1, '(A)', iostat=iostatus, iomsg=emsg) parameter_name

				! Exit the loop when end-of-file is reached
				if (iostatus > 0) then
					print *, trim(emsg)
				else if (iostatus < 0) then
					exit
				end if


				! Split line and get keyword and value
				call split_line(trim(parameter_name), params)
			end do

			close(1)

			print *, 'INPUT PARAMETERS:'
			print *, 'neig: ', params%neig
			print *, 'ndims: ', params%ndims
			print *, 'masses: ', params%masses
			print *, 'grid_sizes: ', params%grid_sizes
			print *, 'min_grid_vals: ', params%min_grid_vals
			print *, 'max_grid_vals: ', params%max_grid_vals

		end subroutine read_inputs

		subroutine split_line(line, params)
			use types
			implicit none
			character(len=*), intent(in) :: line
			type(Parameters), intent(inout) :: params
			character(50) :: keyword
			character(100) :: value_str
			integer :: i, iostat, line_length

			line_length = len_trim(line)
			i = 1
			do while (i <= line_length)

				! Skip whitespaces
				do while (i <= line_length .and. line(i:i) /= '=')
					i = i + 1
				end do

				if (i > line_length) exit  ! End of line reached

				! Check for '=' character
				if (line(i:i) == '=') then
					read(line(:i-1), '(A)', iostat=iostat) keyword
					i = i + 1  ! Skip '=' character

					! Skip whitespaces
					do while (i <= line_length .and. line(i:i) == ' ')
						i = i + 1
					end do

					! Read the value string
					read(line(i:), '(A)', iostat=iostat) value_str
					if (iostat /= 0) exit

					! Process the keyword and assign the value to the corresponding parameter
					call process_keyword(keyword, value_str, params)
				end if

				! Skip to the next keyword
				i = i + len_trim(value_str)
			end do
		end subroutine split_line

		subroutine process_keyword(keyword, value_str, params)
			use types
			implicit none
			character(50), intent(in) :: keyword
			character(100), intent(in) :: value_str
			type(Parameters), intent(inout) :: params

			! Process the keyword and assign the value to the corresponding parameter
			select case (trim(keyword))
			case ('ndims')
				read(value_str, *) params%ndims
			case ('neig')
				read(value_str, *) params%neig
			case ('masses')
				allocate(params%masses(params%ndims))
				read(value_str, *) params%masses
			case ('grid_sizes')
				allocate(params%grid_sizes(params%ndims))
				read(value_str, *) params%grid_sizes
			case ('min_grid_vals')
				allocate(params%min_grid_vals(params%ndims))
				read(value_str, *) params%min_grid_vals
			case ('max_grid_vals')
				allocate(params%max_grid_vals(params%ndims))
				read(value_str, *) params%max_grid_vals
			case ('potential_type')
				read(value_str, *) params%potential
			case ('algorithm')
				read(value_str, *) params%algorithm
			case default
				print *, 'Unknown parameter:', keyword
			end select
		end subroutine process_keyword

end module io