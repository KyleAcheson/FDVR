module dvr_solvers
	use types
	implicit none

	public :: dvr_base, colbert_miller
	private :: kinetic_matrix_nd, init_n_matricies, deallocate_n_matricies, kronecker_product, eye

	contains


	subroutine dvr_base(potential, grids, sizes, nelems, masses, ndims, neig, algorithm)
		use types
		implicit none
		external :: DSYEVR

		abstract interface
			function solver_type(grid, ngrid, mass)
				use types
				real(dp), intent(in) :: grid(:)
				integer, intent(in) :: ngrid
				real(dp), intent(in) :: mass
				real(dp), dimension(ngrid, ngrid) :: solver_type
			end function solver_type
		end interface

		real(dp), intent(in) :: potential(:)
		real(dp), intent(in) :: grids(:)
		real(dp), intent(in) :: masses(:)
		integer, intent(in) :: sizes(:)
		integer, intent(in) :: nelems, ndims, neig, algorithm
		procedure(solver_type), pointer :: algorithm_pointer
		real(dp), allocatable :: H(:, :)
		integer :: i, j, count

		! DSYEVR parameters
		integer :: n, lda, ldz, lwork, liwork, m, info, nb
		real(dp), allocatable :: w(:)
		real(dp), allocatable :: z(:, :)
		real(dp), allocatable :: work(:)
		real(dp) :: abstol, vl, vu
		integer, allocatable :: isuppz(:), iwork(:)
		real(dp) :: dummy(1)
		integer :: idum(1)

		select case(algorithm)
			case (0)
				algorithm_pointer => colbert_miller
			case default
				print *, 'Error: undefined DVR solver'
				stop
		end select

		allocate(H(nelems, nelems))
		H(:, :) = 0.0_dp

		! adds kinetic matrix elems to Hamiltonian
		call kinetic_matrix_nd(algorithm_pointer, grids, masses, sizes, ndims, H)

		count = 1
		do i = 1, nelems
			H(i, i) = H(i, i) + potential(count) ! add potential to Hamiltonian
			count = count + 1
		end do

		! call LAPACK eigensolver
		n = nelems
		nb = 64
		m = neig-1+1
		lda = n
		ldz = n
		lwork = -1
		liwork = -1
		abstol = 0.0_dp
		!vl = 0.0_dp
		!vu = 0.0_dp

		allocate (w(m), z(ldz,m), isuppz(2*m))
		! querry optimal workspace
		call DSYEVR('V', 'I', 'U', n, H, lda, vl, vu, 1, neig, abstol, m, w, z, ldz, isuppz, dummy, lwork, idum, liwork, info)

		lwork = max((nb+6)*n, nint(dummy(1)))
		liwork = max(10*n, idum(1))
		Allocate (work(lwork), iwork(liwork))

		!lwork = nint(dummy(1))
		!liwork = idum(1)
		!allocate(work(lwork), iwork(liwork))

		! do computation
		call DSYEVR('V', 'I', 'U', n, H, lda, vl, vu, 1, neig, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)

		! check for successful computation
		if (info == 0) then
			! Eigenvalues are in w, eigenvectors are in z
			print *, 'SUCCESFULL EIGENVALUE COMPUTATION'
			print *, 'Eigenvalues:', w(:neig)
			!print *, 'Eigenvectors:'
			!do i = 1, neig
			!	print *, 'eigvec', i, ':'
			!	print *, z(:, i)
			!end do
		else
			print *, 'Error:', info
		end if

		! Deallocate workspace
		deallocate(work)
		deallocate(iwork)
		deallocate(isuppz)
		deallocate(w)
		deallocate(z)

	end subroutine dvr_base


	subroutine kinetic_matrix_nd(calculator, grids, masses, sizes, dim, Tmat)
		use types
		implicit none

		abstract interface
			function solver_type(grid, ngrid, mass)
				use types
				real(dp), intent(in) :: grid(:)
				integer, intent(in) :: ngrid
				real(dp), intent(in) :: mass
				real(dp), dimension(ngrid, ngrid) :: solver_type
			end function solver_type
		end interface

		procedure(solver_type), pointer :: calculator
		real(dp), intent(in) :: grids(:)
		real(dp), intent(in) :: masses(:)
		integer, intent(in) :: sizes(:)
		integer, intent(in) :: dim
		real(dp), intent(out) :: Tmat(:,:)
		integer :: i, j, init_size, new_size, total_size, start_ind, end_ind
		real(dp), allocatable :: Tslice(:, :), temp(:, :)
		type(list), allocatable :: matricies_list(:)
		real(dp) :: mass
		integer :: n

		if (dim < 2) then
			print *, "Error: Dimension must be at least 2."
			stop
		end if

		total_size = product(sizes)
		Tmat = 0.0_dp

		init_size = product(sizes(:2))
		allocate(temp(init_size, init_size))
		allocate(Tslice(init_size, init_size))

		start_ind = 1
		end_ind = 0

		do i = 1, dim

			end_ind = end_ind + sizes(i)

			matricies_list = init_n_matricies(sizes, dim)
			mass = masses(i)
			matricies_list(i)%elements(:, :) = calculator(grids(start_ind:end_ind), sizes(i), mass)
			Tslice(:, :) = kronecker_product(matricies_list(1)%elements, matricies_list(2)%elements)

			do j = 3, dim
				new_size = product(sizes(:j))
				temp = kronecker_product(Tslice, matricies_list(j)%elements)
				deallocate(Tslice)
				allocate(Tslice(new_size, new_size))
				Tslice = temp
				deallocate(temp)
			end do

			Tmat(:, :) = Tmat(:, :) + Tslice(:, :)
			deallocate(Tslice)
			!call deallocate_n_matricies(matricies_list, dim)
			allocate(Tslice(init_size, init_size))

			start_ind = start_ind + sizes(i)

		end do

		deallocate(Tslice)

	end subroutine kinetic_matrix_nd


	function init_n_matricies(sizes, dim) result(matricies_list)

		use types
		implicit none
		integer, intent(in) :: dim
		integer, intent(in) :: sizes(dim)
		type(list), allocatable :: matricies_list(:)
		integer :: i, nr

		allocate(matricies_list(dim))
		do i = 1, dim
			nr = sizes(i)
			allocate(matricies_list(i)%elements(nr, nr))
			matricies_list(i)%elements(:, :) = eye(nr)
		end do

	end function init_n_matricies


	subroutine deallocate_n_matricies(matricies_list, dim)
		use types
		implicit none
		integer, intent(in) :: dim
		type(list), allocatable, intent(inout) :: matricies_list(:)
		integer :: i

		do i = 1, dim
			deallocate(matricies_list(i)%elements)
		end do

		deallocate(matricies_list)

	end subroutine deallocate_n_matricies


	function kronecker_product(A, B) result(C)
		use types
		implicit none

		real(dp), intent(in) :: A(:,:)
		real(dp), intent(in) :: B(:,:)
		integer :: m, n, p, q, i, j, k, l
		real(dp), allocatable :: C(:,:)

		m = size(A, 1)
		n = size(A, 2)
		p = size(B, 1)
		q = size(B, 2)
		allocate(C(m * p, n * q))
		C(:, :) = 0.0_dp

		do concurrent(i = 1:m, j = 1:n, k = 1:p, l = 1:q)
			C((i - 1) * p + k, (j - 1) * q + l) = A(i, j) * B(k, l)
		end do

	end function kronecker_product


	function eye(n) result(Imat)
		use types
		implicit none

		integer, intent(in) :: n
		integer :: i, j
		real(dp), allocatable :: Imat(:,:)

		allocate(Imat(n, n))
		Imat = 0.0_dp

		do i = 1, n
			Imat(i, i) = 1.0_dp
		end do

	end function eye


	function colbert_miller(grid, ngrid, mass) result(Tn)
		use types
		implicit none
		integer, intent(in) :: ngrid
		real(dp), intent(in) :: grid(:)
		real(dp), intent(in) :: mass
		real(dp) :: Tn(ngrid, ngrid)
		real(dp) :: dx, pi
		integer :: i, j, hbar

		hbar = 1
		pi = 3.141592653589793_dp
		dx = grid(2) - grid(1)

		do i = 1, ngrid
			do j = 1, ngrid
				if (i == j) then
					Tn(i, j) = ((hbar**2) *  pi**2) / (6 * mass * dx ** 2)
				else if (i /= j) then
					Tn(i, j) = ((hbar**2) * (-1.0_dp)**(i - j)) / (mass * dx**2 * (i - j)**2)
				end if
			end do
		end do

	end function colbert_miller




end module dvr_solvers