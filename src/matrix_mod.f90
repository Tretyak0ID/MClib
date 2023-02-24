module matrix_mod
implicit none

	type, public :: matrix_t

		real(kind=8), allocatable :: content(:, :)
		integer(kind=8)           :: N, M

	contains

		procedure, public :: zeros
		procedure, public :: ones
		procedure, public :: print

		procedure, public :: transpose

	end type matrix_t

contains
	
!------------------------BASE-OPERATIONS------------------------!
	subroutine zeros(this, N, M)
		class(matrix_t), intent(out) :: this
		integer(kind=8), intent(in)  :: N, M
		integer(kind=8)              :: i, j

		allocate(this%content(0 : N - 1, 0: M - 1))

		this%N = N
		this%M = M

		do j = 0, M - 1
		do i = 0, N  - 1
			this%content(i, j) = 0.0_8
		end do
		end do
	end subroutine zeros

	subroutine ones(this, N, M)
		class(matrix_t), intent(out) :: this
		integer(kind=8), intent(in)  :: N, M
		integer(kind=8)              :: i, j

		allocate(this%content(0 : N - 1, 0: M - 1))

		this%N = N
		this%M = M

		do j = 0, M - 1
		do i = 0, N  - 1
			this%content(i, j) = 1.0_8
		end do
		end do
	end subroutine ones

	function copy(m1) result(m)
		class(matrix_t), intent(in)  :: m1
		class(matrix_t), allocatable :: m
		integer(kind=8) :: i, j

		allocate(m)
		call m%zeros(m1%N, m1%M)

		do j = 0, m%M - 1
		do i = 0, m%N  - 1
			m%content(i, j) = m1%content(i, j)
		end do
		end do
	end function copy

	subroutine transpose(this)
		class(matrix_t), intent(inout) :: this
		integer(kind=8) :: i, j
		real(kind=8)    :: buff
	end subroutine transpose

	subroutine print(this)
		class(matrix_t), intent(inout) :: this
		integer(kind=8) :: i
		real(kind=8) :: print_string(this%M)

		do i = 0, this%N - 1
			print_string = this%content(i, :)
			print *, print_string
		end do
	end subroutine print

!------------------------ELEMETRARY-MATRIX-AND-VECTORS-OPERATIONS------------------------!
	function add(m1, m2) result(m)
		class(matrix_t),  intent(in)  :: m1, m2
		class(matrix_t), allocatable  :: m
		integer(kind=8)               :: i, j

		if (m1%N .eq. m2%N .and. m1%M .eq. m2%M) then
			allocate(m)
			call m%zeros(m1%N, m1%M)
			do j = 0, m%M - 1
				do i = 0, m%N  - 1
					m%content(i, j) = m1%content(i, j) + m2%content(i, j)
				end do
			end do
		else
			print *, 'Error: sizes of the matrices must match'
		end if
	end function add

	function mults(scalar1, m1) result(m)
		class(matrix_t), intent(in)  :: m1
		real(kind=8),    intent(in)  :: scalar1
		class(matrix_t), allocatable :: m
		integer(kind=8)              :: i, j

		allocate(m)
		call m%zeros(m1%N, m1%M)
		do j = 0, m%M - 1
			do i = 0, m%N  - 1
				m%content(i, j) = scalar1 * m1%content(i, j)
			end do
		end do
	end function mults

	function mult(m1, m2) result(m)
		!direct matrix multiplication (O(nmr) operations (O(n^3)))
		class(matrix_t),  intent(in)  :: m1, m2
		class(matrix_t), allocatable  :: m
		integer(kind=8)               :: i, j, k

		if(m1%M == m2%N) then
			allocate(m)
			call m%zeros(m1%N, m2%M)
			do j = 0, m%M - 1
				do k = 0, m1%M - 1
					do i = 0, m%N - 1
						m%content(i, j) = m%content(i, j) + m1%content(i, k) * m2%content(k, j)
					end do
				end do
			end do
		else
			print *, 'Error: sizes of the matrices must match'
		end if
	end function mult

	function dot(x, y) result(c)
		!dot-product (O(n) operations)
		real(kind=8),  allocatable :: c
		type(matrix_t), intent(in) :: x, y
		integer(kind=8) :: k

		c = 0.0_8

		if (x%M .eq. 1 .and. y%M .eq. 1) then
			do k = 0, x%N - 1
				c = c + x%content(k, 0) * y%content(k, 0)
			end do
		else
			print *, 'Error: arguments must be vector-columns'
		end if
	end function dot

	function saxpy(alpha, x, y) result(z)
		!saxpy = scalar alpha x plus y (O(n) operations)
		class(matrix_t), allocatable :: z
		type(matrix_t), intent(in)   :: x, y
		real(kind=8),   intent(in)   :: alpha
		integer(kind=8) :: k

		if (x%M .eq. 1 .and. y%M .eq. 1) then
			allocate(z)
			call z%zeros(x%N, x%M)

			do k = 0, x%N - 1
				z%content(k,0) = (alpha * x%content(k, 0)  +  y%content(k, 0))
			end do
		else
			print *, 'Error: arguments must be vector-columns'
		end if
	end function saxpy

	function gaxpy(A, x, y) result(z)
		!gaxpy = eneral Ax plus y (O(mn) operations)
		class(matrix_t), allocatable :: z
		type(matrix_t), intent(in)   :: x, y, A
		integer(kind=8) :: i, j

		if ((A%M .eq. x%N) .and. (x%M .eq. 1) .and. (y%M .eq. 1)) then
			allocate(z)
			z = copy(y)

			do j = 0, x%N - 1
				do i = 0, z%N - 1
					z%content(i, 0) = z%content(i, 0) + A%content(i, j) * x%content(j, 0)
				end do
			end do
		else
			print *, 'Error: inconsistent sizes of input arguments'
		end if
	end function gaxpy

end module matrix_mod