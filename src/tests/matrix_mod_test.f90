program matrix_mod_test
use matrix_mod

	type(matrix_t) :: m1, m2, m3, m4
	type(matrix_t) :: v1, v2, v3, v4
	real(kind=8)   :: a, b, c

	!direct-mult
	call m1%ones(2, 2)
	call m2%ones(2, 2)
	m1%content(1, 0) = 2.0_8
	m1%content(1, 1) = 3.0_8
	m3 = mult(m1, m2)
	call m3%print

	!dot-product
	call v1%zeros(3, 1)
	call v2%zeros(3, 1)
	v1%content(0,0) = 1.0_8
	v1%content(1,0) = 2.0_8
	v1%content(2,0) = 3.0_8
	v2%content(0,0) = 3.0_8
	v2%content(1,0) = 2.0_8
	v2%content(2,0) = 1.0_8
	c = dot(v1, v2)
	print *, c

	!saxpy
	v3 = saxpy(2.0_8, v1, v2)
	call v3%print

	!gaxpy
	call m4%ones(3, 3)
	m4%content(1, 2) = 10.0_8 
	v4 = gaxpy(m4, v1, v2)
	call v4%print

end program matrix_mod_test