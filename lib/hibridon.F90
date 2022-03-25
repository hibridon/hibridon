module mod_hibridon
	use mod_basis, only:ab_basis
	type :: params_type
		class(ab_basis), allocatable  :: basis
	contains
	end type params_type
 end module mod_hibridon