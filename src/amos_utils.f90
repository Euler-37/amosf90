module amos_utils
   implicit none
   real(8),parameter::d1mach1 = 2.d0**(minexponent(1.d0)-1) ! the smallest positive magnitude.
   real(8),parameter::d1mach2 = huge(1.d0)                  ! the largest magnitude.
   real(8),parameter::d1mach3 = 2.d0**(-digits(1.d0))       ! the smallest relative spacing.
   real(8),parameter::d1mach4 = 2.d0**(1-digits(1.d0))      ! the largest relative spacing.
   real(8),parameter::d1mach5 = log10(2.d0)
   interface iszero
      module procedure iszero_real8
      module procedure iszero_complex8
   end interface iszero
contains
   pure logical function iszero_real8(a)result(res)
      real(8),intent(in)::a
      res=transfer(a,1_8)==0
   end function iszero_real8

   pure logical function iszero_complex8(a)result(res)
      complex(8),intent(in)::a
      res=transfer(a%re,1_8)==0 .and. transfer(a%im,1_8)==0
   end function iszero_complex8
end module amos_utils
