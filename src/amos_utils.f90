module amos_utils
   implicit none
   real(8),parameter::d1mach1 = 2.d0**(minexponent(1.d0)-1) ! the smallest positive magnitude.
   real(8),parameter::d1mach2 = huge(1.d0)                  ! the largest magnitude.
   real(8),parameter::d1mach3 = 2.d0**(-digits(1.d0))       ! the smallest relative spacing.
   real(8),parameter::d1mach4 = 2.d0**(1-digits(1.d0))      ! the largest relative spacing.
   real(8),parameter::d1mach5 = log10(2.d0)
   real(8),parameter::rt2=1.41421356237309505d0
   real(8),parameter::dpi=3.14159265358979324d0
   real(8),parameter::rthdpi=1.25331413731550025d0
   real(8),parameter::rtpi=0.159154943091895336d0
   real(8),parameter::spi=1.90985931710274403d0
   real(8),parameter::hpi=1.57079632679489662d0
   real(8),parameter::fpi=1.89769999331517738d0
   real(8),parameter::thpi=4.71238898038468986D+00
   real(8),parameter::tth=6.66666666666666666d-01
   integer,parameter::i1mach1=5
   integer,parameter::i1mach2=6
   integer,parameter::i1mach3=0
   integer,parameter::i1mach4=0
   integer,parameter::i1mach5=bit_size(1)
   integer,parameter::i1mach6=4
   integer,parameter::i1mach7=radix(1)
   integer,parameter::i1mach8=bit_size(1)-1
   integer,parameter::i1mach9=huge(1)
   integer,parameter::i1mach10=radix(1.0)
   integer,parameter::i1mach11=digits(1.0)
   integer,parameter::i1mach12=minexponent(1.0)
   integer,parameter::i1mach13=maxexponent(1.0)
   integer,parameter::i1mach14=digits(1.d0)
   integer,parameter::i1mach15=minexponent(1.d0)
   integer,parameter::i1mach16=maxexponent(1.d0)
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
