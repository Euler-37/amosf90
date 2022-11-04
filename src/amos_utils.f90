module amos_utils
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
