
!// forquill v1.01 beta www.fcode.cn
subroutine zs1s2(zr, s1, s2, nz, ascle, alim, iuf)
   use amos_utils
   implicit none
   !***begin prologue  zs1s2
   !***refer to  zbesk,zairy
   !
   !     zs1s2 tests for a possible underflow resulting from the
   !     addition of the i and k functions in the analytic con-
   !     tinuation formula where s1=k function and s2=i function.
   !     on kode=1 the i and k functions are different orders of
   !     magnitude, but for kode=2 they can be of the same order
   !     of magnitude and the maximum must be at least one
   !     precision above the underflow limit.
   !
   !***end prologue  zs1s2
   !     complex czero,c1,s1,s1d,s2,zr
   complex(8),intent(in)   ::zr
   complex(8),intent(inout)::s1,s2
   real(8),intent(in)      ::ascle,alim
   integer,intent(inout)   ::iuf, nz
   real(8)::aln,as1,as2
   complex(8)::s1d
   nz = 0
   as1=abs(s1)
   as2=abs(s2)
   loop10:do
      if (iszero(s1) ) exit
      if (iszero(as1)) exit
      aln = -2.0d0*zr%re + log(as1)
      s1d = s1
      s1  = 0.0d0
      as1 = 0.0d0
      if (aln<(-alim)) exit
      s1 =exp(log(s1d)-zr-zr)
      as1=abs(s1)
      iuf = iuf + 1
      exit
   end do loop10
   if (max(as1, as2)>ascle) return
   s1 = 0.0d0
   s2 = 0.0d0
   nz = 1
   iuf = 0
end subroutine zs1s2
