module amos
   use amos_utils
   implicit none
contains

   !// forquill v1.01 beta www.fcode.cn
   subroutine zs1s2(zr, s1, s2, nz, ascle, alim, iuf)
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

   !// forquill v1.01 beta www.fcode.cn
   integer function zuchk(y, ascle, tol)result(nz)
      !***begin prologue  zuchk
      !***refer to zseri,zuoik,zunk1,zunk2,zuni1,zuni2,zkscl
      !
      !      y enters as a scaled quantity whose magnitude is greater than
      !      exp(-alim)=ascle=1.0e+3*d1mach(1)/tol. the test is made to see
      !      if the magnitude of the real or imaginary part would underflow
      !      when y is scaled (by tol) to its proper value. y is accepted
      !      if the underflow is at least one precision below the magnitude
      !      of the largest component; otherwise the phase angle does not have
      !      absolute accuracy and an underflow is assumed.
      !
      !***routines called  (none)
      !***end prologue  zuchk
      !
      !     complex y
      complex(8),intent(in)::y
      real(8),intent(in)   ::ascle,tol
      real(8)::st, wr, wi
      nz = 0
      wr = dabs(y%re)
      wi = dabs(y%im)
      st = min(wr, wi)
      if (st>ascle) return
      st = st/tol
      if (max(wr, wi)<st) nz = 1
   end function zuchk
end module amos
