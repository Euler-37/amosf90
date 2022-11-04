!// forquill v1.01 beta www.fcode.cn
integer function f90_zuchk(y, ascle, tol)result(nz)
   implicit none
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
end function f90_zuchk
