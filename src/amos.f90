module amos
   use amos_utils
   implicit none
contains
   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zs1s2(zr, s1, s2, nz, ascle, alim, iuf)
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
   end subroutine f90_zs1s2

   !// forquill v1.01 beta www.fcode.cn
   integer function f90_zuchk(y, ascle, tol)result(nz)
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

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zrati(z, fnu, n, cy, tol)
      !***begin prologue  zrati
      !***refer to  zbesi,zbesk,zbesh
      !
      !     zrati computes ratios of i bessel functions by backward
      !     recurrence.  the starting index is determined by forward
      !     recurrence as described in j. res. of nat. bur. of standards-b,
      !     mathematical sciences, vol 77b, p111-114, september, 1973,
      !     bessel functions i and j of complex argument and integer order,
      !     by d. j. sookne.
      !
      !***routines called  azabs,zdiv
      !***end prologue  zrati
      !     complex z,cy(1),cone,czero,p1,p2,t1,rz,pt,cdfnu
      integer,intent(in)       ::n
      complex(8),intent(in)    ::z
      complex(8),intent(inout) ::cy(n)
      real(8),intent(in)       ::fnu,tol
      real(8) ak, amagz, ap1, ap2, az,dfnu,&
         rak, rap1, rho, test, test1,ptr
      integer i, id, idnu, inu, itime, k, kk, magz
      complex(8)::t1,p1,p2,pt,rz,cdfnu
      az=abs(z)
      inu = int(fnu)
      idnu = inu + n - 1
      magz = int(az)
      amagz = real(magz+1,8)
      id = idnu - magz - 1
      itime = 1
      k = 1
      ptr = 1.0d0/az
      rz = 2*ptr*conjg(z)
      t1 = rz*max(amagz,real(idnu,8))
      p2 = -t1
      p1 = 1.0d0
      t1 = t1 + rz
      if (id>0) id = 0
      ap2=abs(p2)
      ap1=abs(p1)
      !-----------------------------------------------------------------------
      !     the overflow test on k(fnu+i-1,z) before the call to cbknu
      !     guarantees that p2 is on scale. scale test1 and all subsequent
      !     p2 values by ap1 to ensure that an overflow does not occur
      !     prematurely.
      !-----------------------------------------------------------------------
      test1 = sqrt((ap2+ap2)/(ap1*tol))
      test = test1
      rap1 = 1.0d0/ap1
      p1 = p1*rap1
      p2 = p2*rap1
      ap2 = ap2*rap1
      10 continue
      k = k + 1
      ap1 = ap2
      pt = p2
      p2 = p1 - (t1*pt)
      p1 = pt
      t1 = t1 + rz
      ap2=abs(p2)
      if (ap1<=test) goto 10
      if (itime==2) goto 20
      ak = abs(t1)*0.5d0
      rho = min(ap2/ap1, ak + sqrt(ak*ak-1.0d0))
      test = test1*sqrt(rho/(rho*rho-1.0d0))
      itime = 2
      goto 10
      20 continue
      kk = k + 1 - id
      t1 = real(kk,8)
      dfnu = fnu + real(n-1,8)
      p1 = 1.0d0/ap2
      p2 = 0.0d0
      do i = 1, kk
      pt = p1
      p1 = pt*rz*(dfnu + t1%re)+p2
      p2 = pt
      t1 = t1 - 1.0d0
      end do
      if (iszero(p1))then
         p1 = cmplx(tol,tol,8)
      end if
      cy(n)=p2/p1
      if (n==1) return
      k = n - 1
      ak = real(k,8)
      t1= ak
      cdfnu = fnu*rz
      do i = 2, n
      pt = cdfnu + (t1*rz) + cy(k+1)
      ak = abs(pt)
      if (iszero(ak))then
         pt=cmplx(tol,tol,8)
         ak = tol*rt2
      end if
      rak = 1.0d0/ak
      cy(k)=rak*rak*conjg(pt)
      t1 = t1 - 1.0d0
      k = k - 1
      end do
   end subroutine f90_zrati
end module amos
