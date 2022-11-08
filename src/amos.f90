module amos
   use amos_utils
   implicit none
   public::f90_zbiry,f90_zbesi,f90_zbesj,f90_zbesk
   public::f90_zbesh,f90_zbesy
   private
   !f90_zs1s2,f90_zuchk,f90_zrati,f90_zkscl
   !f90_zasyi,f90_zseri,f90_zunik,f90_zbknu
   !f90_zunhj,f90_zuoik,f90_zwrsk,f90_zmlri
   !f90_zacai,f90_zairy,f90_zuni1,f90_zuni2
   !f90_zbuni,f90_zbinu,f90_zacon,f90_zunk1
   !f90_zunk2,f90_zbunk
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
      !DEBUG write(*,*)"call f90_zs1s2"
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
      !DEBUG write(*,*)"call f90_zrati"
      az=abs(z)
      inu = int(fnu)
      idnu = inu + n - 1
      magz = int(az)
      amagz = real(magz+1,8)
      id = idnu - magz - 1
      itime = 1
      k = 1
      ptr = 1.0d0/az
      rz = 2*ptr*ptr*conjg(z)
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
            pt =cmplx(tol,tol,8)
            ak = tol*rt2
         end if
         rak = 1.0d0/ak
         cy(k)=rak*rak*conjg(pt)
         t1 = t1 - 1.0d0
         k = k - 1
      end do
   end subroutine f90_zrati

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zkscl(zr, fnu, n, y, nz, rz, ascle, tol, elim)
      !***begin prologue  zkscl
      !***refer to  zbesk
      !
      !     set k functions to zero on underflow, continue recurrence
      !     on scaled functions until two members come on scale, then
      !     return with min(nz+2,n) values scaled by 1/tol.
      !
      !***routines called  zuchk,azabs,azlog
      !***end prologue  zkscl
      !     complex ck,cs,cy,czero,rz,s1,s2,y,zr,zd,celm
      complex(8),intent(in)::zr
      real(8),intent(in)::fnu,ascle,tol,elim
      integer,intent(inout)::nz
      integer,intent(in)::n
      complex(8),intent(inout)::y(n),rz
      real(8) acs, as ,&
         fn,&
         celmr, elm, helim, alas
      integer i, ic, kk,  nn, nw
      complex(8)::cy(2),s1,s2,cs,ck,zd
      !
      nz = 0
      ic = 0
      nn = min(2, n)
      do 10 i = 1, nn
         s1 = y(i)
         cy(i) = s1
         as = abs(s1)
         acs = -zr%re + log(as)
         nz = nz + 1
         y(i) = 0.0d0
         if (acs<(-elim)) goto 10
         cs=exp(log(s1)-zr)/tol
         nw=f90_zuchk(cs,ascle, tol)
         if (nw/=0) goto 10
         y(i) = cs
         ic = i
         nz = nz - 1
      10 end do
      if (n==1) return
      if (ic>1) goto 20
      y(1) = 0.0d0
      nz = 2
      20 continue
      if (n==2) return
      if (nz==0) return
      fn = fnu + 1.0d0
      ck=fn*rz
      s1 = cy(1)
      s2 = cy(2)
      helim = 0.5d0*elim
      elm = exp(-elim)
      celmr = elm
      zd = zr
      !
      !     find two consecutive y values on scale. scale recurrence if
      !     s2 gets larger than exp(elim/2)
      !
      do 30 i = 3, n
         kk = i
         cs=s2
         s2=ck*cs+s1
         s1=cs
         ck=ck+rz
         as = abs(s2)
         alas = log(as)
         acs = -zd%re + alas
         nz = nz + 1
         y(i) = 0.0d0
         if (acs<(-elim)) goto 25
         cs=exp(log(s2)-zd)/tol
         nw=f90_zuchk(cs, ascle, tol)
         if (nw/=0) goto 25
         y(i) = cs
         nz = nz - 1
         if (ic==kk-1) goto 40
         ic = kk
         goto 30
         25 continue
         if (alas<helim) goto 30
         zd = zd - elim
         s1=s1*celmr
         s2=s2*celmr
      30 end do
      nz = n
      if (ic==n) nz = n - 1
      goto 45
      40 continue
      nz = kk - 2
      45 continue
      do i = 1, nz
         y(i) = 0.0d0
      end do
      return
   end subroutine f90_zkscl


   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zseri(z, fnu, kode, n, y, nz, tol, elim, alim)
      !***begin prologue  zseri
      !***refer to  zbesi,zbesk
      !
      !     zseri computes the i bessel function for real(z)>=0.0 by
      !     means of the power series for large cabs(z) in the
      !     region cabs(z)<=2*sqrt(fnu+1). nz=0 is a normal return.
      !     nz>0 means that the last nz components were set to zero
      !     due to underflow. nz<0 means underflow occurred, but the
      !     condition cabs(z)<=2*sqrt(fnu+1) was violated and the
      !     computation must be completed in another routine with n=n-abs(nz).
      !
      !***routines called  dgamln,d1mach,zuchk,azabs,zdiv,azlog,zmlt
      !***end prologue  zseri
      !     complex ak1,ck,coef,cone,crsc,cscl,cz,czero,hz,rz,s1,s2,y,z
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim
      integer,intent(in)::n,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::y(n)
      real(8) aa, acz, ak, arm, ascle,&
         atol, az, &
         crscr, dfnu, fnup, &
         raz, rs, rtr1, s, ss
      integer::i, ib, iflag, il, k, l, m, nn, nw
      complex(8)::w(2),ak1,ck,coef,cz,hz,rz,s1,s2
      !DEBUG write(*,*)"call f90_zseri"
      nz = 0
      az = abs(z)
      if (iszero(az)) goto 160
      arm = 1.0d+3*d1mach1
      rtr1 = sqrt(arm)
      crscr = 1.0d0
      iflag = 0
      if (az<arm) goto 150
      hz=0.5d0*z
      cz = 0.0d0
      if (az>rtr1) cz=hz*hz
      acz =abs(cz)
      nn = n
      ck=log(hz)
      20 continue
      dfnu = fnu + real(nn-1,8)
      fnup = dfnu + 1.0d0
      !-----------------------------------------------------------------------
      !     underflow test
      !-----------------------------------------------------------------------
      ak1 = (ck*dfnu) - log_gamma(fnup)
      if (kode==2) ak1%re = ak1%re - z%re
      if (ak1%re>(-elim)) goto 40
      30 continue
      nz = nz + 1
      y(nn) = 0.0d0
      if (acz>dfnu)goto 190
      nn = nn - 1
      if (nn==0) return
      goto 20
      40 continue
      if (ak1%re<(-alim))then
         iflag = 1
         ss = 1.0d0/tol
         crscr = tol
         ascle = arm*ss
      end if
      aa = exp(ak1%re)
      if (iflag==1) aa = aa*ss
      coef=exp(ak1)
      atol = tol*acz/fnup
      il = min(2, nn)
      loop90:do i = 1, il
         dfnu = fnu + real(nn-i,8)
         fnup = dfnu + 1.0d0
         s1 = 1.0d0
         if (acz>tol*fnup)then
            ak1 = 1.0d0
            ak = fnup + 2.0d0
            s = fnup
            aa = 2.0d0
            loop60:do
               rs = 1.0d0/s
               ak1=(ak1*cz)*rs
               s1=s1+ak1
               s = s + ak
               ak = ak + 2.0d0
               aa = aa*acz*rs
               if (aa<atol)exit loop60
            end do loop60
         end if
         s2 = s1*coef
         w(i) = s2
         if (iflag/=0)then
            nw=f90_zuchk(s2,ascle,tol)
            if (nw/=0) goto 30
         end if
         m = nn - i + 1
         y(m) = s2*crscr
         if (i==il) cycle loop90
         coef=coef/hz*dfnu
      end do loop90
      if (nn<=2) return
      k = nn - 2
      ak = real(k,8)
      raz = 1.0d0/az
      rz=2*(conjg(z)*raz)*raz
      if (iflag==1)goto 120
      ib = 3
      100 continue
      do i = ib, nn
         y(k) = (ak+fnu)*(rz*y(k+1)) + y(k+2)
         ak = ak - 1.0d0
         k = k - 1
      end do
      return
      !-----------------------------------------------------------------------
      !     RECUR BACKWARD WITH SCALED VALUES
      !-----------------------------------------------------------------------
      120 continue
      !-----------------------------------------------------------------------
      !     exp(-alim)=exp(-elim)/tol=approx. one precision above the
      !     underflow limit = ascle = d1mach(1)*ss*1.0d+3
      !-----------------------------------------------------------------------
      s1 = w(1)
      s2 = w(2)
      do l = 3, nn
         ck = s2
         s2 = s1 + (ak+fnu)*(rz*ck)
         s1 = ck
         ck = s2*crscr
         y(k) = ck
         ak = ak - 1.0d0
         k = k - 1
         if (abs(ck)>ascle) goto 140
      end do
      return
      140 continue

      ib = l + 1
      if (ib>nn) return
      goto 100
      150 continue
      nz=n
      do i = ib, nn
         y(k) = (ak+fnu)*(rz*y(k+1)) + y(k+2)
         ak = ak - 1.0d0
         k = k - 1
      end do
      return
      nz = n
      y = 0.0d0
      if (iszero(fnu)) nz = nz - 1
      160 continue
      y(1) = 0.d0
      if (iszero(fnu))then
         y(1) = 1.d0
      end if
      if (n==1) return
      do i = 2, n
         y(i) = 0.d0
      end do
      return
      !-----------------------------------------------------------------------
      !     return with nz<0 if cabs(z*z/4)>fnu+n-nz-1 complete
      !     the calculation in cbinu with n=n-iabs(nz)
      !-----------------------------------------------------------------------
      190 continue
      nz = -nz
      return
   end subroutine f90_zseri

   !// ForQuill v1.01 Beta www.fcode.cn
   subroutine f90_zasyi(z,  fnu, kode, n, y, nz, rl, tol, elim, alim)
      !***begin prologue  zasyi
      !***refer to  zbesi,zbesk
      !
      !     zasyi computes the i bessel function for real(z)>=0.0 by
      !     means of the asymptotic expansion for large cabs(z) in the
      !     region cabs(z)>max(rl,fnu*fnu/2). nz=0 is a normal return.
      !     nz<0 indicates an overflow on kode=1.
      !
      !***routines called  d1mach,azabs,zdiv,azexp,zmlt,azsqrt
      !***end prologue  zasyi
      !     complex ak1,ck,cone,cs1,cs2,cz,czero,dk,ez,p1,rz,s2,y,z
      integer,intent(in)::n
      integer,intent(inout)::nz
      real(8),intent(in)::fnu,rl,tol,elim,alim
      complex(8),intent(in)::z
      complex(8),intent(inout) ::y(n)
      integer i, ib, il, inu, j, jl, k, kode, koded, m, nn
      !
      complex(8) ak1,ck,cs1,cs2,cz,dk,ez,p1,rz,s2
      real(8)::az,arm,rtr1,dfnu,raz,dnu2,fdn,aez,s,arg,ak,bk
      real(8)::aa,atol,bb,sgn,sqk
      !DEBUG write(*,*)'call f90_zasyi'
      nz = 0
      az = abs(z)
      arm = 1.0d+3*d1mach1
      rtr1 = sqrt(arm)
      il = min(2, n)
      dfnu = fnu + real(n-il,8)
      !-----------------------------------------------------------------------
      !     overflow test
      !-----------------------------------------------------------------------
      raz = 1.0d0/az
      ak1=sqrt(rtpi*conjg(z)*raz*raz)
      cz = z
      if (kode/=2) goto 10
      cz=cmplx(0.0d0,z%im,8)
      10 continue
      if (abs(cz%re)>elim) goto 100
      dnu2 = dfnu + dfnu
      koded = 1
      if ((abs(cz%re)>alim) .and. (n>2)) goto 20
      koded = 0
      ak1=ak1*exp(cz)
      20 continue
      fdn = 0.0d0
      if (dnu2>rtr1) fdn = dnu2*dnu2
      ez = z*8.0d0
      !-----------------------------------------------------------------------
      !     when z is imaginary, the error test must be made relative to the
      !     first reciprocal power since this is the leading term of the
      !     expansion for the imaginary part.
      !-----------------------------------------------------------------------
      aez = 8.0d0*az
      s = tol/aez
      jl = int(rl+rl) + 2
      p1 = 0.0d0
      if (iszero(z%im)) goto 30
      !-----------------------------------------------------------------------
      !     calculate exp(pi*(0.5+fnu+n-il)*i) to minimize losses of
      !     significance when fnu or n is large
      !-----------------------------------------------------------------------
      inu = int(fnu)
      arg = (fnu-real(inu))*dpi
      inu = inu + n - il
      ak = -sin(arg)
      bk =  cos(arg)
      if (z%im<0.0d0) bk = -bk
      p1=cmplx(ak,bk,8)
      if (mod(inu,2)==0) goto 30
      p1=-p1
      30 continue
      do k = 1, il
         sqk = fdn - 1.0d0
         atol = s*dabs(sqk)
         sgn = 1.0d0
         cs1 = 1.0d0
         cs2 = 1.0d0
         ck = 1.0d0
         ak = 0.0d0
         aa = 1.0d0
         bb = aez
         dk = ez
         do j = 1, jl
            ck = (ck/dk)*sqk
            cs2 = cs2 + ck
            sgn = -sgn
            cs1 = cs1 + ck*sgn
            dk = dk + ez
            aa = aa*dabs(sqk)/bb
            bb = bb + aez
            ak = ak + 8.0d0
            sqk = sqk - ak
            if (aa<=atol) goto 50
         end do
         goto 110
         50 continue
         s2 = cs1
         if (z%re+z%re>=elim) goto 60
         s2=s2+exp(-2*z)*p1*cs2
         60 continue
         fdn = fdn + 8.0d0*dfnu + 4.0d0
         p1=-p1
         m = n - il + k
         y(m) = s2*ak1 
      end do
      if (n<=2) return
      nn = n
      k = nn - 2
      ak = real(k,8)
      rz=2*conjg(z)*raz*raz
      ib = 3
      do i = ib, nn
         y(k) = (ak+fnu)*(rz*y(k+1)) + y(k+2)
         ak = ak - 1.0d0
         k = k - 1
      end do
      if (koded==0) return
      ck=exp(cz)
      do i = 1, nn
         y=y*ck
      end do
      return
      100 continue
      nz = -1
      return
      110 continue
      nz = -2
      return
   end subroutine f90_zasyi

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zunik(zr, fnu, ikflg, ipmtr, tol, init, phi, zeta1, zeta2, sum1, cwrk)
      !***begin prologue  zunik
      !***refer to  zbesi,zbesk
      !
      !        zunik computes parameters for the uniform asymptotic
      !        expansions of the i and k functions on ikflg= 1 or 2
      !        respectively by
      !
      !        w(fnu,zr) = phi*exp(zeta)*sum
      !
      !        where       zeta=-zeta1 + zeta2       or
      !                          zeta1 - zeta2
      !
      !        the first call must have init=0. subsequent calls with the
      !        same zr and fnu will return the i or k function on ikflg=
      !        1 or 2 with no change in init. cwrk is a complex work
      !        array. ipmtr=0 computes all parameters. ipmtr=1 computes phi,
      !        zeta1,zeta2.
      !
      !***routines called  zdiv,azlog,azsqrt,d1mach
      !***end prologue  zunik
      !     complex cfn,con,cone,crfn,cwrk,czero,phi,s,sr,sum,t,t2,zeta1,
      !    *zeta2,zn,zr
      complex(8),intent(in)::zr
      integer,intent(in)::ikflg,ipmtr
      integer,intent(inout)::init
      complex(8),intent(inout)::phi,zeta1,zeta2,sum1,cwrk(16)
      real(8),intent(in)::fnu,tol
      complex(8)::crfn,s,sr,t,t2,st
      real(8)::ac,rfn,tr,test
      integer::i,j,l,k
      real(8),parameter::con(2)=[3.98942280401432678D-01,1.25331413731550025D+00]
      real(8),parameter::c(120)=[&
         1.00000000000000000D+00,    -2.08333333333333333D-01,&
         1.25000000000000000D-01,     3.34201388888888889D-01,&
         -4.01041666666666667D-01,     7.03125000000000000D-02,&
         -1.02581259645061728D+00,     1.84646267361111111D+00,&
         -8.91210937500000000D-01,     7.32421875000000000D-02,&
         4.66958442342624743D+00,    -1.12070026162229938D+01,&
         8.78912353515625000D+00,    -2.36408691406250000D+00,&
         1.12152099609375000D-01,    -2.82120725582002449D+01,&
         8.46362176746007346D+01,    -9.18182415432400174D+01,&
         4.25349987453884549D+01,    -7.36879435947963170D+00,&
         2.27108001708984375D-01,     2.12570130039217123D+02,&
         -7.65252468141181642D+02,     1.05999045252799988D+03,&
         -6.99579627376132541D+02,     2.18190511744211590D+02,&
         -2.64914304869515555D+01,     5.72501420974731445D-01,&
         -1.91945766231840700D+03,     8.06172218173730938D+03,&
         -1.35865500064341374D+04,     1.16553933368645332D+04,&
         -5.30564697861340311D+03,     1.20090291321635246D+03,&
         -1.08090919788394656D+02,     1.72772750258445740D+00,&
         2.02042913309661486D+04,    -9.69805983886375135D+04,&
         1.92547001232531532D+05,    -2.03400177280415534D+05,&
         1.22200464983017460D+05,    -4.11926549688975513D+04,&
         7.10951430248936372D+03,    -4.93915304773088012D+02,&
         6.07404200127348304D+00,    -2.42919187900551333D+05,&
         1.31176361466297720D+06,    -2.99801591853810675D+06,&
         3.76327129765640400D+06,    -2.81356322658653411D+06,&
         1.26836527332162478D+06,    -3.31645172484563578D+05,&
         4.52187689813627263D+04,    -2.49983048181120962D+03,&
         2.43805296995560639D+01,     3.28446985307203782D+06,&
         -1.97068191184322269D+07,     5.09526024926646422D+07,&
         -7.41051482115326577D+07,     6.63445122747290267D+07,&
         -3.75671766607633513D+07,     1.32887671664218183D+07,&
         -2.78561812808645469D+06,     3.08186404612662398D+05,&
         -1.38860897537170405D+04,     1.10017140269246738D+02,&
         -4.93292536645099620D+07,     3.25573074185765749D+08,&
         -9.39462359681578403D+08,     1.55359689957058006D+09,&
         -1.62108055210833708D+09,     1.10684281682301447D+09,&
         -4.95889784275030309D+08,     1.42062907797533095D+08,&
         -2.44740627257387285D+07,     2.24376817792244943D+06,&
         -8.40054336030240853D+04,     5.51335896122020586D+02,&
         8.14789096118312115D+08,    -5.86648149205184723D+09,&
         1.86882075092958249D+10,    -3.46320433881587779D+10,&
         4.12801855797539740D+10,    -3.30265997498007231D+10,&
         1.79542137311556001D+10,    -6.56329379261928433D+09,&
         1.55927986487925751D+09,    -2.25105661889415278D+08,&
         1.73951075539781645D+07,    -5.49842327572288687D+05,&
         3.03809051092238427D+03,    -1.46792612476956167D+10,&
         1.14498237732025810D+11,    -3.99096175224466498D+11,&
         8.19218669548577329D+11,    -1.09837515608122331D+12,&
         1.00815810686538209D+12,    -6.45364869245376503D+11,&
         2.87900649906150589D+11,    -8.78670721780232657D+10,&
         1.76347306068349694D+10,    -2.16716498322379509D+09,&
         1.43157876718888981D+08,    -3.87183344257261262D+06,&
         1.82577554742931747D+04,     2.86464035717679043D+11,&
         -2.40629790002850396D+12,     9.10934118523989896D+12,&
         -2.05168994109344374D+13,     3.05651255199353206D+13,&
         -3.16670885847851584D+13,     2.33483640445818409D+13,&
         -1.23204913055982872D+13,     4.61272578084913197D+12,&
         -1.19655288019618160D+12,     2.05914503232410016D+11,&
         -2.18229277575292237D+10,     1.24700929351271032D+09,&
         -2.91883881222208134D+07,     1.18838426256783253D+05 &
         ]
      !DEBUG write(*,*)"call f90_zunik"
      if (init/=0) goto 40
      !-----------------------------------------------------------------------
      !     initialize all variables
      !-----------------------------------------------------------------------
      rfn = 1.0d0/fnu
      !-----------------------------------------------------------------------
      !     overflow test (zr/fnu too small)
      !-----------------------------------------------------------------------
      test = d1mach1*1.0d+3
      ac = fnu*test
      if (abs(zr%re)>ac .or. abs(zr%im)>ac) goto 15
      zeta1 = 2.0d0*dabs(log(test)) + fnu
      zeta2 = fnu
      phi = 1.0d0
      return
      15 continue
      t = zr*rfn
      sr=sqrt(1.0d0+t*t)
      st=log((1+sr)/t)
      zeta1 = fnu*st
      zeta2 = fnu*sr
      t=1/sr
      sr = t*rfn
      cwrk(16)=sqrt(sr)
      phi = cwrk(16)*con(ikflg)
      if (ipmtr/=0) return
      t2=1/s
      cwrk(1) = 1.0d0
      crfn = 1.0d0
      ac = 1.0d0
      l = 1
      do k = 2, 15
         s = 0.0d0
         do j = 1, k
            l = l + 1
            s=t2*s+c(l)
         end do
         crfn=crfn*sr
         cwrk=crfn*s
         ac = ac*rfn
         test = abs(cwrk(k)%re) + dabs(cwrk(k)%im)
         if (ac<tol .and. test<tol) goto 30
      end do
      k = 15
      30 continue
      init = k
      40 continue
      if (ikflg==2) goto 60
      !-----------------------------------------------------------------------
      !     compute sum for the i function
      !-----------------------------------------------------------------------
      s = 0.0d0
      do i = 1, init
         s = s + cwrk(i)
      end do
      sum1 = s
      phi = cwrk(16)*con(1)
      return
      60 continue
      !-----------------------------------------------------------------------
      !     compute sum for the k function
      !-----------------------------------------------------------------------
      s = 0.0d0
      tr = 1.0d0
      do i = 1, init
         s = s + tr*cwrk(i)
         tr = -tr
      end do
      sum1 = s
      phi = cwrk(16)*con(2)
   end subroutine f90_zunik

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbknu(z, fnu, kode, n, y, nz, tol, elim, alim)
      !***begin prologue  zbknu
      !***refer to  zbesi,zbesk,zairy,zbesh
      !
      !     zbknu computes the k bessel function in the right half z plane.
      !
      !***routines called  dgamln,i1mach,d1mach,zkscl,zshch,zuchk,azabs,zdiv,
      !                    azexp,azlog,zmlt,azsqrt
      !***end prologue  zbknu
      !
      complex(8),intent(in)::z
      integer,intent(in)::n,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::y(n)
      real(8),intent(in)::fnu,tol,elim,alim

      real(8) aa, ak, ascle, a1, a2, bb, bk, caz, &
          ckr , crscr, csclr,cbr, dnu, dnu2, etest, fc, fhs,&
          fk, fks,g1, g2, ptr, p1r , p2m, p2r, rak, rcaz, &
         s, smur, str,tm, t1, t2,zr, elm, celmr, as, alas, helim
      integer i, iflag, inu, k, kflag, kk, koded, j, ic, inub, nw
      complex(8):: rz,smu,fmu,f,cz,s1,s2,csh,cch
      complex(8):: ck,p,q,coef,p1,p2,pt,st,cs,zd
      complex(8)::cy(2),cb
      real(8)::bry(3),cssr(3),csrr(3)
      !
      integer,parameter::kmax=30
      real(8),parameter::cc(8)=[&
          5.77215664901532861D-01,    -4.20026350340952355D-02,&
         -4.21977345555443367D-02,     7.21894324666309954D-03,&
         -2.15241674114950973D-04,    -2.01348547807882387D-05,&
          1.13302723198169588D-06,     6.11609510448141582D-09&
         ]
      !DEBUG write (*, *) 'call f90_zbknu'
      caz=abs(z)
      csclr = 1.0d0/tol
      crscr = tol
      cssr(1) = csclr
      cssr(2) = 1.0d0
      cssr(3) = crscr
      csrr(1) = crscr
      csrr(2) = 1.0d0
      csrr(3) = csclr
      bry(1) = 1.0d+3*d1mach1/tol
      bry(2) = 1.0d0/bry(1)
      bry(3) = d1mach2
      nz = 0
      iflag = 0
      koded = kode
      rcaz = 1.0d0/caz
      rz=2*(conjg(z)*rcaz)*rcaz
      inu = nint(fnu)
      dnu = fnu - real(inu,8)
      if (iszero(abs(dnu)-0.5d0)) goto 110
      dnu2 = 0.0d0
      if (abs(dnu)>tol) dnu2 = dnu*dnu
      if (caz>2.0d0) goto 110
      !-----------------------------------------------------------------------
      !     series for cabs(z)<=2.0d0
      !-----------------------------------------------------------------------
      fc = 1.0d0
      smu=log(rz)
      fmu = smu*dnu
      csh=sinh(fmu)
      cch=cosh(fmu)
      if (iszero(dnu)) goto 10
      fc = dnu*dpi
      fc = fc/sin(fc)
      smu=csh/dnu
      10 continue
      a2 = 1.0d0 + dnu
      !-----------------------------------------------------------------------
      !     gam(1-z)*gam(1+z)=pi*z/sin(pi*z), t1=1/gam(1-dnu), t2=1/gam(1+dnu)
      !-----------------------------------------------------------------------
      t2=exp(-log_gamma(a2))
      t1 = 1.0d0/(t2*fc)
      if (abs(dnu)>0.1d0) goto 40
      !-----------------------------------------------------------------------
      !     series for f0 to resolve indeterminacy for small abs(dnu)
      !-----------------------------------------------------------------------
      ak = 1.0d0
      s = cc(1)
      do k = 2, 8
      ak = ak*dnu2
      tm = cc(k)*ak
      s = s + tm
      if (abs(tm)<tol) goto 30
      end do
      30 g1 = -s
      goto 50
      40 continue
      g1 = (t1-t2)/(dnu+dnu)
      50 continue
      g2 = (t1+t2)*0.5d0
      f=cch*g1+smu*g2
      st=exp(fmu)
      p=0.5d0*st/t2
      pt=0.5d0/st
      q=pt/t1
      s1=f
      s2=p
      ak = 1.0d0
      a1 = 1.0d0
      ck=1.0d0
      bk = 1.0d0 - dnu2
      if (inu>0 .or. n>1) goto 80
      !-----------------------------------------------------------------------
      !     generate k(fnu,z), 0.0d0 <= fnu < 0.5d0 and n=1
      !-----------------------------------------------------------------------
      if (caz<tol) goto 70
      cz=zr*zr
      cz = 0.25d0*cz
      t1 = 0.25d0*caz*caz
      60 continue
      f=(f*ak+p+q)/bk
      str = 1.0d0/(ak-dnu)
      p=p*str
      str = 1.0d0/(ak+dnu)
      q = q*str
      rak = 1.0d0/ak
      ck=ck*cz*rak
      s1=ck*f+s1
      a1 = a1*t1*rak
      bk = bk + ak + ak + 1.0d0
      ak = ak + 1.0d0
      if (a1>tol) goto 60
      70 continue
      y(1) = s1
      if (koded==1) return
      st=exp(z)
      y(1)=s1*st
      return
      !-----------------------------------------------------------------------
      !     generate k(dnu,z) and k(dnu+1,z) for forward recurrence
      !-----------------------------------------------------------------------
      80 continue
      if (caz<tol) goto 100
      cz=zr*zr
      cz = 0.25d0*cz
      t1 = 0.25d0*caz*caz
      90 continue
      f=(f*ak+p+q)/bk
      str = 1.0d0/(ak-dnu)
      p = p*str
      str = 1.0d0/(ak+dnu)
      q = q*str
      rak = 1.0d0/ak
      ck=ck*cz*rak
      s1=ck*f+s1
      st=p-f*ak
      s2=ck*st-s2
      a1 = a1*t1*rak
      bk = bk + ak + ak + 1.0d0
      ak = ak + 1.0d0
      if (a1>tol) goto 90
      100 continue
      kflag = 2
      a1 = fnu + 1.0d0
      ak = a1*abs(smur)
      if (ak>alim) kflag = 3
      str = cssr(kflag)
      p2=s2*str
      s2=p2*rz
      s1=s1*str
      if (koded==1) goto 210
      f=exp(z)
      s1=s1*f
      s2=s2*f
      goto 210
      !-----------------------------------------------------------------------
      !     iflag=0 means no underflow occurred
      !     iflag=1 means an underflow occurred- computation proceeds with
      !     koded=2 and a test for on scale values is made during forward
      !     recursion
      !-----------------------------------------------------------------------
      110 continue
      st=sqrt(z)
      coef=rthdpi/st
      kflag = 2
      if (koded==2) goto 120
      if (zr>alim) goto 290
      !     blank line
      st=cssr(kflag)*exp(-z)
      coef=coef*st
      120 continue
      if (iszero(abs(dnu)-0.5d0)) goto 300
      !-----------------------------------------------------------------------
      !     miller algorithm for cabs(z)>2.0d0
      !-----------------------------------------------------------------------
      ak = cos(dpi*dnu)
      ak = abs(ak)
      if (iszero(ak)) goto 300
      fhs = abs(0.25d0-dnu2)
      if (iszero(fhs)) goto 300
      !-----------------------------------------------------------------------
      !     compute r2=f(e). if cabs(z)>=r2, use forward recurrence to
      !     determine the backward index k. r2=f(e) is a straight line on
      !     12<=e<=60. e is computed from 2**(-e)=b**(1-i1mach(14))=
      !     tol where b is the base of the arithmetic.
      !-----------------------------------------------------------------------
      t1 = real(i1mach14-1,8)
      t1 = t1*d1mach5*3.321928094d0
      t1 = max(t1, 12.0d0)
      t1 = min(t1, 60.0d0)
      t2 = tth*t1 - 6.0d0
      if (iszero(z%re)) then
         t1 = hpi
         goto 140
      endif
      t1 = atan(z%im/z%re)
      t1 = abs(t1)
      140 continue
      if (t2>caz) goto 170
      !-----------------------------------------------------------------------
      !     forward recurrence loop when cabs(z)>=r2
      !-----------------------------------------------------------------------
      etest = ak/(dpi*caz*tol)
      fk = 1.0d0
      if (etest<1.0d0) goto 180
      fks = 2.0d0
      ckr = caz + caz + 2.0d0
      p1r = 0.0d0
      p2r = 1.0d0
      do i = 1, kmax
      ak = fhs/fks
      cbr = ckr/(fk+1.0d0)
      ptr = p2r
      p2r = cbr*p2r - p1r*ak
      p1r = ptr
      ckr = ckr + 2.0d0
      fks = fks + fk + fk + 2.0d0
      fhs = fhs + fk + fk
      fk = fk + 1.0d0
      str = abs(p2r)*fk
      if (etest<str) goto 160
      end do
      goto 310
      160 continue
      fk = fk + spi*t1*sqrt(t2/caz)
      fhs = abs(0.25d0-dnu2)
      goto 180
      170 continue
      !-----------------------------------------------------------------------
      !     compute backward index k for cabs(z)<r2
      !-----------------------------------------------------------------------
      a2 = sqrt(caz)
      ak = fpi*ak/(tol*sqrt(a2))
      aa = 3.0d0*t1/(1.0d0+caz)
      bb = 14.7d0*t1/(28.0d0+caz)
      ak = (log(ak)+caz*cos(aa)/(1.0d0+0.008d0*caz))/cos(bb)
      fk = 0.12125d0*ak*ak/caz + 1.5d0
      180 continue
      !-----------------------------------------------------------------------
      !     backward recurrence loop for miller algorithm
      !-----------------------------------------------------------------------
      k = int(fk)
      fk = real(k,8)
      fks = fk*fk
      p1 = 0.0d0
      p2 = tol
      cs = p2
      do i = 1, k
      a1 = fks - fk
      ak = (fks+fk)/(a1+fhs)
      rak = 2.0d0/(fk+1.0d0)
      cb=(fk+z)*rak
      pt=p2
      p2=(pt*cb-p1)*ak
      p1 = pt
      cs = cs + p2
      fks = a1 - fk + 1.0d0
      fk = fk - 1.0d0
      end do
      !-----------------------------------------------------------------------
      !     compute (p2/cs)=(p2/cabs(cs))*(conjg(cs)/cabs(cs)) for better
      !     scaling
      !-----------------------------------------------------------------------
      tm=abs(cs)
      ptr = 1.0d0/tm
      s1 = p2*ptr
      cs = conjg(cs)*ptr
      st = coef * s1
      s1 = st * cs
      if (inu>0 .or. n>1) goto 200
      zd = z
      if (iflag==1) goto 270
      goto 240
      200 continue
      !-----------------------------------------------------------------------
      !     compute p1/p2=(p1/cabs(p2)*conjg(p2)/cabs(p2) for scaling
      !-----------------------------------------------------------------------
      tm = abs(p2)
      ptr = 1.0d0/tm
      p1 = p1*ptr
      p2 = conjg(p2)*ptr
      pt=p1*p2
      st=dnu+0.5d0-pt
      st=st/z
      st = st + 1.0d0
      s2=st*s1
      !-----------------------------------------------------------------------
      !     forward recursion on the three term recursion with relation with
      !     scaling near exponent extremes on kflag=1 or kflag=3
      !-----------------------------------------------------------------------
      210 continue
      str = dnu + 1.0d0
      ck = str*rz
      if (n==1) inu = inu - 1
      if (inu>0) goto 220
      if (n>1) goto 215
      s1 = s2
      215 continue
      zd = z
      if (iflag==1) goto 270
      goto 240
      220 continue
      inub = 1
      if (iflag==1) goto 261
      225 continue
      p1r = csrr(kflag)
      ascle = bry(kflag)
      do i = inub, inu
      st = s2
      s2 = ck*st+s1
      s1 = st
      ck = ck + rz
      if (kflag>=3) goto 230
      p2 = s2*p1r
      p2m = max(abs(p2%re), abs(p2%im))
      if (p2m<=ascle) goto 230
      kflag = kflag + 1
      ascle = bry(kflag)
      s1 = s1*p1r
      s2 = p2
      str = cssr(kflag)
      s1 = s1*str
      s2 = s2*str
      p1r = csrr(kflag)
      230 end do
      if (n/=1) goto 240
      s1 = s2
      240 continue
      str = csrr(kflag)
      y(1) = s1*str
      if (n==1) return
      y(2) = s2*str
      if (n==2) return
      kk = 2
      250 continue
      kk = kk + 1
      if (kk>n) return
      p1r = csrr(kflag)
      ascle = bry(kflag)
      do i = kk, n
      p2 = s2
      s2=ck*p2+s1
      s1 = p2
      ck= ck+ rz
      p2 = s2*p1r
      y(i) = p2
      if (kflag>=3) goto 260
      p2m = max(abs(p2%re), abs(p2%im))
      if (p2m<=ascle) goto 260
      kflag = kflag + 1
      ascle = bry(kflag)
      s1 = s1*p1r
      s2 = p2
      str = cssr(kflag)
      s1 = s1*str
      s2 = s2*str
      p1r = csrr(kflag)
      260 end do
      return
      !-----------------------------------------------------------------------
      !     iflag=1 cases, forward recurrence on scaled values on underflow
      !-----------------------------------------------------------------------
      261 continue
      helim = 0.5d0*elim
      elm = exp(-elim)
      celmr = elm
      ascle = bry(1)
      zd = z
      ic = -1
      j = 2
      do i = 1, inu
      st = s2
      s2=st*ck-s1
      s1 = st
      ck = ck + rz
      as = abs(s2)
      alas = log(as)
      p2r = -zd%re + alas
      if (p2r<(-elim)) goto 263
      st=log(s2)
      p2=-z+st
      p1=exp(p2)/tol
      nw=f90_zuchk(p1,ascle,tol)
      if (nw/=0) goto 263
      j = 3 - j
      cy(j) = p1
      if (ic==(i-1)) goto 264
      ic = i
      goto 262
      263 continue
      if (alas<helim) goto 262
      zd = zd - elim
      s1 = s1*celmr
      s2 = s2*celmr
      262 end do
      if (n/=1) goto 270
      s1 = s2
      goto 270
      264 continue
      kflag = 1
      inub = i + 1
      s2 = cy(j)
      j = 3 - j
      s1 = cy(j)
      if (inub<=inu) goto 225
      if (n/=1) goto 240
      s1 = s2
      goto 240
      270 continue
      y(1) = s1
      if (n==1) goto 280
      y(2) = s2
      280 continue
      ascle = bry(1)
      call f90_zkscl(zd, fnu, n, y, nz, rz, ascle, tol, elim)
      inu = n - nz
      if (inu<=0) return
      kk = nz + 1
      s1 = y(kk)
      y(kk) = s1*csrr(1)
      if (inu==1) return
      kk = nz + 2
      s2 = y(kk)
      y(kk) = s2*csrr(1)
      if (inu==2) return
      t2 = fnu + real(kk-1,8)
      ck = t2*rz
      kflag = 1
      goto 250
      290 continue
      !-----------------------------------------------------------------------
      !     scale by dexp(z), iflag = 1 cases
      !-----------------------------------------------------------------------
      koded = 2
      iflag = 1
      kflag = 2
      goto 120
      !-----------------------------------------------------------------------
      !     fnu=half odd integer case, dnu=-0.5
      !-----------------------------------------------------------------------
      300 continue
      s1 = coef
      s2 = coef
      goto 210
      !
      !
      310 continue
      nz = -2
      return
   end subroutine f90_zbknu


   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zunhj(z,fnu,ipmtr,tol,phi,arg,zeta1,zeta2,asum,bsum)
      !***begin prologue  zunhj
      !***refer to  zbesi,zbesk
      !
      !     references
      !         handbook of mathematical functions by m. abramowitz and i.a.
      !         stegun, ams55, national bureau of standards, 1965, chapter 9.
      !
      !         asymptotics and special functions by f.w.j. olver, academic
      !         press, n.y., 1974, page 420
      !
      !     abstract
      !         zunhj computes parameters for bessel functions c(fnu,z) =
      !         j(fnu,z), y(fnu,z) or h(i,fnu,z) i=1,2 for large orders fnu
      !         by means of the uniform asymptotic expansion
      !
      !         c(fnu,z)=c1*phi*( asum*airy(arg) + c2*bsum*dairy(arg) )
      !
      !         for proper choices of c1, c2, airy and dairy where airy is
      !         an airy function and dairy is its derivative.
      !
      !               (2/3)*fnu*zeta**1.5 = zeta1-zeta2,
      !
      !         zeta1=0.5*fnu*clog((1+w)/(1-w)), zeta2=fnu*w for scaling
      !         purposes in airy functions from cairy or cbiry.
      !
      !         mconj=sign of aimag(z), but is ambiguous when z is real and
      !         must be specified. ipmtr=0 returns all parameters. ipmtr=
      !         1 computes all except asum and bsum.
      !
      !***routines called  azabs,zdiv,azlog,azsqrt,d1mach
      !***end prologue  zunhj
      !     complex arg,asum,bsum,cfnu,cone,cr,czero,dr,p,phi,przth,ptfn,
      !    *rfn13,rtzta,rzth,suma,sumb,tfn,t2,up,w,w2,z,za,zb,zc,zeta,zeta1,
      !    *zeta2,zth
      complex(8),intent(in)::z
      complex(8),intent(inout)::arg,zeta1,zeta2,asum,bsum,phi
      integer,intent(in)::ipmtr
      real(8),intent(in)::tol,fnu
      real(8) :: ang, atol, aw2, azth, btol, fn13, fn23,&
         pp, raw, raw2, razth, rfnu, rfnu2, rfn13, test,ac
      integer ias, ibs, is, j, jr, ju, k, kmax, kp1, ks, &
         l, lr, lrp1, l1, l2, m
      complex(8)::w,zc,zth,zeta,suma,sumb,st,cl,przth,ptfn,rtzt
      complex(8)::rzth,t2,tfn,tza,w2,za,zb
      complex(8)::p(30),up(14),cr(14),dr(14)
      real(8)::ap(30)
      !
      real(8),parameter::ar(14)=[&
         1.00000000000000000D+00,     1.04166666666666667D-01,&
         8.35503472222222222D-02,     1.28226574556327160D-01,&
         2.91849026464140464D-01,     8.81627267443757652D-01,&
         3.32140828186276754D+00,     1.49957629868625547D+01,&
         7.89230130115865181D+01,     4.74451538868264323D+02,&
         3.20749009089066193D+03,     2.40865496408740049D+04,&
         1.98923119169509794D+05,     1.79190200777534383D+06 &
         ]
      real(8),parameter::br(14)=[&
         1.00000000000000000D+00,    -1.45833333333333333D-01,&
         -9.87413194444444444D-02,    -1.43312053915895062D-01,&
         -3.17227202678413548D-01,    -9.42429147957120249D-01,&
         -3.51120304082635426D+00,    -1.57272636203680451D+01,&
         -8.22814390971859444D+01,    -4.92355370523670524D+02,&
         -3.31621856854797251D+03,    -2.48276742452085896D+04,&
         -2.04526587315129788D+05,    -1.83844491706820990D+06 &
         ]
      real(8),parameter::c(105)=[&
         1.00000000000000000D+00,    -2.08333333333333333D-01,&
         1.25000000000000000D-01,     3.34201388888888889D-01,&
         -4.01041666666666667D-01,     7.03125000000000000D-02,&
         -1.02581259645061728D+00,     1.84646267361111111D+00,&
         -8.91210937500000000D-01,     7.32421875000000000D-02,&
         4.66958442342624743D+00,    -1.12070026162229938D+01,&
         8.78912353515625000D+00,    -2.36408691406250000D+00,&
         1.12152099609375000D-01,    -2.82120725582002449D+01,&
         8.46362176746007346D+01,    -9.18182415432400174D+01,&
         4.25349987453884549D+01,    -7.36879435947963170D+00,&
         2.27108001708984375D-01,     2.12570130039217123D+02,&
         -7.65252468141181642D+02,     1.05999045252799988D+03,&
         -6.99579627376132541D+02,     2.18190511744211590D+02,&
         -2.64914304869515555D+01,     5.72501420974731445D-01,&
         -1.91945766231840700D+03,     8.06172218173730938D+03,&
         -1.35865500064341374D+04,     1.16553933368645332D+04,&
         -5.30564697861340311D+03,     1.20090291321635246D+03,&
         -1.08090919788394656D+02,     1.72772750258445740D+00,&
         2.02042913309661486D+04,    -9.69805983886375135D+04,&
         1.92547001232531532D+05,    -2.03400177280415534D+05,&
         1.22200464983017460D+05,    -4.11926549688975513D+04,&
         7.10951430248936372D+03,    -4.93915304773088012D+02,&
         6.07404200127348304D+00,    -2.42919187900551333D+05,&
         1.31176361466297720D+06,    -2.99801591853810675D+06,&
         3.76327129765640400D+06,    -2.81356322658653411D+06,&
         1.26836527332162478D+06,    -3.31645172484563578D+05,&
         4.52187689813627263D+04,    -2.49983048181120962D+03,&
         2.43805296995560639D+01,     3.28446985307203782D+06,&
         -1.97068191184322269D+07,     5.09526024926646422D+07,&
         -7.41051482115326577D+07,     6.63445122747290267D+07,&
         -3.75671766607633513D+07,     1.32887671664218183D+07,&
         -2.78561812808645469D+06,     3.08186404612662398D+05,&
         -1.38860897537170405D+04,     1.10017140269246738D+02,&
         -4.93292536645099620D+07,     3.25573074185765749D+08,&
         -9.39462359681578403D+08,     1.55359689957058006D+09,&
         -1.62108055210833708D+09,     1.10684281682301447D+09,&
         -4.95889784275030309D+08,     1.42062907797533095D+08,&
         -2.44740627257387285D+07,     2.24376817792244943D+06,&
         -8.40054336030240853D+04,     5.51335896122020586D+02,&
         8.14789096118312115D+08,    -5.86648149205184723D+09,&
         1.86882075092958249D+10,    -3.46320433881587779D+10,&
         4.12801855797539740D+10,    -3.30265997498007231D+10,&
         1.79542137311556001D+10,    -6.56329379261928433D+09,&
         1.55927986487925751D+09,    -2.25105661889415278D+08,&
         1.73951075539781645D+07,    -5.49842327572288687D+05,&
         3.03809051092238427D+03,    -1.46792612476956167D+10,&
         1.14498237732025810D+11,    -3.99096175224466498D+11,&
         8.19218669548577329D+11,    -1.09837515608122331D+12,&
         1.00815810686538209D+12,    -6.45364869245376503D+11,&
         2.87900649906150589D+11,    -8.78670721780232657D+10,&
         1.76347306068349694D+10,    -2.16716498322379509D+09,&
         1.43157876718888981D+08,    -3.87183344257261262D+06,&
         1.82577554742931747D+04 ]
      real(8),parameter::alfa(180)=[&
         -4.44444444444444444D-03,    -9.22077922077922078D-04,&
         -8.84892884892884893D-05,     1.65927687832449737D-04,&
         2.46691372741792910D-04,     2.65995589346254780D-04,&
         2.61824297061500945D-04,     2.48730437344655609D-04,&
         2.32721040083232098D-04,     2.16362485712365082D-04,&
         2.00738858762752355D-04,     1.86267636637545172D-04,&
         1.73060775917876493D-04,     1.61091705929015752D-04,&
         1.50274774160908134D-04,     1.40503497391269794D-04,&
         1.31668816545922806D-04,     1.23667445598253261D-04,&
         1.16405271474737902D-04,     1.09798298372713369D-04,&
         1.03772410422992823D-04,     9.82626078369363448D-05,&
         9.32120517249503256D-05,     8.85710852478711718D-05,&
         8.42963105715700223D-05,     8.03497548407791151D-05,&
         7.66981345359207388D-05,     7.33122157481777809D-05,&
         7.01662625163141333D-05,     6.72375633790160292D-05,&
         6.93735541354588974D-04,     2.32241745182921654D-04,&
         -1.41986273556691197D-05,    -1.16444931672048640D-04,&
         -1.50803558053048762D-04,    -1.55121924918096223D-04,&
         -1.46809756646465549D-04,    -1.33815503867491367D-04,&
         -1.19744975684254051D-04,    -1.06184319207974020D-04,&
         -9.37699549891194492D-05,    -8.26923045588193274D-05,&
         -7.29374348155221211D-05,    -6.44042357721016283D-05,&
         -5.69611566009369048D-05,    -5.04731044303561628D-05,&
         -4.48134868008882786D-05,    -3.98688727717598864D-05,&
         -3.55400532972042498D-05,    -3.17414256609022480D-05,&
         -2.83996793904174811D-05,    -2.54522720634870566D-05,&
         -2.28459297164724555D-05,    -2.05352753106480604D-05,&
         -1.84816217627666085D-05,    -1.66519330021393806D-05,&
         -1.50179412980119482D-05,    -1.35554031379040526D-05,&
         -1.22434746473858131D-05,    -1.10641884811308169D-05,&
         -3.54211971457743841D-04,    -1.56161263945159416D-04,&
         3.04465503594936410D-05,     1.30198655773242693D-04,&
         1.67471106699712269D-04,     1.70222587683592569D-04,&
         1.56501427608594704D-04,     1.36339170977445120D-04,&
         1.14886692029825128D-04,     9.45869093034688111D-05,&
         7.64498419250898258D-05,     6.07570334965197354D-05,&
         4.74394299290508799D-05,     3.62757512005344297D-05,&
         2.69939714979224901D-05,     1.93210938247939253D-05,&
         1.30056674793963203D-05,     7.82620866744496661D-06,&
         3.59257485819351583D-06,     1.44040049814251817D-07,&
         -2.65396769697939116D-06,    -4.91346867098485910D-06,&
         -6.72739296091248287D-06,    -8.17269379678657923D-06,&
         -9.31304715093561232D-06,    -1.02011418798016441D-05,&
         -1.08805962510592880D-05,    -1.13875481509603555D-05,&
         -1.17519675674556414D-05,    -1.19987364870944141D-05,&
         3.78194199201772914D-04,     2.02471952761816167D-04,&
         -6.37938506318862408D-05,    -2.38598230603005903D-04,&
         -3.10916256027361568D-04,    -3.13680115247576316D-04,&
         -2.78950273791323387D-04,    -2.28564082619141374D-04,&
         -1.75245280340846749D-04,    -1.25544063060690348D-04,&
         -8.22982872820208365D-05,    -4.62860730588116458D-05,&
         -1.72334302366962267D-05,     5.60690482304602267D-06,&
         2.31395443148286800D-05,     3.62642745856793957D-05,&
         4.58006124490188752D-05,     5.24595294959114050D-05,&
         5.68396208545815266D-05,     5.94349820393104052D-05,&
         6.06478527578421742D-05,     6.08023907788436497D-05,&
         6.01577894539460388D-05,     5.89199657344698500D-05,&
         5.72515823777593053D-05,     5.52804375585852577D-05,&
         5.31063773802880170D-05,     5.08069302012325706D-05,&
         4.84418647620094842D-05,     4.60568581607475370D-05,&
         -6.91141397288294174D-04,    -4.29976633058871912D-04,&
         1.83067735980039018D-04,     6.60088147542014144D-04,&
         8.75964969951185931D-04,     8.77335235958235514D-04,&
         7.49369585378990637D-04,     5.63832329756980918D-04,&
         3.68059319971443156D-04,     1.88464535514455599D-04,&
         3.70663057664904149D-05,    -8.28520220232137023D-05,&
         -1.72751952869172998D-04,    -2.36314873605872983D-04,&
         -2.77966150694906658D-04,    -3.02079514155456919D-04,&
         -3.12594712643820127D-04,    -3.12872558758067163D-04,&
         -3.05678038466324377D-04,    -2.93226470614557331D-04,&
         -2.77255655582934777D-04,    -2.59103928467031709D-04,&
         -2.39784014396480342D-04,    -2.20048260045422848D-04,&
         -2.00443911094971498D-04,    -1.81358692210970687D-04,&
         -1.63057674478657464D-04,    -1.45712672175205844D-04,&
         -1.29425421983924587D-04,    -1.14245691942445952D-04,&
         1.92821964248775885D-03,     1.35592576302022234D-03,&
         -7.17858090421302995D-04,    -2.58084802575270346D-03,&
         -3.49271130826168475D-03,    -3.46986299340960628D-03,&
         -2.82285233351310182D-03,    -1.88103076404891354D-03,&
         -8.89531718383947600D-04,     3.87912102631035228D-06,&
         7.28688540119691412D-04,     1.26566373053457758D-03,&
         1.62518158372674427D-03,     1.83203153216373172D-03,&
         1.91588388990527909D-03,     1.90588846755546138D-03,&
         1.82798982421825727D-03,     1.70389506421121530D-03,&
         1.55097127171097686D-03,     1.38261421852276159D-03,&
         1.20881424230064774D-03,     1.03676532638344962D-03,&
         8.71437918068619115D-04,     7.16080155297701002D-04,&
         5.72637002558129372D-04,     4.42089819465802277D-04,&
         3.24724948503090564D-04,     2.20342042730246599D-04,&
         1.28412898401353882D-04,     4.82005924552095464D-05 &
         ]
      real(8),parameter::beta(210)=[&
         1.79988721413553309D-02,     5.59964911064388073D-03,&
         2.88501402231132779D-03,     1.80096606761053941D-03,&
         1.24753110589199202D-03,     9.22878876572938311D-04,&
         7.14430421727287357D-04,     5.71787281789704872D-04,&
         4.69431007606481533D-04,     3.93232835462916638D-04,&
         3.34818889318297664D-04,     2.88952148495751517D-04,&
         2.52211615549573284D-04,     2.22280580798883327D-04,&
         1.97541838033062524D-04,     1.76836855019718004D-04,&
         1.59316899661821081D-04,     1.44347930197333986D-04,&
         1.31448068119965379D-04,     1.20245444949302884D-04,&
         1.10449144504599392D-04,     1.01828770740567258D-04,&
         9.41998224204237509D-05,     8.74130545753834437D-05,&
         8.13466262162801467D-05,     7.59002269646219339D-05,&
         7.09906300634153481D-05,     6.65482874842468183D-05,&
         6.25146958969275078D-05,     5.88403394426251749D-05,&
         -1.49282953213429172D-03,    -8.78204709546389328D-04,&
         -5.02916549572034614D-04,    -2.94822138512746025D-04,&
         -1.75463996970782828D-04,    -1.04008550460816434D-04,&
         -5.96141953046457895D-05,    -3.12038929076098340D-05,&
         -1.26089735980230047D-05,    -2.42892608575730389D-07,&
         8.05996165414273571D-06,     1.36507009262147391D-05,&
         1.73964125472926261D-05,     1.98672978842133780D-05,&
         2.14463263790822639D-05,     2.23954659232456514D-05,&
         2.28967783814712629D-05,     2.30785389811177817D-05,&
         2.30321976080909144D-05,     2.28236073720348722D-05,&
         2.25005881105292418D-05,     2.20981015361991429D-05,&
         2.16418427448103905D-05,     2.11507649256220843D-05,&
         2.06388749782170737D-05,     2.01165241997081666D-05,&
         1.95913450141179244D-05,     1.90689367910436740D-05,&
         1.85533719641636667D-05,     1.80475722259674218D-05,&
         5.52213076721292790D-04,     4.47932581552384646D-04,&
         2.79520653992020589D-04,     1.52468156198446602D-04,&
         6.93271105657043598D-05,     1.76258683069991397D-05,&
         -1.35744996343269136D-05,    -3.17972413350427135D-05,&
         -4.18861861696693365D-05,    -4.69004889379141029D-05,&
         -4.87665447413787352D-05,    -4.87010031186735069D-05,&
         -4.74755620890086638D-05,    -4.55813058138628452D-05,&
         -4.33309644511266036D-05,    -4.09230193157750364D-05,&
         -3.84822638603221274D-05,    -3.60857167535410501D-05,&
         -3.37793306123367417D-05,    -3.15888560772109621D-05,&
         -2.95269561750807315D-05,    -2.75978914828335759D-05,&
         -2.58006174666883713D-05,    -2.41308356761280200D-05,&
         -2.25823509518346033D-05,    -2.11479656768912971D-05,&
         -1.98200638885294927D-05,    -1.85909870801065077D-05,&
         -1.74532699844210224D-05,    -1.63997823854497997D-05,&
         -4.74617796559959808D-04,    -4.77864567147321487D-04,&
         -3.20390228067037603D-04,    -1.61105016119962282D-04,&
         -4.25778101285435204D-05,     3.44571294294967503D-05,&
         7.97092684075674924D-05,     1.03138236708272200D-04,&
         1.12466775262204158D-04,     1.13103642108481389D-04,&
         1.08651634848774268D-04,     1.01437951597661973D-04,&
         9.29298396593363896D-05,     8.40293133016089978D-05,&
         7.52727991349134062D-05,     6.69632521975730872D-05,&
         5.92564547323194704D-05,     5.22169308826975567D-05,&
         4.58539485165360646D-05,     4.01445513891486808D-05,&
         3.50481730031328081D-05,     3.05157995034346659D-05,&
         2.64956119950516039D-05,     2.29363633690998152D-05,&
         1.97893056664021636D-05,     1.70091984636412623D-05,&
         1.45547428261524004D-05,     1.23886640995878413D-05,&
         1.04775876076583236D-05,     8.79179954978479373D-06,&
         7.36465810572578444D-04,     8.72790805146193976D-04,&
         6.22614862573135066D-04,     2.85998154194304147D-04,&
         3.84737672879366102D-06,    -1.87906003636971558D-04,&
         -2.97603646594554535D-04,    -3.45998126832656348D-04,&
         -3.53382470916037712D-04,    -3.35715635775048757D-04,&
         -3.04321124789039809D-04,    -2.66722723047612821D-04,&
         -2.27654214122819527D-04,    -1.89922611854562356D-04,&
         -1.55058918599093870D-04,    -1.23778240761873630D-04,&
         -9.62926147717644187D-05,    -7.25178327714425337D-05,&
         -5.22070028895633801D-05,    -3.50347750511900522D-05,&
         -2.06489761035551757D-05,    -8.70106096849767054D-06,&
         1.13698686675100290D-06,     9.16426474122778849D-06,&
         1.56477785428872620D-05,     2.08223629482466847D-05,&
         2.48923381004595156D-05,     2.80340509574146325D-05,&
         3.03987774629861915D-05,     3.21156731406700616D-05,&
         -1.80182191963885708D-03,    -2.43402962938042533D-03,&
         -1.83422663549856802D-03,    -7.62204596354009765D-04,&
         2.39079475256927218D-04,     9.49266117176881141D-04,&
         1.34467449701540359D-03,     1.48457495259449178D-03,&
         1.44732339830617591D-03,     1.30268261285657186D-03,&
         1.10351597375642682D-03,     8.86047440419791759D-04,&
         6.73073208165665473D-04,     4.77603872856582378D-04,&
         3.05991926358789362D-04,     1.60315694594721630D-04,&
         4.00749555270613286D-05,    -5.66607461635251611D-05,&
         -1.32506186772982638D-04,    -1.90296187989614057D-04,&
         -2.32811450376937408D-04,    -2.62628811464668841D-04,&
         -2.82050469867598672D-04,    -2.93081563192861167D-04,&
         -2.97435962176316616D-04,    -2.96557334239348078D-04,&
         -2.91647363312090861D-04,    -2.83696203837734166D-04,&
         -2.73512317095673346D-04,    -2.61750155806768580D-04,&
         6.38585891212050914D-03,     9.62374215806377941D-03,&
         7.61878061207001043D-03,     2.83219055545628054D-03,&
         -2.09841352012720090D-03,    -5.73826764216626498D-03,&
         -7.70804244495414620D-03,    -8.21011692264844401D-03,&
         -7.65824520346905413D-03,    -6.47209729391045177D-03,&
         -4.99132412004966473D-03,    -3.45612289713133280D-03,&
         -2.01785580014170775D-03,    -7.59430686781961401D-04,&
         2.84173631523859138D-04,     1.10891667586337403D-03,&
         1.72901493872728771D-03,     2.16812590802684701D-03,&
         2.45357710494539735D-03,     2.61281821058334862D-03,&
         2.67141039656276912D-03,     2.65203073395980430D-03,&
         2.57411652877287315D-03,     2.45389126236094427D-03,&
         2.30460058071795494D-03,     2.13684837686712662D-03,&
         1.95896528478870911D-03,     1.77737008679454412D-03,&
         1.59690280765839059D-03,     1.42111975664438546D-03&
         ]
      real(8),parameter::gama(30)=[&
         6.29960524947436582D-01,     2.51984209978974633D-01,&
         1.54790300415655846D-01,     1.10713062416159013D-01,&
         8.57309395527394825D-02,     6.97161316958684292D-02,&
         5.86085671893713576D-02,     5.04698873536310685D-02,&
         4.42600580689154809D-02,     3.93720661543509966D-02,&
         3.54283195924455368D-02,     3.21818857502098231D-02,&
         2.94646240791157679D-02,     2.71581677112934479D-02,&
         2.51768272973861779D-02,     2.34570755306078891D-02,&
         2.19508390134907203D-02,     2.06210828235646240D-02,&
         1.94388240897880846D-02,     1.83810633800683158D-02,&
         1.74293213231963172D-02,     1.65685837786612353D-02,&
         1.57865285987918445D-02,     1.50729501494095594D-02,&
         1.44193250839954639D-02,     1.38184805735341786D-02,&
         1.32643378994276568D-02,     1.27517121970498651D-02,&
         1.22761545318762767D-02,     1.18338262398482403D-02 &
         ]

      real(8),parameter::ex1 =3.33333333333333333D-01
      real(8),parameter::ex2 =6.66666666666666667D-01
      !DEBUG write (*, *) 'call f90_zunhj'
      rfnu=1.0d0/fnu
      !-----------------------------------------------------------------------
      !     overflow test (z/fnu too small)
      !-----------------------------------------------------------------------
      test=d1mach1*1.0d+3
      ac=fnu*test
      if(abs(z%re)>ac.or.abs(z%im)>ac) goto 15
         zeta1=2.d0*abs(log(test))+fnu
         zeta2=fnu
         phi=1.0d0
         arg=1.0d0
      return
      15 continue
      zb=z*rfnu
      rfnu2 = rfnu*rfnu
      !-----------------------------------------------------------------------
      !     compute in the fourth quadrant
      !-----------------------------------------------------------------------
      fn13=fnu**ex1
      fn23=fn13**2
      rfn13=1.d0/fn13
      w2=1.0d0-zb*zb
      aw2=abs(w2)
      if (aw2>0.25d0) goto 130
      !-----------------------------------------------------------------------
      !     power series for cabs(w2)<=0.25d0
      !-----------------------------------------------------------------------
      k = 1
      p(1) = 1.0d0
      suma = gama(1)
      ap(1) = 1.0d0
      if (aw2<tol) goto 20
      do k = 2, 30
      p(k)=p(k-1)*w2
      suma=suma+gama(k)*p(k)
      ap(k)=ap(k-1)*aw2
      if (ap(k)<tol) goto 20
      end do
      k = 30
      20 continue
      kmax = k
      zeta=w2*suma
      arg=zeta*fn23
      za=sqrt(suma)
      st=sqrt(w2)
      zeta2 = st*fnu
      st = 1.0d0 + ex2*(zeta*za)
      zeta1 = st*zeta2
      za = 2*za
      st = sqrt(za)
      phi =st*rfn13
      if (ipmtr==1) goto 120
      !-----------------------------------------------------------------------
      !     sum series for asum and bsum
      !-----------------------------------------------------------------------
      sumb=0.0d0
      do k=1,kmax
      sumb=sumb+beta(k)*p(k)
      end do
      asum=0.0d0
      bsum=sumb
      l1=0
      l2=30
      btol=tol*(abs(bsum%re)+abs(bsum%im))
      atol=tol
      pp=1.0d0
      ias=0
      ibs=0
      if (rfnu2<tol) goto 110
      do is = 2, 7
      atol = atol/rfnu2
      pp = pp*rfnu2
      if (ias==1) goto 60
      suma = 0.0d0
      do k = 1, kmax
      m = l1 + k
      suma = suma + p(k)*alfa(m)
      if (ap(k)<atol) goto 50
      end do
      50 continue
      asum = asum + suma*pp
      if (pp<tol) ias = 1
      60 continue
      if (ibs==1) goto 90
      sumb = 0.0d0
      do k = 1, kmax
      m = l2 + k
      sumb = sumb + p(k)*beta(m)
      if (ap(k)<atol) goto 80
      end do
      80 continue
      bsum = bsum + sumb*pp
      if (pp<btol) ibs = 1
      90 continue
      if (ias==1 .and. ibs==1) goto 110
      l1 = l1 + 30
      l2 = l2 + 30
      end do
      110 continue
      asum = asum + 1.0d0
      pp = rfnu*rfn13
      bsum = bsum*pp
      120 continue
      return
      !-----------------------------------------------------------------------
      !     cabs(w2)>0.25d0
      !-----------------------------------------------------------------------
      130 continue
      w=sqrt(w2)
      if (w%re<0.0d0) w%re = 0.0d0
      if (w%im<0.0d0) w%im = 0.0d0
      st=1.0d0+w
      za=st/zb
      zc=log(za)
      if (zc%im<0.0d0) zc%im = 0.0d0
      if (zc%im>hpi)   zc%im = hpi
      if (zc%re<0.0d0) zc%re = 0.0d0
      zth=1.5d0*(zc-w)
      zeta1 = zc*fnu
      zeta2 = w*fnu
      azth=abs(zth)
      ang = thpi
      if (zth%re>=0.0d0 .and. zth%im<0.0d0) goto 140
      ang = hpi
      if (iszero(zth%re)) goto 140
      ang = atan(zth%im/zth%re)
      if (zth%re<0.0d0) ang = ang + dpi
      140 continue
      pp = azth**ex2
      ang = ang*ex2
      zeta=cmplx(pp*cos(ang),pp*sin(ang),8)
      if (zeta%im<0.0d0) zeta%im = 0.0d0
      arg = zeta*fn23
      rtzt=zth/zeta
      za=rtzt/w
      tza = za + za
      st=sqrt(tza)
      phi = st*rfn13
      if (ipmtr==1) goto 120
      raw = 1.0d0/sqrt(aw2)
      st=conjg(w)*raw
      tfn = st*rfnu*raw
      razth = 1.0d0/azth
      st = conjg(zth)*razth
      rzth = st*azth*rfnu
      zc = rzth*ar(2)
      raw2 = 1.0d0/aw2
      st = conjg(w2)*raw2
      t2 = st*raw2
      st = t2*c(2) + c(3)
      up(2)=st*tfn
      bsum = up(2) + zc
      asum = 0.0d0
      if (rfnu<tol) goto 220
      przth = rzth
      ptfn = tfn
      up(1) = 1.0d0
      pp = 1.0d0
      btol = tol*(abs(bsum%re)+abs(bsum%im))
      ks = 0
      kp1 = 2
      l = 3
      ias = 0
      ibs = 0
      do lr = 2, 12, 2
      lrp1 = lr + 1
      !-----------------------------------------------------------------------
      !     compute two additional cr, dr, and up for two more terms in
      !     next suma and sumb
      !-----------------------------------------------------------------------
      do k = lr, lrp1
      ks = ks + 1
      kp1 = kp1 + 1
      l = l + 1
      za = c(l)
      do j = 2, kp1
      l = l + 1
      za=za*t2+cl
      end do
      ptfn=ptfn*tfn
      up(kp1)=ptfn*za
      cr(ks) = przth*br(ks+1)
      przth = przth*rzth
      dr(ks) = przth*ar(ks+2)
      end do
      pp = pp*rfnu2
      if (ias==1) goto 180
      suma = up(lrp1)
      ju = lrp1
      do  jr=1,lr
      ju = ju - 1
      suma = suma + cr(jr)*up(ju)
      end do
      asum = asum + suma
      test = abs(suma%re) + abs(suma%im)
      if (pp<tol .and. test<tol) ias = 1
      180 continue
      if (ibs==1) goto 200
      sumb = up(lr+2) + up(lrp1)*zc
      ju = lrp1
      do jr=1,lr
      ju = ju - 1
      sumb = sumb + dr(jr)*up(ju)
      end do
      bsum = bsum + sumb
      test = abs(sumb%re) + abs(sumb%im)
      if (pp<btol .and. test<btol) ibs = 1
      200 continue
      if (ias==1 .and. ibs==1) goto 220
      end do
      220 continue
      asum = asum + 1.0d0
      st = -bsum*rfn13
      bsum=st/rtzt
      goto 120
   end subroutine f90_zunhj

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zuoik(z, fnu, kode, ikflg, n, y, nuf, tol, elim, alim)
      !***begin prologue  zuoik
      !***refer to  zbesi,zbesk,zbesh
      !
      !     zuoik computes the leading terms of the uniform asymptotic
      !     expansions for the i and k functions and compares them
      !     (in logarithmic form) to alim and elim for over and underflow
      !     where alim<elim. if the magnitude, based on the leading
      !     exponential, is less than alim or greater than -alim, then
      !     the result is on scale. if not, then a refined test using other
      !     multipliers (in logarithmic form) is made based on elim. here
      !     exp(-elim)=smallest machine number*1.0e+3 and exp(-alim)=
      !     exp(-elim)/tol
      !
      !     ikflg=1 means the i sequence is tested
      !          =2 means the k sequence is tested
      !     nuf = 0 means the last member of the sequence is on scale
      !         =-1 means an overflow would occur
      !     ikflg=1 and nuf>0 means the last nuf y values were set to zero
      !             the first n-nuf values must be set by another routine
      !     ikflg=2 and nuf==n means all y values were set to zero
      !     ikflg=2 and 0<nuf<n not considered. y must be set by
      !             another routine
      !
      !***routines called  zuchk,zunhj,zunik,d1mach,azabs,azlog
      !***end prologue  zuoik
      !     complex arg,asum,bsum,cwrk,cz,czero,phi,sum,y,z,zb,zeta1,zeta2,zn,
      !    *zr
      integer,intent(in)   ::kode,n
      integer,intent(inout)::nuf
      real(8),intent(in)   ::tol,fnu,elim,alim
      complex(8),intent(in)::z
      complex(8),intent(inout)::y(n)
      real(8)   ::aarg, aphi,ascle, ax, ay, fnn, gnn, gnu, rcz
      integer   ::i, iform, ikflg, init, nn, nw
      complex(8)::arg,asum,bsum,cwrk(16),cz,sum,zb,zeta1,zeta2,zr,phi,zn,st
      real(8),parameter::aic=1.265512123484645396d+00
      !DEBUG write (*, *) 'call f90_zuoik'
      nuf = 0
      nn = n
      zr = z
      if (z%re>=0.0d0) goto 10
      zr = -z
      10 continue
      zb = zr
      ax = abs(z%re)*1.7321d0
      ay = abs(z%im)
      iform = 1
      if (ay>ax) iform = 2
      gnu = max(fnu, 1.0d0)
      if (ikflg==1) goto 20
      fnn = real(nn,8)
      gnn = fnu + fnn - 1.0d0
      gnu = max(gnn, fnn)
      20 continue
      !-----------------------------------------------------------------------
      !     only the magnitude of arg and phi are needed along with the
      !     real parts of zeta1, zeta2 and zb. no attempt is made to get
      !     the sign of the imaginary part correct.
      !-----------------------------------------------------------------------
      if (iform==2) goto 30
      init = 0
      call f90_zunik(zr, gnu, ikflg, 1, tol, init, phi, zeta1, zeta2, sum, cwrk)
      cz=-zeta1+zeta2
      goto 50
      30 continue
      zn = conjg(zr)
      if (z%im>0.0d0) goto 40
      zn%re = -zn%re
      40 continue
      call f90_zunhj(zn, gnu, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
      cz = -zeta1 + zeta2
      aarg = abs(arg)
      50 continue
      if (kode==1) goto 60
      cz = cz - zb
      60 continue
      if (ikflg==1) goto 70
      cz = -cz
      70 continue
      aphi = abs(phi)
      rcz = cz%re
      !-----------------------------------------------------------------------
      !     overflow test
      !-----------------------------------------------------------------------
      if (rcz>elim) goto 210
      if (rcz<alim) goto 80
      rcz = rcz + log(aphi)
      if (iform==2) rcz = rcz - 0.25d0*log(aarg) - aic
      if (rcz>elim) goto 210
      goto 130
      80 continue
      !-----------------------------------------------------------------------
      !     underflow test
      !-----------------------------------------------------------------------
      if (rcz<(-elim)) goto 90
      if (rcz>(-alim)) goto 130
      rcz = rcz + log(aphi)
      if (iform==2) rcz = rcz - 0.25d0*log(aarg) - aic
      if (rcz>(-elim)) goto 110
      90 continue
      do i = 1, nn
      y(i) = 0.0d0
      end do
      nuf = nn
      return
      110 continue
      ascle = 1.0d+3*d1mach1/tol
      st=log(phi)
      cz=cz+st
      if (iform==1) goto 120
      st=log(arg)
      cz=cz-0.25d0*st-aic
      120 continue
      ax = exp(rcz)/tol
      ay = cz%im
      cz%re = ax*cos(ay)
      cz%im = ax*sin(ay)
      nw=f90_zuchk(cz,ascle,tol)
      if (nw/=0) goto 90
      130 continue
      if (ikflg==2) return
      if (n==1) return
      !-----------------------------------------------------------------------
      !     set underflows on i sequence
      !-----------------------------------------------------------------------
      140 continue
      gnu = fnu + real(nn-1,8)
      if (iform==2) goto 150
      init = 0
      call f90_zunik(zr, gnu, ikflg, 1, tol, init, phi, zeta1, zeta2, sum, cwrk)
      cz = -zeta1 + zeta2
      goto 160
      150 continue
      call f90_zunhj(zn, gnu, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
      cz = -zeta1 + zeta2
      aarg = abs(arg)
      160 continue
      if (kode==1) goto 170
      cz = cz - zb
      170 continue
      aphi = abs(phi)
      rcz = cz%re
      if (rcz<(-elim)) goto 180
      if (rcz>(-alim)) return
      rcz = rcz + log(aphi)
      if (iform==2) rcz = rcz - 0.25d0*log(aarg) - aic
      if (rcz>(-elim)) goto 190
      180 continue
      y(nn) = 0.0d0
      nn = nn - 1
      nuf = nuf + 1
      if (nn==0) return
      goto 140
      190 continue
      ascle = 1.0d+3*d1mach1/tol
      st=log(phi)
      cz = cz + st
      if (iform==1) goto 200
      st=log(arg)
      cz = cz - 0.25d0*st - aic
      200 continue
      ax = exp(rcz)/tol
      ay = cz%im
      cz%re = ax*cos(ay)
      cz%im = ax*sin(ay)
      nw=f90_zuchk(cz, ascle, tol)
      if (nw/=0) goto 180
      return
      210 continue
      nuf = -1
      return
   end subroutine f90_zuoik

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zwrsk(zr, fnu, kode, n, y, nz, cw, tol, elim, alim)
      !***begin prologue  zwrsk
      !***refer to  zbesi,zbesk
      !
      !     zwrsk computes the i bessel function for re(z)>=0.0 by
      !     normalizing the i function ratios from zrati by the wronskian
      !
      !***routines called  d1mach,zbknu,zrati,azabs
      !***end prologue  zwrsk
      integer,intent(in)   ::kode,n
      integer,intent(inout)::nz
      real(8),intent(in)   ::tol,fnu,elim,alim
      complex(8),intent(in)::zr
      complex(8),intent(inout)::y(n),cw(2)
      real(8)::act,acw,ascle,csclr,ract
      integer::i,nw
      complex(8)::cinu,ct,c1,c2,st,pt
      !-----------------------------------------------------------------------
      !     i(fnu+i-1,z) by backward recurrence for ratios
      !     y(i)=i(fnu+i,z)/i(fnu+i-1,z) from crati normalized by the
      !     wronskian with k(fnu,z) and k(fnu+1,z) from cbknu.
      !-----------------------------------------------------------------------
      !DEBUG write (*, *) 'call f90_zwrsk'
      nz = 0
      call f90_zbknu(zr, fnu, kode, 2, cw, nw, tol, elim, alim)
      if (nw/=0) goto 50
      call f90_zrati(zr, fnu, n, y, tol)
      !-----------------------------------------------------------------------
      !     recur forward on i(fnu+1,z) = r(fnu,z)*i(fnu,z),
      !     r(fnu+j-1,z)=y(j),  j=1,...,n
      !-----------------------------------------------------------------------
      cinu = 1.0d0
      if (kode==1) goto 10
      cinu%re = cos(zr%im)
      cinu%im = sin(zr%im)
      10 continue
      !-----------------------------------------------------------------------
      !     on low exponent machines the k functions can be close to both
      !     the under and overflow limits and the normalization must be
      !     scaled to prevent over or underflow. cuoik has determined that
      !     the result is on scale.
      !-----------------------------------------------------------------------
      acw = abs(cw(2))
      ascle = 1.0d+3*d1mach1/tol
      csclr = 1.0d0
      if (acw>ascle) goto 20
      csclr = 1.0d0/tol
      goto 30
      20 continue
      ascle = 1.0d0/ascle
      if (acw<ascle) goto 30
      csclr = tol
      30 continue
      c1 = cw(1)*csclr
      c2 = cw(2)*csclr
      st = y(1)
      !-----------------------------------------------------------------------
      !     cinu=cinu*(conjg(ct)/cabs(ct))*(1.0d0/cabs(ct) prevents
      !     under- or overflow prematurely by squaring cabs(ct)
      !-----------------------------------------------------------------------
      pt=st*c1
      pt=pt+c2
      ct=zr*pt
      act = abs(ct)
      ract = 1.0d0/act
      ct=conjg(ct)*ract
      pt = cinu*ract
      cinu = pt*ct
      y(1) = cinu*csclr
      if (n==1) return
      do i = 2, n
      cinu=st*cinu
      st = y(i)
      y(i) = cinu*csclr
      end do
      return
      50 continue
      nz = -1
      if (nw==(-2)) nz = -2
      return
   end subroutine f90_zwrsk


   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zmlri(z, fnu, kode, n, y, nz, tol)
      !***begin prologue  zmlri
      !***refer to  zbesi,zbesk
      !
      !     zmlri computes the i bessel function for re(z)>=0.0 by the
      !     miller algorithm normalized by a neumann series.
      !
      !***routines called  dgamln,d1mach,azabs,azexp,azlog,zmlt
      !***end prologue  zmlri
      !     complex ck,cnorm,cone,ctwo,czero,pt,p1,p2,rz,sum,y,z
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol
      integer,intent(in)::n,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::y(n)
      real(8) ack, ak, ap, at, az, bk, fkap, fkk, flam, fnf, &
         p1r, raz, rho, rho2, scle, tfnf, tst
      integer i, iaz, ifnu, inu, itime, k, kk, km, m
      complex(8)::ck,cnorm,pt,p1,p2,rz,sum,st
      !DEBUG write (*, *) 'call f90_zmlri'
      scle = d1mach1/tol
      nz = 0
      az = abs(z)
      iaz = int(az)
      ifnu = int(fnu)
      inu = ifnu + n - 1
      at = real(iaz,8) + 1.0d0
      raz = 1.0d0/az
      st=conjg(z)*raz
      ck = st*at*raz
      rz = (st+st)*raz
      p1 = 0.0d0
      p2 = 1.0d0
      ack = (at+1.0d0)*raz
      rho = ack + sqrt(ack*ack-1.0d0)
      rho2 = rho*rho
      tst = (rho2+rho2)/((rho2-1.0d0)*(rho-1.0d0))
      tst = tst/tol
      !-----------------------------------------------------------------------
      !     compute relative truncation error index for series
      !-----------------------------------------------------------------------
      ak = at
      do i = 1, 80
         pt=p2
         p2=p1-ck*pt
         p1 = pt
         ck = ck + rz
         ap = abs(p2)
         if (ap>tst*ak*ak) goto 20
         ak = ak + 1.0d0
      end do
      goto 110
      20 continue
      i = i + 1
      k = 0
      if (inu<iaz) goto 40
      !-----------------------------------------------------------------------
      !     compute relative truncation error for ratios
      !-----------------------------------------------------------------------
      p1 = 0.0d0
      p2 = 1.0d0
      at = real(inu,8) + 1.0d0
      st=conjg(z)*raz
      ck = st*at*raz
      ack = at*raz
      tst = sqrt(ack/tol)
      itime = 1
      do k = 1, 80
         pt = p2
         p2 = p1 - ck*pt
         p1 = pt
         ck = ck + rz
         ap = abs(p2)
         if (ap<tst) goto 30
         if (itime==2) goto 40
         ack = abs(ck)
         flam = ack + sqrt(ack*ack-1.0d0)
         fkap = ap/abs(p1)
         rho = min(flam, fkap)
         tst = tst*sqrt(rho/(rho*rho-1.0d0))
         itime = 2
      30 end do
      goto 110
      40 continue
      !-----------------------------------------------------------------------
      !     backward recurrence and sum normalizing relation
      !-----------------------------------------------------------------------
      k = k + 1
      kk = max(i+iaz, k+inu)
      fkk = real(kk,8)
      p1 = 0.0d0
      !-----------------------------------------------------------------------
      !     scale p2 and sum by scle
      !-----------------------------------------------------------------------
      p2 = scle
      fnf = fnu - real(ifnu,8)
      tfnf = fnf + fnf
      bk=log_gamma(fkk+tfnf+1.0d0)-log_gamma(fkk+1.0d0)-log_gamma(tfnf+1.0d0)
      bk = exp(bk)
      sum = 0.0d0
      km = kk - inu
      do i = 1, km
         pt = p2
         p2 = p1 + (fkk+fnf)*(rz*pt)
         p1 = pt
         ak = 1.0d0 - tfnf/(fkk+tfnf)
         ack = bk*ak
         sum = sum + (ack+bk)*p1
         bk = ack
         fkk = fkk - 1.0d0
      end do
      y(n) = p2
      if (n==1) goto 70
      do i = 2, n
         pt = p2
         p2 = p1 + (fkk+fnf)*(rz*pt)
         p1 = pt
         ak = 1.0d0 - tfnf/(fkk+tfnf)
         ack = bk*ak
         sum = sum + (ack+bk)*p1
         bk = ack
         fkk = fkk - 1.0d0
         m = n - i + 1
         y(m) = p2
      end do
      70 continue
      if (ifnu<=0) goto 90
      do i = 1, ifnu
         pt = p2
         p2 = p1 + (fkk+fnf)*(rz*pt)
         p1 = pt
         ak = 1.0d0 - tfnf/(fkk+tfnf)
         ack = bk*ak
         sum = sum + (ack+bk)*p1
         bk = ack
         fkk = fkk - 1.0d0
      end do
      90 continue
      pt = z
      if (kode==2) pt%re = 0.0d0
      p1 = -fnf*log(rz) + pt
      ap=log_gamma(1.0d0+fnf)
      pt = p1 - ap
      !-----------------------------------------------------------------------
      !     the division cexp(pt)/(sum+p2) is altered to avoid overflow
      !     in the denominator by squaring large quantities
      !-----------------------------------------------------------------------
      p2 = p2 + sum
      ap = abs(p2)
      p1r= 1.0d0/ap
      ck = exp(pt)*p1r
      pt = conjg(p2)*p1r
      cnorm=ck*pt
      do i = 1, n
         y(i)=y(i)*cnorm
      end do
      return
      110 continue
      nz = -2
      return
   end subroutine f90_zmlri

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zacai(z, fnu, kode, mr, n, y, nz, rl, tol, elim, alim)
      !***begin prologue  zacai
      !***refer to  zairy
      !
      !     zacai applies the analytic continuation formula
      !
      !         k(fnu,zn*exp(mp))=k(fnu,zn)*exp(-mp*fnu) - mp*i(fnu,zn)
      !                 mp=dpi*mr*cmplx(0.0,1.0)
      !
      !     to continue the k function from the right half to the left
      !     half z plane for use with zairy where fnu=1/3 or 2/3 and n=1.
      !     zacai is the same as zacon with the parts for larger orders and
      !     recurrence removed. a recursive call to zacon can result if zacon
      !     is called from zairy.
      !
      !***routines called  zasyi,zbknu,zmlri,zseri,zs1s2,d1mach,azabs
      !***end prologue  zacai
      !     complex csgn,cspn,c1,c2,y,z,zn,cy
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim,rl
      integer,intent(in)::n,kode,mr
      integer,intent(inout)::nz
      complex(8),intent(inout)::y(n)
      real(8)::arg, ascle, az, dfnu, fmr, sgn, yy
      integer inu, iuf, nn, nw
      complex(8)::cy(2),zn,csgn,cspn,c1,c2
      !DEBUG write(*,*)"call f90_zacai, not test"
      nz = 0
      zn = -z
      az = abs(z)
      nn = n
      dfnu = fnu + real(n-1,8)
      if (az<=2.0d0) goto 10
      if (az*az*0.25d0>dfnu+1.0d0) goto 20
      10 continue
      !-----------------------------------------------------------------------
      !     power series for the i function
      !-----------------------------------------------------------------------
      call f90_zseri(zn, fnu, kode, nn, y, nw, tol, elim, alim)
      goto 40
      20 continue
      if (az<rl) goto 30
      !-----------------------------------------------------------------------
      !     asymptotic expansion for large z for the i function
      !-----------------------------------------------------------------------
      call f90_zasyi(zn, fnu, kode, nn, y, nw, rl, tol, elim, alim)
      if (nw<0) goto 80
      goto 40
      30 continue
      !-----------------------------------------------------------------------
      !     miller algorithm normalized by the series for the i function
      !-----------------------------------------------------------------------
      call f90_zmlri(zn, fnu, kode, nn, y, nw, tol)
      if (nw<0) goto 80
      40 continue
      !-----------------------------------------------------------------------
      !     analytic continuation to the left half plane for the k function
      !-----------------------------------------------------------------------
      call f90_zbknu(zn, fnu, kode, 1, cy, nw, tol, elim, alim)
      if (nw/=0) goto 80
      fmr = real(mr,8)
      sgn = -sign(dpi, fmr)
      csgn=cmplx(0.0d0,sgn,8)
      if (kode==1) goto 50
      yy = -zn%im
      csgn%re = -csgn%im*sin(yy)
      csgn%im =  csgn%im*cos(yy)
      50 continue
      !-----------------------------------------------------------------------
      !     calculate cspn=exp(fnu*dpi*i) to minimize losses of significance
      !     when fnu is large
      !-----------------------------------------------------------------------
      inu = int(fnu)
      arg = (fnu-real(inu,8))*sgn
      cspn%re = cos(arg)
      cspn%im = sin(arg)
      if (mod(inu,2)==0) goto 60
      cspn = -cspn
      60 continue
      c1 = cy(1)
      c2 = y(1)
      if (kode==1) goto 70
      iuf = 0
      ascle = 1.0d+3*d1mach1/tol
      call f90_zs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
      70 continue
      y(1)=cspn*c1+csgn*c2
      return
      80 continue
      nz = -1
      if (nw==(-2)) nz = -2
      return
   end subroutine f90_zacai

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zairy(z, id, kode, ai, nz, ierr)
      !***begin prologue  zairy
      !***date written   830501   (yymmdd)
      !***revision date  890801   (yymmdd)
      !***category no.  b5k
      !***keywords  airy function,bessel functions of order one third
      !***author  amos, donald e., sandia national laboratories
      !***purpose  to compute airy functions ai(z) and dai(z) for complex z
      !***description
      !
      !                      ***a real(8) routine***
      !         on kode=1, zairy computes the complex airy function ai(z) or
      !         its derivative dai(z)/dz on id=0 or id=1 respectively. on
      !         kode=2, a scaling option cexp(zta)*ai(z) or cexp(zta)*
      !         dai(z)/dz is provided to remove the exponential decay in
      !         -pi/3<arg(z)<pi/3 and the exponential growth in
      !         pi/3<abs(arg(z))<pi where zta=(2/3)*z*csqrt(z).
      !
      !         while the airy functions ai(z) and dai(z)/dz are analytic in
      !         the whole z plane, the corresponding scaled functions defined
      !         for kode=2 have a cut along the negative real axis.
      !         defintions and notation are found in the nbs handbook of
      !         mathematical functions (ref. 1).
      !
      !         input      zr,zi are real(8)
      !           zr,zi  - z=cmplx(zr,zi)
      !           id     - order of derivative, id=0 or id=1
      !           kode   - a parameter to indicate the scaling option
      !                    kode= 1  returns
      !                             ai=ai(z)                on id=0 or
      !                             ai=dai(z)/dz            on id=1
      !                        = 2  returns
      !                             ai=cexp(zta)*ai(z)       on id=0 or
      !                             ai=cexp(zta)*dai(z)/dz   on id=1 where
      !                             zta=(2/3)*z*csqrt(z)
      !
      !         output     air,aii are real(8)
      !           air,aii- complex answer depending on the choices for id and
      !                    kode
      !           nz     - underflow indicator
      !                    nz= 0   , normal return
      !                    nz= 1   , ai=cmplx(0.0d0,0.0d0) due to underflow in
      !                              -pi/3<arg(z)<pi/3 on kode=1
      !           ierr   - error flag
      !                    ierr=0, normal return - computation completed
      !                    ierr=1, input error   - no computation
      !                    ierr=2, overflow      - no computation, real(zta)
      !                            too large on kode=1
      !                    ierr=3, cabs(z) large      - computation completed
      !                            losses of signifcance by argument reduction
      !                            produce less than half of machine accuracy
      !                    ierr=4, cabs(z) too large  - no computation
      !                            complete loss of accuracy by argument
      !                            reduction
      !                    ierr=5, error              - no computation,
      !                            algorithm termination condition not met
      !
      !***long description
      !
      !         ai and dai are computed for cabs(z)>1.0 from the k bessel
      !         functions by
      !
      !            ai(z)=c*sqrt(z)*k(1/3,zta) , dai(z)=-c*z*k(2/3,zta)
      !                           c=1.0/(pi*sqrt(3.0))
      !                            zta=(2/3)*z**(3/2)
      !
      !         with the power series for cabs(z)<=1.0.
      !
      !         in most complex variable computation, one must evaluate ele-
      !         mentary functions. when the magnitude of z is large, losses
      !         of significance by argument reduction occur. consequently, if
      !         the magnitude of zeta=(2/3)*z**1.5 exceeds u1=sqrt(0.5/ur),
      !         then losses exceeding half precision are likely and an error
      !         flag ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
      !         real(8) unit roundoff limited to 18 digits precision.
      !         also, if the magnitude of zeta is larger than u2=0.5/ur, then
      !         all significance is lost and ierr=4. in order to use the int
      !         function, zeta must be further restricted not to exceed the
      !         largest integer, u3=i1mach(9). thus, the magnitude of zeta
      !         must be restricted by min(u2,u3). on 32 bit machines, u1,u2,
      !         and u3 are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single
      !         precision arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double
      !         precision arithmetic respectively. this makes u2 and u3 limit-
      !         ing in their respective arithmetics. this means that the mag-
      !         nitude of z cannot exceed 3.1e+4 in single and 2.1e+6 in
      !         real(8) arithmetic. this also means that one can
      !         expect to retain, in the worst cases on 32 bit machines,
      !         no digits in single precision and only 7 digits in double
      !         precision arithmetic. similar considerations hold for other
      !         machines.
      !
      !         the approximate relative error in the magnitude of a complex
      !         bessel function can be expressed by p*10**s where p=max(unit
      !         roundoff,1.0e-18) is the nominal precision and 10**s repre-
      !         sents the increase in error due to argument reduction in the
      !         elementary functions. here, s=max(1,abs(log10(cabs(z))),
      !         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
      !         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
      !         have only absolute accuracy. this is most likely to occur when
      !         one component (in absolute value) is larger than the other by
      !         several orders of magnitude. if one component is 10**k larger
      !         than the other, then one can expect only max(abs(log10(p))-k,
      !         0) significant digits; or, stated another way, when k exceeds
      !         the exponent of p, no significant digits remain in the smaller
      !         component. however, the phase angle retains absolute accuracy
      !         because, in complex arithmetic with precision p, the smaller
      !         component will not (as a rule) decrease below p times the
      !         magnitude of the larger component. in these extreme cases,
      !         the principal phase angle is on the order of +p, -p, pi/2-p,
      !         or -pi/2+p.
      !
      !***references  handbook of mathematical functions by m. abramowitz
      !                 and i. a. stegun, nbs ams series 55, u.s. dept. of
      !                 commerce, 1955.
      !
      !               computation of bessel functions of complex argument
      !                 and large order by d. e. amos, sand83-0643, may, 1983
      !
      !               a subroutine package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, sand85-
      !                 1018, may, 1985
      !
      !               a portable package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, trans.
      !                 math. software, 1986
      !
      !***routines called  zacai,zbknu,azexp,azsqrt,i1mach,d1mach
      !***end prologue  zairy
      !     complex ai,cone,csq,cy,s1,s2,trm1,trm2,z,zta,z3
      complex(8),intent(in)::z
      complex(8),intent(inout)::ai
      integer,intent(in)::id,kode
      integer,intent(inout)::nz,ierr
      real(8)::aa, ad, ak, alim, atrm, az, az3, bk,&
         cc, ck,dig, dk, d1, d2, elim, fid, fnu, rl, sfac,&
         tol, alaz
      integer::iflag, k, mr, nn
      complex(8)::cy(1),csq,s1,s2,trm1,trm2,zta,z3
      real(8),parameter::c1 =3.55028053887817240d-01
      real(8),parameter::c2 =2.58819403792806799d-01
      real(8),parameter::coef=1.83776298473930683d-01
      !***first executable statement  zairy
      !DEBUG write(*,*)"call f90_zairy,not test"
      ierr = 0
      nz = 0
      if (id<0 .or. id>1) ierr = 1
      if (kode<1 .or. kode>2) ierr = 1
      if (ierr/=0) return
      az = abs(z)
      tol = max(d1mach4, 1.0d-18)
      fid = real(id,8)
      if (az>1.0d0) goto 70
      !-----------------------------------------------------------------------
      !     power series for cabs(z)<=1.
      !-----------------------------------------------------------------------
      s1 = 1.d0
      s2 = 1.d0
      if (az<tol) goto 170
      aa = az*az
      if (aa<tol/az) goto 40
      trm1 = 1.d0
      trm2 = 1.d0
      atrm = 1.0d0
      z3=z**3
      az3 = az*aa
      ak = 2.0d0 + fid
      bk = 3.0d0 - fid - fid
      ck = 4.0d0 - fid
      dk = 3.0d0 + fid + fid
      d1 = ak*dk
      d2 = bk*ck
      ad = min(d1, d2)
      ak = 24.0d0 + 9.0d0*fid
      bk = 30.0d0 - 9.0d0*fid
      do k = 1, 25
         trm1=trm1*z3/d1
         s1=s1+trm1
         trm2=trm2*z3/d2
         s2 = s2 + trm2
         atrm = atrm*az3/ad
         d1 = d1 + ak
         d2 = d2 + bk
         ad = min(d1, d2)
         if (atrm<tol*ad) goto 40
         ak = ak + 18.0d0
         bk = bk + 18.0d0
      end do
      40 continue
      if (id==1) goto 50
      ai = s1*c1 - c2*(z*s2)
      if (kode==1) return
      ai=ai*exp(tth*z*sqrt(z))
      return
      50 continue
      ai = -s2*c2
      if (az<=tol) goto 60
      cc = c1/(1.0d0+fid)
      ai = ai + cc*(z*s1*z)
      60 continue
      if (kode==1) return
      ai=ai*exp(tth*(z*sqrt(z)))
      return
      !-----------------------------------------------------------------------
      !     case for cabs(z)>1.0
      !-----------------------------------------------------------------------
      70 continue
      fnu = (1.0d0+fid)/3.0d0
      !-----------------------------------------------------------------------
      !     set parameters related to machine constants.
      !     tol is the approximate unit roundoff limited to 1.0d-18.
      !     elim is the approximate exponential over- and underflow limit.
      !     exp(-elim)<exp(-alim)=exp(-elim)/tol    and
      !     exp(elim)>exp(alim)=exp(elim)*tol       are intervals near
      !     underflow and overflow limits where scaled arithmetic is done.
      !     rl is the lower boundary of the asymptotic expansion for large z.
      !     dig = number of base 10 digits in tol = 10**(-dig).
      !-----------------------------------------------------------------------
      k    = min(abs(i1mach15), abs(i1mach16))
      elim = 2.303d0*(real(k,8)*d1mach5-3.0d0)
      aa   = d1mach5*real(i1mach14 - 1,8)
      dig  = min(aa, 18.0d0)
      aa   = aa*2.303d0
      alim = elim + max(-aa, -41.45d0)
      rl   = 1.2d0*dig + 3.0d0
      alaz = log(az)
      !--------------------------------------------------------------------------
      !     test for proper range
      !-----------------------------------------------------------------------
      aa = 0.5d0/tol
      aa = min(aa,real(i1mach9,8)*0.5d0)
      aa = aa**tth
      if (az>aa) goto 260
      aa = sqrt(aa)
      if (az>aa) ierr = 3
      csq=sqrt(z)
      zta = tth*(z*csq)
      !-----------------------------------------------------------------------
      !     re(zta)<=0 when re(z)<0, especially when im(z) is small
      !-----------------------------------------------------------------------
      iflag = 0
      sfac = 1.0d0
      ak = zta%im
      if (z%re>=0.0d0) goto 80
      bk = zta%re
      ck = -dabs(bk)
      zta=cmplx(ck,ak,8)
      80 continue
      if (.not.iszero(z%im)) goto 90
      if (z%re>0.0d0) goto 90
      zta=cmplx(0.0d0,ak,8)
      90 continue
      aa = zta%re
      if (aa>=0.0d0 .and. z%re>0.0d0) goto 110
      if (kode==2) goto 100
      !-----------------------------------------------------------------------
      !     overflow test
      !-----------------------------------------------------------------------
      if (aa>(-alim)) goto 100
      aa = -aa + 0.25d0*alaz
      iflag = 1
      sfac = tol
      if (aa>elim) goto 270
      100 continue
      !-----------------------------------------------------------------------
      !     cbknu and cacon return exp(zta)*k(fnu,zta) on kode=2
      !-----------------------------------------------------------------------
      mr = 1
      if (z%im<0.0d0) mr = -1
      call f90_zacai(zta, fnu, kode, mr, 1, cy, nn, rl, tol, elim, alim)
      if (nn<0) goto 280
      nz = nz + nn
      goto 130
      110 continue
      if (kode==2) goto 120
      !-----------------------------------------------------------------------
      !     underflow test
      !-----------------------------------------------------------------------
      if (aa<alim) goto 120
      aa = -aa - 0.25d0*alaz
      iflag = 2
      sfac = 1.0d0/tol
      if (aa<(-elim)) goto 210
      120 continue
      call f90_zbknu(zta, fnu, kode, 1, cy, nz, tol, elim, alim)
      130 continue
      s1 = cy(1)*coef
      if (iflag/=0) goto 150
      if (id==1) goto 140
      ai = csq*s1
      return
      140 continue
      ai=-z*s1
      return
      150 continue
      s1 = s1*sfac
      if (id==1) goto 160
      s1=s1*csq
      ai = s1/sfac
      return
      160 continue
      s1=-s1*z
      ai = s1/sfac
      return
      170 continue
      aa = 1.0d+3*d1mach1
      s1 = 0.d0
      if (id==1) goto 190
      if (az<=aa) goto 180
      s1 = c2*z
      180 continue
      ai=c1-s1
      return
      190 continue
      ai = -c2
      aa = sqrt(aa)
      if (az<=aa) goto 200
      s1=0.5d0*z*z
      200 continue
      ai = ai + c1*s1
      return
      210 continue
      nz = 1
      ai = 0.d0
      return
      270 continue
      nz = 0
      ierr = 2
      return
      280 continue
      if (nn==(-1)) goto 270
      nz = 0
      ierr = 5
      return
      260 continue
      ierr = 4
      nz = 0
      return
   end subroutine f90_zairy

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zuni1(z, fnu, kode, n, y, nz, nlast, fnul, tol, elim, alim)
      !***begin prologue  zuni1
      !***refer to  zbesi,zbesk
      !
      !     zuni1 computes i(fnu,z)  by means of the uniform asymptotic
      !     expansion for i(fnu,z) in -pi/3<=arg z<=pi/3.
      !
      !     fnul is the smallest order permitted for the asymptotic
      !     expansion. nlast=0 means all of the y values were set.
      !     nlast/=0 is the number left to be computed by another
      !     formula for orders fnu to fnu+nlast-1 because fnu+nlast-1<fnul.
      !     y(i)=czero for i=nlast+1,n
      !
      !***routines called  zuchk,zunik,zuoik,d1mach,azabs
      !***end prologue  zuni1
      !     complex cfn,cone,crsc,cscl,csr,css,cwrk,czero,c1,c2,phi,rz,sum,s1,
      !    *s2,y,z,zeta1,zeta2
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim,fnul
      integer,intent(in)::n,kode
      integer,intent(inout)::nz,nlast
      complex(8),intent(inout)::y(n)
      real(8)::aphi, ascle, crsc, cscl,c1r, c2m,fn,rast, rs1
      integer::i, iflag, init, k, m, nd, nn, nuf, nw
      complex(8)::cwrk(16),cy(2),s1,s2,c2,phi,rz,st,sum,zeta1,zeta2
      real(8)::bry(3),cssr(3),csrr(3)
      !
      !DEBUG write (*, *) 'call f90_zuni1, not test'
      nz = 0
      nd = n
      nlast = 0
      !-----------------------------------------------------------------------
      !     computed values with exponents between alim and elim in mag-
      !     nitude are scaled to keep intermediate arithmetic on scale,
      !     exp(alim)=exp(elim)*tol
      !-----------------------------------------------------------------------
      cscl = 1.0d0/tol
      crsc = tol
      cssr(1) = cscl
      cssr(2) = 1.d0
      cssr(3) = crsc
      csrr(1) = crsc
      csrr(2) = 1.d0
      csrr(3) = cscl
      bry(1) = 1.0d+3*d1mach1/tol
      !-----------------------------------------------------------------------
      !     check for underflow and overflow on first member
      !-----------------------------------------------------------------------
      fn = max(fnu, 1.0d0)
      init = 0
      call f90_zunik(z, fn, 1, 1, tol, init, phi, zeta1, zeta2, sum, cwrk)
      if (kode==1) goto 10
      st = z + zeta2
      rast = fn/abs(st)
      st=conjg(st)*rast*rast
      s1 = -zeta1 + st
      goto 20
      10 continue
      s1 = -zeta1 + zeta2
      20 continue
      rs1 = s1%re
      if (abs(rs1)>elim) goto 130
      30 continue
      nn = min(2, nd)
      do i = 1, nn
         fn = fnu + real(nd-i,8)
         init = 0
         call f90_zunik(z, fn, 1, 0, tol, init, phi, zeta1, zeta2, sum, cwrk)
         if (kode==1) goto 40
         st = z + zeta2
         rast = fn/abs(st)
         st=conjg(st)*rast*rast
         s1=-zeta1+st+cmplx(0.0d0,z%im,8)
         goto 50
         40 continue
         s1 = -zeta1 + zeta2
         50 continue
         !-----------------------------------------------------------------------
         !     test for underflow and overflow
         !-----------------------------------------------------------------------
         rs1 = s1%re
         if (abs(rs1)>elim) goto 110
         if (i==1) iflag = 2
         if (abs(rs1)<alim) goto 60
         !-----------------------------------------------------------------------
         !     refine  test and scale
         !-----------------------------------------------------------------------
         aphi = abs(phi)
         rs1 = rs1 + log(aphi)
         if (abs(rs1)>elim) goto 110
         if (i==1) iflag = 1
         if (rs1<0.0d0) goto 60
         if (i==1) iflag = 3
         60 continue
         !-----------------------------------------------------------------------
         !     scale s1 if cabs(s1)<ascle
         !-----------------------------------------------------------------------
         s2=phi*sum
         s1=exp(s1)*cssr(iflag)
         s2=s2*s1
         if (iflag/=1) goto 70
         nw=f90_zuchk(s2, bry(1), tol)
         if (nw/=0) goto 110
         70 continue
         cy(i) = s2
         m = nd - i + 1
         y(m) = s2*csrr(iflag)
      end do
      if (nd<=2) goto 100
      rast = 1.0d0/abs(z)
      st=conjg(z)*rast
      rz=(st+st)*rast
      bry(2) = 1.0d0/bry(1)
      bry(3) = d1mach2
      s1 = cy(1)
      s2 = cy(2)
      c1r = csrr(iflag)
      ascle = bry(iflag)
      k = nd - 2
      fn = dble(float(k))
      do i = 3, nd
         c2 = s2
         s2 = s1 + (fnu+fn)*(rz*c2)
         s1 = c2
         c2 = s2*c1r
         y(k) = c2
         k = k - 1
         fn = fn - 1.0d0
         if (iflag>=3) goto 90
         c2m = max(abs(c2%re), abs(c2%im))
         if (c2m<=ascle) goto 90
         iflag = iflag + 1
         ascle = bry(iflag)
         s1 = s1*c1r
         s2 = c2
         s1 = s1*cssr(iflag)
         s2 = s2*cssr(iflag)
         c1r = csrr(iflag)
      90 end do
      100 continue
      return
      !-----------------------------------------------------------------------
      !     set underflow and update parameters
      !-----------------------------------------------------------------------
      110 continue
      if (rs1>0.0d0) goto 120
      y(nd) = 0.0d0
      nz = nz + 1
      nd = nd - 1
      if (nd==0) goto 100
      call f90_zuoik(z, fnu, kode, 1, nd, y, nuf, tol, elim, alim)
      if (nuf<0) goto 120
      nd = nd - nuf
      nz = nz + nuf
      if (nd==0) goto 100
      fn = fnu + real(nd-1,8)
      if (fn>=fnul) goto 30
      nlast = nd
      return
      120 continue
      nz = -1
      return
      130 continue
      if (rs1>0.0d0) goto 120
      nz = n
      do i = 1, n
         y(i) = 0.0d0
      end do
      return
   end subroutine f90_zuni1


   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zuni2(z, fnu, kode, n, y, nz, nlast, fnul, tol, elim, alim)
      !***begin prologue  zuni2
      !***refer to  zbesi,zbesk
      !
      !     zuni2 computes i(fnu,z) in the right half plane by means of
      !     uniform asymptotic expansion for j(fnu,zn) where zn is z*i
      !     or -z*i and zn is in the right half plane also.
      !
      !     fnul is the smallest order permitted for the asymptotic
      !     expansion. nlast=0 means all of the y values were set.
      !     nlast/=0 is the number left to be computed by another
      !     formula for orders fnu to fnu+nlast-1 because fnu+nlast-1<fnul.
      !     y(i)=czero for i=nlast+1,n
      !
      !***routines called  zairy,zuchk,zunhj,zuoik,d1mach,azabs
      !***end prologue  zuni2
      !     complex ai,arg,asum,bsum,cfn,ci,cid,cip,cone,crsc,cscl,csr,css,
      !    *czero,c1,c2,dai,phi,rz,s1,s2,y,z,zb,zeta1,zeta2,zn
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim,fnul
      integer,intent(in)::n,kode
      integer,intent(inout)::nz,nlast
      complex(8),intent(inout)::y(n)
      real(8)::aarg, ang, aphi, ascle, cidi, crsc, cscl, c1r, c2m, fn, rast, raz, rs1
      integer i, iflag, in, inu, j, k, nai, nd, ndai,nn, nuf, nw, idum
      complex(8),parameter::cip(4)=[complex(8)::(1.0d0,0.0d0),(0.0d0,1.0d0),(-1.0d0,0.0d0),(0.0d0,-1.0d0)]
      real(8),parameter::aic=1.265512123484645396d+00
      real(8)::csrr(3),cssr(3),bry(3)
      complex(8)::cy(2),ca,c2,zn,s1,s2,zb,asum,bsum,zeta1,zeta2,st,rz,ai,arg,phi,dai
      !
      !DEBUG write (*, *) 'call f90_zuni2,not test'
      nz = 0
      nd = n
      nlast = 0
      !-----------------------------------------------------------------------
      !     computed values with exponents between alim and elim in mag-
      !     nitude are scaled to keep intermediate arithmetic on scale,
      !     exp(alim)=exp(elim)*tol
      !-----------------------------------------------------------------------
      cscl = 1.0d0/tol
      crsc = tol
      cssr(1) = cscl
      cssr(2) = 1.0d0
      cssr(3) = crsc
      csrr(1) = crsc
      csrr(2) = 1.0d0
      csrr(3) = cscl
      bry(1) = 1.0d+3*d1mach1/tol
      !-----------------------------------------------------------------------
      !     zn is in the right half plane after rotation by ci or -ci
      !-----------------------------------------------------------------------
      zn%re = z%im
      zn%im = -z%re
      zb = z
      cidi = -1.0d0
      inu = int(fnu)
      ang = hpi*(fnu-real(inu,8))
      c2%re = cos(ang)
      c2%im = sin(ang)
      ca = c2
      in = inu + n - 1
      in = mod(in, 4) + 1
      c2=c2*cip(in)
      if (z%im>0.0d0) goto 10
      zn%re = -zn%re
      zb%im = -zb%im
      cidi  = -cidi
      c2%im = -c2%im
      10 continue
      !-----------------------------------------------------------------------
      !     check for underflow and overflow on first member
      !-----------------------------------------------------------------------
      fn = max1(fnu, 1.0d0)
      call f90_zunhj(zn, fn, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
      if (kode==1) goto 20
      st=zb+zeta2
      rast = fn/abs(st)
      st=conjg(st)*rast*rast
      s1 = -zeta1 + st
      goto 30
      20 continue
      s1 = -zeta1 + zeta2
      30 continue
      rs1 = s1%re
      if (abs(rs1)>elim) goto 150
      40 continue
      nn = min(2, nd)
      do i = 1, nn
         fn = fnu + real(nd-i,8)
         call f90_zunhj(zn, fnu, 0, tol, phi, arg, zeta1, zeta2, asum, bsum)
         if (kode==1) goto 50
         st = zb + zeta2
         rast = fn/abs(st)
         st = conjg(st)*rast*rast
         s1 = -zeta1 + st+cmplx(0.0d0,abs(z%im),8)
         goto 60
         50 continue
         s1 = -zeta1 + zeta2
         60 continue
         !-----------------------------------------------------------------------
         !     test for underflow and overflow
         !-----------------------------------------------------------------------
         rs1 = s1%re
         if (abs(rs1)>elim) goto 120
         if (i==1) iflag = 2
         if (abs(rs1)<alim) goto 70
         !-----------------------------------------------------------------------
         !     refine  test and scale
         !-----------------------------------------------------------------------
         !-----------------------------------------------------------------------
         aphi = abs(phi)
         aarg = abs(arg)
         rs1 = rs1 + log(aphi) - 0.25d0*log(aarg) - aic
         if (abs(rs1)>elim) goto 120
         if (i==1) iflag = 1
         if (rs1<0.0d0) goto 70
         if (i==1) iflag = 3
         70 continue
         !-----------------------------------------------------------------------
         !     scale s1 to keep intermediate arithmetic on scale near
         !     exponent extremes
         !-----------------------------------------------------------------------
         call f90_zairy(arg, 0, 2, ai, nai, idum)
         call f90_zairy(arg, 1, 2, dai, ndai, idum)
         st=dai*bsum
         st = st + (ai*asum)
         s2 = phi*st
         s1=exp(s1)*cssr(iflag)
         s2=s2*s1
         if (iflag/=1) goto 80
         nw=f90_zuchk(s2, bry(1), tol)
         if (nw/=0) goto 120
         80 continue
         if (z%im<=0.0d0) s2%im = -s2%im
         s2=s2*c2
         cy(i) = s2
         j = nd - i + 1
         y(j) = s2*csrr(iflag)
         c2=cmplx(-c2%im,c2%re,8)*cidi
      end do
      if (nd<=2) goto 110
      raz = 1.0d0/abs(z)
      st=conjg(z)*raz
      rz = (st+st)*raz
      bry(2) = 1.0d0/bry(1)
      bry(3) = d1mach2
      s1 = cy(1)
      s2 = cy(2)
      c1r = csrr(iflag)
      ascle = bry(iflag)
      k = nd - 2
      fn = real(k,8)
      do i = 3, nd
         c2 = s2
         s2 = s1 + (fnu+fn)*(rz*c2)
         s1 = c2
         c2 = s2*c1r
         y(k) = c2
         k = k - 1
         fn = fn - 1.0d0
         if (iflag>=3) goto 100
         c2m = max(abs(c2%re), abs(c2%im))
         if (c2m<=ascle) goto 100
         iflag = iflag + 1
         ascle = bry(iflag)
         s1 = s1*c1r
         s2 = c2
         s1 = s1*cssr(iflag)
         s2 = s2*cssr(iflag)
         c1r = csrr(iflag)
      100 end do
      110 continue
      return
      120 continue
      if (rs1>0.0d0) goto 140
      !-----------------------------------------------------------------------
      !     set underflow and update parameters
      !-----------------------------------------------------------------------
      y(nd) = 0.0d0
      nz = nz + 1
      nd = nd - 1
      if (nd==0) goto 110
      call f90_zuoik(z, fnu, kode, 1, nd, y, nuf, tol, elim, alim)
      if (nuf<0) goto 140
      nd = nd - nuf
      nz = nz + nuf
      if (nd==0) goto 110
      fn = fnu + real(nd-1,8)
      if (fn<fnul) goto 130
      !      fn = cidi
      !      j = nuf + 1
      !      k = mod(j,4) + 1
      !      s1 = cip(k)
      !      if (fn<0.0d0) s1%im = -s1%im
      !      c2=c2*s1
      in = inu + nd - 1
      in = mod(in, 4) + 1
      c2=ca*cip(in)
      if (z%im<=0.0d0) c2%im = -c2%im
      goto 40
      130 continue
      nlast = nd
      return
      140 continue
      nz = -1
      return
      150 continue
      if (rs1>0.0d0) goto 140
      nz = n
      do i = 1, n
         y(i) = 0.0d0
      end do
      return
   end subroutine f90_zuni2

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbuni(z, fnu, kode, n, y, nz, nui, nlast, fnul, tol, elim, alim)
      !***begin prologue  zbuni
      !***refer to  zbesi,zbesk
      !
      !     zbuni computes the i bessel function for large cabs(z)>
      !     fnul and fnu+n-1<fnul. the order is increased from
      !     fnu+n-1 greater than fnul by adding nui and computing
      !     according to the uniform asymptotic expansion for i(fnu,z)
      !     on iform=1 and the expansion for j(fnu,z) on iform=2
      !
      !***routines called  zuni1,zuni2,azabs,d1mach
      !***end prologue  zbuni
      !     complex cscl,cscr,cy,rz,st,s1,s2,y,z
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim,fnul
      integer,intent(in)::n,kode,nui
      integer,intent(inout)::nz,nlast
      complex(8),intent(inout)::y(n)
      real(8)::ax, ay, csclr, cscrr, dfnu, fnui, gnu, raz, str,ascle, c1m
      integer::i, iflag, iform, k,  nl,nw
      complex(8)::cy(2),st,rz,s1,s2
      real(8)::bry(3)
      !DEBUG write(*,*)"call f90_zbuni, not test"
      nz = 0
      ax = abs(z%re)*1.7321d0
      ay = abs(z%im)
      iform = 1
      if (ay>ax) iform = 2
      if (nui==0) goto 60
      fnui = real(nui,8)
      dfnu = fnu + real(n-1,8)
      gnu = dfnu + fnui
      if (iform==2) goto 10
      !-----------------------------------------------------------------------
      !     asymptotic expansion for i(fnu,z) for large fnu applied in
      !     -pi/3<=arg(z)<=pi/3
      !-----------------------------------------------------------------------
      call f90_zuni1(z, gnu, kode, 2, cy, nw, nlast, fnul, tol, elim, alim)
      goto 20
      10 continue
      !-----------------------------------------------------------------------
      !     asymptotic expansion for j(fnu,z*exp(m*hpi)) for large fnu
      !     applied in pi/3<abs(arg(z))<=pi/2 where m=+i or -i
      !     and hpi=pi/2
      !-----------------------------------------------------------------------
      call f90_zuni2(z, gnu, kode, 2, cy, nw, nlast, fnul, tol, elim, alim)
      20 continue
      if (nw<0) goto 50
      if (nw/=0) goto 90
      str = abs(cy(1))
      !----------------------------------------------------------------------
      !     scale backward recurrence, bry(3) is defined but never used
      !----------------------------------------------------------------------
      bry(1) = 1.0d+3*d1mach1/tol
      bry(2) = 1.0d0/bry(1)
      bry(3) = bry(2)
      iflag = 2
      ascle = bry(2)
      csclr = 1.0d0
      if (str>bry(1)) goto 21
      iflag = 1
      ascle = bry(1)
      csclr = 1.0d0/tol
      goto 25
      21 continue
      if (str<bry(2)) goto 25
      iflag = 3
      ascle = bry(3)
      csclr = tol
      25 continue
      cscrr = 1.0d0/csclr
      s1 = cy(2)*csclr
      s2 = cy(1)*csclr
      raz = 1.0d0/abs(z)
      st = conjg(z)*raz
      rz = (st+st)*raz
      do i = 1, nui
      st = s2
      s2 = (dfnu+fnui)*(rz*st) + s1
      s1 = st
      fnui = fnui - 1.0d0
      if (iflag>=3) goto 30
      st = s2*cscrr
      c1m = max(abs(st%re), abs(st%im))
      if (c1m<=ascle) goto 30
      iflag = iflag + 1
      ascle = bry(iflag)
      s1 = s1*cscrr
      s2 = st
      csclr = csclr*tol
      cscrr = 1.0d0/csclr
      s1 = s1*csclr
      s2 = s2*csclr
      30 end do
      y(n) = s2*cscrr
      if (n==1) return
      nl = n - 1
      fnui = real(nl,8)
      k = nl
      do i = 1, nl
      st = s2
      s2 = (fnu+fnui)*(rz*st) + s1
      s1 = st
      st = s2*cscrr
      y(k) = st
      fnui = fnui - 1.0d0
      k = k - 1
      if (iflag>=3) goto 40
      c1m = max(abs(st%re), abs(st%im))
      if (c1m<=ascle) goto 40
      iflag = iflag + 1
      ascle = bry(iflag)
      s1 = s1*cscrr
      s2 = st
      csclr = csclr*tol
      cscrr = 1.0d0/csclr
      s1 = s1*csclr
      s2 = s2*csclr
      40 end do
      return
      50 continue
      nz = -1
      if (nw==(-2)) nz = -2
      return
      60 continue
      if (iform==2) goto 70
      !-----------------------------------------------------------------------
      !     asymptotic expansion for i(fnu,z) for large fnu applied in
      !     -pi/3<=arg(z)<=pi/3
      !-----------------------------------------------------------------------
      call f90_zuni1(z, fnu, kode, n, y, nw, nlast, fnul, tol, elim, alim)
      goto 80
      70 continue
      !-----------------------------------------------------------------------
      !     asymptotic expansion for j(fnu,z*exp(m*hpi)) for large fnu
      !     applied in pi/3<abs(arg(z))<=pi/2 where m=+i or -i
      !     and hpi=pi/2
      !-----------------------------------------------------------------------
      call f90_zuni2(z, fnu, kode, n, y, nw, nlast, fnul, tol, elim, alim)
      80 continue
      if (nw<0) goto 50
      nz = nw
      return
      90 continue
      nlast = n
      return
   end subroutine f90_zbuni


   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbinu(z, fnu, kode, n, cy, nz, rl, fnul, tol, elim, alim)
      !***begin prologue  zbinu
      !***refer to  zbesh,zbesi,zbesj,zbesk,zairy,zbiry
      !
      !     zbinu computes the i function in the right half z plane
      !
      !***routines called  azabs,zasyi,zbuni,zmlri,zseri,zuoik,zwrsk
      !***end prologue  zbinu
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim,fnul,rl
      integer,intent(in)::n,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::cy(n)
      real(8)::az,dfnu
      integer::i,inw,nlast,nn,nui,nw
      complex(8)::cw(2)
      !DEBUG write (*, *) 'call f90_zbinu'
      !
      nz = 0
      az = abs(z)
      nn = n
      dfnu = fnu + real(n-1,8)
      if (az<=2.0d0) goto 10
      if (az*az*0.25d0>dfnu+1.0d0) goto 20
      10 continue
      !-----------------------------------------------------------------------
      !     power series
      !-----------------------------------------------------------------------
      call f90_zseri(z, fnu, kode, nn, cy, nw, tol, elim, alim)

      inw = abs(nw)
      nz = nz + inw
      nn = nn - inw
      if (nn==0) return
      if (nw>=0) goto 120
      dfnu = fnu + real(nn-1,8)
      20 continue
      if (az<rl) goto 40
      if (dfnu<=1.0d0) goto 30
      if (az+az<dfnu*dfnu) goto 50
      !-----------------------------------------------------------------------
      !     asymptotic expansion for large z
      !-----------------------------------------------------------------------
      30 continue
      call f90_zasyi(z, fnu, kode, nn, cy, nw, rl, tol, elim, alim)
      if (nw<0) goto 130
      goto 120
      40 continue
      if (dfnu<=1.0d0) goto 70
      50 continue
      !-----------------------------------------------------------------------
      !     overflow and underflow test on i sequence for miller algorithm
      !-----------------------------------------------------------------------
      call f90_zuoik(z, fnu, kode, 1, nn, cy, nw, tol, elim, alim)
      if (nw<0) goto 130
      nz = nz + nw
      nn = nn - nw
      if (nn==0) return
      dfnu = fnu + real(nn-1,8)
      if (dfnu>fnul) goto 110
      if (az>fnul) goto 110
      60 continue
      if (az>rl) goto 80
      70 continue
      !-----------------------------------------------------------------------
      !     miller algorithm normalized by the series
      !-----------------------------------------------------------------------
      call f90_zmlri(z, fnu, kode, nn, cy, nw, tol)
      if (nw<0) goto 130
      goto 120
      80 continue
      !-----------------------------------------------------------------------
      !     miller algorithm normalized by the wronskian
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      !     overflow test on k functions used in wronskian
      !-----------------------------------------------------------------------
      call f90_zuoik(z, fnu, kode, 2, 2, cw, nw, tol, elim, alim)
      if (nw>=0) goto 100
      nz = nn
      do i = 1, nn
         cy(i) = 0.0d0
      end do
      return
      100 continue
      if (nw>0) goto 130
      call f90_zwrsk(z, fnu, kode, nn, cy, nw, cw, tol, elim, alim)
      if (nw<0) goto 130
      goto 120
      110 continue
      !-----------------------------------------------------------------------
      !     increment fnu+nn-1 up to fnul, compute and recur backward
      !-----------------------------------------------------------------------
      nui = int(fnul-dfnu) + 1
      nui = max(nui, 0)
      call f90_zbuni(z, fnu, kode, nn, cy, nw, nui, nlast, fnul, tol, elim, alim)
      if (nw<0) goto 130
      nz = nz + nw
      if (nlast==0) goto 120
      nn = nlast
      goto 60
      120 continue
      return
      130 continue
      nz = -1
      if (nw==(-2)) nz = -2
      return
   end subroutine f90_zbinu

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zacon(z, fnu, kode, mr, n, y, nz, rl, fnul, tol, elim, alim)
      !***begin prologue  zacon
      !***refer to  zbesk,zbesh
      !
      !     zacon applies the analytic continuation formula
      !
      !         k(fnu,zn*exp(mp))=k(fnu,zn)*exp(-mp*fnu) - mp*i(fnu,zn)
      !                 mp=dpi*mr*cmplx(0.0,1.0)
      !
      !     to continue the k function from the right half to the left
      !     half z plane
      !
      !***routines called  zbinu,zbknu,zs1s2,d1mach,azabs,zmlt
      !***end prologue  zacon
      !     complex ck,cone,cscl,cscr,csgn,cspn,cy,czero,c1,c2,rz,sc1,sc2,st,
      !    *s1,s2,y,z,zn
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim,fnul,rl
      integer,intent(in)::n,kode,mr
      integer,intent(inout)::nz
      complex(8),intent(inout)::y(n)
      real(8)::arg, ascle, as2, azn, bscle, csr, c1m,  fmr, fn, razn, sgn, yy
      integer i, inu, iuf, kflag, nn, nw
      real(8)::bry(3),cssr(3),csrr(3),cscl,cscr
      complex(8)::ck,csgn,cspn,c1,c2,rz,sc1,sc2,st,s1,s2,zn,cp,cy(2)
      !DEBUG write (*, *) 'call f90_zacon,not test'
      nz = 0
      zn = -z
      nn = n
      call f90_zbinu(zn, fnu, kode, nn, y, nw, rl, fnul, tol, elim, alim)
      if (nw<0) goto 90
      !-----------------------------------------------------------------------
      !     analytic continuation to the left half plane for the k function
      !-----------------------------------------------------------------------
      nn = min(2, n)
      call f90_zbknu(zn, fnu, kode, nn, cy, nw, tol, elim, alim)
      if (nw/=0) goto 90
      s1 = cy(1)
      fmr = real(mr,8)
      sgn = -sign(dpi, fmr)
      csgn=cmplx(0.0d0,sgn,8)
      if (kode==1) goto 10
      yy = -zn%im
      cp=cmplx(cos(yy),sin(yy),8)
      csgn=csgn*cp
      10 continue
      !-----------------------------------------------------------------------
      !     calculate cspn=exp(fnu*dpi*i) to minimize losses of significance
      !     when fnu is large
      !-----------------------------------------------------------------------
      inu = int(fnu)
      arg = (fnu-real(inu,8))*sgn
      cp=cmplx(cos(arg),sin(arg),8)
      cspn=cp
      if (mod(inu,2)==0) goto 20
      cspn = -cspn
      20 continue
      iuf = 0
      c1 = s1
      c2 = y(1)
      ascle = 1.0d+3*d1mach1/tol
      if (kode==1) goto 30
      call f90_zs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
      sc1 = c1
      30 continue
      y(1) = (cspn*c1) + (csgn*c2)
      if (n==1) return
      cspn = -cspn
      s2 = cy(2)
      c1 = s2
      c2 = y(2)
      if (kode==1) goto 40
      call f90_zs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
      sc2 = c1
      40 continue
      y(2) = (cspn*c1) + (csgn*c2)
      if (n==2) return
      cspn = -cspn
      azn = abs(zn)
      razn = 1.0d0/azn
      st = conjg(zn)*razn
      rz = (st+st)*razn
      fn = fnu + 1.0d0
      ck = fn*rz
      !-----------------------------------------------------------------------
      !     scale near exponent extremes during recurrence on k functions
      !-----------------------------------------------------------------------
      cscl = 1.0d0/tol
      cscr = tol
      cssr(1) = cscl
      cssr(2) = 1.0d0
      cssr(3) = cscr
      csrr(1) = cscr
      csrr(2) = 1.0d0
      csrr(3) = cscl
      bry(1) = ascle
      bry(2) = 1.0d0/ascle
      bry(3) = d1mach2
      as2 = abs(s2)
      kflag = 2
      if (as2>bry(1)) goto 50
      kflag = 1
      goto 60
      50 continue
      if (as2<bry(2)) goto 60
      kflag = 3
      60 continue
      bscle = bry(kflag)
      s1 = s1*cssr(kflag)
      s2 = s2*cssr(kflag)
      csr = csrr(kflag)
      do i = 3, n
      st = s2
      s2 = ck*st + s1
      s1 = st
      c1 = s2*csr
      st = c1
      c2 = y(i)
      if (kode==1) goto 70
      if (iuf<0) goto 70
      call f90_zs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
      sc1 = sc2
      sc2 = c1
      if (iuf/=3) goto 70
      iuf = -4
      s1 = sc1*cssr(kflag)
      s2 = sc2*cssr(kflag)
      st = sc2
      70 continue
      y(i) = (cspn*c1) + csgn*c2
      ck = ck + rz
      cspn = -cspn
      if (kflag>=3) goto 80
      c1m = max(abs(c1%re), abs(c1%im))
      if (c1m<=bscle) goto 80
      kflag = kflag + 1
      bscle = bry(kflag)
      s1 = s1*csr
      s2 = st
      s1 = s1*cssr(kflag)
      s2 = s2*cssr(kflag)
      csr = csrr(kflag)
      80 end do
      return
      90 continue
      nz = -1
      if (nw==(-2)) nz = -2
      return
   end subroutine f90_zacon

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbiry(z, id, kode, bi, ierr)
      !***begin prologue  zbiry
      !***date written   830501   (yymmdd)
      !***revision date  890801   (yymmdd)
      !***category no.  b5k
      !***keywords  airy function,bessel functions of order one third
      !***author  amos, donald e., sandia national laboratories
      !***purpose  to compute airy functions bi(z) and dbi(z) for complex z
      !***description
      !
      !                      ***a real(8) routine***
      !         on kode=1, cbiry computes the complex airy function bi(z) or
      !         its derivative dbi(z)/dz on id=0 or id=1 respectively. on
      !         kode=2, a scaling option cexp(-axzta)*bi(z) or cexp(-axzta)*
      !         dbi(z)/dz is provided to remove the exponential behavior in
      !         both the left and right half planes where
      !         zta=(2/3)*z*csqrt(z)=cmplx(xzta,yzta) and axzta=abs(xzta).
      !         defintions and notation are found in the nbs handbook of
      !         mathematical functions (ref. 1).
      !
      !         input      zr,zi are real(8)
      !           zr,zi  - z=cmplx(zr,zi)
      !           id     - order of derivative, id=0 or id=1
      !           kode   - a parameter to indicate the scaling option
      !                    kode= 1  returns
      !                             bi=bi(z)                 on id=0 or
      !                             bi=dbi(z)/dz             on id=1
      !                        = 2  returns
      !                             bi=cexp(-axzta)*bi(z)     on id=0 or
      !                             bi=cexp(-axzta)*dbi(z)/dz on id=1 where
      !                             zta=(2/3)*z*csqrt(z)=cmplx(xzta,yzta)
      !                             and axzta=abs(xzta)
      !
      !         output     bir,bii are real(8)
      !           bir,bii- complex answer depending on the choices for id and
      !                    kode
      !           ierr   - error flag
      !                    ierr=0, normal return - computation completed
      !                    ierr=1, input error   - no computation
      !                    ierr=2, overflow      - no computation, real(z)
      !                            too large on kode=1
      !                    ierr=3, cabs(z) large      - computation completed
      !                            losses of signifcance by argument reduction
      !                            produce less than half of machine accuracy
      !                    ierr=4, cabs(z) too large  - no computation
      !                            complete loss of accuracy by argument
      !                            reduction
      !                    ierr=5, error              - no computation,
      !                            algorithm termination condition not met
      !
      !***long description
      !
      !         bi and dbi are computed for cabs(z)>1.0 from the i bessel
      !         functions by
      !
      !                bi(z)=c*sqrt(z)*( i(-1/3,zta) + i(1/3,zta) )
      !               dbi(z)=c *  z  * ( i(-2/3,zta) + i(2/3,zta) )
      !                               c=1.0/sqrt(3.0)
      !                             zta=(2/3)*z**(3/2)
      !
      !         with the power series for cabs(z)<=1.0.
      !
      !         in most complex variable computation, one must evaluate ele-
      !         mentary functions. when the magnitude of z is large, losses
      !         of significance by argument reduction occur. consequently, if
      !         the magnitude of zeta=(2/3)*z**1.5 exceeds u1=sqrt(0.5/ur),
      !         then losses exceeding half precision are likely and an error
      !         flag ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
      !         real(8) unit roundoff limited to 18 digits precision.
      !         also, if the magnitude of zeta is larger than u2=0.5/ur, then
      !         all significance is lost and ierr=4. in order to use the int
      !         function, zeta must be further restricted not to exceed the
      !         largest integer, u3=i1mach(9). thus, the magnitude of zeta
      !         must be restricted by min(u2,u3). on 32 bit machines, u1,u2,
      !         and u3 are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single
      !         precision arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double
      !         precision arithmetic respectively. this makes u2 and u3 limit-
      !         ing in their respective arithmetics. this means that the mag-
      !         nitude of z cannot exceed 3.1e+4 in single and 2.1e+6 in
      !         real(8) arithmetic. this also means that one can
      !         expect to retain, in the worst cases on 32 bit machines,
      !         no digits in single precision and only 7 digits in double
      !         precision arithmetic. similar considerations hold for other
      !         machines.
      !
      !         the approximate relative error in the magnitude of a complex
      !         bessel function can be expressed by p*10**s where p=max(unit
      !         roundoff,1.0e-18) is the nominal precision and 10**s repre-
      !         sents the increase in error due to argument reduction in the
      !         elementary functions. here, s=max(1,abs(log10(cabs(z))),
      !         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
      !         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
      !         have only absolute accuracy. this is most likely to occur when
      !         one component (in absolute value) is larger than the other by
      !         several orders of magnitude. if one component is 10**k larger
      !         than the other, then one can expect only max(abs(log10(p))-k,
      !         0) significant digits; or, stated another way, when k exceeds
      !         the exponent of p, no significant digits remain in the smaller
      !         component. however, the phase angle retains absolute accuracy
      !         because, in complex arithmetic with precision p, the smaller
      !         component will not (as a rule) decrease below p times the
      !         magnitude of the larger component. in these extreme cases,
      !         the principal phase angle is on the order of +p, -p, dpi/2-p,
      !         or -dpi/2+p.
      !
      !***references  handbook of mathematical functions by m. abramowitz
      !                 and i. a. stegun, nbs ams series 55, u.s. dept. of
      !                 commerce, 1955.
      !
      !               computation of bessel functions of complex argument
      !                 and large order by d. e. amos, sand83-0643, may, 1983
      !
      !               a subroutine package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, sand85-
      !                 1018, may, 1985
      !
      !               a portable package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, trans.
      !                 math. software, 1986
      !
      !***routines called  zbinu,azabs,zdiv,azsqrt,d1mach,i1mach
      !***end prologue  zbiry
      !     complex bi,cone,csq,cy,s1,s2,trm1,trm2,z,zta,z3
      complex(8),intent(in)::z
      integer,intent(in)::id,kode
      complex(8),intent(inout)::bi
      integer,intent(inout)::ierr
      real(8)::aa, ad, ak, alim, atrm, az, az3, bb, &
         bk, cc, ck, &
         dig, dk, d1, d2, eaa, elim, fid, fmr,&
         fnu, fnul, rl, sfac, tol, z3r
      integer ::k, nz
      complex(8)::csq,cy(2),s1,s2,trm1,trm2,zta,z3,st
      real(8),parameter::c1  =6.14926627446000736d-01
      real(8),parameter::c2  =4.48288357353826359d-01 
      real(8),parameter::coef=5.77350269189625765d-01
      !***first executable statement  zbiry
      !DEBUG write (*, *) 'call f90_zbiry.not test'
      ierr = 0
      nz = 0
      if (id<0 .or. id>1) ierr = 1
      if (kode<1 .or. kode>2) ierr = 1
      if (ierr/=0) return
      az = abs(z)
      tol = max(d1mach4, 1.0d-18)
      fid = real(id,8)
      if (az>1.0d0) goto 70
      !-----------------------------------------------------------------------
      !     power series for cabs(z)<=1.
      !-----------------------------------------------------------------------
      s1 = 1.0d0
      s2 = 1.0d0
      if (az<tol) goto 130
      aa = az*az
      if (aa<tol/az) goto 40
      trm1 = 1.0d0
      trm2 = 1.0d0
      atrm = 1.0d0
      z3=z**3
      az3 = az*aa
      ak = 2.0d0 + fid
      bk = 3.0d0 - fid - fid
      ck = 4.0d0 - fid
      dk = 3.0d0 + fid + fid
      d1 = ak*dk
      d2 = bk*ck
      ad = min(d1, d2)
      ak = 24.0d0 + 9.0d0*fid
      bk = 30.0d0 - 9.0d0*fid
      do k = 1, 25
      trm1=trm1*z3/d1
      s1 = s1 + trm1
      trm2=trm2*z3/d2
      s2 = s2 + trm2
      atrm = atrm*az3/ad
      d1 = d1 + ak
      d2 = d2 + bk
      ad = min(d1, d2)
      if (atrm<tol*ad) goto 40
      ak = ak + 18.0d0
      bk = bk + 18.0d0
      end do
      40 continue
      if (id==1) goto 50
      bi = c1*s1 + c2*(z*s2)
      if (kode==1) return
      st=sqrt(z)
      zta = tth*(z*st)
      aa = zta%re
      aa = -abs(aa)
      eaa = exp(aa)
      bi = bi*eaa
      return
      50 continue
      bi = s2*c2
      if (az<=tol) goto 60
      cc = c1/(1.0d0+fid)
      st = s1*z
      bi = bi + cc*(st*z)
      60 continue
      if (kode==1) return
      st=sqrt(z)
      zta = tth*(z*st)
      aa = zta%re
      aa = -abs(aa)
      eaa = exp(aa)
      bi = bi*eaa
      return
      !-----------------------------------------------------------------------
      !     case for cabs(z)>1.0
      !-----------------------------------------------------------------------
      70 continue
      fnu = (1.0d0+fid)/3.0d0
      !-----------------------------------------------------------------------
      !     set parameters related to machine constants.
      !     tol is the approximate unit roundoff limited to 1.0e-18.
      !     elim is the approximate exponential over- and underflow limit.
      !     exp(-elim)<exp(-alim)=exp(-elim)/tol    and
      !     exp(elim)>exp(alim)=exp(elim)*tol       are intervals near
      !     underflow and overflow limits where scaled arithmetic is done.
      !     rl is the lower boundary of the asymptotic expansion for large z.
      !     dig = number of base 10 digits in tol = 10**(-dig).
      !     fnul is the lower boundary of the asymptotic series for large fnu.
      !-----------------------------------------------------------------------
      k = min(abs(i1mach15), (i1mach16))
      elim = 2.303d0*(real(k,8)*d1mach5-3.0d0)
      aa = d1mach5*real(i1mach14 - 1,8)
      dig = dmin1(aa, 18.0d0)
      aa = aa*2.303d0
      alim = elim + max(-aa, -41.45d0)
      rl = 1.2d0*dig + 3.0d0
      fnul = 10.0d0 + 6.0d0*(dig-3.0d0)
      !-----------------------------------------------------------------------
      !     test for range
      !-----------------------------------------------------------------------
      aa = 0.5d0/tol
      bb = real(i1mach9,8)*0.5d0
      aa = min(aa, bb)
      aa = aa**tth
      if (az>aa) goto 260
      aa = sqrt(aa)
      if (az>aa) ierr = 3
      csq=sqrt(z)
      zta = tth*z*csq
      !-----------------------------------------------------------------------
      !     re(zta)<=0 when re(z)<0, especially when im(z) is small
      !-----------------------------------------------------------------------
      sfac = 1.0d0
      ak = zta%im
      if (z%re>=0.0d0) goto 80
      bk = zta%re
      ck = -abs(bk)
      zta=cmplx(ck,ak,8)
      80 continue
      if (.not.iszero(z%im) .or. z%re>0.0d0) goto 90
      zta=cmplx(0.0d0,ak,8)
      90 continue
      aa = zta%re
      if (kode==2) goto 100
      !-----------------------------------------------------------------------
      !     overflow test
      !-----------------------------------------------------------------------
      bb = abs(aa)
      if (bb<alim) goto 100
      bb = bb + 0.25d0*log(az)
      sfac = tol
      if (bb>elim) goto 190
      100 continue
      fmr = 0.0d0
      if (aa>=0.0d0 .and. z%re>0.0d0) goto 110
      fmr = dpi
      if (z%im<0.0d0) fmr = -dpi
      zta = -zta
      110 continue
      !-----------------------------------------------------------------------
      !     aa=factor for analytic continuation of i(fnu,zta)
      !     kode=2 returns exp(-abs(xzta))*i(fnu,zta) from cbesi
      !-----------------------------------------------------------------------
      call f90_zbinu(zta, fnu, kode, 1, cy, nz, rl, fnul, tol, elim, alim)
      if (nz<0) goto 200
      aa = fmr*fnu
      z3r = sfac
      st = cmplx(cos(aa),sin(aa),8)
      s1 = (st*cy(1))*z3r
      fnu = (2.0d0-fid)/3.0d0
      call f90_zbinu(zta, fnu, kode, 2, cy, nz, rl, fnul, tol, elim, alim)
      cy(1) = cy(1)*z3r
      cy(2) = cy(2)*z3r
      !-----------------------------------------------------------------------
      !     backward recur one step for orders -1/3 or -2/3
      !-----------------------------------------------------------------------
      st=cy(1)/zta
      s2 = (fnu+fnu)*st + cy(2)
      aa = fmr*(fnu-1.0d0)
      st=cmplx(cos(aa),sin(aa),8)
      s1 = coef*(s1+s2*st)
      if (id==1) goto 120
      s1=csq*s1
      bi = s1/sfac
      return
      120 continue
      s1=z*s1
      bi = s1/sfac
      return
      130 continue
      aa = c1*(1.0d0-fid) + fid*c2
      bi = aa
      return
      190 continue
      ierr = 2
      nz = 0
      return
      200 continue
      if (nz==(-1)) goto 190
      nz = 0
      ierr = 5
      return
      260 continue
      ierr = 4
      nz = 0
      return
   end subroutine f90_zbiry

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zunk1(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
      !***begin prologue  zunk1
      !***refer to  zbesk
      !
      !     zunk1 computes k(fnu,z) and its analytic continuation from the
      !     right half plane to the left half plane by means of the
      !     uniform asymptotic expansion.
      !     mr indicates the direction of rotation for analytic continuation.
      !     nz=-1 means an overflow will occur
      !
      !***routines called  zkscl,zs1s2,zuchk,zunik,d1mach,azabs
      !***end prologue  zunk1
      !     complex cfn,ck,cone,crsc,cs,cscl,csgn,cspn,csr,css,cwrk,cy,czero,
      !    *c1,c2,phi,phid,rz,sum,sumd,s1,s2,y,z,zeta1,zeta1d,zeta2,zeta2d,zr
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim
      integer,intent(in)::n,kode,mr
      integer,intent(inout)::nz
      complex(8),intent(inout)::y(n)
      real(8)::ang, aphi, asc, ascle, csgni,&
          c1r, c2m,fmr, fn, fnf, rast, razr, rs1,sgn,cscl,crsc
      integer i, ib, iflag, ifn, il, inu, iuf, k, kdflg, kflag, kk, nw, initd, ic, ipard, j,m
      real(8)::bry(3), cssr(3),csrr(3)
      integer::init(2)
      complex(8)::sum(2),zeta1(2),zeta2(2),cy(2),cwrk(16,3),phi(2)
      complex(8)::ck,cspn,csr,c1,c2,phid,rz,sumd,s1,s2,zr
      complex(8)::st,zet1d,zet2d
      !
      !DEBUG write (*, *) 'call f90_zunk1,not test'
      kdflg = 1
      nz = 0
      !-----------------------------------------------------------------------
      !     exp(-alim)=exp(-elim)/tol=approx. one precision greater than
      !     the underflow limit
      !-----------------------------------------------------------------------
      cscl = 1.0d0/tol
      crsc = tol
      cssr(1) = cscl
      cssr(2) = 1.0d0
      cssr(3) = crsc
      csrr(1) = crsc
      csrr(2) = 1.0d0
      csrr(3) = cscl
      bry(1) = 1.0d+3*d1mach1/tol
      bry(2) = 1.0d0/bry(1)
      bry(3) = d1mach2
      zr = z
      if (z%re>=0.0d0) goto 10
      zr = -z
      10 continue
      j = 2
      do i = 1, n
         !-----------------------------------------------------------------------
         !     j flip flops between 1 and 2 in j = 3 - j
         !-----------------------------------------------------------------------
         j = 3 - j
         fn = fnu + real(i-1,8)
         init(j) = 0
         call f90_zunik(zr, fn, 2, 0, tol, init(j), phi(j), zeta1(j), zeta2(j), sum(j), cwrk(:,j))
         if (kode==1) goto 20
         st = zr + zeta2(j)
         rast = fn/abs(st)
         st=conjg(st)*rast*rast
         s1 = zeta1(j) - st
         goto 30
         20 continue
         s1 = zeta1(j) - zeta2(j)
         30 continue
         rs1 = s1%re
         !-----------------------------------------------------------------------
         !     test for underflow and overflow
         !-----------------------------------------------------------------------
         if (abs(rs1)>elim) goto 60
         if (kdflg==1) kflag = 2
         if (abs(rs1)<alim) goto 40
         !-----------------------------------------------------------------------
         !     refine  test and scale
         !-----------------------------------------------------------------------
         aphi = abs(phi(j))
         rs1 = rs1 + log(aphi)
         if (abs(rs1)>elim) goto 60
         if (kdflg==1) kflag = 1
         if (rs1<0.0d0) goto 40
         if (kdflg==1) kflag = 3
         40 continue
         !-----------------------------------------------------------------------
         !     scale s1 to keep intermediate arithmetic on scale near
         !     exponent extremes
         !-----------------------------------------------------------------------
         s2 = phi(j)*sum(j)
         s1=cssr(kflag)*exp(s1)
         s2=s2*s1
         if (kflag/=1) goto 50
         nw=f90_zuchk(s2,bry(1),tol)
         if (nw/=0) goto 60
         50 continue
         cy(kdflg) = s2
         y(i) = s2*csrr(kflag)
         if (kdflg==2) goto 75
         kdflg = 2
         goto 70
         60 continue
         if (rs1>0.0d0) goto 300
         !-----------------------------------------------------------------------
         !     for zr<0.0, the i function to be added will overflow
         !-----------------------------------------------------------------------
         if (z%re<0.0d0) goto 300
         kdflg = 1
         y(i) = 0.0d0
         nz = nz + 1
         if (i==1) goto 70
         if (iszero(y(i-1))) goto 70
         y(i-1) = 0.0d0
         nz = nz + 1
      70 end do
      i = n
      75 continue
      razr = 1.0d0/abs(zr)
      st=conjg(zr)*razr
      rz = (st+st)*razr
      ck = fn*rz
      ib = i + 1
      if (n<ib) goto 160
      !-----------------------------------------------------------------------
      !     test last member for underflow and overflow. set sequence to zero
      !     on underflow.
      !-----------------------------------------------------------------------
      fn = fnu + real(n-1,8)
      ipard = 1
      if (mr/=0) ipard = 0
      initd = 0
      call f90_zunik(zr, fn, 2, ipard, tol, initd, phid, zet1d, zet2d, sumd, cwrk(:,3))
      if (kode==1) goto 80
      st = zr + zet2d
      rast = fn/abs(st)
      st=conjg(st)*rast*rast
      s1 = zet1d - st
      goto 90
      80 continue
      s1 = zet1d - zet2d
      90 continue
      rs1 = s1%re
      if (abs(rs1)>elim) goto 95
      if (abs(rs1)<alim) goto 100
      !----------------------------------------------------------------------------
      !     refine estimate and test
      !-------------------------------------------------------------------------
      aphi = abs(phid)
      rs1 = rs1 + log(aphi)
      if (abs(rs1)<elim) goto 100
      95 continue
      if (abs(rs1)>0.0d0) goto 300
      !-----------------------------------------------------------------------
      !     for zr<0.0, the i function to be added will overflow
      !-----------------------------------------------------------------------
      if (z%re<0.0d0) goto 300
      nz = n
      do i = 1, n
         y(i) = 0.0d0
      end do
      return
      !---------------------------------------------------------------------------
      !     forward recur for remainder of the sequence
      !----------------------------------------------------------------------------
      100 continue
      s1 = cy(1)
      s2 = cy(2)
      c1r = csrr(kflag)
      ascle = bry(kflag)
      do i = ib, n
         c2 = s2
         s2 = ck*c2+ s1
         s1 = c2
         ck = ck + rz
         c2 = s2*c1r
         y(i) = c2
         if (kflag>=3) goto 120
         c2m = max(abs(c2%re), abs(c2%im))
         if (c2m<=ascle) goto 120
         kflag = kflag + 1
         ascle = bry(kflag)
         s1 = s1*c1r
         s2 = c2
         s1 = s1*cssr(kflag)
         s2 = s2*cssr(kflag)
         c1r = csrr(kflag)
      120 end do
      160 continue
      if (mr==0) return
      !-----------------------------------------------------------------------
      !     analytic continuation for re(z)<0.0d0
      !-----------------------------------------------------------------------
      nz = 0
      fmr = real(mr,8)
      sgn = -dsign(dpi, fmr)
      !-----------------------------------------------------------------------
      !     cspn and csgn are coeff of k and i functions resp.
      !-----------------------------------------------------------------------
      csgni = sgn
      inu = int(fnu)
      fnf = fnu - real(inu,8)
      ifn = inu + n - 1
      ang = fnf*sgn
      cspn = cmplx(cos(ang),sin(ang),8)
      if (mod(ifn,2)==0) goto 170
      cspn = -cspn
      170 continue
      asc = bry(1)
      iuf = 0
      kk = n
      kdflg = 1
      ib = ib - 1
      ic = ib - 1
      do k = 1, n
         fn = fnu + real(kk-1,8)
         !-----------------------------------------------------------------------
         !     logic to sort out cases whose parameters were set for the k
         !     function above
         !-----------------------------------------------------------------------
         m = 3
         if (n>2) goto 175
         172 continue
         initd = init(j)
         phid = phi(j)
         zet1d = zeta1(j)
         zet2d = zeta2(j)
         sumd = sum(j)
         sumd = sum(j)
         m = j
         j = 3 - j
         goto 180
         175 continue
         if ((kk==n) .and. (ib<n)) goto 180
         if ((kk==ib) .or. (kk==ic)) goto 172
         initd = 0
         180 continue
         call f90_zunik(zr, fn, 1, 0, tol, initd, phid, zet1d, zet2d, sumd, cwrk(:,m))
         if (kode==1) goto 200
         st = zr + zet2d
         rast = fn/abs(st)
         st=conjg(st)*rast*rast
         s1 = -zet1d + st
         goto 210
         200 continue
         s1 = -zet1d + zet2d
         210 continue
         !-----------------------------------------------------------------------
         !     test for underflow and overflow
         !-----------------------------------------------------------------------
         rs1 = s1%re
         if (abs(rs1)>elim) goto 260
         if (kdflg==1) iflag = 2
         if (abs(rs1)<alim) goto 220
         !-----------------------------------------------------------------------
         !     refine  test and scale
         !-----------------------------------------------------------------------
         aphi = abs(phid)
         rs1 = rs1 + log(aphi)
         if (abs(rs1)>elim) goto 260
         if (kdflg==1) iflag = 1
         if (rs1<0.0d0) goto 220
         if (kdflg==1) iflag = 3
         220 continue
         st = phid*sumd
         s2=csgni*cmplx(-st%im,st%re,8)
         s1=exp(s1)*cssr(iflag)
         s2=s2*s1
         if (iflag/=1) goto 230
         nw=f90_zuchk(s2, bry(1), tol)
         if (nw==0) goto 230
         s2 = 0.0d0
         230 continue
         cy(kdflg) = s2
         c2 = s2
         s2 = s2*csrr(iflag)
         !-----------------------------------------------------------------------
         !     add i and k functions, k sequence in y(i), i=1,n
         !-----------------------------------------------------------------------
         s1 = y(kk)
         if (kode==1) goto 250
         call f90_zs1s2(zr, s1, s2, nw, asc, alim, iuf)
         nz = nz + nw
         250 continue
         y(kk) = s1*cspn + s2
         kk = kk - 1
         cspn = -cspn
         if (.not.iszero(c2)) goto 255
         kdflg = 1
         goto 270
         255 continue
         if (kdflg==2) goto 275
         kdflg = 2
         goto 270
         260 continue
         if (rs1>0.0d0) goto 300
         s2 = 0.0d0
         goto 230
      270 end do
      k = n
      275 continue
      il = n - k
      if (il==0) return
      !-----------------------------------------------------------------------
      !     recur backward for remainder of i sequence and add in the
      !     k functions, scaling the i sequence during recurrence to keep
      !     intermediate arithmetic on scale near exponent extremes.
      !-----------------------------------------------------------------------
      s1 = cy(1)
      s2 = cy(2)
      csr = csrr(iflag)
      ascle = bry(iflag)
      fn = real(inu+il,8)
      do i = 1, il
         c2 = s2
         s2 = s1 + (fn+fnf)*(rz*c2)
         s1 = c2
         fn = fn - 1.0d0
         c2 = s2*csr
         ck = c2
         c1 = y(kk)
         if (kode==1) goto 280
         call f90_zs1s2(zr, c1, c2, nw, asc, alim, iuf)
         nz = nz + nw
         280 continue
         y(kk) = c1*cspn + c2
         kk = kk - 1
         cspn = -cspn
         if (iflag>=3) goto 290
         c2m = max(abs(ck%re), abs(ck%im))
         if (c2m<=ascle) goto 290
         iflag = iflag + 1
         ascle = bry(iflag)
         s1 = s1*csr
         s2 = ck
         s1 = s1*cssr(iflag)
         s2 = s2*cssr(iflag)
         csr = csrr(iflag)
      290 end do
      return
      300 continue
      nz = -1
      return
   end subroutine f90_zunk1

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zunk2(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
      !***begin prologue  zunk2
      !***refer to  zbesk
      !
      !     zunk2 computes k(fnu,z) and its analytic continuation from the
      !     right half plane to the left half plane by means of the
      !     uniform asymptotic expansions for h(kind,fnu,zn) and j(fnu,zn)
      !     where zn is in the right half plane, kind=(3-mr)/2, mr=+1 or
      !     -1. here zn=zr*i or -zr*i where zr=z if z is in the right
      !     half plane or zr=-z if z is in the left half plane. mr indic-
      !     ates the direction of rotation for analytic continuation.
      !     nz=-1 means an overflow will occur
      !
      !***routines called  zairy,zkscl,zs1s2,zuchk,zunhj,d1mach,azabs
      !***end prologue  zunk2
      !     complex ai,arg,argd,asum,asumd,bsum,bsumd,cfn,ci,cip,ck,cone,crsc,
      !    *cr1,cr2,cs,cscl,csgn,cspn,csr,css,cy,czero,c1,c2,dai,phi,phid,rz,
      !    *s1,s2,y,z,zb,zeta1,zeta1d,zeta2,zeta2d,zn,zr
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim
      integer,intent(in)::n,kode,mr
      integer,intent(inout)::nz
      complex(8),intent(inout)::y(n)
      real(8)::aarg, ang, aphi, asc, ascle,&
          crsc, cscl, csgni,  csr, c1r,  c2m, fmr, fn, fnf,&
          rast, razr, rs1,  sgn, yy
      integer::i, ib, iflag, ifn, il, in, inu, iuf, k, kdflg, kflag, kk, nai, ndai, nw, idum, j, ipard, ic
      real(8)::bry(3),cssr(3),csrr(3)
      complex(8)::asum(2),bsum(2),phi(2),arg(2),zeta1(2),zeta2(2),cy(2)
      complex(8)::ca,zn,cs,s1,s2,zr,c2,ck,ai,asumd,bsumd,argd,c1,cspn,dai,phid,pt,rz,st,zb,zet1d,zet2d
      complex(8),parameter::cr1=(1.0d0, 1.73205080756887729d0)
      complex(8),parameter::cr2=(-0.5d0, -8.66025403784438647d-01)
      real(8),parameter::aic=1.26551212348464539d+00
      complex(8),parameter::cip(4)=[complex(8)::&
        ( 1.0d0,  0.0d0), ( 0.0d0, -1.0d0), (-1.0d0,0.0d0),( 0.0d0,  1.0d0)]
      !
      kdflg = 1
      nz = 0
      !-----------------------------------------------------------------------
      !     exp(-alim)=exp(-elim)/tol=approx. one precision greater than
      !     the underflow limit
      !-----------------------------------------------------------------------
      cscl = 1.0d0/tol
      crsc = tol
      cssr(1) = cscl
      cssr(2) = 1.0d0
      cssr(3) = crsc
      csrr(1) = crsc
      csrr(2) = 1.0d0
      csrr(3) = cscl
      bry(1) = 1.0d+3*d1mach1/tol
      bry(2) = 1.0d0/bry(1)
      bry(3) = d1mach2
      zr = z
      if (z%re>=0.0d0) goto 10
      zr = -z
      10 continue
      yy = zr%im
      zn%re = zr%im
      zn%im = -zr%re
      zb = zr
      inu = int(fnu)
      fnf = fnu - real(inu,8)
      ang = -hpi*fnf
      ca=cmplx(cos(ang),sin(ang),8)
      c2=hpi*conjg(ca)
      kk = mod(inu, 4) + 1
      st = c2*cip(kk) 
      cs = cr1*st
      if (yy>0.0d0) goto 20
      zn = -zn
      20 continue
      !-----------------------------------------------------------------------
      !     k(fnu,z) is computed from h(2,fnu,-i*z) where z is in the first
      !     quadrant. fourth quadrant values (yy<=0.0e0) are computed by
      !     conjugation since the k function is real on the positive real axis
      !-----------------------------------------------------------------------
      j = 2
      do i = 1, n
         !-----------------------------------------------------------------------
         !     j flip flops between 1 and 2 in j = 3 - j
         !-----------------------------------------------------------------------
         j = 3 - j
         fn = fnu + real(i-1,8)
         call f90_zunhj(zn, fn, 0, tol, phi(j), arg(j), zeta1(j), zeta2(j), asum(j), bsum(j))
         if (kode==1) goto 30
         st = zb + zeta2(j)
         rast = fn/abs(st)
         st=conjg(st)*rast*rast
         s1 = zeta1(j) - st
         goto 40
         30 continue
         s1 = zeta1(j) - zeta2(j)
         40 continue
         !-----------------------------------------------------------------------
         !     test for underflow and overflow
         !-----------------------------------------------------------------------
         rs1 = s1%re
         if (abs(rs1)>elim) goto 70
         if (kdflg==1) kflag = 2
         if (abs(rs1)<alim) goto 50
         !-----------------------------------------------------------------------
         !     refine  test and scale
         !-----------------------------------------------------------------------
         aphi = abs(phi(j))
         aarg = abs(arg(j))
         rs1 = rs1 + log(aphi) - 0.25d0*log(aarg) - aic
         if (abs(rs1)>elim) goto 70
         if (kdflg==1) kflag = 1
         if (rs1<0.0d0) goto 50
         if (kdflg==1) kflag = 3
         50 continue
         !-----------------------------------------------------------------------
         !     scale s1 to keep intermediate arithmetic on scale near
         !     exponent extremes
         !-----------------------------------------------------------------------
         c2 = arg(j)*cr2
         call f90_zairy(c2, 0, 2, ai, nai, idum)
         call f90_zairy(c2, 1, 2, dai, ndai, idum)
         st = dai*bsum(j)
         pt = st*cr2
         st = pt + (ai*asum(j))
         pt = st*phi(j)
         s2 = pt*cs
         s1=exp(s1)*cssr(kflag)
         s2=s1*s2
         if (kflag/=1) goto 60
         nw=f90_zuchk(s2, bry(1), tol)
         if (nw/=0) goto 70
         60 continue
         if (yy<=0.0d0) s2%im = -s2%im
         cy(kdflg) = s2
         y(i) = s2*csrr(kflag)
         cs=cmplx(cs%im,-cs%re,8)
         if (kdflg==2) goto 85
         kdflg = 2
         goto 80
         70 continue
         if (rs1>0.0d0) goto 320
         !-----------------------------------------------------------------------
         !     for zr<0.0, the i function to be added will overflow
         !-----------------------------------------------------------------------
         if (z%re<0.0d0) goto 320
         kdflg = 1
         y(i) = 0.0d0
         nz = nz + 1
         cs=cmplx(cs%im,-cs%re,8)
         if (i==1) goto 80
         if (iszero(y(i-1))) goto 80
         y(i-1) = 0.0d0
         nz = nz + 1
      80 end do
      i = n
      85 continue
      razr = 1.0d0/abs(zr)
      st=conjg(zr)*razr
      rz =(st+st)*razr
      ck = fn*rz
      ib = i + 1
      if (n<ib) goto 180
      !-----------------------------------------------------------------------
      !     test last member for underflow and overflow. set sequence to zero
      !     on underflow.
      !-----------------------------------------------------------------------
      fn = fnu + real(n-1,8)
      ipard = 1
      if (mr/=0) ipard = 0
      call f90_zunhj(zn, fn, ipard, tol, phid, argd, zet1d, zet2d, asumd, bsumd)
      if (kode==1) goto 90
      st = zb + zet2d
      rast = fn/abs(st)
      st=conjg(st)*rast*rast
      s1 = zet1d - st
      goto 100
      90 continue
      s1 = zet1d - zet2d
      100 continue
      rs1 = s1%re
      if (abs(rs1)>elim) goto 105
      if (abs(rs1)<alim) goto 120
      !----------------------------------------------------------------------------
      !     refine estimate and test
      !-------------------------------------------------------------------------
      aphi = abs(phid)
      rs1 = rs1 + log(aphi)
      if (abs(rs1)<elim) goto 120
      105 continue
      if (rs1>0.0d0) goto 320
      !-----------------------------------------------------------------------
      !     for zr<0.0, the i function to be added will overflow
      !-----------------------------------------------------------------------
      if (z%re<0.0d0) goto 320
      nz = n
      do i = 1, n
         y(i) = 0.0d0
      end do
      return
      120 continue
      s1 = cy(1)
      s2 = cy(2)
      c1r = csrr(kflag)
      ascle = bry(kflag)
      do i = ib, n
         c2 = s2
         s2 = ck*c2  + s1
         s1 = c2
         ck = ck + rz
         c2 = s2*c1r
         y(i) = c2
         if (kflag>=3) goto 130
         c2m = max(abs(c2%re), abs(c2%im))
         if (c2m<=ascle) goto 130
         kflag = kflag + 1
         ascle = bry(kflag)
         s1 = s1*c1r
         s2 = c2
         s1 = s1*cssr(kflag)
         s2 = s2*cssr(kflag)
         c1r = csrr(kflag)
      130 end do
      180 continue
      if (mr==0) return
      !-----------------------------------------------------------------------
      !     analytic continuation for re(z)<0.0d0
      !-----------------------------------------------------------------------
      nz = 0
      fmr = real(mr,8)
      sgn = -sign(dpi, fmr)
      !-----------------------------------------------------------------------
      !     cspn and csgn are coeff of k and i funcions resp.
      !-----------------------------------------------------------------------
      csgni = sgn
      if (yy<=0.0d0) csgni = -csgni
      ifn = inu + n - 1
      ang = fnf*sgn
      cspn=cmplx(cos(ang),sin(ang),8)
      if (mod(ifn,2)==0) goto 190
      cspn = -cspn
      190 continue
      !-----------------------------------------------------------------------
      !     cs=coeff of the j function to get the i function. i(fnu,z) is
      !     computed from exp(i*fnu*hpi)*j(fnu,-i*z) where z is in the first
      !     quadrant. fourth quadrant values (yy<=0.0e0) are computed by
      !     conjugation since the i function is real on the positive real axis
      !-----------------------------------------------------------------------
      cs=csgni*cmplx(ca%im,ca%re,8)
      in = mod(ifn, 4) + 1
      c2 = cip(in)
      cs=cs*conjg(c2)
      asc = bry(1)
      iuf = 0
      kk = n
      kdflg = 1
      ib = ib - 1
      ic = ib - 1
      do k = 1, n
         fn = fnu + real(kk-1,8)
         !-----------------------------------------------------------------------
         !     logic to sort out cases whose parameters were set for the k
         !     function above
         !-----------------------------------------------------------------------
         if (n>2) goto 175
         172 continue
         phid  = phi(j)
         argd  = arg(j)
         zet1d = zeta1(j)
         zet2d = zeta2(j)
         asumd = asum(j)
         bsumd = bsum(j)
         j = 3 - j
         goto 210
         175 continue
         if ((kk==n) .and. (ib<n)) goto 210
         if ((kk==ib) .or. (kk==ic)) goto 172
         call f90_zunhj(zn, fn, 0, tol, phid, argd, zet1d, zet2d, asumd, bsumd)
         210 continue
         if (kode==1) goto 220
         st = zb + zet2d
         rast = fn/abs(st)
         st=conjg(st)*rast*rast
         s1 = -zet1d + st
         goto 230
         220 continue
         s1 = -zet1d + zet2d
         230 continue
         !-----------------------------------------------------------------------
         !     test for underflow and overflow
         !-----------------------------------------------------------------------
         rs1 = s1%re
         if (abs(rs1)>elim) goto 280
         if (kdflg==1) iflag = 2
         if (abs(rs1)<alim) goto 240
         !-----------------------------------------------------------------------
         !     refine  test and scale
         !-----------------------------------------------------------------------
         aphi = abs(phid)
         aarg = abs(argd)
         rs1 = rs1 + log(aphi) - 0.25d0*log(aarg) - aic
         if (abs(rs1)>elim) goto 280
         if (kdflg==1) iflag = 1
         if (rs1<0.0d0) goto 240
         if (kdflg==1) iflag = 3
         240 continue
         call f90_zairy(argd, 0, 2, ai, nai, idum)
         call f90_zairy(argd, 1, 2, dai, ndai, idum)
         st = dai*bsumd
         st = st + (ai*asumd)
         pt = st*phid
         s2 = pt*cs
         s1=exp(s1)*cssr(iflag)
         s2=s2*s1
         if (iflag/=1) goto 250
         nw=f90_zuchk(s2, bry(1), tol)
         if (nw==0) goto 250
         s2 = 0.0d0
         250 continue
         if (yy<=0.0d0) s2%im = -s2%im
         cy(kdflg) = s2
         c2 = s2
         s2 = s2*csrr(iflag)
         !-----------------------------------------------------------------------
         !     add i and k functions, k sequence in y(i), i=1,n
         !-----------------------------------------------------------------------
         s1 = y(kk)
         if (kode==1) goto 270
         call f90_zs1s2(zr, s1, s2, nw, asc, alim, iuf)
         nz = nz + nw
         270 continue
         y(kk) = s1*cspn+ s2
         kk = kk - 1
         cspn = -cspn
         cs=cmplx(cs%im,-cs%re,8)
         if (.not.(iszero(c2))) goto 255
         kdflg = 1
         goto 290
         255 continue
         if (kdflg==2) goto 295
         kdflg = 2
         goto 290
         280 continue
         if (rs1>0.0d0) goto 320
         s2 = 0.0d0
         goto 250
      290 end do
      k = n
      295 continue
      il = n - k
      if (il==0) return
      !-----------------------------------------------------------------------
      !     recur backward for remainder of i sequence and add in the
      !     k functions, scaling the i sequence during recurrence to keep
      !     intermediate arithmetic on scale near exponent extremes.
      !-----------------------------------------------------------------------
      s1 = cy(1)
      s2 = cy(2)
      csr = csrr(iflag)
      ascle = bry(iflag)
      fn = real(inu+il,8)
      do i = 1, il
         c2 = s2
         s2 = s1 + (fn+fnf)*(rz*c2)
         s1 = c2
         fn = fn - 1.0d0
         c2 = s2*csr
         ck = c2
         c1 = y(kk)
         if (kode==1) goto 300
         call f90_zs1s2(zr, c1, c2, nw, asc, alim, iuf)
         nz = nz + nw
         300 continue
         y(kk) = c1*cspn+ c2
         kk = kk - 1
         cspn = -cspn
         if (iflag>=3) goto 310
         c2m = max(abs(ck%re),abs(ck%im))
         if (c2m<=ascle) goto 310
         iflag = iflag + 1
         ascle = bry(iflag)
         s1 = s1*csr
         s2 = ck
         s1 = s1*cssr(iflag)
         s2 = s2*cssr(iflag)
         csr = csrr(iflag)
      310 end do
      return
      320 continue
      nz = -1
      return
   end subroutine f90_zunk2

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbunk(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
      !***begin prologue  zbunk
      !***refer to  zbesk,zbesh
      !
      !     zbunk computes the k bessel function for fnu>fnul.
      !     according to the uniform asymptotic expansion for k(fnu,z)
      !     in zunk1 and the expansion for h(2,fnu,z) in zunk2
      !
      !***routines called  zunk1,zunk2
      !***end prologue  zbunk
      !     complex y,z
      complex(8),intent(in)::z
      real(8),intent(in)::fnu,tol,elim,alim
      integer,intent(in)::n,kode,mr
      integer,intent(inout)::nz
      complex(8),intent(inout)::y(n)
      real(8)::ax,ay
      !DEBUG write (*, *) 'call f90_zbunk, not test'
      nz = 0
      ax = abs(z%re)*1.7321d0
      ay = abs(z%im)
      if (ay>ax) goto 10
      !-----------------------------------------------------------------------
      !     asymptotic expansion for k(fnu,z) for large fnu applied in
      !     -pi/3<=arg(z)<=pi/3
      !-----------------------------------------------------------------------
      call f90_zunk1(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
      goto 20
      10 continue
      !-----------------------------------------------------------------------
      !     asymptotic expansion for h(2,fnu,z*exp(m*hpi)) for large fnu
      !     applied in pi/3<abs(arg(z))<=pi/2 where m=+i or -i
      !     and hpi=pi/2
      !-----------------------------------------------------------------------
      call f90_zunk2(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
      20 continue
      return
   end subroutine f90_zbunk

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbesi(z, fnu, kode, n, cy, nz, ierr)
      !***begin prologue  zbesi
      !***date written   830501   (yymmdd)
      !***revision date  890801   (yymmdd)
      !***category no.  b5k
      !***keywords  i-bessel function,complex bessel function,
      !             modified bessel function of the first kind
      !***author  amos, donald e., sandia national laboratories
      !***purpose  to compute i-bessel functions of complex argument
      !***description
      !
      !                    ***a real(8) routine***
      !         on kode=1, zbesi computes an n member sequence of complex
      !         bessel functions cy(j)=i(fnu+j-1,z) for real, nonnegative
      !         orders fnu+j-1, j=1,...,n and complex z in the cut plane
      !         -dpi<arg(z)<=dpi. on kode=2, zbesi returns the scaled
      !         functions
      !
      !         cy(j)=exp(-abs(x))*i(fnu+j-1,z)   j = 1,...,n , x=real(z)
      !
      !         with the exponential growth removed in both the left and
      !         right half planes for z to infinity. definitions and notation
      !         are found in the nbs handbook of mathematical functions
      !         (ref. 1).
      !
      !         input      zr,zi,fnu are real(8)
      !           zr,zi  - z=cmplx(zr,zi),  -dpi<arg(z)<=dpi
      !           fnu    - order of initial i function, fnu>=0.0d0
      !           kode   - a parameter to indicate the scaling option
      !                    kode= 1  returns
      !                             cy(j)=i(fnu+j-1,z), j=1,...,n
      !                        = 2  returns
      !                             cy(j)=i(fnu+j-1,z)*exp(-abs(x)), j=1,...,n
      !           n      - number of members of the sequence, n>=1
      !
      !         output     cyr,cyi are real(8)
      !           cyr,cyi- real(8) vectors whose first n components
      !                    contain real and imaginary parts for the sequence
      !                    cy(j)=i(fnu+j-1,z)  or
      !                    cy(j)=i(fnu+j-1,z)*exp(-abs(x))  j=1,...,n
      !                    depending on kode, x=real(z)
      !           nz     - number of components set to zero due to underflow,
      !                    nz= 0   , normal return
      !                    nz>0 , last nz components of cy set to zero
      !                              to underflow, cy(j)=cmplx(0.0d0,0.0d0)
      !                              j = n-nz+1,...,n
      !           ierr   - error flag
      !                    ierr=0, normal return - computation completed
      !                    ierr=1, input error   - no computation
      !                    ierr=2, overflow      - no computation, real(z) too
      !                            large on kode=1
      !                    ierr=3, cabs(z) or fnu+n-1 large - computation done
      !                            but losses of signifcance by argument
      !                            reduction produce less than half of machine
      !                            accuracy
      !                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
      !                            tion because of complete losses of signifi-
      !                            cance by argument reduction
      !                    ierr=5, error              - no computation,
      !                            algorithm termination condition not met
      !
      !***long description
      !
      !         the computation is carried out by the power series for
      !         small cabs(z), the asymptotic expansion for large cabs(z),
      !         the miller algorithm normalized by the wronskian and a
      !         neumann series for imtermediate magnitudes, and the
      !         uniform asymptotic expansions for i(fnu,z) and j(fnu,z)
      !         for large orders. backward recurrence is used to generate
      !         sequences or reduce orders when necessary.
      !
      !         the calculations above are done in the right half plane and
      !         continued into the left half plane by the formula
      !
      !         i(fnu,z*exp(m*dpi)) = exp(m*dpi*fnu)*i(fnu,z)  real(z)>0.0
      !                       m = +i or -i,  i**2=-1
      !
      !         for negative orders,the formula
      !
      !              i(-fnu,z) = i(fnu,z) + (2/dpi)*sin(dpi*fnu)*k(fnu,z)
      !
      !         can be used. however,for large orders close to integers, the
      !         the function changes radically. when fnu is a large positive
      !         integer,the magnitude of i(-fnu,z)=i(fnu,z) is a large
      !         negative power of ten. but when fnu is not an integer,
      !         k(fnu,z) dominates in magnitude with a large positive power of
      !         ten and the most that the second term can be reduced is by
      !         unit roundoff from the coefficient. thus, wide changes can
      !         occur within unit roundoff of a large integer for fnu. here,
      !         large means fnu>cabs(z).
      !
      !         in most complex variable computation, one must evaluate ele-
      !         mentary functions. when the magnitude of z or fnu+n-1 is
      !         large, losses of significance by argument reduction occur.
      !         consequently, if either one exceeds u1=sqrt(0.5/ur), then
      !         losses exceeding half precision are likely and an error flag
      !         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
      !         real(8) unit roundoff limited to 18 digits precision.
      !         if either is larger than u2=0.5/ur, then all significance is
      !         lost and ierr=4. in order to use the int function, arguments
      !         must be further restricted not to exceed the largest machine
      !         integer, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
      !         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
      !         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
      !         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in real(8)
      !         arithmetic respectively. this makes u2 and u3 limiting in
      !         their respective arithmetics. this means that one can expect
      !         to retain, in the worst cases on 32 bit machines, no digits
      !         in single and only 7 digits in real(8) arithmetic.
      !         similar considerations hold for other machines.
      !
      !         the approximate relative error in the magnitude of a complex
      !         bessel function can be expressed by p*10**s where p=max(unit
      !         roundoff,1.0e-18) is the nominal precision and 10**s repre-
      !         sents the increase in error due to argument reduction in the
      !         elementary functions. here, s=max(1,abs(log10(cabs(z))),
      !         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
      !         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
      !         have only absolute accuracy. this is most likely to occur when
      !         one component (in absolute value) is larger than the other by
      !         several orders of magnitude. if one component is 10**k larger
      !         than the other, then one can expect only max(abs(log10(p))-k,
      !         0) significant digits; or, stated another way, when k exceeds
      !         the exponent of p, no significant digits remain in the smaller
      !         component. however, the phase angle retains absolute accuracy
      !         because, in complex arithmetic with precision p, the smaller
      !         component will not (as a rule) decrease below p times the
      !         magnitude of the larger component. in these extreme cases,
      !         the principal phase angle is on the order of +p, -p, dpi/2-p,
      !         or -dpi/2+p.
      !
      !***references  handbook of mathematical functions by m. abramowitz
      !                 and i. a. stegun, nbs ams series 55, u.s. dept. of
      !                 commerce, 1955.
      !
      !               computation of bessel functions of complex argument
      !                 by d. e. amos, sand83-0083, may, 1983.
      !
      !               computation of bessel functions of complex argument
      !                 and large order by d. e. amos, sand83-0643, may, 1983
      !
      !               a subroutine package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, sand85-
      !                 1018, may, 1985
      !
      !               a portable package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, trans.
      !                 math. software, 1986
      !
      !***routines called  zbinu,i1mach,d1mach
      !***end prologue  zbesi
      !     complex cone,csgn,cw,cy,czero,z,zn
      complex(8),intent(in)::z
      real(8),intent(in)::fnu
      integer,intent(in)::n,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::cy(n)
      real(8)::aa, alim, arg,dig,&
         elim, fnul, rl ,tol,&
         az, bb, fn, ascle, rtol, atol
      integer::i, ierr, inu, k, nn
      complex(8)::zn,csgn
      !
      !***first executable statement  zbesi
      !DEBUG write(*,*)"call f90_zbesi,not test"
      ierr = 0
      nz = 0
      if (fnu<0.0d0) ierr = 1
      if (kode<1 .or. kode>2) ierr = 1
      if (n<1) ierr = 1
      if (ierr/=0) return
      !-----------------------------------------------------------------------
      !     set parameters related to machine constants.
      !     tol is the approximate unit roundoff limited to 1.0e-18.
      !     elim is the approximate exponential over- and underflow limit.
      !     exp(-elim)<exp(-alim)=exp(-elim)/tol    and
      !     exp(elim)>exp(alim)=exp(elim)*tol       are intervals near
      !     underflow and overflow limits where scaled arithmetic is done.
      !     rl is the lower boundary of the asymptotic expansion for large z.
      !     dig = number of base 10 digits in tol = 10**(-dig).
      !     fnul is the lower boundary of the asymptotic series for large fnu.
      !-----------------------------------------------------------------------
      tol = max(d1mach4, 1.0d-18)
      k = min(abs(i1mach15), abs(i1mach16))
      elim = 2.303d0*(real(k,8)*d1mach5-3.0d0)
      aa = d1mach5*real(i1mach14 - 1,8)
      dig = min(aa, 18.0d0)
      aa = aa*2.303d0
      alim = elim + max(-aa, -41.45d0)
      rl = 1.2d0*dig + 3.0d0
      fnul = 10.0d0 + 6.0d0*(dig-3.0d0)
      !-----------------------------------------------------------------------------
      !     test for proper range
      !-----------------------------------------------------------------------
      az = abs(z)
      fn = fnu + real(n-1,8)
      aa = 0.5d0/tol
      bb = real(i1mach9,8)*0.5d0
      aa = min(aa, bb)
      if (az>aa) goto 260
      if (fn>aa) goto 260
      aa = sqrt(aa)
      if (az>aa) ierr = 3
      if (fn>aa) ierr = 3
      zn = z
      csgn = 1.0d0
      if (z%re>=0.0d0) goto 40
      zn = -z
      !-----------------------------------------------------------------------
      !     calculate csgn=exp(fnu*dpi*i) to minimize losses of significance
      !     when fnu is large
      !-----------------------------------------------------------------------
      inu = int(fnu)
      arg = (fnu-real(inu,8))*dpi
      if (z%im<0.0d0) arg = -arg
      csgn=cmplx(cos(arg),sin(arg),8)
      if (mod(inu,2)==0) goto 40
      csgn = -csgn
      40 continue
      call f90_zbinu(zn, fnu, kode, n, cy, nz, rl, fnul, tol, elim, alim)
      if (nz<0) goto 120
      if (z%re>=0.0d0) return
      !-----------------------------------------------------------------------
      !     analytic continuation to the left half plane
      !-----------------------------------------------------------------------
      nn = n - nz
      if (nn==0) return
      rtol = 1.0d0/tol
      ascle = d1mach1*rtol*1.0d+3
      do i = 1, nn
         !       str = cyr(i)*csgnr - cyi(i)*csgni
         !       cyi(i) = cyr(i)*csgni + cyi(i)*csgnr
         !       cyr(i) = str
         atol = 1.0d0
         if (max(abs(cy(i)%re),abs(cy(i)%im))>ascle) goto 55
         cy(i)=cy(i)*rtol
         atol = tol
         55 continue
         cy(i)=cy(i)*csgn*atol
         csgn = -csgn
      end do
      return
      120 continue
      if (nz==(-2)) goto 130
      nz = 0
      ierr = 2
      return
      130 continue
      nz = 0
      ierr = 5
      return
      260 continue
      nz = 0
      ierr = 4
      return
   end subroutine f90_zbesi

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbesj(z, fnu, kode, n, cy, nz, ierr)
      !***begin prologue  zbesj
      !***date written   830501   (yymmdd)
      !***revision date  890801   (yymmdd)
      !***category no.  b5k
      !***keywords  j-bessel function,bessel function of complex argument,
      !             bessel function of first kind
      !***author  amos, donald e., sandia national laboratories
      !***purpose  to compute the j-bessel function of a complex argument
      !***description
      !
      !                      ***a real(8) routine***
      !         on kode=1, cbesj computes an n member  sequence of complex
      !         bessel functions cy(i)=j(fnu+i-1,z) for real, nonnegative
      !         orders fnu+i-1, i=1,...,n and complex z in the cut plane
      !         -pi<arg(z)<=pi. on kode=2, cbesj returns the scaled
      !         functions
      !
      !         cy(i)=exp(-abs(y))*j(fnu+i-1,z)   i = 1,...,n , y=aimag(z)
      !
      !         which remove the exponential growth in both the upper and
      !         lower half planes for z to infinity. definitions and notation
      !         are found in the nbs handbook of mathematical functions
      !         (ref. 1).
      !
      !         input      zr,zi,fnu are real(8)
      !           zr,zi  - z=cmplx(zr,zi),  -pi<arg(z)<=pi
      !           fnu    - order of initial j function, fnu>=0.0d0
      !           kode   - a parameter to indicate the scaling option
      !                    kode= 1  returns
      !                             cy(i)=j(fnu+i-1,z), i=1,...,n
      !                        = 2  returns
      !                             cy(i)=j(fnu+i-1,z)exp(-abs(y)), i=1,...,n
      !           n      - number of members of the sequence, n>=1
      !
      !         output     cyr,cyi are real(8)
      !           cyr,cyi- real(8) vectors whose first n components
      !                    contain real and imaginary parts for the sequence
      !                    cy(i)=j(fnu+i-1,z)  or
      !                    cy(i)=j(fnu+i-1,z)exp(-abs(y))  i=1,...,n
      !                    depending on kode, y=aimag(z).
      !           nz     - number of components set to zero due to underflow,
      !                    nz= 0   , normal return
      !                    nz>0 , last nz components of cy set  zero due
      !                              to underflow, cy(i)=cmplx(0.0d0,0.0d0),
      !                              i = n-nz+1,...,n
      !           ierr   - error flag
      !                    ierr=0, normal return - computation completed
      !                    ierr=1, input error   - no computation
      !                    ierr=2, overflow      - no computation, aimag(z)
      !                            too large on kode=1
      !                    ierr=3, cabs(z) or fnu+n-1 large - computation done
      !                            but losses of signifcance by argument
      !                            reduction produce less than half of machine
      !                            accuracy
      !                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
      !                            tion because of complete losses of signifi-
      !                            cance by argument reduction
      !                    ierr=5, error              - no computation,
      !                            algorithm termination condition not met
      !
      !***long description
      !
      !         the computation is carried out by the formula
      !
      !         j(fnu,z)=exp( fnu*pi*i/2)*i(fnu,-i*z)    aimag(z)>=0.0
      !
      !         j(fnu,z)=exp(-fnu*pi*i/2)*i(fnu, i*z)    aimag(z)<0.0
      !
      !         where i**2 = -1 and i(fnu,z) is the i bessel function.
      !
      !         for negative orders,the formula
      !
      !              j(-fnu,z) = j(fnu,z)*cos(pi*fnu) - y(fnu,z)*sin(pi*fnu)
      !
      !         can be used. however,for large orders close to integers, the
      !         the function changes radically. when fnu is a large positive
      !         integer,the magnitude of j(-fnu,z)=j(fnu,z)*cos(pi*fnu) is a
      !         large negative power of ten. but when fnu is not an integer,
      !         y(fnu,z) dominates in magnitude with a large positive power of
      !         ten and the most that the second term can be reduced is by
      !         unit roundoff from the coefficient. thus, wide changes can
      !         occur within unit roundoff of a large integer for fnu. here,
      !         large means fnu>cabs(z).
      !
      !         in most complex variable computation, one must evaluate ele-
      !         mentary functions. when the magnitude of z or fnu+n-1 is
      !         large, losses of significance by argument reduction occur.
      !         consequently, if either one exceeds u1=sqrt(0.5/ur), then
      !         losses exceeding half precision are likely and an error flag
      !         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
      !         real(8) unit roundoff limited to 18 digits precision.
      !         if either is larger than u2=0.5/ur, then all significance is
      !         lost and ierr=4. in order to use the int function, arguments
      !         must be further restricted not to exceed the largest machine
      !         integer, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
      !         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
      !         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
      !         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in real(8)
      !         arithmetic respectively. this makes u2 and u3 limiting in
      !         their respective arithmetics. this means that one can expect
      !         to retain, in the worst cases on 32 bit machines, no digits
      !         in single and only 7 digits in real(8) arithmetic.
      !         similar considerations hold for other machines.
      !
      !         the approximate relative error in the magnitude of a complex
      !         bessel function can be expressed by p*10**s where p=max(unit
      !         roundoff,1.0e-18) is the nominal precision and 10**s repre-
      !         sents the increase in error due to argument reduction in the
      !         elementary functions. here, s=max(1,abs(log10(cabs(z))),
      !         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
      !         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
      !         have only absolute accuracy. this is most likely to occur when
      !         one component (in absolute value) is larger than the other by
      !         several orders of magnitude. if one component is 10**k larger
      !         than the other, then one can expect only max(abs(log10(p))-k,
      !         0) significant digits; or, stated another way, when k exceeds
      !         the exponent of p, no significant digits remain in the smaller
      !         component. however, the phase angle retains absolute accuracy
      !         because, in complex arithmetic with precision p, the smaller
      !         component will not (as a rule) decrease below p times the
      !         magnitude of the larger component. in these extreme cases,
      !         the principal phase angle is on the order of +p, -p, pi/2-p,
      !         or -pi/2+p.
      !
      !***references  handbook of mathematical functions by m. abramowitz
      !                 and i. a. stegun, nbs ams series 55, u.s. dept. of
      !                 commerce, 1955.
      !
      !               computation of bessel functions of complex argument
      !                 by d. e. amos, sand83-0083, may, 1983.
      !
      !               computation of bessel functions of complex argument
      !                 and large order by d. e. amos, sand83-0643, may, 1983
      !
      !               a subroutine package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, sand85-
      !                 1018, may, 1985
      !
      !               a portable package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, trans.
      !                 math. software, 1986
      !
      !***routines called  zbinu,i1mach,d1mach
      !***end prologue  zbesj
      !
      !     complex ci,csgn,cy,z,zn
      complex(8),intent(in)::z
      real(8),intent(in)::fnu
      integer,intent(in)::n,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::cy(n)
      real(8)::aa, alim, arg, cii, csgni, dig, elim, fnul, rl, tol, &
         bb, fn, az, ascle, rtol, atol
      integer::i, ierr, inu, inuh, ir, k, nl
      complex(8)::zn,csgn
      !DEBUG write(*,*)"call f90_zbesj,not test"
      !
      !***first executable statement  zbesj
      ierr = 0
      nz = 0
      if (fnu<0.0d0) ierr = 1
      if (kode<1 .or. kode>2) ierr = 1
      if (n<1) ierr = 1
      if (ierr/=0) return
      !-----------------------------------------------------------------------
      !     set parameters related to machine constants.
      !     tol is the approximate unit roundoff limited to 1.0e-18.
      !     elim is the approximate exponential over- and underflow limit.
      !     exp(-elim)<exp(-alim)=exp(-elim)/tol    and
      !     exp(elim)>exp(alim)=exp(elim)*tol       are intervals near
      !     underflow and overflow limits where scaled arithmetic is done.
      !     rl is the lower boundary of the asymptotic expansion for large z.
      !     dig = number of base 10 digits in tol = 10**(-dig).
      !     fnul is the lower boundary of the asymptotic series for large fnu.
      !-----------------------------------------------------------------------
      tol = max(d1mach4, 1.0d-18)
      k = min(abs(i1mach15), abs(i1mach16))
      elim = 2.303d0*(real(k,8)*d1mach5-3.0d0)
      aa = d1mach5*real(i1mach14 - 1,8)
      dig = min(aa, 18.0d0)
      aa = aa*2.303d0
      alim = elim + max(-aa, -41.45d0)
      rl = 1.2d0*dig + 3.0d0
      fnul = 10.0d0 + 6.0d0*(dig-3.0d0)
      !-----------------------------------------------------------------------
      !     test for proper range
      !-----------------------------------------------------------------------
      az = abs(z)
      fn = fnu + real(n-1,8)
      aa = 0.5d0/tol
      bb = real(i1mach9,8)*0.5d0
      aa = min(aa, bb)
      if (az>aa) goto 260
      if (fn>aa) goto 260
      aa = sqrt(aa)
      if (az>aa) ierr = 3
      if (fn>aa) ierr = 3
      !-----------------------------------------------------------------------
      !     calculate csgn=exp(fnu*hpi*i) to minimize losses of significance
      !     when fnu is large
      !-----------------------------------------------------------------------
      cii = 1.0d0
      inu = int(fnu)
      inuh = inu/2
      ir = inu - 2*inuh
      arg = (fnu-real(inu-ir,8))*hpi
      csgn=cmplx(cos(arg),sin(arg),8)
      if (mod(inuh,2)==0) goto 40
      csgn = -csgn
      40 continue
      !-----------------------------------------------------------------------
      !     zn is in the right half plane
      !-----------------------------------------------------------------------
      zn=cmplx(z%im,-z%re,8)
      if (z%im>=0.0d0) goto 50
      zn = -zn
      csgni = -csgni
      cii = -cii
      50 continue
      call f90_zbinu(zn, fnu, kode, n, cy, nz, rl, fnul, tol, elim, alim)
      if (nz<0) goto 130
      nl = n - nz
      if (nl==0) return
      rtol = 1.0d0/tol
      ascle = d1mach1*rtol*1.0d+3
      do i = 1, nl
         !       str = cyr(i)*csgnr - cyi(i)*csgni
         !       cyi(i) = cyr(i)*csgni + cyi(i)*csgnr
         !       cyr(i) = str
         atol = 1.0d0
         if (max(abs(cy(i)%re),abs(cy(i)%im))>ascle) goto 55
         cy(i)=cy(i)*rtol
         atol = tol
         55 continue
         cy(i)=cy(i)*csgn*atol
         csgn = cii*cmplx(-csgn%im,csgn%re,8)
      end do
      return
      130 continue
      if (nz==(-2)) goto 140
      nz = 0
      ierr = 2
      return
      140 continue
      nz = 0
      ierr = 5
      return
      260 continue
      nz = 0
      ierr = 4
      return
   end subroutine f90_zbesj

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbesk(z, fnu, kode, n, cy, nz, ierr)
      !***begin prologue  zbesk
      !***date written   830501   (yymmdd)
      !***revision date  890801   (yymmdd)
      !***category no.  b5k
      !***keywords  k-bessel function,complex bessel function,
      !             modified bessel function of the second kind,
      !             bessel function of the third kind
      !***author  amos, donald e., sandia national laboratories
      !***purpose  to compute k-bessel functions of complex argument
      !***description
      !
      !                      ***a real(8) routine***
      !
      !         on kode=1, cbesk computes an n member sequence of complex
      !         bessel functions cy(j)=k(fnu+j-1,z) for real, nonnegative
      !         orders fnu+j-1, j=1,...,n and complex z/=cmplx(0.0,0.0)
      !         in the cut plane -pi<arg(z)<=pi. on kode=2, cbesk
      !         returns the scaled k functions,
      !
      !         cy(j)=exp(z)*k(fnu+j-1,z) , j=1,...,n,
      !
      !         which remove the exponential behavior in both the left and
      !         right half planes for z to infinity. definitions and
      !         notation are found in the nbs handbook of mathematical
      !         functions (ref. 1).
      !
      !         input      zr,zi,fnu are real(8)
      !           zr,zi  - z=cmplx(zr,zi), z/=cmplx(0.0d0,0.0d0),
      !                    -pi<arg(z)<=pi
      !           fnu    - order of initial k function, fnu>=0.0d0
      !           n      - number of members of the sequence, n>=1
      !           kode   - a parameter to indicate the scaling option
      !                    kode= 1  returns
      !                             cy(i)=k(fnu+i-1,z), i=1,...,n
      !                        = 2  returns
      !                             cy(i)=k(fnu+i-1,z)*exp(z), i=1,...,n
      !
      !         output     cyr,cyi are real(8)
      !           cyr,cyi- real(8) vectors whose first n components
      !                    contain real and imaginary parts for the sequence
      !                    cy(i)=k(fnu+i-1,z), i=1,...,n or
      !                    cy(i)=k(fnu+i-1,z)*exp(z), i=1,...,n
      !                    depending on kode
      !           nz     - number of components set to zero due to underflow.
      !                    nz= 0   , normal return
      !                    nz>0 , first nz components of cy set to zero due
      !                              to underflow, cy(i)=cmplx(0.0d0,0.0d0),
      !                              i=1,...,n when x>=0.0. when x<0.0
      !                              nz states only the number of underflows
      !                              in the sequence.
      !
      !           ierr   - error flag
      !                    ierr=0, normal return - computation completed
      !                    ierr=1, input error   - no computation
      !                    ierr=2, overflow      - no computation, fnu is
      !                            too large or cabs(z) is too small or both
      !                    ierr=3, cabs(z) or fnu+n-1 large - computation done
      !                            but losses of signifcance by argument
      !                            reduction produce less than half of machine
      !                            accuracy
      !                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
      !                            tion because of complete losses of signifi-
      !                            cance by argument reduction
      !                    ierr=5, error              - no computation,
      !                            algorithm termination condition not met
      !
      !***long description
      !
      !         equations of the reference are implemented for small orders
      !         dnu and dnu+1.0 in the right half plane x>=0.0. forward
      !         recurrence generates higher orders. k is continued to the left
      !         half plane by the relation
      !
      !         k(fnu,z*exp(mp)) = exp(-mp*fnu)*k(fnu,z)-mp*i(fnu,z)
      !         mp=mr*pi*i, mr=+1 or -1, re(z)>0, i**2=-1
      !
      !         where i(fnu,z) is the i bessel function.
      !
      !         for large orders, fnu>fnul, the k function is computed
      !         by means of its uniform asymptotic expansions.
      !
      !         for negative orders, the formula
      !
      !                       k(-fnu,z) = k(fnu,z)
      !
      !         can be used.
      !
      !         cbesk assumes that a significant digit sinh(x) function is
      !         available.
      !
      !         in most complex variable computation, one must evaluate ele-
      !         mentary functions. when the magnitude of z or fnu+n-1 is
      !         large, losses of significance by argument reduction occur.
      !         consequently, if either one exceeds u1=sqrt(0.5/ur), then
      !         losses exceeding half precision are likely and an error flag
      !         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
      !         real(8) unit roundoff limited to 18 digits precision.
      !         if either is larger than u2=0.5/ur, then all significance is
      !         lost and ierr=4. in order to use the int function, arguments
      !         must be further restricted not to exceed the largest machine
      !         integer, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
      !         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
      !         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
      !         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in real(8)
      !         arithmetic respectively. this makes u2 and u3 limiting in
      !         their respective arithmetics. this means that one can expect
      !         to retain, in the worst cases on 32 bit machines, no digits
      !         in single and only 7 digits in real(8) arithmetic.
      !         similar considerations hold for other machines.
      !
      !         the approximate relative error in the magnitude of a complex
      !         bessel function can be expressed by p*10**s where p=max(unit
      !         roundoff,1.0e-18) is the nominal precision and 10**s repre-
      !         sents the increase in error due to argument reduction in the
      !         elementary functions. here, s=max(1,abs(log10(cabs(z))),
      !         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
      !         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
      !         have only absolute accuracy. this is most likely to occur when
      !         one component (in absolute value) is larger than the other by
      !         several orders of magnitude. if one component is 10**k larger
      !         than the other, then one can expect only max(abs(log10(p))-k,
      !         0) significant digits; or, stated another way, when k exceeds
      !         the exponent of p, no significant digits remain in the smaller
      !         component. however, the phase angle retains absolute accuracy
      !         because, in complex arithmetic with precision p, the smaller
      !         component will not (as a rule) decrease below p times the
      !         magnitude of the larger component. in these extreme cases,
      !         the principal phase angle is on the order of +p, -p, pi/2-p,
      !         or -pi/2+p.
      !
      !***references  handbook of mathematical functions by m. abramowitz
      !                 and i. a. stegun, nbs ams series 55, u.s. dept. of
      !                 commerce, 1955.
      !
      !               computation of bessel functions of complex argument
      !                 by d. e. amos, sand83-0083, may, 1983.
      !
      !               computation of bessel functions of complex argument
      !                 and large order by d. e. amos, sand83-0643, may, 1983.
      !
      !               a subroutine package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, sand85-
      !                 1018, may, 1985
      !
      !               a portable package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, trans.
      !                 math. software, 1986
      !
      !***routines called  zacon,zbknu,zbunk,zuoik,azabs,i1mach,d1mach
      !***end prologue  zbesk
      !
      !     complex cy,z
      complex(8),intent(in)::z
      real(8),intent(in)::fnu
      integer,intent(in)::n,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::cy(n)
      real(8)::aa, alim, aln, arg, az, dig, elim, fn, fnul, rl, tol, ufl, bb
      integer::ierr, k, mr, nn, nuf, nw
      !***first executable statement  zbesk
      !DEBUG write(*,*)"call f90_zbesk,not test"
      ierr = 0
      nz = 0
      if (iszero(z)) ierr = 1
      if (fnu<0.0d0) ierr = 1
      if (kode<1 .or. kode>2) ierr = 1
      if (n<1) ierr = 1
      if (ierr/=0) return
      nn = n
      !-----------------------------------------------------------------------
      !     set parameters related to machine constants.
      !     tol is the approximate unit roundoff limited to 1.0e-18.
      !     elim is the approximate exponential over- and underflow limit.
      !     exp(-elim)<exp(-alim)=exp(-elim)/tol    and
      !     exp(elim)>exp(alim)=exp(elim)*tol       are intervals near
      !     underflow and overflow limits where scaled arithmetic is done.
      !     rl is the lower boundary of the asymptotic expansion for large z.
      !     dig = number of base 10 digits in tol = 10**(-dig).
      !     fnul is the lower boundary of the asymptotic series for large fnu
      !-----------------------------------------------------------------------
      tol = max(d1mach4, 1.0d-18)
      k = min(abs(i1mach15), abs(i1mach16))
      elim = 2.303d0*(real(k,8)*d1mach5-3.0d0)
      aa = d1mach5*real(i1mach14 - 1,8)
      dig = min(aa, 18.0d0)
      aa = aa*2.303d0
      alim = elim + max(-aa, -41.45d0)
      rl = 1.2d0*dig + 3.0d0
      fnul = 10.0d0 + 6.0d0*(dig-3.0d0)
      !-----------------------------------------------------------------------------
      !     test for proper range
      !-----------------------------------------------------------------------
      az = abs(z)
      fn = fnu + real(n-1,8)
      aa = 0.5d0/tol
      bb = real(i1mach9,8)*0.5d0
      aa = min(aa, bb)
      if (az>aa) goto 260
      if (fn>aa) goto 260
      aa = sqrt(aa)
      if (az>aa) ierr = 3
      if (fn>aa) ierr = 3
      !-----------------------------------------------------------------------
      !     overflow test on the last member of the sequence
      !-----------------------------------------------------------------------
      !     ufl = dexp(-elim)
      ufl = d1mach1*1.0d+3
      if (az<ufl) goto 180
      if (fnu>fnul) goto 80
      if (fn<=1.0d0) goto 60
      if (fn>2.0d0) goto 50
      if (az>tol) goto 60
      arg = 0.5d0*az
      aln = -fn*log(arg)
      if (aln>elim) goto 180
      goto 60
      50 continue
      call f90_zuoik(z, fnu, kode, 2, nn, cy, nuf, tol, elim, alim)
      if (nuf<0) goto 180
      nz = nz + nuf
      nn = nn - nuf
      !-----------------------------------------------------------------------
      !     here nn=n or nn=0 since nuf=0,nn, or -1 on return from cuoik
      !     if nuf=nn, then cy(i)=czero for all i
      !-----------------------------------------------------------------------
      if (nn==0) goto 100
      60 continue
      if (z%re<0.0d0) goto 70
      !-----------------------------------------------------------------------
      !     right half plane computation, real(z)>=0.
      !-----------------------------------------------------------------------
      call f90_zbknu(z, fnu, kode, nn, cy, nw, tol, elim, alim)
      if (nw<0) goto 200
      nz = nw
      return
      !-----------------------------------------------------------------------
      !     left half plane computation
      !     pi/2<arg(z)<=pi and -pi<arg(z)<-pi/2.
      !-----------------------------------------------------------------------
      70 continue
      if (nz/=0) goto 180
      mr = 1
      if (z%im<0.0d0) mr = -1
      call f90_zacon(z, fnu, kode, mr, nn, cy, nw, rl, fnul, tol, elim, alim)
      if (nw<0) goto 200
      nz = nw
      return
      !-----------------------------------------------------------------------
      !     uniform asymptotic expansions for fnu>fnul
      !-----------------------------------------------------------------------
      80 continue
      mr = 0
      if (z%re>=0.0d0) goto 90
      mr = 1
      if (z%im<0.0d0) mr = -1
      90 continue
      call f90_zbunk(z, fnu, kode, mr, nn, cy, nw, tol, elim, alim)
      if (nw<0) goto 200
      nz = nz + nw
      return
      100 continue
      if (z%re<0.0d0) goto 180
      return
      180 continue
      nz = 0
      ierr = 2
      return
      200 continue
      if (nw==(-1)) goto 180
      nz = 0
      ierr = 5
      return
      260 continue
      nz = 0
      ierr = 4
      return
   end subroutine f90_zbesk

   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbesh(z, fnu, kode, m, n, cy, nz, ierr)
      !***begin prologue  zbesh
      !***date written   830501   (yymmdd)
      !***revision date  890801   (yymmdd)
      !***category no.  b5k
      !***keywords  h-bessel functions,bessel functions of complex argument,
      !             bessel functions of third kind,hankel functions
      !***author  amos, donald e., sandia national laboratories
      !***purpose  to compute the h-bessel functions of a complex argument
      !***description
      !
      !                      ***a real(8) routine***
      !         on kode=1, zbesh computes an n member sequence of complex
      !         hankel (bessel) functions cy(j)=h(m,fnu+j-1,z) for kinds m=1
      !         or 2, real, nonnegative orders fnu+j-1, j=1,...,n, and complex
      !         z/=cmplx(0.0,0.0) in the cut plane -pi<arg(z)<=pi.
      !         on kode=2, zbesh returns the scaled hankel functions
      !
      !         cy(i)=exp(-mm*z*i)*h(m,fnu+j-1,z)       mm=3-2*m,   i**2=-1.
      !
      !         which removes the exponential behavior in both the upper and
      !         lower half planes. definitions and notation are found in the
      !         nbs handbook of mathematical functions (ref. 1).
      !
      !         input      zr,zi,fnu are real(8)
      !           zr,zi  - z=cmplx(zr,zi), z/=cmplx(0.0d0,0.0d0),
      !                    -pt<arg(z)<=pi
      !           fnu    - order of initial h function, fnu>=0.0d0
      !           kode   - a parameter to indicate the scaling option
      !                    kode= 1  returns
      !                             cy(j)=h(m,fnu+j-1,z),   j=1,...,n
      !                        = 2  returns
      !                             cy(j)=h(m,fnu+j-1,z)*exp(-i*z*(3-2m))
      !                                  j=1,...,n  ,  i**2=-1
      !           m      - kind of hankel function, m=1 or 2
      !           n      - number of members in the sequence, n>=1
      !
      !         output     cyr,cyi are real(8)
      !           cyr,cyi- real(8) vectors whose first n components
      !                    contain real and imaginary parts for the sequence
      !                    cy(j)=h(m,fnu+j-1,z)  or
      !                    cy(j)=h(m,fnu+j-1,z)*exp(-i*z*(3-2m))  j=1,...,n
      !                    depending on kode, i**2=-1.
      !           nz     - number of components set to zero due to underflow,
      !                    nz= 0   , normal return
      !                    nz>0 , first nz components of cy set to zero due
      !                              to underflow, cy(j)=cmplx(0.0d0,0.0d0)
      !                              j=1,...,nz when y>0.0 and m=1 or
      !                              y<0.0 and m=2. for the complmentary
      !                              half planes, nz states only the number
      !                              of underflows.
      !           ierr   - error flag
      !                    ierr=0, normal return - computation completed
      !                    ierr=1, input error   - no computation
      !                    ierr=2, overflow      - no computation, fnu too
      !                            large or cabs(z) too small or both
      !                    ierr=3, cabs(z) or fnu+n-1 large - computation done
      !                            but losses of signifcance by argument
      !                            reduction produce less than half of machine
      !                            accuracy
      !                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
      !                            tion because of complete losses of signifi-
      !                            cance by argument reduction
      !                    ierr=5, error              - no computation,
      !                            algorithm termination condition not met
      !
      !***long description
      !
      !         the computation is carried out by the relation
      !
      !         h(m,fnu,z)=(1/mp)*exp(-mp*fnu)*k(fnu,z*exp(-mp))
      !             mp=mm*hpi*i,  mm=3-2*m,  hpi=pi/2,  i**2=-1
      !
      !         for m=1 or 2 where the k bessel function is computed for the
      !         right half plane re(z)>=0.0. the k function is continued
      !         to the left half plane by the relation
      !
      !         k(fnu,z*exp(mp)) = exp(-mp*fnu)*k(fnu,z)-mp*i(fnu,z)
      !         mp=mr*pi*i, mr=+1 or -1, re(z)>0, i**2=-1
      !
      !         where i(fnu,z) is the i bessel function.
      !
      !         exponential decay of h(m,fnu,z) occurs in the upper half z
      !         plane for m=1 and the lower half z plane for m=2.  exponential
      !         growth occurs in the complementary half planes.  scaling
      !         by exp(-mm*z*i) removes the exponential behavior in the
      !         whole z plane for z to infinity.
      !
      !         for negative orders,the formulae
      !
      !               h(1,-fnu,z) = h(1,fnu,z)*cexp( pi*fnu*i)
      !               h(2,-fnu,z) = h(2,fnu,z)*cexp(-pi*fnu*i)
      !                         i**2=-1
      !
      !         can be used.
      !
      !         in most complex variable computation, one must evaluate ele-
      !         mentary functions. when the magnitude of z or fnu+n-1 is
      !         large, losses of significance by argument reduction occur.
      !         consequently, if either one exceeds u1=sqrt(0.5/ur), then
      !         losses exceeding half precision are likely and an error flag
      !         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
      !         real(8) unit roundoff limited to 18 digits precision.
      !         if either is larger than u2=0.5/ur, then all significance is
      !         lost and ierr=4. in order to use the int function, arguments
      !         must be further restricted not to exceed the largest machine
      !         integer, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
      !         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
      !         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
      !         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in real(8)
      !         arithmetic respectively. this makes u2 and u3 limiting in
      !         their respective arithmetics. this means that one can expect
      !         to retain, in the worst cases on 32 bit machines, no digits
      !         in single and only 7 digits in real(8) arithmetic.
      !         similar considerations hold for other machines.
      !
      !         the approximate relative error in the magnitude of a complex
      !         bessel function can be expressed by p*10**s where p=max(unit
      !         roundoff,1.0d-18) is the nominal precision and 10**s repre-
      !         sents the increase in error due to argument reduction in the
      !         elementary functions. here, s=max(1,abs(log10(cabs(z))),
      !         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
      !         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
      !         have only absolute accuracy. this is most likely to occur when
      !         one component (in absolute value) is larger than the other by
      !         several orders of magnitude. if one component is 10**k larger
      !         than the other, then one can expect only max(abs(log10(p))-k,
      !         0) significant digits; or, stated another way, when k exceeds
      !         the exponent of p, no significant digits remain in the smaller
      !         component. however, the phase angle retains absolute accuracy
      !         because, in complex arithmetic with precision p, the smaller
      !         component will not (as a rule) decrease below p times the
      !         magnitude of the larger component. in these extreme cases,
      !         the principal phase angle is on the order of +p, -p, pi/2-p,
      !         or -pi/2+p.
      !
      !***references  handbook of mathematical functions by m. abramowitz
      !                 and i. a. stegun, nbs ams series 55, u.s. dept. of
      !                 commerce, 1955.
      !
      !               computation of bessel functions of complex argument
      !                 by d. e. amos, sand83-0083, may, 1983.
      !
      !               computation of bessel functions of complex argument
      !                 and large order by d. e. amos, sand83-0643, may, 1983
      !
      !               a subroutine package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, sand85-
      !                 1018, may, 1985
      !
      !               a portable package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, trans.
      !                 math. software, 1986
      !
      !***routines called  zacon,zbknu,zbunk,zuoik,azabs,i1mach,d1mach
      !***end prologue  zbesh
      !
      !     complex cy,z,zn,zt,csgn
      complex(8),intent(in)::z
      real(8),intent(in)::fnu
      integer,intent(in)::n,m,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::cy(n)
      real(8)::aa, alim, aln, arg, az, dig, elim, fmm, &
         fn, fnul, rhpi, rl, sgn, tol, ufl,&
         zti, bb, ascle, rtol, atol
      integer::i, ierr, inu, inuh, ir, k, mm, mr, nn, nuf, nw 
      complex(8)::zn,csgn
      !
      !
      !***first executable statement  zbesh
      !DEBUG write(*,*)"call f90_zbesh,not test"
      ierr = 0
      nz = 0
      if (iszero(z)) ierr = 1
      if (fnu<0.0d0) ierr = 1
      if (m<1 .or. m>2) ierr = 1
      if (kode<1 .or. kode>2) ierr = 1
      if (n<1) ierr = 1
      if (ierr/=0) return
      nn = n
      !-----------------------------------------------------------------------
      !     set parameters related to machine constants.
      !     tol is the approximate unit roundoff limited to 1.0e-18.
      !     elim is the approximate exponential over- and underflow limit.
      !     exp(-elim)<exp(-alim)=exp(-elim)/tol    and
      !     exp(elim)>exp(alim)=exp(elim)*tol       are intervals near
      !     underflow and overflow limits where scaled arithmetic is done.
      !     rl is the lower boundary of the asymptotic expansion for large z.
      !     dig = number of base 10 digits in tol = 10**(-dig).
      !     fnul is the lower boundary of the asymptotic series for large fnu
      !-----------------------------------------------------------------------
      tol = max(d1mach4, 1.0d-18)
      k = min(abs(i1mach15), abs(i1mach16))
      elim = 2.303d0*(real(k,8)*d1mach5-3.0d0)
      aa = d1mach5*real(i1mach14 - 1,8)
      dig = min(aa, 18.0d0)
      aa = aa*2.303d0
      alim = elim + max(-aa, -41.45d0)
      rl = 1.2d0*dig + 3.0d0
      fnul = 10.0d0 + 6.0d0*(dig-3.0d0)

      fn = fnu + real(nn-1,8)
      mm = 3 - m - m
      fmm = real(mm,8)
      zn=fmm*cmplx(z%im,-z%re,8)
      !-----------------------------------------------------------------------
      !     test for proper range
      !-----------------------------------------------------------------------
      az = abs(z)
      fn = fnu + real(n-1,8)
      aa = 0.5d0/tol
      bb = real(i1mach9,8)*0.5d0
      aa = min(aa, bb)
      if (az>aa) goto 260
      if (fn>aa) goto 260
      aa = sqrt(aa)
      if (az>aa) ierr = 3
      if (fn>aa) ierr = 3
      !-----------------------------------------------------------------------
      !     overflow test on the last member of the sequence
      !-----------------------------------------------------------------------
      ufl = d1mach1*1.0d+3
      if (az<ufl) goto 230
      if (fnu>fnul) goto 90
      if (fn<=1.0d0) goto 70
      if (fn>2.0d0) goto 60
      if (az>tol) goto 70
      arg = 0.5d0*az
      aln = -fn*log(arg)
      if (aln>elim) goto 230
      goto 70
      60 continue
      call f90_zuoik(zn, fnu, kode, 2, nn, cy, nuf, tol, elim, alim)
      if (nuf<0) goto 230
      nz = nz + nuf
      nn = nn - nuf
      !-----------------------------------------------------------------------
      !     here nn=n or nn=0 since nuf=0,nn, or -1 on return from cuoik
      !     if nuf=nn, then cy(i)=czero for all i
      !-----------------------------------------------------------------------
      if (nn==0) goto 140
      70 continue
      if ((zn%re<0.0d0) .or. (iszero(zn%re) .and. zn%im<0.0d0 .and. m==2)) goto 80
      !-----------------------------------------------------------------------
      !     right half plane computation, xn>=0. .and. (xn/=0. .or.
      !     yn>=0. .or. m=1)
      !-----------------------------------------------------------------------
      call f90_zbknu(zn, fnu, kode, nn, cy, nz, tol, elim, alim)
      goto 110
      !-----------------------------------------------------------------------
      !     left half plane computation
      !-----------------------------------------------------------------------
      80 continue
      mr = -mm
      call f90_zacon(zn, fnu, kode, mr, nn, cy, nw, rl, fnul, tol, elim, alim)
      if (nw<0) goto 240
      nz = nw
      goto 110
      90 continue
      !-----------------------------------------------------------------------
      !     uniform asymptotic expansions for fnu>fnul
      !-----------------------------------------------------------------------
      mr = 0
      if ((zn%re>=0.0d0) .and. .not.iszero(zn%re).or. zn%im>=0.0d0 .or. m/=2) goto 100
      mr = -mm
      if (.not.iszero(zn%re) .or. zn%im>=0.0d0) goto 100
      zn = -zn
      100 continue
      call f90_zbunk(zn, fnu, kode, mr, nn, cy, nw, tol, elim, alim)
      if (nw<0) goto 240
      nz = nz + nw
      110 continue
      !-----------------------------------------------------------------------
      !     h(m,fnu,z) = -fmm*(i/hpi)*(zt**fnu)*k(fnu,-z*zt)
      !
      !     zt=exp(-fmm*hpi*i) = cmplx(0.0,-fmm), fmm=3-2*m, m=1,2
      !-----------------------------------------------------------------------
      sgn = sign(hpi, -fmm)
      !-----------------------------------------------------------------------
      !     calculate exp(fnu*hpi*i) to minimize losses of significance
      !     when fnu is large
      !-----------------------------------------------------------------------
      inu = int(fnu)
      inuh = inu/2
      ir = inu - 2*inuh
      arg = (fnu-real(inu-ir,8))*sgn
      rhpi = 1.0d0/sgn
      !     zni = rhpi*dcos(arg)
      !     znr = -rhpi*dsin(arg)
      csgn=rhpi*cmplx(-sin(arg),cos(arg),8)
      if (mod(inuh,2)==0) goto 120
      !     znr = -znr
      !     zni = -zni
      csgn = -csgn
      120 continue
      zti = -fmm
      rtol = 1.0d0/tol
      ascle = ufl*rtol
      do i = 1, nn
         !       str = cyr(i)*znr - cyi(i)*zni
         !       cyi(i) = cyr(i)*zni + cyi(i)*znr
         !       cyr(i) = str
         !       str = -zni*zti
         !       zni = znr*zti
         !       znr = str
         atol = 1.0d0
         if (max(abs(cy(i)%re),abs(cy(i)%im))>ascle) goto 135
         cy(i)=cy(i)*rtol
         atol = tol
         135 continue
         cy(i)=cy(i)*csgn*atol
         csgn = zti*cmplx(-csgn%im,csgn%re,8)
      end do
      return
      140 continue
      if (zn%re<0.0d0) goto 230
      return
      230 continue
      nz = 0
      ierr = 2
      return
      240 continue
      if (nw==(-1)) goto 230
      nz = 0
      ierr = 5
      return
      260 continue
      nz = 0
      ierr = 4
      return
   end subroutine f90_zbesh


   !// forquill v1.01 beta www.fcode.cn
   subroutine f90_zbesy(z, fnu, kode, n, cy, nz, cwrk, ierr)
      !***begin prologue  zbesy
      !***date written   830501   (yymmdd)
      !***revision date  890801   (yymmdd)
      !***category no.  b5k
      !***keywords  y-bessel function,bessel function of complex argument,
      !             bessel function of second kind
      !***author  amos, donald e., sandia national laboratories
      !***purpose  to compute the y-bessel function of a complex argument
      !***description
      !
      !                      ***a real(8) routine***
      !
      !         on kode=1, cbesy computes an n member sequence of complex
      !         bessel functions cy(i)=y(fnu+i-1,z) for real, nonnegative
      !         orders fnu+i-1, i=1,...,n and complex z in the cut plane
      !         -pi<arg(z)<=pi. on kode=2, cbesy returns the scaled
      !         functions
      !
      !         cy(i)=exp(-abs(y))*y(fnu+i-1,z)   i = 1,...,n , y=aimag(z)
      !
      !         which remove the exponential growth in both the upper and
      !         lower half planes for z to infinity. definitions and notation
      !         are found in the nbs handbook of mathematical functions
      !         (ref. 1).
      !
      !         input      zr,zi,fnu are real(8)
      !           zr,zi  - z=cmplx(zr,zi), z/=cmplx(0.0d0,0.0d0),
      !                    -pi<arg(z)<=pi
      !           fnu    - order of initial y function, fnu>=0.0d0
      !           kode   - a parameter to indicate the scaling option
      !                    kode= 1  returns
      !                             cy(i)=y(fnu+i-1,z), i=1,...,n
      !                        = 2  returns
      !                             cy(i)=y(fnu+i-1,z)*exp(-abs(y)), i=1,...,n
      !                             where y=aimag(z)
      !           n      - number of members of the sequence, n>=1
      !           cwrkr, - real(8) work vectors of dimension at
      !           cwrki    at least n
      !
      !         output     cyr,cyi are real(8)
      !           cyr,cyi- real(8) vectors whose first n components
      !                    contain real and imaginary parts for the sequence
      !                    cy(i)=y(fnu+i-1,z)  or
      !                    cy(i)=y(fnu+i-1,z)*exp(-abs(y))  i=1,...,n
      !                    depending on kode.
      !           nz     - nz=0 , a normal return
      !                    nz>0 , nz components of cy set to zero due to
      !                    underflow (generally on kode=2)
      !           ierr   - error flag
      !                    ierr=0, normal return - computation completed
      !                    ierr=1, input error   - no computation
      !                    ierr=2, overflow      - no computation, fnu is
      !                            too large or cabs(z) is too small or both
      !                    ierr=3, cabs(z) or fnu+n-1 large - computation done
      !                            but losses of signifcance by argument
      !                            reduction produce less than half of machine
      !                            accuracy
      !                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
      !                            tion because of complete losses of signifi-
      !                            cance by argument reduction
      !                    ierr=5, error              - no computation,
      !                            algorithm termination condition not met
      !
      !***long description
      !
      !         the computation is carried out by the formula
      !
      !         y(fnu,z)=0.5*(h(1,fnu,z)-h(2,fnu,z))/i
      !
      !         where i**2 = -1 and the hankel bessel functions h(1,fnu,z)
      !         and h(2,fnu,z) are calculated in cbesh.
      !
      !         for negative orders,the formula
      !
      !              y(-fnu,z) = y(fnu,z)*cos(pi*fnu) + j(fnu,z)*sin(pi*fnu)
      !
      !         can be used. however,for large orders close to half odd
      !         integers the function changes radically. when fnu is a large
      !         positive half odd integer,the magnitude of y(-fnu,z)=j(fnu,z)*
      !         sin(pi*fnu) is a large negative power of ten. but when fnu is
      !         not a half odd integer, y(fnu,z) dominates in magnitude with a
      !         large positive power of ten and the most that the second term
      !         can be reduced is by unit roundoff from the coefficient. thus,
      !         wide changes can occur within unit roundoff of a large half
      !         odd integer. here, large means fnu>cabs(z).
      !
      !         in most complex variable computation, one must evaluate ele-
      !         mentary functions. when the magnitude of z or fnu+n-1 is
      !         large, losses of significance by argument reduction occur.
      !         consequently, if either one exceeds u1=sqrt(0.5/ur), then
      !         losses exceeding half precision are likely and an error flag
      !         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
      !         real(8) unit roundoff limited to 18 digits precision.
      !         if either is larger than u2=0.5/ur, then all significance is
      !         lost and ierr=4. in order to use the int function, arguments
      !         must be further restricted not to exceed the largest machine
      !         integer, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
      !         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
      !         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
      !         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in real(8)
      !         arithmetic respectively. this makes u2 and u3 limiting in
      !         their respective arithmetics. this means that one can expect
      !         to retain, in the worst cases on 32 bit machines, no digits
      !         in single and only 7 digits in real(8) arithmetic.
      !         similar considerations hold for other machines.
      !
      !         the approximate relative error in the magnitude of a complex
      !         bessel function can be expressed by p*10**s where p=max(unit
      !         roundoff,1.0e-18) is the nominal precision and 10**s repre-
      !         sents the increase in error due to argument reduction in the
      !         elementary functions. here, s=max(1,abs(log10(cabs(z))),
      !         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
      !         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
      !         have only absolute accuracy. this is most likely to occur when
      !         one component (in absolute value) is larger than the other by
      !         several orders of magnitude. if one component is 10**k larger
      !         than the other, then one can expect only max(abs(log10(p))-k,
      !         0) significant digits; or, stated another way, when k exceeds
      !         the exponent of p, no significant digits remain in the smaller
      !         component. however, the phase angle retains absolute accuracy
      !         because, in complex arithmetic with precision p, the smaller
      !         component will not (as a rule) decrease below p times the
      !         magnitude of the larger component. in these extreme cases,
      !         the principal phase angle is on the order of +p, -p, pi/2-p,
      !         or -pi/2+p.
      !
      !***references  handbook of mathematical functions by m. abramowitz
      !                 and i. a. stegun, nbs ams series 55, u.s. dept. of
      !                 commerce, 1955.
      !
      !               computation of bessel functions of complex argument
      !                 by d. e. amos, sand83-0083, may, 1983.
      !
      !               computation of bessel functions of complex argument
      !                 and large order by d. e. amos, sand83-0643, may, 1983
      !
      !               a subroutine package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, sand85-
      !                 1018, may, 1985
      !
      !               a portable package for bessel functions of a complex
      !                 argument and nonnegative order by d. e. amos, trans.
      !                 math. software, 1986
      !
      !***routines called  zbesh,i1mach,d1mach
      !***end prologue  zbesy
      !
      !     complex cwrk,cy,c1,c2,ex,hci,z,zu,zv
      complex(8),intent(in)::z
      real(8),intent(in)::fnu
      integer,intent(in)::n,kode
      integer,intent(inout)::nz
      complex(8),intent(inout)::cy(n),cwrk(n)
      real(8) ::elim, ey, hcii, tay,&
         ascle, rtol, atol, aa, bb, tol
      integer i, ierr, k, nz1, nz2
      complex(8)::st,ex,c1,c2
      !***first executable statement  zbesy
      ierr = 0
      nz = 0
      if (iszero(z)) ierr = 1
      if (fnu<0.0d0) ierr = 1
      if (kode<1 .or. kode>2) ierr = 1
      if (n<1) ierr = 1
      if (ierr/=0) return
      hcii = 0.5d0
      call f90_zbesh(z, fnu, kode, 1, n, cy, nz1, ierr)
      if (ierr/=0 .and. ierr/=3) goto 170
      call f90_zbesh(z, fnu, kode, 2, n, cwrk, nz2, ierr)
      if (ierr/=0 .and. ierr/=3) goto 170
      nz = min(nz1, nz2)
      if (kode==2) goto 60
      do i = 1, n
         st = cwrk(i) - cy(i)
         cy(i)=cmplx(-st%im,st%re,8)*hcii
      end do
      return
      60 continue
      tol = max(d1mach4, 1.0d-18)
      k = min(abs(i1mach15), abs(i1mach16))
      !-----------------------------------------------------------------------
      !     elim is the approximate exponential under- and overflow limit
      !-----------------------------------------------------------------------
      elim = 2.303d0*(real(k,8)*d1mach5-3.0d0)
      ex%re = cos(z%re)
      ex%im = sin(z%re)
      ey = 0.0d0
      tay = abs(z%im+z%im)
      if (tay<elim) ey = dexp(-tay)
      if (z%im<0.0d0) goto 90
      c1 = ex*ey
      c2 = conjg(ex)
      70 continue
      nz = 0
      rtol = 1.0d0/tol
      ascle = d1mach1*rtol*1.0d+3
      do i = 1, n
         !       str = c1r*cyr(i) - c1i*cyi(i)
         !       sti = c1r*cyi(i) + c1i*cyr(i)
         !       str = -str + c2r*cwrkr(i) - c2i*cwrki(i)
         !       sti = -sti + c2r*cwrki(i) + c2i*cwrkr(i)
         !       cyr(i) = -sti*hcii
         !       cyi(i) = str*hcii
         aa = cwrk(i)%re
         bb = cwrk(i)%im
         atol = 1.0d0
         if (max(abs(aa),abs(bb))>ascle) goto 75
         aa = aa*rtol
         bb = bb*rtol
         atol = tol
         75 continue
         st%re = (aa*c2%re-bb*c2%im)*atol
         st%im = (aa*c2%im+bb*c2%re)*atol

         aa = cy(i)%re
         bb = cy(i)%im
         atol = 1.0d0
         if (max(abs(aa),abs(bb))>ascle) goto 85
         aa = aa*rtol
         bb = bb*rtol
         atol = tol
         85 continue
         st%re = st%re - (aa*c1%re-bb*c1%im)*atol
         st%im = st%im - (aa*c1%im+bb*c1%re)*atol
         cy(i)%re = -st%im*hcii
         cy(i)%im = st%re*hcii
         if (iszero(st).and. iszero(ey)) nz = nz + 1
      end do
      return
      90 continue
      c1 = ex
      c2 = conjg(ex)*ey
      goto 70
      170 continue
      nz = 0
      return
   end subroutine f90_zbesy
end module amos
