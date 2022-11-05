module amos
   use amos_utils
   implicit none
   private
   public::f90_zs1s2,f90_zuchk,f90_zrati
   public::f90_zasyi,f90_zseri,f90_zunik
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
         acs = -zr%re + dlog(as)
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
      !     zasyi computes the i bessel function for real(z).ge.0.0 by
      !     means of the asymptotic expansion for large cabs(z) in the
      !     region cabs(z).gt.max(rl,fnu*fnu/2). nz=0 is a normal return.
      !     nz.lt.0 indicates an overflow on kode=1.
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
      write(*,*)'call f90_zasyi'
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
      complex(8)::crfn,s,sr,t,t2,zn,st
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
      write(*,*)"call f90_zunik"
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
      zeta1 = 2.0d0*dabs(dlog(test)) + fnu
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
end module amos
