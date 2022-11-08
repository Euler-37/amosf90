program check
   use amos
   implicit none
   integer,parameter::n=10
   integer,parameter::step=200
   integer::i,ierr,nz,j
   complex(8)::z,rc(n),cwork(n)
   call random_seed(put=[(i,i=1,33)])
   tmpx=1.d+1
   tmpy=1.d+1
   do i=1,step
       call random_number(z%re)
       call random_number(z%im)
       call f90_zbiry(z, 0, 2, zz, ierr)
       call f90_zbesi(z, 1.4d0, 2, n, rc, nz, ierr)
       call f90_zbesj(z, 0.1d0, 1, n, rc, nz, ierr)
       call f90_zbesk(z, 0.0d0, 2, n, rc, nz, ierr)
       call f90_zbesh(z, 1.4d0, 2, 1, n , rc, nz   ,ierr)
       call f90_zbesy(z, 0.1d0, 2, n, rc, nz, cwork,ierr)
   end do
end program check
