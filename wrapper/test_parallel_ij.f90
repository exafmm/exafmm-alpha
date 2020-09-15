program main
  implicit none
  include 'mpif.h'
  integer i,ni,nj,images,ierr,mpisize,mpirank
  integer,dimension (128) :: iseed
  real(8),parameter :: pi=3.14159265358979312d0
  real(8) diff,norm
  real(8),allocatable,dimension(:) :: xi,ui,ud,xj,gj

  real(8),parameter :: kappa = 1.003078798D-7
  real(8),parameter :: kappa4pi=7.982247452D-9
  integer npl
  real(8) xi_min,xi_max,xj_min,xj_max,gj_min,gj_max
  integer nir,njr

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,mpisize,ierr)
  call mpi_comm_rank(mpi_comm_world,mpirank,ierr)

write(*,*) "mpisize ", mpisize

!  open(70,file='fort.70')
!  open(85,file='fort.85')
  open(85,file='fort.50')
  read(85,*) ni

!  open(71,file='fort.71')
!  open(86,file='fort.86')
  open(86,file='fort.51')
  read(86,*) nj, npl
  nj = nj*npl

!  ni = 10000
!  nj = 20000
!  ni = 10000
!  nj = 20000
!  ni = 6
!  nj = 6
  nir = ni
  njr = nj
if( mpisize.ne.1 ) then
  ni = ni/mpisize
  nj = nj/mpisize
  nir = ni*mpisize
  njr = nj*mpisize
end if
  write(*,'(a,i3,4i7)') "[before partition] mpirank, nir, njr, ni, nj= ", mpirank, nir, njr, ni, nj

!!  allocate( xi(3*ni),ui(3*ni),ud(3*ni),xj(3*nj),gj(3*nj) )
  allocate( xi(mpisize*3*ni),ui(mpisize*3*ni),ud(mpisize*3*ni),xj(mpisize*3*nj),gj(mpisize*3*nj) )
!  do i = 1,128
!     iseed(i) = mpirank
!     iseed(i) = 0
!  enddo
!  call random_seed(put=iseed)
!  call random_number(xi)
!  call random_number(xj)
!  call random_number(gj)

  do i = 1,njr
!     xj(3*i-2) = (xj(3*i-2) - 0.5) * pi
!     xj(3*i-1) = (xj(3*i-1) - 0.5) * pi
!     xj(3*i-0) = (xj(3*i-0) - 0.5) * pi
!     gj(3*i-2) = (gj(3*i-2) - 0.5) / ni
!     gj(3*i-1) = (gj(3*i-1) - 0.5) / ni
!     gj(3*i-0) = (gj(3*i-0) - 0.5) / ni
!     write(71,*) xj(3*i-2), gj(3*i-2)
!     write(71,*) xj(3*i-1), gj(3*i-1)
!     write(71,*) xj(3*i-0), gj(3*i-0)
     read(86,*) xj(3*i-2), gj(3*i-2)
     read(86,*) xj(3*i-1), gj(3*i-1)
     read(86,*) xj(3*i-0), gj(3*i-0)
  enddo
  close(86)

  do i = 1,nir
!     xi(3*i-2) = (xi(3*i-2) - 0.5) * pi
!     xi(3*i-1) = (xi(3*i-1) - 0.5) * pi
!     xi(3*i-0) = (xi(3*i-0) - 0.5) * pi
     ui(3*i-2) = 0
     ui(3*i-1) = 0
     ui(3*i-0) = 0
     ud(3*i-2) = 0
     ud(3*i-1) = 0
     ud(3*i-0) = 0
!     write(70,*) xi(3*i-2)
!     write(70,*) xi(3*i-1)
!     write(70,*) xi(3*i-0)
     read(85,*) xi(3*i-2)
     read(85,*) xi(3*i-1)
     read(85,*) xi(3*i-0)
  enddo
  close(85)



! if( mpirank==0 ) then
!  do i = 1,nir
!     write(*,*) mpirank, i, xi(3*i-2)
!     write(*,*) mpirank, i, xi(3*i-1)
!     write(*,*) mpirank, i, xi(3*i-0)
!  end do
!     write(*,*)
!  do i = 1,njr
!     write(*,*) mpirank, i, xj(3*i-2), gj(3*i-2)
!     write(*,*) mpirank, i, xj(3*i-1), gj(3*i-1)
!     write(*,*) mpirank, i, xj(3*i-0), gj(3*i-0)
!  end do
! end if
!     write(*,*)
! if( mpirank==1 ) then
!  do i = 1,nir
!     write(*,*) mpirank, i, xi(3*i-2)
!     write(*,*) mpirank, i, xi(3*i-1)
!     write(*,*) mpirank, i, xi(3*i-0)
!  end do
!     write(*,*)
!  do i = 1,njr
!     write(*,*) mpirank, i, xj(3*i-2), gj(3*i-2)
!     write(*,*) mpirank, i, xj(3*i-1), gj(3*i-1)
!     write(*,*) mpirank, i, xj(3*i-0), gj(3*i-0)
!  end do
! end if
!     write(*,*)
!
 if( mpirank==1 ) then
  do i = 1,ni
     xi(3*i-2) = xi(3*i-2+3*ni)
     xi(3*i-1) = xi(3*i-1+3*ni)
     xi(3*i-0) = xi(3*i-0+3*ni)
!     write(*,*) mpirank, i, xi(3*i-2)
!     write(*,*) mpirank, i, xi(3*i-1)
!     write(*,*) mpirank, i, xi(3*i-0)
  enddo
!!     write(*,*)
  do i = 1,nj
     xj(3*i-2) = xj(3*i-2+3*nj)
     xj(3*i-1) = xj(3*i-1+3*nj)
     xj(3*i-0) = xj(3*i-0+3*nj)
     gj(3*i-2) = gj(3*i-2+3*nj)
     gj(3*i-1) = gj(3*i-1+3*nj)
     gj(3*i-0) = gj(3*i-0+3*nj)
!     write(*,*) mpirank, i, xj(3*i-2), gj(3*i-2)
!     write(*,*) mpirank, i, xj(3*i-1), gj(3*i-1)
!     write(*,*) mpirank, i, xj(3*i-0), gj(3*i-0)
  enddo
 end if

  xi_min = 1.0D50
  xi_max =-1.0D50
  do i = 1,ni*3
    xi_min = min(xi(i),xi_min)
    xi_max = max(xi(i),xi_max)
  enddo
  xj_min = 1.0D50
  gj_min = 1.0D50
  xj_max =-1.0D50
  gj_max =-1.0D50
  do i = 1,nj*3
    xj_min = min(xj(i),xj_min)
    gj_min = min(gj(i),gj_min)
    xj_max = max(xj(i),xj_max)
    gj_max = max(gj(i),gj_max)
  enddo
if(mpirank==0) then
  write(*,*) "xi= ",xi_min,xi_max
  write(*,*) "xj= ",xj_min,xj_max
  write(*,*) "gj= ",gj_min,gj_max
end if

  do i = 1,ni*3
   xi(i) = xi(i) * pi*1.0D3
!   xi(i) = xi(i) * 1.0D3
  enddo
  do i = 1,nj*3
   xj(i) = xj(i) * pi*1.0D3
   gj(i) = gj(i) * pi*1.0D3
!   xj(i) = xj(i) * 1.0D3
!   gj(i) = gj(i) * 1.0D3
  enddo

  xi_min = 1.0D50
  xi_max =-1.0D50
  do i = 1,ni*3
    xi_min = min(xi(i),xi_min)
    xi_max = max(xi(i),xi_max)
  enddo
  xj_min = 1.0D50
  gj_min = 1.0D50
  xj_max =-1.0D50
  gj_max =-1.0D50
  do i = 1,nj*3
    xj_min = min(xj(i),xj_min)
    gj_min = min(gj(i),gj_min)
    xj_max = max(xj(i),xj_max)
    gj_max = max(gj(i),gj_max)
  enddo
if(mpirank==0) then
  write(*,*) "xi= ",xi_min,xi_max
  write(*,*) "xj= ",xj_min,xj_max
  write(*,*) "gj= ",gj_min,gj_max
end if

!--
!  call mpi_init(ierr)
!  call mpi_comm_size(mpi_comm_world,mpisize,ierr)
!  call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
!--

  images = 3
  call fmm_init(images)
  call fmm_partition(ni,xi,nj,xj,gj)
  call fmm_biot_savart(ni,xi,ui,nj,xj,gj)
  call direct_biot_savart(ni,xi,ud,nj,xj,gj)

  do i = 1,ni*3
     ui(i)=ui(i)*pi*1.0D3
     ud(i)=ud(i)*pi*1.0D3
  enddo


 if( mpirank==0 ) then
 open(80,file='fort.80')
  do i = 1,ni
     write(80,*) xi(3*i-2)
     write(80,*) xi(3*i-1)
     write(80,*) xi(3*i-0)
  end do
 close(80)
 end if
 if( mpirank==1 ) then
 open(81,file='fort.81')
  do i = 1,ni
     write(81,*) xi(3*i-2)
     write(81,*) xi(3*i-1)
     write(81,*) xi(3*i-0)
  end do
 close(81)
 end if
 if( mpirank==0 ) then
 open(82,file='fort.82')
  do i = 1,nj
     write(82,*) xj(3*i-2), gj(3*i-2)
     write(82,*) xj(3*i-1), gj(3*i-1)
     write(82,*) xj(3*i-0), gj(3*i-0)
  end do
 close(82)
 end if
 if( mpirank==1 ) then
 open(83,file='fort.83')
  do i = 1,nj
     write(83,*) xj(3*i-2), gj(3*i-2)
     write(83,*) xj(3*i-1), gj(3*i-1)
     write(83,*) xj(3*i-0), gj(3*i-0)
  end do
 close(83)
 end if


if(mpirank==0) then
  write(*,'(a,i3,4i7)') "mpirank, nir, njr, ni, nj= ", mpirank, nir, njr, ni, nj
end if

  open(73,file='u-mpirank0-images0')
!  open(73,file='u-mpirank0-images1')
!  open(73,file='u-mpirank0-images3')
!   do i = 1,nir*3
!!     write(73,*) ui(i)
!     write(73,*) ud(i)
!     write(*,*) i, ud(i)
!   enddo
!   do i = 1,nir*3
!     read(73,*) ud(i)
!   enddo
  close(73)

!if(mpirank==0) then
!   do i = 1,nir
!    print '(a,i5,3es15.6)',"i, ud ",i,ud(3*i-2)*kappa,ud(3*i-1)*kappa,ud(3*i-0)*kappa
!   enddo
!end if
if(mpirank==1) then
  write(*,'(a,i3,4i7)') "mpirank, nir, njr, ni, nj= ", mpirank, nir, njr, ni, nj
!  open(74,file='u-mpirank1-images0')
!  open(74,file='u-mpirank1-images1')
!  open(74,file='u-mpirank1-images3')
!   do i = 1,nir*3
!     write(74,*) ui(i)
!   enddo
!   do i = 1,nir*3
!     read(74,*) ud(i)
!   enddo
!  close(74)
end if

  diff = 0
  norm = 0
  do i = 1,ni
     diff = diff + (ui(3*i-2) - ud(3*i-2)) ** 2
     diff = diff + (ui(3*i-1) - ud(3*i-1)) ** 2
     diff = diff + (ui(3*i-0) - ud(3*i-0)) ** 2
     norm = norm + ud(3*i-2) ** 2
     norm = norm + ud(3*i-1) ** 2
     norm = norm + ud(3*i-0) ** 2
  enddo
!  call mpi_allreduce(diff,diff,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!  call mpi_allreduce(norm,norm,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
if(mpirank==0) then
  write(*,'(a,i3,4i7)') "mpirank, nir, njr, ni, nj= ", mpirank, nir, njr, ni, nj
  print '(a,es12.5)',"error         :",sqrt(diff/norm)
!  do i=1, ni
  do i=1, 3
    print '(a,i5,3es15.6)',"i, ui ",i,ui(3*i-2)*kappa,ui(3*i-1)*kappa,ui(3*i-0)*kappa
    print '(a,i5,3es15.6)',"i, ud ",i,ud(3*i-2)*kappa,ud(3*i-1)*kappa,ud(3*i-0)*kappa
  end do
end if
if(mpirank==1) then
  write(*,'(a,i3,4i7)') "mpirank, nir, njr, ni, nj= ", mpirank, nir, njr, ni, nj
  print '(a,es12.5)',"error         :",sqrt(diff/norm)
!  do i=1, ni
  do i=1, 3
    print '(a,i5,3es15.6)',"i, ui ",i,ui(3*i-2)*kappa,ui(3*i-1)*kappa,ui(3*i-0)*kappa
    print '(a,i5,3es15.6)',"i, ud ",i,ud(3*i-2)*kappa,ud(3*i-1)*kappa,ud(3*i-0)*kappa
  end do
end if
  call mpi_finalize(ierr)
  call fmm_finalize()
  deallocate( xi,ui,ud,xj,gj )
end program main
