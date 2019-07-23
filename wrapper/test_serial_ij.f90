program main
  implicit none
  integer i,ni,nj
  integer,dimension (128) :: iseed
  real(8) diff,norm
  real(8),allocatable,dimension(:) :: xi,ui,ud,xj,gj
  real(8) kappa4pi, kappa, pi
  integer npl

  pi = 4.0D0*atan(1.0D0)
  kappa = 1.003078798D-7
  kappa4pi=7.982247452D-9

  open(50,file='fort.50')
  read(50,*) ni

  open(51,file='fort.51')
  read(51,*) nj, npl
  nj = nj*npl

!  ni = 1000
!  nj = 20000
!  nj = 80000
!  nj = 640*8
!  nj = 8
!  nj = 8*2
!  nj = 8*3
  allocate( xi(3*ni),ui(3*ni),ud(3*ni),xj(3*nj),gj(3*nj) )
!  do i = 1,128
!     iseed(i) = 0
!  enddo
!  call random_seed(put=iseed)
!  call random_number(xi)
!  call random_number(xj)
!  call random_number(gj)

  do i = 1,nj
!     gj(3*i-2) = (gj(3*i-2) - 0.5) / ni
!     gj(3*i-1) = (gj(3*i-1) - 0.5) / ni
!     gj(3*i-0) = (gj(3*i-0) - 0.5) / ni
     read(51,*) xj(3*i-2), gj(3*i-2)
     read(51,*) xj(3*i-1), gj(3*i-1)
     read(51,*) xj(3*i-0), gj(3*i-0)
  enddo
  close(51)

  do i = 1,ni
     ui(3*i-2) = 0
     ui(3*i-1) = 0
     ui(3*i-0) = 0
     ud(3*i-2) = 0
     ud(3*i-1) = 0
     ud(3*i-0) = 0
     read(50,*) xi(3*i-2)
     read(50,*) xi(3*i-1)
     read(50,*) xi(3*i-0)
  enddo
  close(50)

  call fmm_init()
  call fmm_biot_savart(ni,xi,ui,nj,xj,gj)
  call direct_biot_savart(ni,xi,ud,nj,xj,gj)
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
  print '(a,es12.5)',"error         :",sqrt(diff/norm)
!  print '(a,3es15.6)',"ui(1:3) ",ui(1),ui(2),ui(3)
!  print '(a,3es15.6)',"ud(1:3) ",ud(1),ud(2),ud(3)
!  print '(a,3es15.6)',"ui(1:3) ",ui(1)*4*pi,ui(2)*4*pi,ui(3)*4*pi
!  print '(a,3es15.6)',"ud(1:3) ",ud(1)*4*pi,ud(2)*4*pi,ud(3)*4*pi
  print '(a,3es15.6)',"ui(1:3) ",ui(1)*kappa,ui(2)*kappa,ui(3)*kappa
  print '(a,3es15.6)',"ud(1:3) ",ud(1)*kappa,ud(2)*kappa,ud(3)*kappa
  call fmm_finalize()
  deallocate( xi,ui,ud,xj,gj )
end program main
