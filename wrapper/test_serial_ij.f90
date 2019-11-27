program main
  implicit none
  integer i,ni,nj,images
  integer,dimension (128) :: iseed
  real(8) diff,norm
  real(8),allocatable,dimension(:) :: xi,ui,ud,xj,gj
  real(8) kappa4pi, kappa, pi
  integer npl

  pi = 4.0D0*atan(1.0D0)
  kappa = 1.003078798D-7
  kappa4pi=7.982247452D-9

!  open(50,file='fort.50')
!  read(50,*) ni

!  open(51,file='fort.51')
!  read(51,*) nj, npl
!  nj = nj*npl

  ni = 10000
  nj = 20000
  allocate( xi(3*ni),ui(3*ni),ud(3*ni),xj(3*nj),gj(3*nj) )
  do i = 1,128
     iseed(i) = 0
  enddo
  call random_seed(put=iseed)
  call random_number(xi)
  call random_number(xj)
  call random_number(gj)

  do i = 1,nj
     xj(3*i-2) = (xj(3*i-2) - 0.5) * pi
     xj(3*i-1) = (xj(3*i-1) - 0.5) * pi
     xj(3*i-0) = (xj(3*i-0) - 0.5) * pi
     gj(3*i-2) = (gj(3*i-2) - 0.5) / ni
     gj(3*i-1) = (gj(3*i-1) - 0.5) / ni
     gj(3*i-0) = (gj(3*i-0) - 0.5) / ni
!     read(51,*) xj(3*i-2), gj(3*i-2)
!     read(51,*) xj(3*i-1), gj(3*i-1)
!     read(51,*) xj(3*i-0), gj(3*i-0)
  enddo
!  close(51)

  do i = 1,ni
     xi(3*i-2) = (xi(3*i-2) - 0.5) * pi
     xi(3*i-1) = (xi(3*i-1) - 0.5) * pi
     xi(3*i-0) = (xi(3*i-0) - 0.5) * pi
     ui(3*i-2) = 0
     ui(3*i-1) = 0
     ui(3*i-0) = 0
     ud(3*i-2) = 0
     ud(3*i-1) = 0
     ud(3*i-0) = 0
!     read(50,*) xi(3*i-2)
!     read(50,*) xi(3*i-1)
!     read(50,*) xi(3*i-0)
  enddo
!  close(50)

  images = 3
  call fmm_init(images)
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
  call fmm_finalize()
  deallocate( xi,ui,ud,xj,gj )
end program main
