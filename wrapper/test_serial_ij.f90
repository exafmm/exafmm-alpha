program main
  implicit none
  integer i,ni,nj,images
  integer,dimension (128) :: iseed
  real(8) diff,norm
  real(8),allocatable,dimension(:) :: xi,ui,ud,xj,gj
  real(8) kappa4pi, kappa, pi
  integer npl

  real(8) xi_min,xi_max,xj_min,xj_max,gj_min,gj_max

  pi = 4.0D0*atan(1.0D0)
  kappa = 1.003078798D-7
  kappa4pi=7.982247452D-9

!t  open(50,file='fort.50')
!t  read(50,*) ni

!t  open(51,file='fort.51')
!t  read(51,*) nj, npl
!t  nj = nj*npl

  ni = 10000
  nj = 20000
!  nj = 80000
!  nj = 640*8
!  nj = 8
!  nj = 8*2
!  nj = 8*3
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
!t     read(51,*) xj(3*i-2), gj(3*i-2)
!t     read(51,*) xj(3*i-1), gj(3*i-1)
!t     read(51,*) xj(3*i-0), gj(3*i-0)
  enddo
!t  close(51)

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
!t     read(50,*) xi(3*i-2)
!t     read(50,*) xi(3*i-1)
!t     read(50,*) xi(3*i-0)
  enddo
!t  close(50)

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
  write(*,*) "xi= ",xi_min,xi_max
  write(*,*) "xj= ",xj_min,xj_max
  write(*,*) "gj= ",gj_min,gj_max

!t  do i = 1,ni*3
!t   xi(i) = xi(i) * pi*1.0D3
!t!   xi(i) = xi(i) * 1.0D3
!t  enddo
!t  do i = 1,nj*3
!t   xj(i) = xj(i) * pi*1.0D3
!t   gj(i) = gj(i) * pi*1.0D3
!t!   xj(i) = xj(i) * 1.0D3
!t!   gj(i) = gj(i) * 1.0D3
!t  enddo

!t  xi_min = 1.0D50
!t  xi_max =-1.0D50
!t  do i = 1,ni*3
!t    xi_min = min(xi(i),xi_min)
!t    xi_max = max(xi(i),xi_max)
!t  enddo
!t  xj_min = 1.0D50
!t  gj_min = 1.0D50
!t  xj_max =-1.0D50
!t  gj_max =-1.0D50
!t  do i = 1,nj*3
!t    xj_min = min(xj(i),xj_min)
!t    gj_min = min(gj(i),gj_min)
!t    xj_max = max(xj(i),xj_max)
!t    gj_max = max(gj(i),gj_max)
!t  enddo
!t  write(*,*) "xi= ",xi_min,xi_max
!t  write(*,*) "xj= ",xj_min,xj_max
!t  write(*,*) "gj= ",gj_min,gj_max

  images = 1
  call fmm_init(images)
  call fmm_biot_savart(ni,xi,ui,nj,xj,gj)
  call direct_biot_savart(ni,xi,ud,nj,xj,gj)

!t  do i = 1,ni*3
!t     ui(i)=ui(i)*pi*1.0D3
!t     ud(i)=ud(i)*pi*1.0D3
!t  enddo

!  open(60,file='ud-images0-ring')
!  open(60,file='ud-images1-ring')
!  open(60,file='ud-images3-ring')
!   do i = 1,ni*3
!     write(60,*) ud(i)
!   enddo
!  close(60)
!  open(60,file='ud-images0-ring')
!t  open(60,file='ud-images1-ring')
!  open(60,file='ud-images3-ring')
!t   do i = 1,ni*3
!t     read(60,*) ui(i)
!t   enddo
!t  close(60)

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
  print '(a,3es15.6)',"ui(1:3) ",ui(1),ui(2),ui(3)
  print '(a,3es15.6)',"ud(1:3) ",ud(1),ud(2),ud(3)
!  print '(a,3es15.6)',"ui(1:3) ",ui(1)*4*pi,ui(2)*4*pi,ui(3)*4*pi
!  print '(a,3es15.6)',"ud(1:3) ",ud(1)*4*pi,ud(2)*4*pi,ud(3)*4*pi
!t  print '(a,3es15.6)',"ui(1:3) ",ui(1)*kappa,ui(2)*kappa,ui(3)*kappa
!t  print '(a,3es15.6)',"ud(1:3) ",ud(1)*kappa,ud(2)*kappa,ud(3)*kappa
!t  print '(a,3es15.6)',"ui(4:6) ",ui(4)*kappa,ui(5)*kappa,ui(6)*kappa
!t  print '(a,3es15.6)',"ud(4:6) ",ud(4)*kappa,ud(5)*kappa,ud(6)*kappa
  call fmm_finalize()
  deallocate( xi,ui,ud,xj,gj )
end program main
