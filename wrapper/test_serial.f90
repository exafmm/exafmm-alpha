program main
  implicit none
  integer i,ni,nj
  integer,dimension (128) :: iseed
  real(8) diff,norm
  real(8),allocatable,dimension(:) :: xi,ui,ud,xj,gj
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
     gj(3*i-2) = (gj(3*i-2) - 0.5) / ni
     gj(3*i-1) = (gj(3*i-1) - 0.5) / ni
     gj(3*i-0) = (gj(3*i-0) - 0.5) / ni
  enddo
  do i = 1,ni
     ui(3*i-2) = 0
     ui(3*i-1) = 0
     ui(3*i-0) = 0
     ud(3*i-2) = 0
     ud(3*i-1) = 0
     ud(3*i-0) = 0
  enddo
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
  call fmm_finalize()
  deallocate( xi,ui,ud,xj,gj )
end program main
