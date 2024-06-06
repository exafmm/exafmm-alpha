program main
  implicit none
  include 'mpif.h'
  integer i,ni,nj,images,ierr,mpisize,mpirank
  integer,dimension (3) :: image_dim
  integer,dimension (128) :: iseed
  real(8) diff,norm
  real(8),parameter :: pi=3.14159265358979312d0
  real(8),allocatable,dimension(:) :: xi,ui,ud,xj,gj
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,mpisize,ierr)
  call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
  ni = 10000
  nj = 20000
  allocate( xi(6*ni),ui(6*ni),ud(6*ni),xj(6*nj),gj(6*nj) )
  do i = 1,128
     iseed(i) = mpirank+1
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
  enddo
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
  enddo
  images = 0
  image_dim(1) = 1
  image_dim(2) = 1
  image_dim(3) = 1
  call fmm_init(images,image_dim)
  call fmm_partition(ni,xi,nj,xj,gj)
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
  call mpi_allreduce(diff,diff,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(norm,norm,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  if (mpirank.eq.0) print '(a,es12.5)',"error         :",sqrt(diff/norm)
  call mpi_finalize(ierr)
  call fmm_finalize()
  deallocate( xi,ui,ud,xj,gj )
end program main
