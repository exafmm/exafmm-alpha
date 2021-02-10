program main
  implicit none
  include 'mpif.h'
  integer i,n,images,ierr,mpisize,mpirank
  integer,dimension (128) :: iseed
  real(8) diff,norm
  real(8),parameter :: pi=3.14159265358979312d0
  real(8),allocatable,dimension(:) :: x,g,u,ud
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,mpisize,ierr)
  call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
  n = 10000
  allocate( x(3*n),g(3*n),u(3*n),ud(3*n) )
  do i = 1,128
     iseed(i) = mpirank
  enddo
  call random_seed(put=iseed)
  call random_number(x)
  call random_number(g)
  do i = 1,n
     x(3*i-2) = (x(3*i-2) - 0.5) * pi
     x(3*i-1) = (x(3*i-1) - 0.5) * pi
     x(3*i-0) = (x(3*i-0) - 0.5) * pi
     g(3*i-2) = (g(3*i-2) - 0.5) / n
     g(3*i-1) = (g(3*i-1) - 0.5) / n
     g(3*i-0) = (g(3*i-0) - 0.5) / n
     u(3*i-2) = 0
     u(3*i-1) = 0
     u(3*i-0) = 0
     ud(3*i-2) = 0
     ud(3*i-1) = 0
     ud(3*i-0) = 0
  enddo
  images = 0
  call fmm_init(images)
  call fmm_partition(n,x,g,u)
  call fmm_biot_savart(n,x,g,u)
  call direct_biot_savart(n,x,g,ud)
  diff = 0
  norm = 0
  do i = 1,n
     diff = diff + (u(3*i-2) - ud(3*i-2)) ** 2
     diff = diff + (u(3*i-1) - ud(3*i-1)) ** 2
     diff = diff + (u(3*i-0) - ud(3*i-0)) ** 2
     norm = norm + ud(3*i-2) ** 2
     norm = norm + ud(3*i-1) ** 2
     norm = norm + ud(3*i-0) ** 2
  enddo
  call mpi_allreduce(diff,diff,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(norm,norm,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  if (mpirank.eq.0) print '(a,es12.5)',"error         :",sqrt(diff/norm)
  call mpi_finalize(ierr)
  call fmm_finalize()
  deallocate( x,g,u,ud )
end program main
