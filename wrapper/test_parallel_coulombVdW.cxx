#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
extern "C" {
#include "md.h"
#include "vtgrapeproto.h"
}
//#define DATAFILE

extern "C" void FMMcalccoulomb_ij(int ni, double* xi, double* qi, double* fi,
  int nj, double* xj, double* qj, double rscale, int tblno, double size, int periodicflag);

extern "C" void FMMcalcvdw_ij(int ni, double* xi, int* atypei, double* fi,
  int nj, double* xj, int* atypej, int nat, double* gscale, double* rscale,
  int tblno, double size, int periodicflag);

void MPI_ShiftI(int *var, int n, int mpisize, int mpirank) {
  int *buf = new int [n];
  const int isend = (mpirank + 1          ) % mpisize;
  const int irecv = (mpirank - 1 + mpisize) % mpisize;
  MPI_Request sreq, rreq;
  MPI_Isend(var,n,MPI_INT,irecv,1,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(buf,n,MPI_INT,isend,1,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);
  int i;
  for( i=0; i<n; i++ ) {
    var[i] = buf[i];
  }
  delete[] buf;
}

void MPI_Shift(double *var, int n, int mpisize, int mpirank) {
  double *buf = new double [n];
  const int isend = (mpirank + 1          ) % mpisize;
  const int irecv = (mpirank - 1 + mpisize) % mpisize;
  MPI_Request sreq, rreq;
  MPI_Isend(var,n,MPI_DOUBLE,irecv,1,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(buf,n,MPI_DOUBLE,isend,1,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);
  int i;
  for( i=0; i<n; i++ ) {
    var[i] = buf[i];
  }
  delete[] buf;
}

int main(int argc, char **argv) {
  MPI_Init(&argc,&argv);
  int mpisize, mpirank;
  MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
#ifdef DATAFILE
  const int N = 78537;
  const int nat = 22;
  const double size = 90;
#else
  const int N = 10000;
  const int nat = 16;
  const double size = 2;
#endif
  double *xi     = new double [3*N];
  double *qi     = new double [N];
  double *pi     = new double [3*N];
  double *fi     = new double [3*N];
  double *pd     = new double [3*N];
  double *fd     = new double [3*N];
  double *xj     = new double [3*N];
  double *qj     = new double [N];
  int *atypei    = new int [N];
  int *atypej    = new int [N];
  double *rscale = new double [64*64];
  double *gscale = new double [64*64];
  double *frscale = new double [64*64];
  double *fgscale = new double [64*64];
  int *numex  = new int [N];
  int *natex  = new int [N];
  MR3init();

#ifdef DATAFILE
  std::ifstream fid("datafile",std::ios::in);
  for( int i=0; i<nat*nat; i++ ) {
    fid >> rscale[i] >> gscale[i] >> frscale[i] >> fgscale[i];
  }
  int ic;
  for( int i=0; i<N; i++ ) {
    fid >> ic >> atypei[i] >> xi[3*i+0] >> xi[3*i+1] >> xi[3*i+2] >> qi[i];
  }
  fid.close();
  for( int i=0; i<N; i++ ) {
    xj[3*i+0] = xi[3*i+0];
    xj[3*i+1] = xi[3*i+1];
    xj[3*i+2] = xi[3*i+2];
    qj[i] = qi[i];
    atypei[i]--;
    atypej[i] = atypei[i];
  }
#else
  srand48(mpirank);
  float average = 0;
  for( int i=0; i<N; i++ ) {
    xi[3*i+0] = drand48() * size - size/2;
    xi[3*i+1] = drand48() * size - size/2;
    xi[3*i+2] = drand48() * size - size/2;
    qi[i] = drand48()*2.0-1.0;
    average += qi[i];
    atypei[i] = drand48() * nat;
  }
  average /= N;
  for( int i=0; i<N; i++ ) {
    qi[i] -= average;
  }
  average = 0;
  for( int i=0; i<N; i++ ) {
    xj[3*i+0] = drand48() * size - size/2;
    xj[3*i+1] = drand48() * size - size/2;
    xj[3*i+2] = drand48() * size - size/2;
    qj[i] = drand48()*2.0-1.0;
    atypej[i] = drand48() * nat;
  }
  average /= N;
  for( int i=0; i<N; i++ ) {
    qj[i] -= average;
  }

  if( mpirank == 0 ) {
    for( int i=0; i<nat; i++ ) {
      gscale[i*nat+i] = drand48();
      drand48();
      rscale[i*nat+i] = drand48();
      rscale[i*nat+i] = 1;
    }
    for( int i=0; i<nat; i++ ) {
      for( int j=0; j<nat; j++ ) {
        if( i != j ) {
          gscale[i*nat+j] = sqrt(gscale[i*nat+i]*gscale[j*nat+j]);
          rscale[i*nat+j] = (sqrt(rscale[i*nat+i]) + sqrt(rscale[j*nat+j])) * 0.5;
          rscale[i*nat+j] *= rscale[i*nat+j];
        }
      }
    }
  }
  MPI_Bcast(rscale,nat*nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(gscale,nat*nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
  for( int i=0; i<nat; i++ ) {
    for( int j=0; j<nat; j++ ) {
      frscale[i*nat+j] = rscale[i*nat+j];
      fgscale[i*nat+j] = gscale[i*nat+j];
    }
  }
#endif
  for( int i=0; i<N; i++ ) {
    pi[3*i+0] = pi[3*i+1] = pi[3*i+2] = 0;
    fi[3*i+0] = fi[3*i+1] = fi[3*i+2] = 0;
    pd[3*i+0] = pd[3*i+1] = pd[3*i+2] = 0;
    fd[3*i+0] = fd[3*i+1] = fd[3*i+2] = 0;
    numex[i] = 0;
    natex[i] = 0;
  }

  FMMcalccoulomb_ij(N, xi, qi, pi, N, xj, qj, 0.0, 1, size, 0);
  FMMcalccoulomb_ij(N, xi, qi, fi, N, xj, qj, 0.0, 0, size, 0);
  for( int irank=0; irank<mpisize; irank++ ) {
    MPI_Shift(xj,3*N,mpisize,mpirank);
    MPI_Shift(qj,N,mpisize,mpirank);
#if 1
    MR3calccoulomb_ij_host(N, xi, qi, pd, N, xj, qj, 1.0, 1, size, 2);
    MR3calccoulomb_ij_host(N, xi, qi, fd, N, xj, qj, 1.0, 0, size, 2);
#else
    for( int i=0; i<N; i++ ) {
      double P = 0, Fx = 0, Fy = 0, Fz = 0;
      for( int j=0; j<N; j++ ) {
        double dx = xi[3*i+0] - xj[3*j+0];
        double dy = xi[3*i+1] - xj[3*j+1];
        double dz = xi[3*i+2] - xj[3*j+2];
        double R2 = dx * dx + dy * dy + dz * dz;
        double invR = 1 / std::sqrt(R2);
        if( R2 == 0 ) invR = 0;
        double invR3 = qj[j] * invR * invR * invR;
        P += qj[j] * invR;
        Fx += dx * invR3;
        Fy += dy * invR3;
        Fz += dz * invR3;
      }
      pd[3*i+0] += P;
      fd[3*i+0] += Fx;
      fd[3*i+1] += Fy;
      fd[3*i+2] += Fz;
    }
#endif
  }
  double fc[3];
  for( int d=0; d<3; ++d ) fc[d]=0;
  for( int i=0; i<N; ++i ) {
    for( int d=0; d<3; ++d ) {
      fc[d] += qi[i] * xi[3*i+d];
    }
  }
  for( int i=0; i<N; ++i ) {
    pd[3*i+0] += 2.0 * M_PI / (3.0 * size * size * size)
              * (fc[0] * fc[0] + fc[1] * fc[1] + fc[2] * fc[2]) / N;
  }
  for( int i=0; i<N; ++i ) {
    for( int d=0; d<3; ++d ) {
      fd[3*i+d] -= 4.0 * M_PI * qi[i] * fc[d] / (3.0 * size * size * size);
    }
  }
  double Pd = 0, Pn = 0, Fd = 0, Fn = 0;
  for( int i=0; i<N; i++ ) {
    pd[3*i+0] *= 0.5;
    Pd += (pi[3*i+0] - pd[3*i+0]) * (pi[3*i+0] - pd[3*i+0]);
    Pn += pd[3*i+0] * pd[3*i+0];
    Fd += (fi[3*i+0] - fd[3*i+0]) * (fi[3*i+0] - fd[3*i+0])
        + (fi[3*i+1] - fd[3*i+1]) * (fi[3*i+1] - fd[3*i+1])
        + (fi[3*i+2] - fd[3*i+2]) * (fi[3*i+2] - fd[3*i+2]);
    Fn += fd[3*i+0] * fd[3*i+0] + fd[3*i+1] * fd[3*i+1] + fd[3*i+2] * fd[3*i+2];
    pi[3*i+0] = pi[3*i+1] = pi[3*i+2] = 0;
    fi[3*i+0] = fi[3*i+1] = fi[3*i+2] = 0;
    pd[3*i+0] = pd[3*i+1] = pd[3*i+2] = 0;
    fd[3*i+0] = fd[3*i+1] = fd[3*i+2] = 0;
  }
  std::cout << "Coulomb       potential @ rank " << mpirank << " : " << sqrtf(Pd/Pn) << std::endl;
  std::cout << "Coulomb       force     @ rank " << mpirank << " : " << sqrtf(Fd/Fn) << std::endl;

  FMMcalcvdw_ij(N,xi,atypei,pi,N,xj,atypej,nat,gscale,rscale,3,size,0);
  FMMcalcvdw_ij(N,xi,atypei,fi,N,xj,atypej,nat,fgscale,frscale,2,size,0);
  for( int irank=0; irank<mpisize; irank++ ) {
    MPI_Shift(xj,3*N,mpisize,mpirank);
    MPI_ShiftI(atypej,N,mpisize,mpirank);
#if 1
//    MR3calcvdw_ij_exlist(N,xi,atypei,pd,N,xj,atypej,nat,gscale,rscale,3,size,0,numex,natex);
//    MR3calcvdw_ij_exlist(N,xi,atypei,fd,N,xj,atypej,nat,fgscale,frscale,2,size,0,numex,natex);
    MR3calcvdw_ij_host(N,xi,atypei,pd,N,xj,atypej,nat,gscale,rscale,3,size,0);
    MR3calcvdw_ij_host(N,xi,atypei,fd,N,xj,atypej,nat,fgscale,frscale,2,size,0);
#else
    for( int i=0; i<N; i++ ) {
      double P = 0, Fx = 0, Fy = 0, Fz = 0;
      for( int j=0; j<N; j++ ) {
        double dx = xi[3*i+0] - xj[3*j+0];
        double dy = xi[3*i+1] - xj[3*j+1];
        double dz = xi[3*i+2] - xj[3*j+2];
        double R2 = dx * dx + dy * dy + dz * dz;
        if( R2 != 0 ) {
          double rs = rscale[atypei[i]*nat+atypej[j]];
          double gs = gscale[atypei[i]*nat+atypej[j]];
          double R2s = R2 * rs;
          if( R2MIN <= R2 && R2 < R2MAX ) {
            double invR2 = 1.0 / R2s;
            double invR6 = invR2 * invR2 * invR2;
            double dtmp = gs * invR6 * invR2 * (2.0 * invR6 - 1.0);
            P += gs * invR6 * (invR6 - 1.0);
            Fx += dx * dtmp;
            Fy += dy * dtmp;
            Fz += dz * dtmp;
          }
        }
      }
      pd[3*i+0] += P;
      fd[3*i+0] += Fx;
      fd[3*i+1] += Fy;
      fd[3*i+2] += Fz;
    }
#endif
  }
  Pd = Pn = Fd = Fn = 0;
  for( int i=0; i<N; i++ ) {
    Pd += (pi[3*i+0] - pd[3*i+0]) * (pi[3*i+0] - pd[3*i+0]);
    Pn += pd[3*i+0] * pd[3*i+0];
    Fd += (fi[3*i+0] - fd[3*i+0]) * (fi[3*i+0] - fd[3*i+0])
        + (fi[3*i+1] - fd[3*i+1]) * (fi[3*i+1] - fd[3*i+1])
        + (fi[3*i+2] - fd[3*i+2]) * (fi[3*i+2] - fd[3*i+2]);
    Fn += fd[3*i+0] * fd[3*i+0] + fd[3*i+1] * fd[3*i+1] + fd[3*i+2] * fd[3*i+2];
  }
  std::cout << "Van der Waals potential @ rank " << mpirank << " : " << sqrtf(Pd/Pn) << std::endl;
  std::cout << "Van der Waals force     @ rank " << mpirank << " : " << sqrtf(Fd/Fn) << std::endl;

  delete[] xi;
  delete[] qi;
  delete[] pi;
  delete[] fi;
  delete[] pd;
  delete[] fd;
  delete[] xj;
  delete[] qj;
  delete[] atypei;
  delete[] atypej;
  delete[] rscale;
  delete[] gscale;
  delete[] frscale;
  delete[] fgscale;
  delete[] numex;
  delete[] natex;

  MPI_Finalize();
}
