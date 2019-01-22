#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>

#include "vtgrape.h"
#include "md.h"            // application specific include file
#include "vtgrapeproto.h"
//#include "md.c"

void lj_mother(void *);
void ljpot_mother(void *);
void grav_mother(void *);
void gravpot_mother(void *);
void real_mother(void *);
void realpot_mother(void *);
void wave_mother(void *);
void wavepot_mother(void *);
void nacl_mother(void *);
void r1_mother(void *);
void rsqrt_mother(void *);
void mul_mother(void *);
void gravfb_mother(void *);
void gravpotfb_mother(void *);
void ljfb_mother(void *);
void ljpotfb_mother(void *);
void realfb_mother(void *);
void realpotfb_mother(void *);
#ifdef MD_CELLINDEX
void ljfbci_mother(void *);
void ljpotfbci_mother(void *);
void realci_mother(void *);
void realpotci_mother(void *);
void realfbci_mother(void *);
void realpotfbci_mother(void *);
void realljfbci_mother(void *);
void realljfbci2_mother(void *);
void realljpotfbci_nojdiv_mother(void *);
void realljpotfbci_mother(void *);
#endif


typedef union {
  int i;
  float f;
} FI;

typedef union {
  struct{
    FI fi0;
    FI fi1;
  } fi2;
  double d;
} DI2;

#if MD_PERIODIC_FIXED==1
#define COPY_POS(vgvec,pos,volume_1) \
        /*di2.d=(pos)*xmax_1+0x180000;*/\
        di2.d=(pos)*(volume_1)+0x180000;\
        (vgvec)=di2.fi2.fi0.f;\
        /*printf("in COPY_POS : pos=%e xmax_1=%e vgvec=%08x\n",pos,xmax_1,vgvec);*/
#define COPY_POS_INT(vgvec,pos,volume_1) \
        di2.d=(pos)*(volume_1)+0x180000;\
        (vgvec)=di2.fi2.fi0.i;\
        /*printf("in COPY_POS_INT : pos=%e xmax_1=%e vgvec=%08x\n",pos,xmax_1,vgvec);*/
#define COPY_POS_DEFINE \
        /*double xmax_1=1.0/unit->xmax;*/\
        double volume_1[3]={1.0/unit->volume[0],1.0/unit->volume[1],1.0/unit->volume[2]};\
        DI2 di2;
#else
#define COPY_POS(vgvec,pos,dum) \
        (vgvec)=(float)(pos);
#define COPY_POS_INT(vgvec,pos,dum) \
        (vgvec)=(float)(pos);
#define COPY_POS_DEFINE \
        double volume_1[3]={1.0/unit->volume[0],1.0/unit->volume[1],1.0/unit->volume[2]};
#endif


static double  Coulomb_vdw_factor=1.0;
static VG_UNIT *MR3_unit;
#ifdef MD_CELLINDEX
static int Ldim[3]={1,1,1};
#endif

#ifdef MD_GENERATE_EMUTABLE
static float *Emutable_r1=NULL;
static float *Emutable_rsqrt=NULL;
#endif
#ifdef MD_MEASURE_TIME
double Time[MD_MEASURE_TIME][3];
                  /* [][0]-laptime, [][1]-sprittime, [][2]-sum
		      0 -  9 : 
		     10 - 19 :
                     20 - 29 : MR3calccoulomb_vdw_nlist_ij_emu_work
                                 depending on tblno
		                 0 - coulomb force
				 1 - coulomb potential
				 6 - real-space force
				 7 - real-space potential
                     30 - 39 : MR3calccoulomb_vdw_ij_ci_virial
		                 0 - total
                                 1 - before assign cell
                                 2 - assign cell
                                 3 - copy positions to cell
				 4 - other setup on cell
				 5 - copy j data to GPU
				 6 - copy i data to GPU
				 7 - starting kernel on GPU
				 8 - wait before copying result
				 9 - copy result
                     40 - 49 : 
                                 0 - total mr3calccoulomb_vdw_ij_ci_exlist_
                                 1 - total MR3calccoulomb_vdw_ij_ci_exlist
				 2 -         MR3calccoulomb_vdw_ij_ci_virial
				 3 -         MR3calccoulomb_vdw_nlist_ij_emu
                                 4 - total mr3calccoulomb_vdw_exlist_
                                 5 -         before GPU call
                                 6 -         calling MR3calccoulomb_vdw_exlist
                                 7 -                  
				 8 -         MR3calccoulomb_vdw_ij_ci
                     50 - 59 : 
                                 0 - total MR3calccoulomb_vdw_ij_exlist
				 1 -         nonoverlap thread : MR3calccoulomb_vdw_ij_exlist_core
				 2 -         overlap thread : total
                                 3 -     potential emu : MR3calccoulomb_vdw_nlist_ij_emu in MR3calccoulomb_vdw_ij_exlist_core
				 4 -     potentail GPU : this part is divided into 60-69
				 5 -     force     GPU : this part is divided into 60-69
                                 6 -     force     emu : 
				 7 -     force copy    : 
                     60 - 69 :
                                 0 - total MR3calccoulomb_vdw_ij_ci_old_09
                                 1 - before assign cell
				 2 - assign cell
				 3 - copy positions to cell
				 4 - other setup on cell
				 5 - copy j data to GPU
				 6 - copy i data to GPU
				 7 - starting kernel on GPU
				 8 - wait before copying result
				 9 - copy result
		  */
#endif


static void vg_get_cputime(double *laptime, double *sprittime)
{
  struct timeval tv;
  struct timezone tz;
  double sec,microsec;

  gettimeofday(&tv, &tz);
  sec=tv.tv_sec;
  microsec=tv.tv_usec;

  *sprittime = sec + microsec * 1e-6 - *laptime;
  *laptime = sec + microsec * 1e-6;
  //  printf("    vg_get_cputime is called: ltime=%e stime=%e\n",*laptime,*sprittime);
}


#ifdef MD_MEASURE_TIME
static void vg_start_timer(int timerindex)
{
  vg_get_cputime(&Time[timerindex][0], &Time[timerindex][1]);
}


static void vg_stop_and_accumulate_timer(int timerindex)
{
  vg_get_cputime(&Time[timerindex][0], &Time[timerindex][1]);
  Time[timerindex][2]+=Time[timerindex][1];
}
#endif


void *MR3_my_malloc2(size_t n, char *s)
{
  void *ret;

#if 0
  {
    static int count=0;
    if(count<34){
      printf("** count=%d doubling the allocation **\n",count);
      n*=2;
    }
    count++;
  }
#endif
  printf("Allocating %d bytes for %s\n",(int)n,s);
  if(n==0){
    fprintf(stderr,"** warning : mallocing %d **\n",(int)n);
    //    sleep(1000);
    MR3_exit(1);
  }
  ret=malloc(n);
  printf("Allocated pointer=%016llx which is %s\n",
	 (unsigned long long)ret,s);

  return ret;
}


void MR3_my_free2(void **p, char *s)
{
  printf("Freeing pointer=%016llx which is %s\n",
	 (unsigned long long)(*p),s);
  if(*p==NULL){
    printf("** warning : pointer=%s is NULL before free **\n",s);
  }
  free(*p);
  *p=NULL;
}


VG_UNIT *m3_get_unit(void)
{
  return MR3_unit;
}


#if 1
void vg_exit(int ret)
{
  exit(ret);
}
#else
int vg_exit(int ret)
{
  return ret;
}
#endif


VG_UNIT *vg_allocate_unit(void (*function)(void *))
{
  int i;
  VG_UNIT *unit;
  char *s;

#ifdef MD_QAUNION_ATYPEBIT
  if((1<<MD_QAUNION_ATYPEBIT)*(1<<MD_QAUNION_ATYPEBIT)
     <VG_MINIMUM_ATYPE_BLOCK){
    fprintf(stderr,"** error : (2^MD_QAUNION_ATYPEBIT)^2 must not be lower than VG_MINIMUM_ATYPE_BLOCK **\n");
    vg_exit(1);
  }
#endif
  if((unit=(VG_UNIT *)MR3_malloc_pointer(sizeof(VG_UNIT),"*vg_allocate_unit"))==NULL){
    //    fprintf(stderr,"** error : can't malloc unit **\n")
    return NULL;
  }
  unit->deviceid=0;
  s=getenv("VG_DEVICEID");
  if(s!=NULL){
    sscanf(s,"%d",&(unit->deviceid));
    printf("VG_DEVICEID is set %d\n",unit->deviceid);
  }
  unit->nj=unit->ni=unit->nati=unit->natj=unit->nf=0;
  unit->jvectors=unit->ivectors=unit->matrices=NULL;
  unit->scalers=unit->fvectors=unit->pscalers=NULL;
  unit->calculate[0]=function;
  for(i=1;i<VG_MAX_FUNCTIONS;i++) unit->calculate[i]=NULL;
  unit->r1=NULL;
  unit->rsqrt=NULL;
  //  unit->fvec=NULL;
  unit->fthread=NULL;
  unit->thread=0;
  unit->ni_overlap=0;
  unit->potc=0.0;
  unit->potv=0.0;
  unit->rcut2=0.0;
  unit->cuton=CUTON_DEFAULT;
  unit->debug_flag=0;
  unit->jsize=unit->isize=unit->fsize=0;

  // for GPU overlap
  unit->gpuoverlapflag=unit->function_index=0;
#ifdef MD_QAUNION_ATYPEBIT
  //  unit->scaleqi_1=unit->scaleqj_1=0.0;
#endif

  if(unit->scalers==NULL){
    unit->ssize=sizeof(VG_SCALER);
    if((unit->scalers=(void *)MR3_malloc_pointer(unit->ssize,"vg_allocate_unit"))==NULL){
      fprintf(stderr,"** error : can't malloc scalers in vg_allocate_unit **\n");
      vg_exit(1);
    }
  }

  return unit;
}


void vg_free_unit(VG_UNIT *unit)
{
  MR3_free_pointer(unit->jvectors,"jvectors in vg_free_unit");
  MR3_free_pointer(unit->ivectors,"ivectors in vg_free_unit");
  MR3_free_pointer(unit->matrices,"matrices in vg_free_unit");
  MR3_free_pointer(unit->scalers,"scalers in vg_free_unit");
  MR3_free_pointer(unit->fvectors,"fvectoers in vg_free_unit");
  MR3_free_pointer(unit->pscalers,"pscalers in vg_free_unit");
  MR3_free_pointer(unit->r1,"r1 in vg_free_unit");
  MR3_free_pointer(unit->rsqrt,"rsqrt in vg_free_unit");
  //  MR3_free_pointer(unit->fvec,"fvec in vg_free_unit");
  MR3_free_pointer(unit->fthread,"fthread in vg_free_unit");
  MR3_free_pointer(unit,"unit in vg_free_unit");
}


void vg_set_function(VG_UNIT *unit, int function_index, void (*function)(void *))
{
  if(function_index>=VG_MAX_FUNCTIONS){
    fprintf(stderr,"** error : function_index=%d is out of range **\n",
	    function_index);
    vg_exit(1);
  }
  unit->calculate[function_index]=function;
}


void vg_set_vectors(VG_UNIT *unit, int nj, void *jvectors)
{
  int size;

  printf("** warning: free and malloc of jvectors should be modified **\n");
  size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
  unit->nj=nj;
  //  fprintf(stderr,"in vg_set_vectors: unit->jvectors=%x before free\n",unit->jvectors);
  MR3_free_pointer(unit->jvectors,"vg_set_vectors");
  if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors"))==NULL){
    fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors **\n");
    vg_exit(1);
  }
  //  fprintf(stderr,"in vg_set_vectors: unit->jvectors=%x after malloc %d\n",unit->jvectors,size);
  memcpy(unit->jvectors,jvectors,size);
}


void vg_set_matrices(VG_UNIT *unit, int nati, int natj, void *matrices)
{
  int size;

  if(nati*natj>VG_MINIMUM_ATYPE_BLOCK){
    fprintf(stderr,"** error : nati*natj>VG_MINIMUM_ATYPE_BLOCK(%d) **\n",
	    VG_MINIMUM_ATYPE_BLOCK);
    vg_exit(1);
  }
  size=sizeof(VG_MATRIX)*MATRIX_ROUNDUP(nati*natj);
  unit->nati=nati;
  unit->natj=natj;
  MR3_free_pointer(unit->matrices,"vg_set_matrices");
  if((unit->matrices=(void *)MR3_malloc_pointer(size,"vg_set_matrices"))==NULL){
    fprintf(stderr,"** error : can't malloc matrices in vg_set_vectors **\n");
    vg_exit(1);
  }
  memcpy(unit->matrices,matrices,size);
}


void vg_set_scalers(VG_UNIT *unit, void *scalers)
{
  unit->ssize=sizeof(VG_SCALER);
  MR3_free_pointer(unit->scalers,"vg_set_scalers");
  if((unit->scalers=(void *)MR3_malloc_pointer(unit->ssize,"vg_set_scalers"))==NULL){
    fprintf(stderr,"** error : can't malloc scalers in vg_set_vectors **\n");
    vg_exit(1);
  }
  memcpy(unit->scalers,scalers,unit->ssize);
}


void vg_set_pipeline_vectors(VG_UNIT *unit, int ni, void *ivectors)
{
  int size;

  printf("** warning: free and malloc of ivectors should be modified **\n");
  size=sizeof(VG_IVEC)*NI_ROUNDUP(ni);
  //  if(ni==0){fprintf(stderr,"** error : ni=0 in vg_set_pipeline_vectors **\n");exit(1);}
  unit->ni=ni;
  MR3_free_pointer(unit->ivectors,"vg_set_pipeline_vectors");
  if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors"))==NULL){
    fprintf(stderr,"** error : can't malloc ivectors in vg_set_vectors **\n");
    vg_exit(1);
  }
  memcpy(unit->ivectors,ivectors,size);
}


void vg_calculate(VG_UNIT *unit, int function_index, 
		  int nf, void *fvectors, void *pscalers)
{
  int fsize,psize;

  fsize=sizeof(VG_FVEC)*NF_ROUNDUP(nf);
  unit->nf=nf;
  MR3_free_pointer(unit->fvectors,"vg_calculate");
  if((unit->fvectors=(void *)MR3_malloc_pointer(fsize,"vg_calculate"))==NULL){
    fprintf(stderr,"** error : can't malloc fvectors in vg_calculate **\n");
    vg_exit(1);
  }
  bzero(unit->fvectors,fsize);
  psize=sizeof(VG_PSCALER);
#ifdef MD_CELLINDEX
  if(unit->pscalers==NULL){
    if((unit->pscalers=(VG_PSCALER *)MR3_malloc_pointer(sizeof(VG_PSCALER),"vg_calculate"))==NULL){
      fprintf(stderr,"** error : can't malloc unit->pscalers in vg_calculate **\n");
      vg_exit(1);
    }
  }
#else
  MR3_free_pointer(unit->pscalers,"vg_calculate");
  if((unit->pscalers=(void *)MR3_malloc_pointer(psize,"vg_calculate"))==NULL){
    fprintf(stderr,"** error : can't malloc pscalers in vg_calculate **\n");
    vg_exit(1);
  }
  bzero(unit->pscalers,psize);
#endif
  unit->calculate[function_index]((void *)unit);
  memcpy(fvectors,unit->fvectors,fsize);
  memcpy(pscalers,unit->pscalers,psize);
}


void vg_set_vectors_pos_charge_atype(VG_UNIT *unit, int nj,
				     double pos[][3], double charge[],
				     int atype[])
{
#if defined(MD_USE_QAUNION) && defined(MD_QAUNION_ATYPEBIT)
  int i,j,size;
  VG_JVEC *jvec;
  double scaleq;
  float scaleq_1;
  COPY_POS_DEFINE;

  size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
  //  size*=2;printf("** size of j-particles is doubled in vg_set_vectors_pos_charge_atype\n");

  unit->nj=nj;
  if(unit->jvectors==NULL){
    if(size<unit->jsize) size=unit->jsize;
    unit->jsize=size;
#ifdef MD_PRINT_WARN
    printf("allocating %d jvectors\n",size);
#endif
    if((unit->jvectors=(void *)MR3_malloc_pointer(size,"void vg_set_vectors_pos_charge_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge_atype **\n");
      vg_exit(1);
    }
  }
  else if(size>unit->jsize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d jvectors\n",size);
#endif
    MR3_free_pointer(unit->jvectors,"vg_set_vectors_pos_charge_atype");
    if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_charge_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge_atype **\n");
      vg_exit(1);
    }
    unit->jsize=size;
  }
  vg_set_scaleq(nj,charge,&scaleq,&scaleq_1);
  if(unit->scalers!=NULL){
    ((VG_SCALER *)(unit->scalers))->scaleqj_1=scaleq_1;
  }
  else{
    fprintf(stderr,"** error : unit->scalers=NULL **\n");
    vg_exit(1);
  }
  for(i=0;i<nj;i++){
    jvec=(VG_JVEC *)(unit->jvectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(jvec->r[j],pos[i][j],volume_1[j]);
    }
    jvec->qatype.atype=Q_DOUBLE_TO_INT(charge[i],scaleq);
    jvec->qatype.atype|=atype[i] & MASK(MD_QAUNION_ATYPEBIT);
  }
#elif defined(MD_NACL)
  int i,j,size;
  VG_JVEC *jvec;
  COPY_POS_DEFINE;

  size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
  unit->nj=nj;
  if(unit->jvectors==NULL){
    if(size<unit->jsize) size=unit->jsize;
    unit->jsize=size;
#ifdef MD_PRINT_WARN
    printf("allocating %d jvectors\n",size);
#endif
    if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_charge_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge_atype **\n");
      vg_exit(1);
    }
  }
  else if(size>unit->jsize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d jvectors\n",size);
#endif
    MR3_free_pointer(unit->jvectors,"vg_set_vectors_pos_charge");
    if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_charge_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge_atype **\n");
      vg_exit(1);
    }
  }
  for(i=0;i<nj;i++){
    jvec=(VG_JVEC *)(unit->jvectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(jvec->r[j],pos[i][j],volume_1[j]);
    }
    jvec->q=(float)(charge[i]);
#if MD_USE_LARGE_VAL==2
    jvec->q*=(float)(LOWER_VAL_FACTOR_ORG);
#endif
    jvec->atype=atype[i];
  }
#else
  fprintf(stderr,"** error : MD_QAUNION_ATYPEBIT must be defined for vg_set_vectors_pos_charge_atype **\n");
  vg_exit(1);
#endif
}


#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
void vg_set_vectors_pos_halfsigma_sqrtepsilon(VG_UNIT *unit, int nj,
					      double pos[][3],
					      double halfsigma[], 
					      double sqrtepsilon[])
{
  int i,j,size;
  VG_JVEC *jvec;
  COPY_POS_DEFINE;

  printf("** warning: free and malloc of jvectors should be modified **\n");
  size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
  unit->nj=nj;
  MR3_free_pointer(unit->jvectors,"vg_set_vectors_pos_halfsigma_sqrtepsilon");
  if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_halfsigma_sqrtepsilon"))==NULL){
    fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge_atype **\n");
    vg_exit(1);
  }
  for(i=0;i<nj;i++){
    jvec=(VG_JVEC *)(unit->jvectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(jvec->r[j],pos[i][j],volume_1[j]);
    }
    jvec->halfsigma=(float)(halfsigma[i]);
    jvec->sqrtepsilon=(float)(sqrtepsilon[i]);
  }
}
#endif


#ifdef MD_USE_QAUNION
void vg_set_vectors_pos_charge(VG_UNIT *unit, int nj,
			       double pos[][3], double charge[])
{
  int i,j,size;
  VG_JVEC *jvec;
#ifdef MD_QAUNION_ATYPEBIT
  double scaleq;
  float scaleq_1;
#endif
  COPY_POS_DEFINE;

  size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
  //size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj)*2;printf("size is doubled in vg_set_vectors_pos_charge\n");
  //  size=sizeof(VG_JVEC)*(NJ_ROUNDUP(nj)+VG_MINIMUM_PARTICLE_BLOCK_J);printf("size is modified in vg_set_vectors_pos_charge\n");

  unit->nj=nj;
  //  fprintf(stderr,"in vg_set_vectors_pos_charge: unit->jvectors=%x before free\n",unit->jvectors);
#if 0
  MR3_free_pointer(unit->jvectors,"vg_set_vectors_pos_charges");
  if((unit->jvectors=(void *)MR3_malloc_pointer(size,"unit->jvectors in vg_set_vectors_pos_charge"))==NULL){
    fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge **\n");
    vg_exit(1);
  }
#endif
#if 0
  if(unit->jvectors==NULL && (unit->jvectors=(void *)MR3_malloc_pointer(size,"unit->jvectors in vg_set_vectors_pos_charge"))==NULL){
    fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge **\n");
    vg_exit(1);
  }
  fprintf(stderr,"** warning : malloc of unit->jvectors is skipped in vg_set_vectors_pos_charge **\n");
#endif
#if 1
  if(unit->jvectors==NULL){
    if(size<unit->jsize) size=unit->jsize;
    unit->jsize=size;
#ifdef MD_PRINT_WARN
    printf("allocating %d jvectors\n",size);
#endif
    if((unit->jvectors=(void *)MR3_malloc_pointer(size,"void vg_set_vectors_pos_charge"))==NULL){
      fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge **\n");
      vg_exit(1);
    }
  }
  else if(size>unit->jsize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d jvectors\n",size);
#endif
    MR3_free_pointer(unit->jvectors,"vg_set_vectors_pos_charge");
    if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_charge"))==NULL){
      fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_charge **\n");
      vg_exit(1);
    }
    unit->jsize=size;
  }
#endif
  //  fprintf(stderr,"in vg_set_vectors_pos_charge: unit->jvectors=%x after malloc %d\n",unit->jvectors,size);
#ifdef MD_QAUNION_ATYPEBIT
  vg_set_scaleq(nj,charge,&scaleq,&scaleq_1);
  if(unit->scalers!=NULL){
    ((VG_SCALER *)(unit->scalers))->scaleqj_1=scaleq_1;
#if MD_USE_LARGE_VAL==2 || MD_USE_LARGE_VAL==20
    //    ((VG_SCALER *)(unit->scalers))->scaleqj_1*=LOWER_VAL_FACTOR_ORG;
#endif
  }
  else{
    fprintf(stderr,"** error : unit->scalers=NULL **\n");
    vg_exit(1);
  }
#endif
  for(i=0;i<nj;i++){
    jvec=(VG_JVEC *)(unit->jvectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(jvec->r[j],pos[i][j],volume_1[j]);
    }
#ifdef MD_QAUNION_ATYPEBIT
    jvec->qatype.atype=Q_DOUBLE_TO_INT(charge[i],scaleq) | \
      (MASK(MD_QAUNION_ATYPEBIT) & jvec->qatype.atype);
#else
    jvec->qatype.q=(float)(charge[i]);
#if MD_USE_LARGE_VAL==2
    jvec->qatype.q*=(float)(LOWER_VAL_FACTOR_ORG);
#endif
#endif
  }
}


void vg_set_vectors_pos_atype(VG_UNIT *unit, int nj,
			      double pos[][3], int atype[])
{
  int i,j,size;
  VG_JVEC *jvec;
  COPY_POS_DEFINE;

  size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
  //  size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj)*2;printf("size is doubled in vg_set_vectors_pos_charge\n");
  //  size=sizeof(VG_JVEC)*(NJ_ROUNDUP(nj)+VG_MINIMUM_PARTICLE_BLOCK_J);printf("size is modified in vg_set_vectors_pos_charge\n");

  unit->nj=nj;
  //  fprintf(stderr,"in vg_set_vectors_pos_atype: unit->jvectors=%x before free\n",unit->jvectors);
#if 0
  MR3_free_pointer(unit->jvectors,"vg_set_vectors_pos_atype");
  if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_atype"))==NULL){
    fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_atype **\n");
    vg_exit(1);
  }
#endif
#if 0
  if(unit->jvectors==NULL && (unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_atype"))==NULL){
    fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_atype **\n");
    vg_exit(1);
  }
  fprintf(stderr,"** warning : malloc of unit->jvectors is skipped in vg_set_vectors_pos_atype **\n");
#endif
#if 1
  if(unit->jvectors==NULL){
    if(size<unit->jsize) size=unit->jsize;
    unit->jsize=size;
#ifdef MD_PRINT_WARN
    printf("allocating %d jvectors\n",size);
#endif
    if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_atype **\n");
      vg_exit(1);
    }
  }
  else if(size>unit->jsize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d jvectors\n",size);
#endif
    MR3_free_pointer(unit->jvectors,"vg_set_vectors_pos_atype");
    if((unit->jvectors=(void *)MR3_malloc_pointer(size,"vg_set_vectors_pos_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc jvectors in vg_set_vectors_pos_atype **\n");
      vg_exit(1);
    }
    unit->jsize=size;
  }
#endif
  //  fprintf(stderr,"in vg_set_vectors_pos_atype: unit->jvectors=%x after malloc %d\n",unit->jvectors,size);
  for(i=0;i<nj;i++){
    jvec=(VG_JVEC *)(unit->jvectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(jvec->r[j],pos[i][j],volume_1[j]);
      //      printf("  j=%d pos[%d]=%e jvec=%08x\n",i,j,pos[i][j],jvec->r[j]);
    }
#ifdef MD_QAUNION_ATYPEBIT
    jvec->qatype.atype=atype[i] | \
      (~(MASK(MD_QAUNION_ATYPEBIT)) & jvec->qatype.atype);
#else
    jvec->qatype.atype=atype[i];
#endif
  }
}
#endif


void vg_set_matrices_gscales_rscales(VG_UNIT *unit, int nati_org, int natj_org, 
				     double gscales[], double rscales[])
{
  int i,j,size;
  VG_MATRIX *mat;
  int nati=nati_org,natj=natj_org;

#if VDW_SHIFT>=1
  double rcut2;
  if(unit->rcut2==0.0){
#if 0
    printf("** warning : unit->rcut2 is not set. using default %f **\n",
	   MD_LJ_R2MAX);
    unit->rcut2=MD_LJ_R2MAX;
#else
    fprintf(stderr,"** error : unit->rcut2 must be set **\n");
    vg_exit(1);
#endif
  }
  rcut2=unit->rcut2;
#endif
#if defined(MD_CELLINDEX) && 0
  nati++;
  natj++;
#endif
  if(nati*natj>VG_MINIMUM_ATYPE_BLOCK){
    fprintf(stderr,"** error : nati*natj>VG_MINIMUM_ATYPE_BLOCK(%d) **\n",
	    VG_MINIMUM_ATYPE_BLOCK);
    vg_exit(1);
  }
  size=sizeof(VG_MATRIX)*MATRIX_ROUNDUP(nati*natj);
  //  size=sizeof(VG_MATRIX)*MATRIX_ROUNDUP(nati*natj)*2;printf("size is doubled in vg_set_matrices_gscales_rscales\n");

  unit->nati=nati;
  unit->natj=natj;
  MR3_free_pointer(unit->matrices,"vg_set_matrices_gscales_rscales");
  if((unit->matrices=(void *)MR3_malloc_pointer(size,"vg_set_matrices_gscales_rscales"))==NULL){
    fprintf(stderr,"** error : can't malloc matrices in vg_set_vectors **\n");
    vg_exit(1);
  }
#if defined(MD_CELLINDEX) && 0
#if 0 // atype=nati_org or atype=natj_org is dummy
  for(i=0;i<nati;i++){
    for(j=0;j<natj;j++){
      mat=(VG_MATRIX *)(unit->matrices)+i*natj+j;
      if(i==nati_org || j==natj_org){
	mat->gscale=0.0f;
	mat->rscale=1.0f;
#if VDW_SHIFT>=1
	mat->shift=1.0f;
#endif
      }
      else{
	mat->gscale=(float)(gscales[i*natj_org+j]);
#if MD_USE_LARGE_VAL==2 || MD_USE_LARGE_VAL==20
	mat->gscale*=(float)(LOWER_VAL_FACTOR_ORG); // Narumi's method v2
#endif
	mat->rscale=(float)(rscales[i*natj_org+j]);
#if VDW_SHIFT==1
	mat->shift=(float)(pow(rscales[i*natj_org+j]*rcut2,-3.0));
#elif VDW_SHIFT==2
	mat->shift=(float)((unit->cuton)*(unit->cuton));
#endif
      }
    }
  }
#else // atype=0 is dummy
  for(i=0;i<nati;i++){
    for(j=0;j<natj;j++){
      mat=(VG_MATRIX *)(unit->matrices)+i*natj+j;
      if(i==0 || j==0){
	mat->gscale=0.0f;
	mat->rscale=1.0f;
#if VDW_SHIFT>=1
	mat->shift=1.0f;
#endif
      }
      else{
	mat->gscale=(float)(gscales[(i-1)*natj_org+j-1]);
#if MD_USE_LARGE_VAL==2 || MD_USE_LARGE_VAL==20
	mat->gscale*=(float)(LOWER_VAL_FACTOR_ORG); // Narumi's method v2
#endif
	mat->rscale=(float)(rscales[(i-1)*natj_org+j-1]);
#if VDW_SHIFT==1
	mat->shift=(float)(pow(rscales[(i-1)*natj_org+j-1]*rcut2,-3.0));
#elif VDW_SHIFT==2
	mat->shift=(float)((unit->cuton)*(unit->cuton));
#endif
      }
    }
  }
#endif
#else // else of MD_CELLINDEX
  for(i=0;i<nati_org;i++){
    for(j=0;j<natj_org;j++){
      mat=(VG_MATRIX *)(unit->matrices)+i*natj+j;
      mat->gscale=(float)(gscales[i*natj_org+j]);
#if MD_USE_LARGE_VAL==2 || MD_USE_LARGE_VAL==20
      mat->gscale*=(float)(LOWER_VAL_FACTOR_ORG); // Narumi's method v2
#endif
      mat->rscale=(float)(rscales[i*natj_org+j]);
#if VDW_SHIFT==1
      mat->shift=(float)(pow(rscales[i*natj_org+j]*rcut2,-3.0));
#elif VDW_SHIFT==2
      mat->shift=(float)((unit->cuton)*(unit->cuton));
#endif
    }
  }
#endif // end of MD_CELLINDEX
}


void vg_set_pol_sigm_ipotro_pc_pd_zz(VG_UNIT *unit, int nati, int natj,
				     double pol[], double sigm[], double ipotro[],
				     double pc[], double pd[], double zz[])
{
#ifdef MD_NACL
  int i,j,size;
  VG_MATRIX *mat;

  if(nati*natj>VG_MINIMUM_ATYPE_BLOCK){
    fprintf(stderr,"** error : nati*natj>VG_MINIMUM_ATYPE_BLOCK(%d) **\n",
	    VG_MINIMUM_ATYPE_BLOCK);
    vg_exit(1);
  }
  size=sizeof(VG_MATRIX)*MATRIX_ROUNDUP(nati*natj);
  unit->nati=nati;
  unit->natj=natj;
  MR3_free_pointer(unit->matrices,"vg_set_pol_sigm_ipotro_pc_pd_zz");
  if((unit->matrices=(void *)MR3_malloc_pointer(size,"vg_set_pol_sigm_ipotro_pc_pd_zz"))==NULL){
    fprintf(stderr,"** error : can't malloc matrices in vg_set_vectors **\n");
    vg_exit(1);
  }
  for(i=0;i<nati;i++){
    for(j=0;j<natj;j++){
      mat=(VG_MATRIX *)(unit->matrices)+i*natj+j;
      //      mat->gscale=(float)(gscales[i*natj+j]);
      //      mat->rscale=(float)(rscales[i*natj+j]);
      mat->pol=(float)(pol[i*natj+j]);
      mat->sigm=(float)(sigm[i*natj+j]);
      mat->ipotro=(float)(ipotro[i*natj+j]);
      mat->pc=(float)(pc[i*natj+j]);
      mat->pd=(float)(pd[i*natj+j]);
      mat->zz=(float)(zz[i*natj+j]);
    }
  }
#else
  fprintf(stderr,"**  error : MD_NACL should be defined in md.h **\n");
  vg_exit(1);
#endif
}


void vg_set_scalers_volume_alpha(VG_UNIT *unit, double volume[3], double alpha)
{
  int i;
  VG_SCALER *scal;

  if(unit->scalers==NULL){
    unit->ssize=sizeof(VG_SCALER);
    //    unit->ssize=sizeof(VG_SCALER)*2;printf("size is doubled in vg_set_scalers_volume_alpha\n");

    //    MR3_free_pointer(unit->scalers,"vg_set_scalers_volume_alpha");
    if((unit->scalers=(void *)MR3_malloc_pointer(unit->ssize,"vg_set_scalers_volume_alpha"))==NULL){
      fprintf(stderr,"** error : can't malloc scalers in vg_set_vectors **\n");
      vg_exit(1);
    }
  }
  scal=(VG_SCALER *)(unit->scalers);
  for(i=0;i<3;i++) scal->volume[i]=(float)(volume[i]);
  scal->alpha=(float)alpha;
  scal->alphafac=alpha*Coulomb_vdw_factor;
  scal->alpha3fac=alpha*alpha*alpha*Coulomb_vdw_factor;
#if defined(COULOMB_SHIFT) || 0
  scal->rcut21=(float)(1.0/unit->rcut2);
#endif
#if 1
  {
    static int ini=0;
    if(ini==0){
#ifdef MD_PRINT_WARN
      printf("eps2 is set zero in vg_set_scalers_volume_alpha in vtgrape.c **\n");
#endif
      ini=1;
    }
  }
  //  scal->eps2=0.0f;
#endif
}


#ifdef MD_MATRIX_IN_SCALER
void vg_set_scalers_volume_alpha_gscales_rscales(VG_UNIT *unit, double volume[3], 
						 double alpha, int nat,
						 double gscalesf[],
						 double gscalesp[], double rscales[])
{
  int i;
  VG_SCALER *scal;

  if(unit->scalers==NULL){
    unit->ssize=sizeof(VG_SCALER);
    //    MR3_free_pointer(unit->scalers,"vg_set_scalers_volume_alpha_gscales_rscales");
    if((unit->scalers=(void *)MR3_malloc_pointer(unit->ssize,"unit->scalers in vg_set_scalers_volume_alpha_gscales_rscales"))==NULL){
      fprintf(stderr,"** error : can't malloc scalers in vg_set_vectors **\n");
      vg_exit(1);
    }
  }
  scal=(VG_SCALER *)(unit->scalers);
  for(i=0;i<3;i++) scal->volume[i]=(float)(volume[i]);
  scal->alpha=(float)alpha;
  for(i=0;i<nat*nat;i++){
    scal->gscalesf[i]=(float)(gscalesf[i]);
    scal->gscalesp[i]=(float)(gscalesp[i]);
    scal->rscales[i]=(float)(rscales[i]);
  }
#if 1
  {
    static int ini=0;
    if(ini==0){
#ifdef MD_PRINT_WARN
      printf("eps2 is set zero in vg_set_scalers_volume_alpha in vtgrape.c **\n");
#endif
      ini=1;
    }
  }
  //  scal->eps2=0.0f;
#endif
}
#endif


void vg_set_pipeline_vectors_pos_charge_atype(VG_UNIT *unit, int ni, 
					      double pos[][3], double charge[],
					      int atype[])
{
#if defined(MD_USE_QAUNION) && defined(MD_QAUNION_ATYPEBIT)
  int i,j,size;
  VG_IVEC *ivec;
  double scaleq;
  float scaleq_1;
  COPY_POS_DEFINE;

  size=sizeof(VG_IVEC)*NI_ROUNDUP(ni);

  unit->ni=ni;
  //  printf("unit->ni=%d ni=%d in set_pipeline_vectors_pos_charge_atype\n",unit->ni,ni);
  if(unit->ivectors==NULL){
    if(size<unit->isize) size=unit->isize;
    unit->isize=size;
#ifdef MD_PRINT_WARN
    printf("allocating %d ivectors\n",size);
#endif
    if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_charge_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc ivectors in vg_set_pipeline_vectors_pos_charge_atype **\n");
      vg_exit(1);
    }
  }
  else if(size>unit->isize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d ivectors\n",size);
#endif
    MR3_free_pointer(unit->ivectors,"vg_set_pipeline_vectors_pos_charge_atype");
    if((unit->ivectors=(void *)MR3_malloc_pointer(size,"void vg_set_pipeline_vectors_pos_charge_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc ivectors in vg_set_pipeline_vectors_pos_charge_atype **\n");
      vg_exit(1);
    }
    unit->isize=size;
  }
  vg_set_scaleq(ni,charge,&scaleq,&scaleq_1);
  if(unit->scalers!=NULL){
    ((VG_SCALER *)(unit->scalers))->scaleqi_1=scaleq_1;
#if MD_USE_LARGE_VAL==2 || MD_USE_LARGE_VAL==20
    //    ((VG_SCALER *)(unit->scalers))->scaleqi_1*=LOWER_VAL_FACTOR_ORG;
#endif
  }
  else{
    fprintf(stderr,"** error : unit->scalers=NULL **\n");
    vg_exit(1);
  }
  for(i=0;i<ni;i++){
    ivec=(VG_IVEC *)(unit->ivectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(ivec->r[j],pos[i][j],volume_1[j]);
    }
    ivec->qatype.atype=Q_DOUBLE_TO_INT(charge[i],scaleq);
    ivec->qatype.atype|=atype[i] & MASK(MD_QAUNION_ATYPEBIT);
  }
#elif defined(MD_NACL)
  int i,j,size;
  VG_IVEC *ivec;
  COPY_POS_DEFINE;

  size=sizeof(VG_IVEC)*NI_ROUNDUP(ni);
  unit->ni=ni;

  if(unit->ivectors==NULL){
    if(size<unit->isize) size=unit->isize;
    unit->isize=size;
#ifdef MD_PRINT_WARN
    printf("allocating %d ivectors\n",size);
#endif
    if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_charge_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc ivectors in vg_set_vectors_pos_charge_atype **\n");
      vg_exit(1);
    }
  }
  else if(size>unit->isize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d ivectors\n",size);
#endif
    MR3_free_pointer(unit->ivectors,"vg_set_pipeline_vectors_pos_charge");
    if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_charge_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc ivectors in vg_set_vectors_pos_charge_atype **\n");
      vg_exit(1);
    }
  }
  for(i=0;i<ni;i++){
    ivec=(VG_IVEC *)(unit->ivectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(ivec->r[j],pos[i][j],volume_1[j]);
    }
    ivec->q=(float)(charge[i]);
    ivec->atype=atype[i];
  }
#else
  fprintf(stderr,"** error : MD_QAUNION_ATYPEBIT must be defined for vg_set_pipeline_vectors_pos_charge_atype **\n");
  vg_exit(1);
#endif
}


#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
void vg_set_pipeline_vectors_pos_halfsigma_sqrtepsilon(VG_UNIT *unit, int ni, 
						       double pos[][3], 
						       double halfsigma[],
						       double sqrtepsilon[])
{
  int i,j,size;
  VG_IVEC *ivec;
  COPY_POS_DEFINE;

  printf("** warning: free and malloc of ivectors should be modified **\n");
  size=sizeof(VG_IVEC)*NI_ROUNDUP(ni);
  unit->ni=ni;
  MR3_free_pointer(unit->ivectors,"vg_set_pipeline_vectors_pos_halfsigma_sqrtepsilon");
  if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_halfsigma_sqrtepsilon"))==NULL){
    fprintf(stderr,"** error : can't malloc ivectors in vg_set_vectors_pos_charge_atype **\n");
    vg_exit(1);
  }
  for(i=0;i<ni;i++){
    ivec=(VG_IVEC *)(unit->ivectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(ivec->r[j],pos[i][j],volume_1[j]);
    }
    ivec->halfsigma=(float)(halfsigma[i]);
    ivec->sqrtepsilon=(float)(sqrtepsilon[i]);
  }
}
#endif


#ifdef MD_USE_QAUNION
void vg_set_pipeline_vectors_pos_charge(VG_UNIT *unit, int ni, 
					double pos[][3], double charge[])
{
  int i,j,size;
  VG_IVEC *ivec;
#ifdef MD_QAUNION_ATYPEBIT
  double scaleq;
  float scaleq_1;
#endif
  COPY_POS_DEFINE;

  size=sizeof(VG_IVEC)*NI_ROUNDUP(ni);
  //  size=sizeof(VG_IVEC)*NI_ROUNDUP(ni)*2;printf("size is doubled in vg_set_pipeline_vectors_pos_charge\n");

  unit->ni=ni;
#if 0
  MR3_free_pointer(unit->ivectors,"vg_set_pipeline_vectors_pos_charge");
  if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_charge"))==NULL){
    fprintf(stderr,"** error : can't malloc ivectors in vg_set_vectors_pos_charge **\n");
    vg_exit(1);
  }
#else
  if(unit->ivectors==NULL){
    if(size<unit->isize) size=unit->isize;
    unit->isize=size;
#ifdef MD_PRINT_WARN
    printf("allocating %d ivectors\n",size);
#endif
    if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_charge"))==NULL){
      fprintf(stderr,"** error : can't malloc ivectors in vg_set_pipeline_vectors_pos_charge **\n");
      vg_exit(1);
    }
  }
  else if(size>unit->isize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d ivectors\n",size);
#endif
    MR3_free_pointer(unit->ivectors,"vg_set_pipeline_vectors_pos_charge");
    if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_charge"))==NULL){
      fprintf(stderr,"** error : can't malloc ivectors in vg_set_pipeline_vectors_pos_charge **\n");
      vg_exit(1);
    }
    unit->isize=size;
  }
#endif
#ifdef MD_QAUNION_ATYPEBIT
  vg_set_scaleq(ni,charge,&scaleq,&scaleq_1);
  if(unit->scalers!=NULL){
    ((VG_SCALER *)(unit->scalers))->scaleqi_1=scaleq_1;
  }
  else{
    fprintf(stderr,"** error : unit->scalers=NULL **\n");
    vg_exit(1);
  }
#endif
  for(i=0;i<ni;i++){
    ivec=(VG_IVEC *)(unit->ivectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(ivec->r[j],pos[i][j],volume_1[j]);
    }
#ifdef MD_QAUNION_ATYPEBIT
    ivec->qatype.atype=Q_DOUBLE_TO_INT(charge[i],scaleq) | \
      (MASK(MD_QAUNION_ATYPEBIT) & ivec->qatype.atype);
#else
    ivec->qatype.q=(float)(charge[i]);
#endif
  }
}


void vg_set_pipeline_vectors_pos_atype(VG_UNIT *unit, int ni, 
				       double pos[][3], int atype[])
{
  int i,j,size;
  VG_IVEC *ivec;
  COPY_POS_DEFINE;

  size=sizeof(VG_IVEC)*NI_ROUNDUP(ni);
  //  size=sizeof(VG_IVEC)*NI_ROUNDUP(ni);printf("size is doubled in vg_set_pipeline_vectors_pos_atype\n");

  unit->ni=ni;
#if 0
  MR3_free_pointer(unit->ivectors,"vg_set_pipeline_vectors_pos_atype");
  if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_atype"))==NULL){
    fprintf(stderr,"** error : can't malloc ivectors in vg_set_vectors_pos_atype **\n");
    vg_exit(1);
  }
#else
  if(unit->ivectors==NULL){
    if(size<unit->isize) size=unit->isize;
    unit->isize=size;
#ifdef MD_PRINT_WARN
    printf("allocating %d ivectors\n",size);
#endif
    if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc ivectors in vg_set_pipeline_vectors_pos_atype **\n");
      vg_exit(1);
    }
  }
  else if(size>unit->isize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d ivectors\n",size);
#endif
    MR3_free_pointer(unit->ivectors,"vg_set_pipeline_vectors_pos_atype");
    if((unit->ivectors=(void *)MR3_malloc_pointer(size,"vg_set_pipeline_vectors_pos_atype"))==NULL){
      fprintf(stderr,"** error : can't malloc ivectors in vg_set_pipeline_vectors_pos_atype **\n");
      vg_exit(1);
    }
    unit->isize=size;
  }
#endif
  for(i=0;i<ni;i++){
    ivec=(VG_IVEC *)(unit->ivectors)+i;
    for(j=0;j<3;j++){
      COPY_POS(ivec->r[j],pos[i][j],volume_1[j]);
    }
#ifdef MD_QAUNION_ATYPEBIT
    ivec->qatype.atype=atype[i] | \
      (~(MASK(MD_QAUNION_ATYPEBIT)) & ivec->qatype.atype);
#else
    ivec->qatype.atype=atype[i];
#endif
  }
}
#endif


void vg_convert_forces(VG_UNIT *unit, int nf, int function_index, double fi[][3])
{
  int i,j;
  VG_FVEC *fvec;

  for(i=0;i<nf;i++){
    fvec=(VG_FVEC *)(unit->fvectors)+i;
#ifdef ACCUMULATE_PAIR_FLOAT
#if ACCUMULATE_PAIR_FLOAT==3
    for(j=0;j<3;j++) fi[i][j]=((double *)(fvec->fi))[j];
#else // else of ACCUMULATE_PAIR_FLOAT==3
#if MD_USE_LARGE_VAL==1
    for(j=0;j<3;j++) fi[i][j]=(double)(fvec->fi[j])-(double)LARGE_VAL+(double)(fvec->fi2[j]);
#elif MD_USE_LARGE_VAL==2
    for(j=0;j<3;j++){
      fi[i][j]=(double)(fvec->fi[j])-(double)LARGE_VAL
	+((int *)(fvec->fi2))[j]*((double)LOWER_VAL_FACTOR_1);
      fi[i][j]*=(double)LOWER_VAL_FACTOR_1_ORG;// Narumi's method v2
    }
#elif MD_USE_LARGE_VAL==20
    if(function_index==2 || function_index==3
#ifdef MD_CELLINDEX
       || function_index==12 || function_index==13
#endif
#ifdef MD_QAUNION_ATYPEBIT
       || function_index==36 || function_index==37 
       || function_index==46 || function_index==47
#endif
       ){
      for(j=0;j<3;j++){
	fi[i][j]=(double)(fvec->fi[j])-(double)LARGE_VAL
	  +((int *)(fvec->fi2))[j]*((double)LOWER_VAL_FACTOR_1);
	fi[i][j]*=(double)LOWER_VAL_FACTOR_1_ORG;// Narumi's method v2
      }
    }
    else{
      for(j=0;j<3;j++) fi[i][j]=(double)(fvec->fi[j])+(double)(fvec->fi2[j]);
    }
#else
    for(j=0;j<3;j++) fi[i][j]=(double)(fvec->fi[j])+(double)(fvec->fi2[j]);
    //    if(i<3) printf("i=%d fi=%20.12e %20.12e %20.12e fi2=%20.12e %20.12e %20.12e\n",i,fvec->fi[0],fvec->fi[1],fvec->fi[2],fvec->fi2[0],fvec->fi2[1],fvec->fi2[2]);
#endif
#endif // end of ACCUMULATE_PAIR_FLOAT==3
#else
    for(j=0;j<3;j++) fi[i][j]=(double)(fvec->fi[j]);
#endif
  }
}


void vg_calculate_forces_potentials(VG_UNIT *unit, int function_index,
				    int nf, double fi[][3],
				    double *potc, double *potv)
{
  int fsize,psize;
  int i,j;
  VG_FVEC *fvec;
  VG_PSCALER *pscal;

  //  printf("in vg_calculate_forces_potentials : tblno=%d\n",function_index);
  fsize=sizeof(VG_FVEC)*NF_ROUNDUP(nf);
  //  fsize=sizeof(VG_FVEC)*NF_ROUNDUP(nf)*2;printf("size is doubled in vg_calculate_forces_potentials\n");
  unit->nf=nf;
#if 0
  //  printf("unit->fvectors=%x before free\n",unit->fvectors);
  MR3_free_pointer(unit->fvectors,"vg_calculate_forces_potentials");
  if((unit->fvectors=(void *)MR3_malloc_pointer(fsize,"vg_calculate_forces_potentials"))==NULL){
    fprintf(stderr,"** error : can't malloc fvectors in vg_calculate **\n");
    vg_exit(1);
  }
  //  printf("unit->fvectors=%x after malloc\n",unit->fvectors);
#else
  //  printf("unit->fvectors=%x in vg_calculate_potentials\n",unit->fvectors);
  if(unit->fvectors==NULL){
    //    printf("fsize=%d unit->fsize=%d nf=%d(*VG_FVEC=%d)\n",fsize,unit->fsize,nf,nf*sizeof(VG_FVEC));
    if(fsize<unit->fsize) fsize=unit->fsize;
    unit->fsize=fsize;
#ifdef MD_PRINT_WARN
    printf("allocating %d fvectors\n",fsize);
#endif
    if((unit->fvectors=(void *)MR3_malloc_pointer(fsize,"vg_calculate_forces_potentials"))==NULL){
      fprintf(stderr,"** error : can't malloc fvectors in vg_calculate_forces_potentials **\n");
      vg_exit(1);
    }
  }
  else if(fsize>unit->fsize){
#ifdef MD_PRINT_WARN
    printf("reallocating %d fvectors\n",fsize);
#endif
    MR3_free_pointer(unit->fvectors,"vg_calculate_forces_potentials");
    if((unit->fvectors=(void *)MR3_malloc_pointer(fsize,"vg_calculate_forces_potentials"))==NULL){
      fprintf(stderr,"** error : can't malloc fvectors in vg_calculate_forces_potentials **\n");
      vg_exit(1);
    }
    unit->fsize=fsize;
  }
#endif
  bzero(unit->fvectors,fsize);
  psize=sizeof(VG_PSCALER);
  //  psize=sizeof(VG_PSCALER)*2;printf("size is doubled in vg_calculate_forces_potentials\n");

#ifdef MD_CELLINDEX
  if(unit->pscalers==NULL){
    if((unit->pscalers=(VG_PSCALER *)MR3_malloc_pointer(sizeof(VG_PSCALER),"vg_calculate_forces_potentials"))==NULL){
      fprintf(stderr,"** error : can't malloc unit->pscalers in vg_calculate **\n");
      vg_exit(1);
    }
  }
#else
  MR3_free_pointer(unit->pscalers,"vg_calculate_forces_potentials");
  if((unit->pscalers=(void *)MR3_malloc_pointer(psize,"vg_calculate_forces_potentials"))==NULL){
    fprintf(stderr,"** error : can't malloc pscalers in vg_calculate **\n");
    vg_exit(1);
  }
  bzero(unit->pscalers,psize);
#endif
  //  if(unit->debug_flag==0) 
  //  printf("in vg_calculate_forces_potentials: tblno is added by 10\n",function_index+=10;
  //  printf("in vg_calculate_forces_potentials: function_index=%d\n",function_index); 
  unit->calculate[function_index]((void *)unit);
  //  else debug_mother((void *)unit);

#if defined(MD_PRINT_WARN) && 0
  printf("  unit->gpuoverlapflag=%d in vg_calculate_forces_potentials\n",unit->gpuoverlapflag);
#endif
  if(unit->gpuoverlapflag==0){
    vg_convert_forces(unit,nf,function_index,fi);
  }
  else{
    unit->function_index=function_index;
    unit->fthread=(double *)fi;
  }
  pscal=(VG_PSCALER *)(unit->pscalers);
  *potc=(double)(pscal->potc);
  *potv=(double)(pscal->potv);
}


static void vg_get_funcresult(int n_org, float r[], float r1[], int function_index)
{
#ifdef MD_USE_QAUNION
  int i,j,ni,nj,n,l,nii;
  VG_UNION_FI fi;
  double vol[3]={256.0,256.0,256.0},rscale=1.0;
  double potc,potv;
  double *xi;
  double *xj;
  double *force;
  double *qi;
  double *qj;
  VG_UNIT *unit=MR3_unit;

  if(unit==NULL){
    fprintf(stderr,"** error : unit is null **\n");
    exit(1);
  }
  if(n_org<0) n=-n_org;
  else        n=n_org;
  
  if((xi=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"void vg_get_funcresult"))==NULL){
    fprintf(stderr,"** error : can't malloc xi **\n");
    exit(1);
  }
  if((xj=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"void vg_get_funcresult"))==NULL){
    fprintf(stderr,"** error : can't malloc xj **\n");
    exit(1);
  }
  if((force=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"void vg_get_funcresult"))==NULL){
    fprintf(stderr,"** error : can't malloc force **\n");
    exit(1);
  }
  if((qi=(double *)MR3_malloc_pointer(sizeof(double)*n,"void vg_get_funcresult"))==NULL){
    fprintf(stderr,"** error : can't malloc qi **\n");
    exit(1);
  }
  if((qj=(double *)MR3_malloc_pointer(sizeof(double)*n,"void vg_get_funcresult"))==NULL){
    fprintf(stderr,"** error : can't malloc qj **\n");
    exit(1);
  }

  ni=n;
#if 1
  nj=1;
#else
  nj=VG_MINIMUM_PARTICLE_BLOCK_I;
  printf("**** nj is set %d\n",nj);
#endif

  if(n_org>0){
    for(i=0;i<n;i++){
      if(i<(1<<23)){
	fi.f=1.0f;
	fi.i=(fi.i & 0xff800000) | i;
	//    printf("i=%d fi.i=%08x\n",i,fi.i);
      }
      else{
	fi.f=2.0f;
	fi.i=(fi.i & 0xff800000) | (i-(1<<23));
      }
      r[i]=fi.f;
    }
  }

  for(i=0;i<ni;i++){
    xi[i*3+0]=r[i];
    xi[i*3+1]=xi[i*3+2]=0.0;
    qi[i]=1.0;
    force[i*3]=force[i*3+1]=force[i*3+2]=0.0;
  }
  for(i=0;i<VG_MINIMUM_PARTICLE_BLOCK_I;i++){
    xj[i*3]=xj[i*3+1]=xj[i*3+2]=0.0;
    qj[i]=0.0;
  }
  qj[0]=1.0;

  for(l=0;l<(ni+(1<<23)-1)/(1<<23);l++){
    int ii,niii;
    nii=l*(1<<23)+(1<<23)>ni ? ni-l*(1<<23):(1<<23);
    vg_set_vectors_pos_charge(unit,nj,(double (*)[3])xj,qj);
    vg_set_scalers_volume_alpha(unit,vol,rscale);
    for(ii=0;ii<nii;ii+=MD_MAX_I_PER_KERNEL){
      niii=ii+MD_MAX_I_PER_KERNEL>nii ? nii-ii:MD_MAX_I_PER_KERNEL;
      vg_set_pipeline_vectors_pos_charge(unit,niii,((double (*)[3])xi)+l*(1<<23)+ii,qi);
      vg_calculate_forces_potentials(unit,function_index,niii,(double (*)[3])force+ii,&potc,&potv);
    }
    for(i=0;i<nii;i++){
      r1[i+l*(1<<23)]=force[i*3];
    }
  }
  /*
  for(l=0;l<(ni+(1<<23)-1)/(1<<23);l++){
    nii=l*(1<<23)+(1<<23)>ni ? ni-l*(1<<23):(1<<23);
    vg_set_vectors_pos_charge(unit,nj,(double (*)[3])xj,qj);
    vg_set_scalers_volume_alpha(unit,vol,rscale);
    vg_set_pipeline_vectors_pos_charge(unit,nii,((double (*)[3])xi)+l*(1<<23),qi);
    vg_calculate_forces_potentials(unit,function_index,nii,(double (*)[3])force,&potc,&potv);
    for(i=0;i<nii;i++){
      r1[i+l*(1<<23)]=force[i*3];
    }
    }*/
  //  for(i=0;i<3;i++) printf("i=%d r1=%e\n",i,r1[i]);
  MR3_free_pointer(xi,"void vg_get_funcresult");
  MR3_free_pointer(xj,"void vg_get_funcresult");
  MR3_free_pointer(qi,"void vg_get_funcresult");
  MR3_free_pointer(qj,"void vg_get_funcresult");
  MR3_free_pointer(force,"void vg_get_funcresult");
#else
  fprintf(stderr,"** MD_USE_QAUNION must be defined in vg_get_r1result **\n");
#endif
}


void vg_get_r1result(int n_org, float r[], float r1[])
{
#if 1
  vg_get_funcresult(n_org,r,r1,20);
#else
  vg_get_funcresult(n_org,r,r1,0);
  printf("************** function_index is set 0\n");
#endif
}


void vg_get_rsqrtresult(int n_org, float r[], float r1[])
{
  vg_get_funcresult(n_org,r,r1,21);
  //  printf("************** function_index is set 20\n");
}


void vg_get_mulresult(int n_org, float r[], float q[], float rq[])
{
#ifdef MD_USE_QAUNION
  int function_index=22;
  int i,j,ni,nj,n,l,nii;
  VG_UNION_FI fi;
  double vol[3]={256.0,256.0,256.0},rscale=1.0;
  double potc,potv;
  double *xi;
  double *xj;
  double *force;
  double *qi;
  double *qj;
  VG_UNIT *unit=MR3_unit;

  if(n_org<0) n=-n_org;
  else        n=n_org;
  
  if((xi=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"vg_get_mulresult"))==NULL){
    fprintf(stderr,"** error : can't malloc xi **\n");
    exit(1);
  }
  if((xj=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"vg_get_mulresult"))==NULL){
    fprintf(stderr,"** error : can't malloc xj **\n");
    exit(1);
  }
  if((force=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"vg_get_mulresult"))==NULL){
    fprintf(stderr,"** error : can't malloc force **\n");
    exit(1);
  }
  if((qi=(double *)MR3_malloc_pointer(sizeof(double)*n,"vg_get_mulresult"))==NULL){
    fprintf(stderr,"** error : can't malloc qi **\n");
    exit(1);
  }
  if((qj=(double *)MR3_malloc_pointer(sizeof(double)*n,"vg_get_mulresult"))==NULL){
    fprintf(stderr,"** error : can't malloc qj **\n");
    exit(1);
  }

  ni=n;
#if 1
  nj=1;
#else
  nj=VG_MINIMUM_PARTICLE_BLOCK_I;
  printf("**** nj is set %d\n",nj);
#endif

  if(n_org>0){
    printf("n_org>0 is not supported\n");
    vg_exit(1);
  }

  for(i=0;i<ni;i++){
    xi[i*3+0]=r[i];
    xi[i*3+1]=xi[i*3+2]=0.0;
    qi[i]=1.0;
    force[i*3]=force[i*3+1]=force[i*3+2]=0.0;
  }
  for(i=0;i<nj;i++){
    xj[i*3]=xj[i*3+1]=xj[i*3+2]=0.0;
    qj[i]=q[i];
  }
  //  for(i=0;i<VG_MINIMUM_PARTICLE_BLOCK_I;i++){
  //    xj[i*3]=xj[i*3+1]=xj[i*3+2]=0.0;
  //    qj[i]=q[i];
  //  }
  //  qj[0]=1.0;

  for(l=0;l<(ni+(1<<23)-1)/(1<<23);l++){
    int ii,niii;
    nii=l*(1<<23)+(1<<23)>ni ? ni-l*(1<<23):(1<<23);
    vg_set_vectors_pos_charge(unit,nj,(double (*)[3])xj,qj);
    vg_set_scalers_volume_alpha(unit,vol,rscale);
    for(ii=0;ii<nii;ii+=MD_MAX_I_PER_KERNEL){
      niii=ii+MD_MAX_I_PER_KERNEL>nii ? nii-ii:MD_MAX_I_PER_KERNEL;
      vg_set_pipeline_vectors_pos_charge(unit,niii,((double (*)[3])xi)+l*(1<<23)+ii,qi);
      vg_calculate_forces_potentials(unit,function_index,niii,(double (*)[3])force+ii,&potc,&potv);
    }
    for(i=0;i<nii;i++){
      rq[i+l*(1<<23)]=force[i*3];
    }
  }
  /*
  for(l=0;l<(ni+(1<<23)-1)/(1<<23);l++){
    nii=l*(1<<23)+(1<<23)>ni ? ni-l*(1<<23):(1<<23);
    vg_set_vectors_pos_charge(unit,nj,(double (*)[3])xj,qj);
    vg_set_scalers_volume_alpha(unit,vol,rscale);
    vg_set_pipeline_vectors_pos_charge(unit,nii,((double (*)[3])xi)+l*(1<<23),qi);
    vg_calculate_forces_potentials(unit,function_index,nii,(double (*)[3])force,&potc,&potv);
    for(i=0;i<nii;i++){
      r1[i+l*(1<<23)]=force[i*3];
    }
    }*/
  //  for(i=0;i<3;i++) printf("i=%d r1=%e\n",i,r1[i]);
  MR3_free_pointer(xi,"vg_get_mulresult");
  MR3_free_pointer(xj,"vg_get_mulresult");
  MR3_free_pointer(qi,"vg_get_mulresult");
  MR3_free_pointer(qj,"vg_get_mulresult");
  MR3_free_pointer(force,"vg_get_mulresult");
#else
  fprintf(stderr,"** MD_USE_QAUNION must be defined in vg_get_r1result **\n");
#endif
}


int vg_write_emufile(int flag)
{
  /* 
     flag 0 -- write emu file
          1 -- do not write emu file, but compare MD5SUM
   */
  float *r;
  int i;
  VG_UNIT *unit=m3_get_unit();
  FILE *fp;
  char fname[1024];

  if(unit==NULL){
    fprintf(stderr,"** error : unit is NULL **\n");
    return 1;
  }
#ifdef MD_GENERATE_EMUTABLE
  if(Emutable_r1!=NULL && Emutable_rsqrt!=NULL){
    if(unit->r1==NULL){
      if((unit->r1=(float *)MR3_malloc_pointer(sizeof(float)*(1<<23),"vg_write_emufile"))==NULL){
	fprintf(stderr,"** error : can't malloc r1 **\n");
	vg_exit(1);
      }
    }
    memcpy(unit->r1,Emutable_r1,sizeof(float)*(1<<23));
    if(unit->rsqrt==NULL){
      if((unit->rsqrt=(float *)MR3_malloc_pointer(sizeof(float)*(1<<24),"vg_write_emufile"))==NULL){
	fprintf(stderr,"** error : can't malloc rsqrt **\n");
	vg_exit(1);
      }
    }
    memcpy(unit->rsqrt,Emutable_rsqrt,sizeof(float)*(1<<24));
  }
  else if(Emutable_r1==NULL && Emutable_rsqrt==NULL){
#endif
  if((r=(float *)MR3_malloc_pointer(sizeof(float)*(1<<24),"vg_write_emufile"))==NULL){
    fprintf(stderr,"** error : can't malloc r **\n");
    vg_exit(1);
  }
  if(unit->r1==NULL){
    if((unit->r1=(float *)MR3_malloc_pointer(sizeof(float)*(1<<23),"vg_write_emufile"))==NULL){
      fprintf(stderr,"** error : can't malloc r1 **\n");
      vg_exit(1);
    }
  }
  if(unit->rsqrt==NULL){
    if((unit->rsqrt=(float *)MR3_malloc_pointer(sizeof(float)*(1<<24),"vg_write_emufile"))==NULL){
      fprintf(stderr,"** error : can't malloc rsqrt **\n");
      vg_exit(1);
    }
  }
  vg_get_r1result((1<<23),r,unit->r1);
#if 1
  vg_get_rsqrtresult((1<<24),r,unit->rsqrt);
#else
  printf("** vg_set_rsqrtemu is skipped in MR3init_emu\n");
#endif
  MR3_free_pointer(r,"vg_write_emufile");
#ifdef MD_GENERATE_EMUTABLE
  if((Emutable_r1=(float *)MR3_malloc_pointer(sizeof(float)*(1<<23),"vg_write_emufile"))==NULL){
    fprintf(stderr,"** error : can't malloc Emutable_r1 **\n");
    vg_exit(1);
  }
  if((Emutable_rsqrt=(float *)MR3_malloc_pointer(sizeof(float)*(1<<24),"vg_write_emufile"))==NULL){
    fprintf(stderr,"** error : can't malloc Emutable_rsqrt **\n");
    vg_exit(1);
  }
  memcpy(Emutable_r1,unit->r1,sizeof(float)*(1<<23));
  memcpy(Emutable_rsqrt,unit->rsqrt,sizeof(float)*(1<<24));
  }
  else{
    fprintf(stderr,"** error : Emutable_r1 or Emutable_sqrt is strange **\n");
    vg_exit(1);
  }
#endif

  if(flag==0){ // write to file
    sprintf(fname,"%s",MD_R1_EMU_FNAME);
    if((fp=fopen(fname,"w"))==NULL){
      fprintf(stderr,"** error : can't open %s **\n",fname);
      vg_exit(1);
    }
    printf("writing %s ... ",fname);
    fwrite(unit->r1,sizeof(float),(1<<23),fp);
    fclose(fp);
    printf("finished.\n");
    sprintf(fname,"%s",MD_RSQRT_EMU_FNAME);
    if((fp=fopen(fname,"w"))==NULL){
      fprintf(stderr,"** error : can't open %s **\n",fname);
      vg_exit(1);
    }
    printf("writing %s ... ",fname);
    fwrite(unit->rsqrt,sizeof(float),(1<<24),fp);
    fclose(fp);
    printf("finished.\n");
  }
  else{ // calc md5sum and compare
    unsigned char sum[16];
    char s[1024];
    md5sum((unsigned char *)(unit->r1),sizeof(float)*(1<<23),sum);
    for(i=0;i<16;i++){
      //      printf("%02x",sum[i]);
      sprintf(s+i*2,"%02x",sum[i]);
    }
#ifdef MD_PRINT_WARN      
    //    printf("\n");
    printf("%s\n",s);
#endif
    if(strcmp(s,MD_R1_EMU_MD5SUM)!=0){
      fprintf(stderr,"** error : generating r1 table for GPU emulator failed **\n");
      //      printf("** ignore error **\n");return 0;
      return 1;
    }

    md5sum((unsigned char *)(unit->rsqrt),sizeof(float)*(1<<24),sum);
    for(i=0;i<16;i++){
      //      printf("%02x",sum[i]);
      sprintf(s+i*2,"%02x",sum[i]);
    }
#ifdef MD_PRINT_WARN      
    //    printf("\n");
    printf("%s\n",s);
#endif
    if(strcmp(s,MD_RSQRT_EMU_MD5SUM)!=0){
      fprintf(stderr,"** error : generating rsqrt table for GPU emulator failed **\n");
      //      printf("** ignore error **\n");return 0;
      return 1;
    }
  }
  return 0;
}


void vg_initialize_emu(void)
{
#if 0
  printf("MR3init_emu is skipped\n");
  return;
#endif
  FILE *fp;
  char fname[1024],emudir[1024],*emuvar;

  VG_UNIT *unit=m3_get_unit();
  if(unit->r1==NULL){
    if((unit->r1=(float *)MR3_malloc_pointer(sizeof(float)*(1<<23),"vg_initialize_emu"))==NULL){
      fprintf(stderr,"** error : can't malloc r1 **\n");
      vg_exit(1);
    }
  }
  if(unit->rsqrt==NULL){
    if((unit->rsqrt=(float *)MR3_malloc_pointer(sizeof(float)*(1<<24),"vg_initialize_emu"))==NULL){
      fprintf(stderr,"** error : can't malloc rsqrt **\n");
      vg_exit(1);
    }
  }
  emuvar=getenv("VG_EMUDIR");
#ifdef MD_GENERATE_EMUTABLE
  if(emuvar!=NULL){
    //  printf("vg_write_emufile(1) is called\n");vg_write_emufile(1);
#endif
  if(emuvar==NULL) sprintf(emudir,".");
  else             sprintf(emudir,"%s",emuvar);
  sprintf(fname,"%s/%s",emudir,MD_R1_EMU_FNAME);
#ifdef MD_GENERATE_EMUTABLE
  printf("VG_EMUDIR is set %s, ",emuvar);
#endif
#if defined(MD_GENERATE_EMUTABLE) || defined(MD_PRINT_WARN)
  printf("reading emu table file %s ... ",fname);
#endif
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"** error : can't open %s **\n",fname);
    vg_exit(1);
  }
  fread(unit->r1,sizeof(float),(1<<23),fp);
  fclose(fp);
#if defined(MD_GENERATE_EMUTABLE) || defined(MD_PRINT_WARN)
  printf("finished.\n");
#endif
  sprintf(fname,"%s/%s",emudir,MD_RSQRT_EMU_FNAME);
#ifdef MD_GENERATE_EMUTABLE
  printf("VG_EMUDIR is set %s, ",emuvar);
#endif
#if defined(MD_GENERATE_EMUTABLE) || defined(MD_PRINT_WARN)
  printf("reading emu table file %s ... ",fname);
#endif
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"** error : can't open %s **\n",fname);
    vg_exit(1);
  }
  fread(unit->rsqrt,sizeof(float),(1<<24),fp);
  fclose(fp);
#if defined(MD_GENERATE_EMUTABLE) || defined(MD_PRINT_WARN)
  printf("finished.\n");
#endif
#ifdef MD_GENERATE_EMUTABLE
  }
  else{ // generate emu file on the fly
#if 1
    if(vg_write_emufile(1)!=0){
      fprintf(stderr,"** error : generating emu table failed **\n");
      exit(1); 
    }
#endif
  }
#endif
}


float vg_r1emu(float x)
{
  int index;
  float exp,ret;
  VG_UNION_FI fi;
  VG_UNIT *unit=MR3_unit;

  if(unit->r1==NULL) vg_initialize_emu();
  fi.f=x;
  index=(fi.i & 0x7fffff);
  fi.i&=0xff800000;
  exp=fi.f;
  ret=unit->r1[index]/exp;

  return ret;
}


#ifdef VG_GCC
v4sf vg_r1emu_vec(v4sf xv)
{
  v4sf retv;
  int i;
  
  for(i=0;i<4;i++){
    ((float *)&retv)[i]=vg_r1emu(((float *)&xv)[i]);
  }

  return retv;
}
#endif


float vg_rsqrtemu(float x)
{
  int index;
  float exp,ret;
  VG_UNION_FI fi;
  VG_UNIT *unit=MR3_unit;

  if(unit->rsqrt==NULL) vg_initialize_emu();
  fi.f=x;
  index=(fi.i & 0x7fffff);
  fi.i&=0xff800000;
  exp=fi.f;
  if((*((unsigned int *)&exp)>>23) % 2==0){
    exp*=0.5f;
    index+=(1<<23);
  }
  ret=unit->rsqrt[index]/sqrtf(exp);

  return ret;
}


//#define MR3_malloc_pointer(size,string) malloc(size)
//#define MR3_free_pointer(p,string)      free(p)


#ifdef MD_QAUNION_ATYPEBIT
void vg_set_scaleq(int n, double q[], 
		   double *scaleq, float *scaleq_1)
{
  int i;
  double maxq=0.0;

  for(i=0;i<n;i++){
    if(fabs(q[i])>maxq){
      maxq=fabs(q[i]);
    }
  }
  *scaleq=0x7fffff00/maxq;
  *scaleq_1=1.0/(*scaleq);
  //  printf("maxq=%e scaleq=%e scaleq_1=%e\n",maxq,*scaleq,*scaleq_1);
}
#endif


void m3_set_softening(VG_UNIT *unit, double eps)
{
  VG_SCALER *scal;
#if 0
  static int ini=0;
  if(ini==0){
    printf("*** MR3_unit is used instead of unit in m3_set_softening in vtgrape.c **\n");
    ini=1;
  }
  scal=(VG_SCALER *)(MR3_unit->scalers);
#else
  scal=(VG_SCALER *)(unit->scalers);
#endif
  if(scal==NULL){
    unit->ssize=sizeof(VG_SCALER);
    if((unit->scalers=(void *)MR3_malloc_pointer(unit->ssize,"m3_set_softening"))==NULL){
      fprintf(stderr,"** error : can't malloc scalers in m3_set_softening **\n");
      vg_exit(1);
    }
    scal=unit->scalers;
  }
  scal->eps2=(float)(eps*eps);
}


void MR3_exit(int ret)
{
  vg_exit(ret);
}


void MR3init(void)
{
  VG_UNIT *unit;

#ifdef MD_PRINT_WARN
  printf("** MD_PRINT_WARN is defined. Undef if you want **\n");
#endif
#ifdef MD_MEASURE_TIME
  {
    int i;
    for(i=0;i<MD_MEASURE_TIME;i++){
      Time[i][0]=Time[i][1]=Time[i][2]=0.0;
    }
  }
#endif
#if 0 // print the size of struct
  printf("in MR3init : size of VG_UNIT=%d, M3_CELL=%d, VG_JVEC=%d VG_IVEC=%d, VG_MATRIX=%d, VG_SCALER=%d, VG_FVEC=%d, VG_PSCALER=%d, VG_UNION_FI=%d\n",
	 sizeof(VG_UNIT),sizeof(M3_CELL),sizeof(VG_JVEC),sizeof(VG_IVEC),
	 sizeof(VG_MATRIX),sizeof(VG_SCALER),sizeof(VG_FVEC),
	 sizeof(VG_PSCALER),sizeof(VG_UNION_FI));
  debug_mother(unit);
#endif  

#if VG_MINIMUM_PARTICLE_BLOCK_I==256
  printf("** warning: VG_MINIMUM_PARTICLE_BLOCK_I=256 sometimes slow **\n");
#endif
  if(MD_LJ_R2MIN!=0.25f){
#ifdef MD_PRINT_WARN
    printf("** warning: MD_LJ_R2MIN is set %f **\n",MD_LJ_R2MIN);
#endif
  }
  if(MD_LJ_R2MAX!=64.0f){
#ifdef MD_PRINT_WARN
    printf("** warning: MD_LJ_R2MAX is set %f **\n",MD_LJ_R2MAX);
#endif
  }
  while((unit=vg_allocate_unit(grav_mother))==NULL){
    sleep(1);
    printf(".");
  }
  MR3_unit=unit;
#ifdef MD_GENERATE_EMUTABLE
  vg_set_function(unit,20,r1_mother);
  vg_set_function(unit,21,rsqrt_mother);
  vg_set_function(unit,22,mul_mother);
  vg_initialize_emu();
#ifdef MD_PRINT_WARN
  printf("** warning: vg_initialize_emu is called in MR3init **\n");
#endif
  vg_free_unit(MR3_unit);
  while((unit=vg_allocate_unit(grav_mother))==NULL){
    sleep(1);
    printf(".");
  }
  MR3_unit=unit;
#endif
#ifdef MD_PRINT_WARN
  printf("MR3init succeeded in vtgrape.c\n");
#endif

#if 1 // GPU

#if 1 // free boundary optimized routine is used
  vg_set_function(unit,0,gravfb_mother);
  vg_set_function(unit,1,gravpotfb_mother);
#if defined(MD_CELLINDEX) && 0
  vg_set_function(unit,2,ljfbci_mother);
  vg_set_function(unit,3,ljpotfbci_mother);
#else
  vg_set_function(unit,2,ljfb_mother);
  vg_set_function(unit,3,ljpotfb_mother);
  //  vg_set_function(unit,3,ljpot_mother);
#endif
  //  vg_set_function(unit,6,realfb_mother);
  //  vg_set_function(unit,7,realpotfb_mother);
  vg_set_function(unit,6,real_mother);
  vg_set_function(unit,7,realpot_mother);
#else
  printf("** non optimized routine is used **\n");
  vg_set_function(unit,0,grav_mother);
  vg_set_function(unit,1,gravpot_mother);
  vg_set_function(unit,2,lj_mother);
  vg_set_function(unit,3,ljpot_mother);
  vg_set_function(unit,6,real_mother);
  vg_set_function(unit,7,realpot_mother);
#endif
  vg_set_function(unit,4,wave_mother);
  vg_set_function(unit,5,wavepot_mother);
#ifdef MD_CELLINDEX
  //  vg_set_function(unit,12,ljfb_mother);
  //  vg_set_function(unit,13,ljpotfb_mother);
  vg_set_function(unit,12,ljfbci_mother);
  vg_set_function(unit,13,ljpotfbci_mother);
  vg_set_function(unit,16,realci_mother);
  vg_set_function(unit,17,realpotci_mother);
  //  vg_set_function(unit,16,realfbci_mother);
  //  vg_set_function(unit,17,realpotfbci_mother);
#endif
  vg_set_function(unit,20,r1_mother);
  vg_set_function(unit,21,rsqrt_mother);
  vg_set_function(unit,22,mul_mother);
#if defined(MD_QAUNION_ATYPEBIT) && 1
  vg_set_function(unit,36,realljfbci_mother);
  vg_set_function(unit,37,realljpotfbci_nojdiv_mother); // potential of real and lj
                                                 // real pot is force[i*3], lj pot is force[i*3+1]
  vg_set_function(unit,46,realljfbci2_mother);
  vg_set_function(unit,47,realljpotfbci_mother); // potential of real and lj
                                                 // real pot is force[i*3], lj pot is force[i*3+1]
#endif
  /*
  vg_set_function(unit,100,grav_fb_eps_mother);
  vg_set_function(unit,110,grav_fb_isrzero_mother);
  vg_set_function(unit,120,grav_periodic_nrint_mother);
  vg_set_function(unit,130,grav_periodic_fixed_mother);
  vg_set_function(unit,140,grav_periodic_fixed_pairfloat_simple_mother);
  vg_set_function(unit,150,grav_periodic_fixed_pairfloat_mid_mother);
  vg_set_function(unit,160,grav_periodic_fixed_pairfloat_full_mother);
  */
  /*
  vg_set_function(unit,102,lj_fb_iszero_mother);
  vg_set_function(unit,112,lj_fb_atomtype_constant_mother);
  vg_set_function(unit,122,lj_fb_atomtype_shared_mother);
  vg_set_function(unit,132,lj_fb_halfsigmasqrtepsilon_mother);
  vg_set_function(unit,142,lj_periodic_nrint_atomtype_constant_mother);
  vg_set_function(unit,152,lj_periodic_fixed_atomtype_constant_mother);
  vg_set_function(unit,162,lj_periodic_fixed_atomtype_constant_pairfloat_simple_mother);
  vg_set_function(unit,172,lj_periodic_fixed_atomtype_constant_pairfloat_mid_mother);
  vg_set_function(unit,182,lj_periodic_fixed_atomtype_constant_pairfloat_full_mother);
  */
#ifdef MD_NACL
  printf("** MD_NACL is defined\n");
  vg_set_function(unit,8,nacl_mother);
#endif
#else // Host
  vg_set_function(unit,0,md_grav);
  vg_set_function(unit,1,md_gravpot);
  vg_set_function(unit,2,md_lj);
  vg_set_function(unit,3,md_ljpot);
  vg_set_function(unit,4,md_wave);
  vg_set_function(unit,5,md_wavepot);
  vg_set_function(unit,6,md_real);
  vg_set_function(unit,7,md_realpot);
#endif

#if MD_LJ_05SIGMA_8SIGMA==2
  vg_set_rcut2(unit,MD_LJ_R2MAX);
#elif MD_LJ_05SIGMA_8SIGMA==0
  printf("** LJ range is not limited to [0.5 sigma, 8 sigma) or some cutoff\n");
#endif
#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
  printf("** MD_SIGMA_EPSION_IN_VGSTRUCT is defined\n");
#endif
#if VDW_SHIFT>=1 && MD_LJ_05SIGMA_8SIGMA!=2
  fprintf(stderr,"** error : MD_LJ_05SIGMA_8SIGMA must be 2 when VDW_SHIFT is defined **\n");
  vg_exit(1);
#endif
#if MD_LJ_05SIGMA_8SIGMA==2 && (defined(MD_SIGMA_EPSION_IN_VGSTRUCT) \
    || defined(MD_MATRIX_IN_SCALER) || defined(MD_MATRIX_COPY_TO_SHARED))
  fprintf(stderr,"** error : illeagal combination for MD_LJ_05SIGMA_8SIGMA=2 **\n");
  vg_exit(1);
#endif

#ifdef MD_CELLINDEX
  {
    char *s;
    if((s=getenv("VG_LDIMX"))!=NULL){
      sscanf(s,"%d",&Ldim[0]);
      printf("VG_LDIMX is set %d\n",Ldim[0]);
    }
    if((s=getenv("VG_LDIMY"))!=NULL){
      sscanf(s,"%d",&Ldim[1]);
      printf("VG_LDIMY is set %d\n",Ldim[1]);
    }
    if((s=getenv("VG_LDIMZ"))!=NULL){
      sscanf(s,"%d",&Ldim[2]);
      printf("VG_LDIMZ is set %d\n",Ldim[2]);
    }
  }
#endif
#if defined(MD_GENERATE_EMUTABLE) && 0
    if(vg_write_emufile(1)!=0){
      fprintf(stderr,"** error : generating emu table failed **\n");
      exit(1); 
    }
#endif
}


void mr3init_(void)
{
  MR3init();
}


void mr3init__(void)
{
  mr3init_();
}


void MR3free(void)
{
  vg_free_unit(MR3_unit);
#ifdef MD_MEASURE_TIME
  {
    int i,j,k;
    printf("* Timing report  ");
    for(j=0;j<10;j++) printf("%8d ",j);
    printf("*\n");
    for(i=0;i<MD_MEASURE_TIME;i+=10){
      k=0;
      for(j=0;j<10;j++) if(Time[i+j][2]!=0.0) k=1;
      if(k==1){
	printf("         Time[%3d] ",i);
	for(j=0;j<10;j++) printf("%8.3f ",Time[i+j][2]);
	printf("\n");
      }
    }
  }
#endif
}


void mr3free_(void)
{
  MR3free();
}


void mr3free__(void)
{
  mr3free_();
}


void MR3calccoulomb(double x[], int n, double q[], double rscale,
		    int tblno, double xmax, int periodicflag,
		    int natchangeflag, double force[])
{
  int i,j,*atype;
  VG_UNIT *unit=MR3_unit;
  double vol[3],potc,potv,(*ftmp)[3],alpha3;

#ifndef MD_USE_QAUNION
  if((atype=(int *)MR3_malloc_pointer(sizeof(int)*n,"MR3calccoulomb"))==NULL){
    fprintf(stderr,"** error : can't malloc atype in MR3calccoulomb **\n");
    vg_exit(1);
  }
  for(i=0;i<n;i++) atype[i]=0;
#endif
  if((ftmp=(double (*)[3])MR3_malloc_pointer(sizeof(double)*n*3,"MR3calccoulomb"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp in MR3calccoulomb **\n");
    vg_exit(1);
  }
  if((periodicflag & 1)==0) xmax*=2;
  for(i=0;i<3;i++) vol[i]=xmax;
#if MD_PERIODIC_FIXED==1
  //  unit->xmax=xmax;
  for(i=0;i<3;i++) unit->volume[i]=xmax;
#endif

#ifdef MD_USE_QAUNION
  vg_set_vectors_pos_charge(unit,n,(double (*)[3])x,q);
#else
  vg_set_vectors_pos_charge_atype(unit,n,(double (*)[3])x,q,atype);
#endif
  vg_set_scalers_volume_alpha(unit,vol,rscale);
#ifdef MD_USE_QAUNION
  vg_set_pipeline_vectors_pos_charge(unit,n,(double (*)[3])x,q);
#else
  vg_set_pipeline_vectors_pos_charge_atype(unit,n,(double (*)[3])x,q,atype);
#endif
  vg_calculate_forces_potentials(unit,tblno,n,ftmp,&potc,&potv);
  switch(tblno){
  case 0:
  case 1:
    if((periodicflag & 2)==0){
      for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j];
    }
    else{
      for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j]*q[i];
    }
    break;
  case 6:
    alpha3=rscale*rscale*rscale;
    for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j]*q[i]*alpha3;
    break;
  case 7:
    for(i=0;i<n;i++) force[i*3]+=ftmp[i][0]*q[i]*rscale;
    break;
  default:
    fprintf(stderr,"** error : not supported mode=%d **\n",tblno);
    vg_exit(1);
    break;
  }

#ifndef MD_USE_QAUNION
  MR3_free_pointer(atypMR3inite);
#endif
  MR3_free_pointer(ftmp,"MR3calccoulomb");
}


void mr3calccoulomb_(double x[], int *n, double q[], double *rscale,
		     int *tblno, double *xmax, int *periodicflag,
		     int *natchangeflag, double force[])
{
  MR3calccoulomb(x,*n,q,*rscale,*tblno,*xmax,*periodicflag,
		 *natchangeflag,force);
}


void mr3calccoulomb__(double x[], int *n, double q[], double *rscale,
		      int *tblno, double *xmax, int *periodicflag,
		      int *natchangeflag, double force[])
{
  mr3calccoulomb_(x,n,q,rscale,
		  tblno,xmax,periodicflag,
		  natchangeflag,force);
}


void MR3calccoulomb_ij(int ni, double xi[], double qi[], double force[],
		       int nj, double xj[], double qj[],
		       double rscale, 
		       int tblno, double xmax, int periodicflag)
{
  /*
      periodicflag bit 0 ---- 0 : non periodic
                              1 : periodic
                   bit 1 ---- 0 : do not multiply charge nor rscale
                              1 : multiply charge and rscale
   */

  int i,j,*atypei,*atypej;
  VG_UNIT *unit=MR3_unit;
  double vol[3],potc,potv,(*ftmp)[3],alpha3;
  int multiplyq=0;

#if 0
  printf("in MR3calccoulomb_ij: Coulomb_vdw_factor=%e tblno=%d\n",Coulomb_vdw_factor,tblno);
  for(i=0;i<3;i++) printf("  before force_host[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
  MR3calccoulomb_ij_host(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,periodicflag);
  for(i=0;i<3;i++) printf("  after force_host[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
#endif
#ifndef MD_USE_QAUNION
  if((atypei=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calccoulomb_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc atypei in MR3calccoulomb_ij **\n");
    vg_exit(1);
  }
  if((atypej=(int *)MR3_malloc_pointer(sizeof(int)*nj,"MR3calccoulomb_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc atypej in MR3calccoulomb_ij **\n");
    vg_exit(1);
  }
  for(i=0;i<ni;i++) atypei[i]=0;
  for(i=0;i<nj;i++) atypej[i]=0;
#endif
  if((ftmp=(double (*)[3])MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp in MR3calccoulomb_ij **\n");
    vg_exit(1);
  }
  if((periodicflag & 2)!=0) multiplyq=1;
  if((periodicflag & 1)==0) xmax*=2;
  for(i=0;i<3;i++) vol[i]=xmax;
#if MD_PERIODIC_FIXED==1
  //  unit->xmax=xmax;
  for(i=0;i<3;i++) unit->volume[i]=xmax;
#endif

#ifdef MD_USE_QAUNION
  vg_set_vectors_pos_charge(unit,nj,(double (*)[3])xj,qj);
#else
  vg_set_vectors_pos_charge_atype(unit,nj,(double (*)[3])xj,qj,atypej);
#endif
  vg_set_scalers_volume_alpha(unit,vol,rscale);
#ifdef MD_USE_QAUNION
  vg_set_pipeline_vectors_pos_charge(unit,ni,(double (*)[3])xi,qi);
#else
  vg_set_pipeline_vectors_pos_charge_atype(unit,ni,(double (*)[3])xi,qi,atypei);
#endif
  //  if(unit->debug_flag==0) 
  vg_calculate_forces_potentials(unit,tblno,ni,ftmp,&potc,&potv);
  //  else debug_mother((void *)unit);
  if(tblno>=100) tblno=(tblno-100) % 10;
  switch(tblno){
  case 0:
  case 1:
    if(multiplyq){
      for(i=0;i<ni;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j]*qi[i]*Coulomb_vdw_factor;
    }
    else{
      for(i=0;i<ni;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j];
    }
    break;
  case 6:
#ifdef MD_CELLINDEX
  case 16:
#endif
    alpha3=rscale*rscale*rscale;
    for(i=0;i<ni;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j]*qi[i]*alpha3*Coulomb_vdw_factor;
    break;
  case 7:
#ifdef MD_CELLINDEX
  case 17:
#endif
    for(i=0;i<ni;i++) force[i*3]+=ftmp[i][0]*qi[i]*rscale*Coulomb_vdw_factor;
    break;
  default:
    fprintf(stderr,"** error : not supported mode=%d **\n",tblno);
    vg_exit(1);
    break;
  }
  //  for(i=0;i<3;i++) printf("  force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);

#ifndef MD_USE_QAUNION
  MR3_free_pointer(atypei,"MR3calccoulomb_ij");
  MR3_free_pointer(atypej,"MR3calccoulomb_ij");
#endif
  MR3_free_pointer(ftmp,"MR3calccoulomb_ij");
}


void mr3calccoulomb_ij_(int *ni, double xi[], double qi[], double force[],
			int *nj, double xj[], double qj[],
			double *rscale, 
			int *tblno, double *xmax, int *periodicflag)
{
  MR3calccoulomb_ij(*ni,xi,qi,force,
		    *nj,xj,qj,
		    *rscale, 
		    *tblno,*xmax,*periodicflag);
}


void mr3calccoulomb_ij__(int *ni, double xi[], double qi[], double force[],
			 int *nj, double xj[], double qj[],
			 double *rscale, 
			 int *tblno, double *xmax, int *periodicflag)
{
  mr3calccoulomb_ij_(ni,xi,qi,force,
		     nj,xj,qj,
		     rscale, 
		     tblno,xmax,periodicflag);
}


void MR3_sort_with_index_int(int data[], int index[], int left, int right)
{
  int tmp,middle,data_middle,index_middle,i,j;
  if(left<right){
    i=left-1;
    j=right+1;
    middle=(left+right)/2;
    data_middle=data[middle];
    index_middle=index[middle]; 
    while(1){
      while(data[++i]<data_middle || 
	     (data[i]==data_middle && index[i]<index_middle));
      while(data[--j]>data_middle ||
	     (data[j]==data_middle && index[j]>index_middle));
      if(i>=j) break;
      tmp=data[i];  data[i]=data[j];   data[j]=tmp;
      tmp=index[i]; index[i]=index[j]; index[j]=tmp;
    }
    MR3_sort_with_index_int(data,index,left,i-1);
    MR3_sort_with_index_int(data,index,j+1,right);
  }
}


void MR3calcvdw_ij_core(int ni, double xi[], int atypei[], double force[],
			int nj, double xj[], int atypej[],
			int nat, double gscale[], double rscale[],
			int tblno, double xmax, int periodicflag)
{
  int i,j;
  VG_UNIT *unit=MR3_unit;
  double vol[3],potc,potv,(*ftmp)[3],*qi,*qj,alpha=1.0,*sqrtrscale;
#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
  double *halfsigmai,*halfsigmaj,*sqrtepsiloni,*sqrtepsilonj;
#endif
#ifdef MD_SORT_ATYPEI
  double *xi2;
  int *atypei2,*index;
#endif

  //  printf("in MR3calcvdw_ij_core: ni=%d nj=%d nat=%d tblno=%d xmax=%e periodicflag=%d\n",ni,nj,nat,tblno,xmax,periodicflag);

#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
  if((sqrtrscale=(double *)MR3_malloc_pointer(sizeof(double)*nat,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc sqrtrscale **\n");
    exit(1);
  }
  if((halfsigmai=(double *)MR3_malloc_pointer(sizeof(double)*ni,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc halfsigmai **\n");
    exit(1);
  }
  if((halfsigmaj=(double *)MR3_malloc_pointer(sizeof(double)*nj,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc halfsigmaj **\n");
    exit(1);
  }
  if((sqrtepsiloni=(double *)MR3_malloc_pointer(sizeof(double)*ni,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc sqrtepsiloni **\n");
    exit(1);
  }
  if((sqrtepsilonj=(double *)MR3_malloc_pointer(sizeof(double)*nj,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc sqrtepsilonj **\n");
    exit(1);
  }
  for(i=0;i<nat;i++) sqrtrscale[i]=sqrt(rscale[i*nat+i]);
  for(i=0;i<ni;i++){
    halfsigmai[i]=0.5*sqrtrscale[atypei[i]];
    //    halfsigmai[i]=0.5*sqrt(rscale[atypei[i]*nat+atypei[i]]);
    sqrtepsiloni[i]=sqrt(gscale[atypei[i]*nat+atypei[i]]);
  }
  for(i=0;i<nj;i++){
    halfsigmaj[i]=0.5*sqrtrscale[atypej[i]];
    //    halfsigmaj[i]=0.5*sqrt(rscale[atypej[i]*nat+atypej[i]]);
    sqrtepsilonj[i]=sqrt(gscale[atypej[i]*nat+atypej[i]]);
  }
#endif

#ifndef MD_USE_QAUNION
  if((qi=(double *)MR3_malloc_pointer(sizeof(double)*ni,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc qi in MR3calcvdw **\n");
    vg_exit(1);
  }
  if((qj=(double *)MR3_malloc_pointer(sizeof(double)*nj,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc qj in MR3calcvdw **\n");
    vg_exit(1);
  }
  bzero(qi,sizeof(double)*ni);
  bzero(qj,sizeof(double)*nj);
  //  for(i=0;i<nj;i++) q[i]=0.0;
#endif
#ifdef MD_SORT_ATYPEI
  if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc xi2 **\n");
    vg_exit(1);
  }
  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc atypei2 **\n");
    vg_exit(1);
  }
  if((index=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc index **\n");
    vg_exit(1);
  }
  for(i=0;i<ni;i++) atypei2[i]=atypei[i];
  for(i=0;i<ni;i++) index[i]=i;
  MR3_sort_with_index_int(atypei2,index,0,ni-1);
  for(i=0;i<ni;i++) for(j=0;j<3;j++) xi2[i*3+j]=xi[index[i]*3+j];
  for(i=0;i<ni;i+=VG_MINIMUM_PARTICLE_BLOCK_I){
    if(i+VG_MINIMUM_PARTICLE_BLOCK_I<ni){
      int flag=0;
      for(j=1;j<VG_MINIMUM_PARTICLE_BLOCK_I;j++){
	if(atypei2[i]!=atypei2[i+j]){
	  flag=1;
	  break;
	}
      }
      if(flag==0){
	//if(1){
      //if(0){
	for(j=0;j<VG_MINIMUM_PARTICLE_BLOCK_I;j++){
	  atypei2[i+j]+=256;
	  //	  printf("atypei2[%d]=%d\n",i,atypei2[i]);
	}
      }
    }
  }
#if 0
  for(i=0;i<ni;i++){
    printf("i=%d atypei=%d atypei2=%d index=%d\n",i,atypei[i],atypei2[i],index[i]);
  }
#endif
#endif
  if((ftmp=(double (*)[3])MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calcvdw_ij_core"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp in MR3calcvdw **\n");
    vg_exit(1);
  }
  if((periodicflag & 1)==0) xmax*=2;
  for(i=0;i<3;i++) vol[i]=xmax;
#if MD_PERIODIC_FIXED==1
  //  unit->xmax=xmax;
  for(i=0;i<3;i++) unit->volume[i]=xmax;
#endif

#if defined(MD_SIGMA_EPSION_IN_VGSTRUCT)
  vg_set_vectors_pos_halfsigma_sqrtepsilon(unit,nj,(double (*)[3])xj,
					   halfsigmaj,sqrtepsilonj);
#elif defined(MD_USE_QAUNION)
  vg_set_vectors_pos_atype(unit,nj,(double (*)[3])xj,atypej);
#else
  vg_set_vectors_pos_charge_atype(unit,nj,(double (*)[3])xj,qj,atypej);
#endif
#ifdef MD_MATRIX_IN_SCALER
  vg_set_scalers_volume_alpha_gscales_rscales(unit,vol,alpha,
					      nat,gscale,gscale,rscale);
#else
  vg_set_matrices_gscales_rscales(unit,nat,nat,gscale,rscale);
  vg_set_scalers_volume_alpha(unit,vol,alpha);
#endif
#if defined(MD_SIGMA_EPSION_IN_VGSTRUCT)
  vg_set_pipeline_vectors_pos_halfsigma_sqrtepsilon(unit,ni,(double (*)[3])xi,
						    halfsigmai,sqrtepsiloni);
#elif defined(MD_USE_QAUNION)
#ifdef MD_SORT_ATYPEI
  vg_set_pipeline_vectors_pos_atype(unit,ni,(double (*)[3])xi2,atypei2);
#else
  vg_set_pipeline_vectors_pos_atype(unit,ni,(double (*)[3])xi,atypei);
#endif
#else
  vg_set_pipeline_vectors_pos_charge_atype(unit,ni,(double (*)[3])xi,qi,atypei);
#endif
  vg_calculate_forces_potentials(unit,tblno,ni,ftmp,&potc,&potv);
#ifdef MD_SORT_ATYPEI
  for(i=0;i<ni;i++) for(j=0;j<3;j++) force[index[i]*3+j]+=ftmp[i][j];
#else
  for(i=0;i<ni;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j];
#endif

#ifndef MD_USE_QAUNION
  MR3_free_pointer(qi,"MR3calcvdw_ij_core");
  MR3_free_pointer(qj,"MR3calcvdw_ij_core");
#endif
  MR3_free_pointer(ftmp,"MR3calcvdw_ij_core");
#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
  MR3_free_pointer(sqrtrscale,"MR3calcvdw_ij_core");
  MR3_free_pointer(halfsigmai,"MR3calcvdw_ij_core");
  MR3_free_pointer(halfsigmaj,"MR3calcvdw_ij_core");
  MR3_free_pointer(sqrtepsiloni,"MR3calcvdw_ij_core");
  MR3_free_pointer(sqrtepsilonj,"MR3calcvdw_ij_core");
#endif
#ifdef MD_SORT_ATYPEI
  MR3_free_pointer(xi2,"MR3calcvdw_ij_core");
  MR3_free_pointer(atypei2,"MR3calcvdw_ij_core");
  MR3_free_pointer(index,"MR3calcvdw_ij_core");
#endif
}


void MR3calcvdw_ij(int ni, double xi[], int atypei[], double force[],
		   int nj, double xj[], int atypej[],
		   int nat, double gscale[], double rscale[],
		   int tblno, double xmax, int periodicflag)
{
  int i,j;
  int *delatype;
  int ni2=0,nj2=0,*atypei2,*atypej2;
  double *xi2,*xj2,*fi2;
  int *new2orgi,*new2orgj;

#if 0
  for(i=0;i<nat;i++){
    for(j=0;j<nat;j++){
      printf("i=%d j=%d gscale=%e rscale=%e\n",i,j,gscale[i*nat+j],rscale[i*nat+j]);
    }
  }
#endif
  if((delatype=(int *)MR3_malloc_pointer(sizeof(int)*nat,"delatype in MR3calcvdw_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc delatype **\n");
    MR3_exit(1);
  }
  if((new2orgi=(int *)MR3_malloc_pointer(sizeof(int)*ni,"new2orgi in MR3calcvdw_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc new2orgi **\n");
    MR3_exit(1);
  }
  if((new2orgj=(int *)MR3_malloc_pointer(sizeof(int)*nj,"new2orgj in MR3calcvdw_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc new2orgj **\n");
    MR3_exit(1);
  }
#if 1
  for(i=0;i<nat;i++){
    delatype[i]=1;
    for(j=0;j<nat;j++){
      if(gscale[i*nat+j]!=0.0) delatype[i]=0;
    }
  }
#else
  for(i=0;i<nat;i++){
    if(gscale[i*nat+i]==0.0){
      delatype[i]=1;
    }
    else{
      delatype[i]=0;
    }
  }
#endif

  for(i=0;i<ni;i++){
    if(delatype[atypei[i]]==0){
      new2orgi[ni2++]=i;
    }
  }
  for(i=0;i<nj;i++){
    if(delatype[atypej[i]]==0){
      new2orgj[nj2++]=i;
    }
  }

  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*ni2,"atypei2 in MR3calcvdw_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc atypei2 **\n");
    MR3_exit(1);
  }
  if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"xi2 in MR3calcvdw_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc xi2 **\n");
    MR3_exit(1);
  }
  if((fi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"fi2 in MR3calcvdw_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc fi2 **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*nj2,"atypej2 in MR3calcvdw_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc atypej2 **\n");
    MR3_exit(1);
  }
  if((xj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2*3,"xj2 in MR3calcvdw_ij"))==NULL){
    fprintf(stderr,"** error : can't malloc xj2 **\n");
    MR3_exit(1);
  }

  for(i=0;i<ni2;i++){
    for(j=0;j<3;j++) xi2[i*3+j]=xi[new2orgi[i]*3+j];
    atypei2[i]=atypei[new2orgi[i]];
    for(j=0;j<3;j++) fi2[i*3+j]=0.0;
  }
  for(i=0;i<nj2;i++){
    for(j=0;j<3;j++) xj2[i*3+j]=xj[new2orgj[i]*3+j];
    atypej2[i]=atypej[new2orgj[i]];
  }

  {
    int nj3,njstart=0,i2;
    do{
      nj3=(nj2-njstart>MD_MAX_J_PER_RUN ? MD_MAX_J_PER_RUN:nj2-njstart);
      MR3calcvdw_ij_core(ni2,xi2,atypei2,fi2,
			 nj3,xj2+njstart*3,atypej2+njstart,
			 nat,gscale,rscale,tblno,xmax,periodicflag);
      njstart+=MD_MAX_J_PER_RUN;
      
    }
    while(nj2-(njstart-MD_MAX_J_PER_RUN)>MD_MAX_J_PER_RUN);
  }

  for(i=0;i<ni2;i++){
    for(j=0;j<3;j++){
      force[new2orgi[i]*3+j]+=fi2[i*3+j];
    }
  }

  MR3_free_pointer(delatype,"delatype in MR3calcvdw_ij");
  MR3_free_pointer(new2orgi,"new2orgi in MR3calcvdw_ij");
  MR3_free_pointer(new2orgj,"new2orgj in MR3calcvdw_ij");
  MR3_free_pointer(atypei2,"atypei2 in MR3calcvdw_ij");
  MR3_free_pointer(xi2,"xi2 in MR3calcvdw_ij");
  MR3_free_pointer(fi2,"fi2 in MR3calcvdw_ij");
  MR3_free_pointer(atypej2,"atypej2 in MR3calcvdw_ij");
  MR3_free_pointer(xj2,"xj2 in MR3calcvdw_ij");
}


void mr3calcvdw_ij_(int *ni, double xi[], int atypei[], double force[],
		    int *nj, double xj[], int atypej[],
		    int *nat, double gscale[], double rscale[],
		    int *tblno, double *xmax, int *periodicflag)
{
    int *atypei2,*atypej2;
    int i;

    if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"atypei2 in mr3calcvdw_ij_"))==NULL){
        fprintf(stderr,
                "** error at malloc atypei2 in mr3calcvdw_ij_ **\n");
        MR3_exit(1);
    }
    if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"atypej2 in mr3calcvdw_ij_"))==NULL){
        fprintf(stderr,
                "** error at malloc atypej2 in mr3calcvdw_ij_ **\n");
        MR3_exit(1);
    }
#if 1
    for(i=0;i<*ni;i++){
      if(atypei[i]>0){
	atypei2[i]=atypei[i]-1;
      }
      else{
	printf("  warning : atypei[%d]=%d should be positive\n",
	       i,atypei[i]);
	atypei2[i]=0;
      }
    }
    for(i=0;i<*nj;i++){
      if(atypej[i]>0){
	atypej2[i]=atypej[i]-1;
      }
      else{
	printf("  warning : atypej[%d]=%d should be positive\n",
	       i,atypej[i]);
      }
    }
#else
    for(i=0;i<*ni;i++) atypei2[i]=atypei[i]-1;
    for(i=0;i<*nj;i++) atypej2[i]=atypej[i]-1;
#endif
    
    MR3calcvdw_ij(*ni,xi,atypei2,force,
		  *nj,xj,atypej2,
		  *nat,gscale,rscale,
		  *tblno,*xmax,*periodicflag);
    MR3_free_pointer(atypei2,"atypei2 in mr3calcvdw_ij_");
    MR3_free_pointer(atypej2,"atypej2 in mr3calcvdw_ij_");
}


void mr3calcvdw_ij__(int *ni, double xi[], int atypei[], double force[],
		     int *nj, double xj[], int atypej[],
		     int *nat, double gscale[], double rscale[],
		     int *tblno, double *xmax, int *periodicflag)
{
  mr3calcvdw_ij_(ni,xi,atypei,force,
		 nj,xj,atypej,
		 nat,gscale,rscale,
		 tblno,xmax,periodicflag);
}


void MR3calcvdw(double x[], int n, int atype[], int nat,
		double gscale[], double rscale[],
		int tblno, double xmax, int periodicflag,
		int natchangeflag, double force[])
{
#if 1
  MR3calcvdw_ij(n,x,atype,force,
		n,x,atype,
		nat,gscale,rscale,tblno,xmax,periodicflag);
#else
  int i,j;
  VG_UNIT *unit=MR3_unit;
  double vol[3],potc,potv,(*ftmp)[3],*q,alpha=1.0;

#ifndef MD_USE_QAUNION
  if((q=(double *)MR3_malloc_pointer(sizeof(double)*n,"MR3calcvdw"))==NULL){
    fprintf(stderr,"** error : can't malloc q in MR3calcvdw **\n");
    vg_exit(1);
  }
  for(i=0;i<n;i++) q[i]=0.0;
#endif
  if((ftmp=(double (*)[3])MR3_malloc_pointer(sizeof(double)*n*3,"MR3calcvdw"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp in MR3calcvdw **\n");
    vg_exit(1);
  }
  if((periodicflag & 1)==0) xmax*=2;
  for(i=0;i<3;i++) vol[i]=xmax;
#if MD_PERIODIC_FIXED==1
  //  unit->xmax=xmax;
  for(i=0;i<3;i++) unit->volume[i]=xmax;
#endif

#ifdef MD_USE_QAUNION
  vg_set_vectors_pos_atype(unit,n,(double (*)[3])x,atype);
#else
  vg_set_vectors_pos_charge_atype(unit,n,(double (*)[3])x,q,atype);
#endif
  vg_set_matrices_gscales_rscales(unit,nat,nat,gscale,rscale);
  vg_set_scalers_volume_alpha(unit,vol,alpha);
#ifdef MD_USE_QAUNION
  vg_set_pipeline_vectors_pos_atype(unit,n,(double (*)[3])x,atype);
#else
  vg_set_pipeline_vectors_pos_charge_atype(unit,n,(double (*)[3])x,q,atype);
#endif
  vg_calculate_forces_potentials(unit,tblno,n,ftmp,&potc,&potv);
  for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j];

#ifndef MD_USE_QAUNION
  MR3_free_pointer(q,"MR3calcvdw");
#endif
  MR3_free_pointer(ftmp,"MR3calcvdw");
#endif
}


void mr3calcvdw_(double x[], int *n, int atype[], int *nat,
		 double gscale[], double rscale[],
		 int *tblno, double *xmax, int *periodicflag,
		 int *natchangeflag, double force[])
{
  int *atype2;
  int i;

  if((atype2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"mr3calcvdw_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atype2 in mr3calcvdw_ **\n");
    MR3_exit(1);
  }
  for(i=0;i<*n;i++){
    if(atype[i]>0){
      atype2[i]=atype[i]-1;
    }
    else{
      printf("  warning : atype[%d]=%d should be positive\n",
	     i,atype[i]);
      atype2[i]=0;
    }
  }
  MR3calcvdw(x,*n,atype2,*nat,
	     gscale,rscale,
	     *tblno,*xmax,*periodicflag,
	     *natchangeflag,force);
  MR3_free_pointer(atype2,"mr3calcvdw_");
}


void mr3calcvdw__(double x[], int *n, int atype[], int *nat,
		  double gscale[], double rscale[],
		  int *tblno, double *xmax, int *periodicflag,
		  int *natchangeflag, double force[])
{
  mr3calcvdw_(x,n,atype,nat,
	      gscale,rscale,
	      tblno,xmax,periodicflag,
	      natchangeflag,force);
}


void MR3calcewald(int *k, int knum_org, double *x, int n, double *chg,
		  double alpha, double epsilon, double cell[3][3],
		  double *force, double *tpot, double stress[3][3])
{
}


void mr3calcewald_(int *k, int *knum, double *x, int *n,
		   double *chg, double *alpha, double *epsilon,
		   double cell[3][3], double *force,
		   double *tpot, double stress[3][3])
{
    MR3calcewald(k,*knum,x,*n,chg,*alpha,*epsilon,cell,
		 force,tpot,stress);
}


void mr3calcewald__(int *k, int *knum, double *x, int *n,
		    double *chg, double *alpha, double *epsilon,
		    double cell[3][3], double *force,
		    double *tpot, double stress[3][3])
{
    mr3calcewald_(k,knum,x,n,chg,alpha,epsilon,cell,
		  force,tpot,stress);
}


void MR3calcnacl(double x[], int n, int atype[], int nat, 
		 double pol[], double sigm[], double ipotro[], 
		 double pc[], double pd[], double zz[],
		 int tblno, double xmax, int periodicflag,
		 double force[])
{
  int i,j;
  VG_UNIT *unit=MR3_unit;
  double vol[3],potc,potv,(*ftmp)[3],*q,alpha=1.0;

  if((q=(double *)MR3_malloc_pointer(sizeof(double)*n,"MR3calcnacl"))==NULL){
    fprintf(stderr,"** error : can't malloc q in MR3calcvdw **\n");
    vg_exit(1);
  }
  if((ftmp=(double (*)[3])MR3_malloc_pointer(sizeof(double)*n*3,"MR3calcnacl"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp in MR3calcvdw **\n");
    vg_exit(1);
  }
  for(i=0;i<n;i++) q[i]=0.0;
  if((periodicflag & 1)==0) xmax*=2;
  for(i=0;i<3;i++) vol[i]=xmax;
#if MD_PERIODIC_FIXED==1
  //  unit->xmax=xmax;
  for(i=0;i<3;i++) unit->volume[i]=xmax;
#endif

  vg_set_vectors_pos_charge_atype(unit,n,(double (*)[3])x,q,atype);
#if 1
  vg_set_pol_sigm_ipotro_pc_pd_zz(unit,nat,nat,pol,sigm,ipotro,pc,pd,zz);
#else // speed test with vdw
  double gscale[4]={1.0,1.0,1.0,1.0},rscale[4]={1.0,1.0,1.0,1.0};
  vg_set_matrices_gscales_rscales(unit,nat,nat,gscale,rscale);
#endif
  vg_set_scalers_volume_alpha(unit,vol,alpha);
  vg_set_pipeline_vectors_pos_charge_atype(unit,n,(double (*)[3])x,q,atype);
  vg_calculate_forces_potentials(unit,tblno,n,ftmp,&potc,&potv);

  for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]+=ftmp[i][j];
  printf("potc=%e potv=%e force=%e %e %e\n",
	 potc,potv,force[0],force[1],force[2]);

  MR3_free_pointer(q,"MR3calcnacl");
  MR3_free_pointer(ftmp,"MR3calcnacl");
}


void vtgrape_force(double xj[][3], // position of j-th particles
		   double mj[],    // mass of j-th particles
		   double xi[][3], // position of i-th particles
		   double eps2,    // softening parameter
		   double a[][3], // force of i-th particles
		   int ni,         // number of i-th particles
		   int nj)
{
  double epsinv;
  int i, j, k;
  double cycle;
  double *p;
  long tsc_end, tsc_start;
  double eps;
  static int ini=0;
  VG_UNIT *unit;
  double size=10.0;

  if(ini==0){
    MR3init();
    ini=1;
  }
  eps=sqrt(eps2);
  //  eps=0.0f;
  unit=m3_get_unit();
  m3_set_softening(unit,eps);
#if 1
  // mi is not come from upper routine
  MR3calccoulomb_ij(ni,(double *)xi,mj,(double *)a,
		    nj,(double *)xj,mj,
		    1.0,0,size,0);
#else  
  {
    static int ini=0;
    if(ini==0){
      printf("** VTGRAPE does not support MR3calccoulomb_ij **\n");
      ini=1;
    }
  }
  MR3calccoulomb((double *)xi,ni,mj,1.0,
		 0,size,0,0,(double *)a);
#endif
  for(i=0;i<ni;i++) for(j=0;j<3;j++) a[i][j]=-a[i][j];
  //  for(i=0;i<ni;i++) for(j=0;j<3;j++) a[i][j]*=-mi[i];
}


#if MD_LJ_05SIGMA_8SIGMA==2 || 1
void vg_set_rcut2(VG_UNIT *unit, double rcut2)
{
  //  printf("vg_set_rcut2 is called. rcut2=%f\n",rcut2);
  VG_SCALER *scal=unit->scalers;
  unit->rcut2=rcut2;
  scal->r2max=(float)rcut2;
#ifdef MD_PRINT_WARN
  {
    static int ini=0;
    if(ini==0){
      if(rcut2!=MD_LJ_R2MAX) printf("rcut2 is modified to %e in vg_set_rcut2\n",rcut2);
      ini=1;
    }
  }
#endif
#if defined(COULOMB_SHIFT) || 1
  if(scal==NULL){
    fprintf(stderr,"** error : scal must be set **\n");
    vg_exit(1);
  }
  scal->rcut21=(float)(1.0/unit->rcut2);
#endif
}
#endif


void mr3_set_rcut2_(double *rcut2)
{
  VG_UNIT *unit=m3_get_unit();
  if(unit!=NULL){
    vg_set_rcut2(unit,*rcut2);
  }
  else{
    fprintf(stderr,"** error : unit is NULL in mr3_set_rcut2_ **\n");
    MR3_exit(1);
  }
}


void vg_set_cuton(VG_UNIT *unit, double cuton)
{
  unit->cuton=cuton;
}


void mr3_set_cuton_(double *cuton)
{
  VG_UNIT *unit=m3_get_unit();
  if(unit!=NULL){
    vg_set_cuton(unit,*cuton);
  }
  else{
    fprintf(stderr,"** error : unit is NULL in mr3_set_cuton_ **\n");
    MR3_exit(1);
  }
}


#include "mddriver.c"
#include "vtgrape_mr3.c"
#include "vtgrape_cellindex.c"
#include "mr3_host.c"
#include "vtgrape_dummy.c"
#include "sock_lib.c"

#undef LARGE_VAL
#undef LOWER_VAL_FACTOR
#undef LOWER_VAL_FACTOR_1
#define LARGE_VAL              LARGE_VAL_ORG
#define LOWER_VAL_FACTOR       LOWER_VAL_FACTOR_ORG
#define LOWER_VAL_FACTOR_1     LOWER_VAL_FACTOR_1_ORG
#include "vtgrape_mixed.c"

