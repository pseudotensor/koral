#include "ko.h"

int 
mpi_exchangedata()
{
#ifdef MPI
  MPI_Request reqs[12];
  int nreqs=0;
  mpi_senddata(reqs,&nreqs);
  mpi_recvdata(reqs,&nreqs);
  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
  mpi_savedata();

  /*
  int i;
  printf("sent:\n");
  for(i=0;i<10;i++)
    printf("%e.",msgbufs[0][i]);
  printf("\n");

  
  printf("recvd:\n");
  for(i=0;i<10;i++)
    printf("%e.",msgbufs[6][i]);
  printf("\n");
  */

  MPI_Barrier(MPI_COMM_WORLD);

  
#endif
  return 0;
}

#ifdef MPI //thick bracket comments out recv and send so that non-mpi compiler can compile without MPI_Request
int
mpi_recvdata(MPI_Request *reqs, int *nreqs)
{
  
  int i,j,k,iv;
  MPI_Status status;
  double temp;
  int verbose=0;
  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      MPI_Irecv(msgbufs[6], NY*NZ*NV*NG, MPI_DOUBLE,
	       mpi_tile2procid(TI+1,TJ,TK), MPI_MSG_XLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_XLO from %d\n",PROCID,mpi_tile2procid(TI+1,TJ,TK));
   }
  
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      MPI_Irecv(msgbufs[7], NG*NY*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI-1,TJ,TK), MPI_MSG_XHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_XHI from %d\n",PROCID,mpi_tile2procid(TI-1,TJ,TK));
    } 

  //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      MPI_Irecv(msgbufs[8], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(TI,TJ+1,TK), MPI_MSG_YLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_YLO from %d\n",PROCID,mpi_tile2procid(TI,TJ+1,TK));
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      MPI_Irecv(msgbufs[9], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(TI,TJ-1,TK), MPI_MSG_YHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_YHI from %d\n",PROCID,mpi_tile2procid(TI,TJ-1,TK));
    }
  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      MPI_Irecv(msgbufs[10], NX*NY*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(TI,TJ,TK+1), MPI_MSG_ZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      MPI_Irecv(msgbufs[11], NX*NY*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(TI,TJ,TK-1), MPI_MSG_ZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_ZHI from %d\n",PROCID,mpi_tile2procid(TI,TJ,TK-1));
    }
 
}


int
mpi_savedata()
{
  
  int i,j,k,iv;
  MPI_Status status;
  double temp;
  int verbose=0;
  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[6][(i-NX)*NY*NZ*NV + j*NZ*NV + k*NV + iv]);
	      //printf("r %d > %d %d > %e\n",PROCID,i,j,msgbufs[6][(i-NX)*NY*NZ*NV + j*NZ*NV + k*NV + 0]);
	    }
    }
  
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[7][(i+NG)*NY*NZ*NV + j*NZ*NV + k*NV + iv]);
    } 

 //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[8][i*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=-NG;j<0;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[9][i*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
    }

  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[1][i*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]);
    }
   //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[11][i*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]);
    }
 
}

int
mpi_senddata(MPI_Request *reqs, int *nreqs)
{  
  int i,j,k,iv;
  int verbose=0;
  ldouble temp;

  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      for(i=0;i<NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[0][(i)*NY*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      //printf("s %d > %d %d > %e %e\n",PROCID,i,j,msgbufs[0][(i)*NY*NZ*NV + j*NZ*NV + k*NV + 0],get_u(p,RHO,i,j,k));
	    }
      MPI_Isend(msgbufs[0], NY*NZ*NV*NG, MPI_DOUBLE,
		mpi_tile2procid(TI-1,TJ,TK), MPI_MSG_XLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d sent MPI_MSG_XLO to %d\n",PROCID,mpi_tile2procid(TI-1,TJ,TK));
    }

  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[1][(i-NX+NG)*NY*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[1],NG*NY*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI+1,TJ,TK), MPI_MSG_XHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d sent MPI_MSG_XHI to %d\n",PROCID,mpi_tile2procid(TI+1,TJ,TK));
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NG;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[2][i*NG*NZ*NV + (j)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[2], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(TI,TJ-1,TK), MPI_MSG_YLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }

  //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[3][i*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[3], NX*NG*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ+1,TK), MPI_MSG_YHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NG;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[4][i*NY*NG*NV + j*NG*NV + (k)*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[4], NX*NY*NG*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK-1), MPI_MSG_ZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[5][i*NY*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[5], NX*NY*NG*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK+1), MPI_MSG_ZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }

  return 0;
}
#endif

//verify if there is real BC at all, if set_bc() needed
int
mpi_hasBC()
{
#ifndef MPI
  return 1;
#else
  if(TI==0 || TI==NTX-1 || TJ==0 || TJ==NTY-1 || TK==0 || TK==NTZ-1)
    return 1; //this cell has some real BC
  else
    return 0;
#endif
}

//verify if given cell from set_bc() falls into a real BC
int
mpi_isitBC(int BCtype)
{
#ifndef MPI
  return 1;
#else //check here if real BC
  if(BCtype==XBCLO && TI==0)
    return 1;
  if(BCtype==XBCHI && TI==NTX-1)
    return 1;
  if(BCtype==YBCLO && TJ==0)
    return 1;
  if(BCtype==YBCHI && TJ==NTY-1)
    return 1;
  if(BCtype==ZBCLO && TK==0)
    return 1;
  if(BCtype==ZBCHI && TK==NTZ-1)
    return 1;

  return 0; 
#endif
}

void
mpi_synchtiming(ldouble *time)
{
#ifdef MPI

  MPI_Barrier(MPI_COMM_WORLD);
  
  global_time=*time;
  
  MPI_Allreduce(&tstepdenmax, &global_tstepdenmax, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);

  MPI_Bcast(&global_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

  *time=global_time;
  tstepdenmax=global_tstepdenmax;
 #endif
}

void
mpi_myinit(int argc, char *argv[])
{
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &PROCID);
  MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);

  if(NPROCS!=NTX*NTY*NTZ)
    {
      printf("Wrong number of processes. Problem set up for: %d x %d x %d = %d processes.\n",NTX,NTY,NTZ,NTX*NTY*NTZ);
      exit(-1);
    }
  
  mpi_procid2tile(PROCID,&TI,&TJ,&TK);
  mpi_tileorigin(TI,TJ,TK,&TOI,&TOJ,&TOK);

  printf("pid: %d/%d; tot.res: %dx%dx%d; tile.res:  %dx%dx%d\n"
	 "tile: %d,%d,%d; tile orig.: %d,%d,%d\n",PROCID,NPROCS,TNX,TNY,TNZ,NX,NY,NZ,TI,TJ,TK,TOI,TOJ,TOK);

#endif
}

void
mpi_myfinalize()
{
#ifdef MPI
  MPI_Finalize();
#endif
}

void
mpi_tileorigin(int ti, int tj, int tk, int* toi, int* toj, int* tok)
{
  *toi = ti * NX;
  *toj = tj * NY;
  *tok = tk * NZ;
}

void
mpi_procid2tile(int procid, int* tilei, int* tilej, int* tilek)
{
  *tilek = floor(procid / (NTX * NTY));
  *tilej = floor((procid - (*tilek) * NTX * NTY) / NTX);
  *tilei = procid - NTX * NTY * (*tilek) - NTX * (*tilej);
}

int
mpi_tile2procid(int tilei, int tilej, int tilek)
{
  return tilek * NTX * NTY + tilej * NTX + tilei;
}