#include "ko.h"

int
mpi_recvdata()
{
  
#ifdef MPI
  int i,j,k,iv;
  MPI_Status status;
  double temp;
  int verbose=1;
  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      MPI_Recv(msgbuf1, NY*NZ*NV*NG, MPI_DOUBLE,
	       mpi_tile2procid(TI+1,TJ,TK), MPI_MSG_XLO, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_XLO from %d\n",PROCID,mpi_tile2procid(TI+1,TJ,TK));
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbuf1[(i+NG)*NY*NZ*NV + j*NZ*NV + k*NV + iv]);
    }

  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      MPI_Recv(msgbuf1, NY*NZ*NV*NG, MPI_DOUBLE,
	       mpi_tile2procid(TI-1,TJ,TK), MPI_MSG_XHI, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_XHI from %d\n",PROCID,mpi_tile2procid(TI-1,TJ,TK));
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbuf1[(i-NX)*NY*NZ*NV + j*NZ*NV + k*NV + iv]);
    } 
  return 0;
 //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      MPI_Recv(msgbuf1, NX*NG*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ+1,TK), MPI_MSG_YLO, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_YLO from %d\n",PROCID,mpi_tile2procid(TI,TJ+1,TK));
      for(i=0;i<NX;i++)
	for(j=-NG;j<0;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbuf1[i*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);

    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      MPI_Recv(msgbuf1, NX*NG*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ-1,TK), MPI_MSG_YHI, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_YHI from %d\n",PROCID,mpi_tile2procid(TI,TJ-1,TK));
      for(i=0;i<NX;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbuf1[i*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);

    }
 //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      MPI_Recv(msgbuf1, NX*NY*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK+1), MPI_MSG_ZLO, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_ZLO from %d\n",PROCID,mpi_tile2procid(TI,TJ,TK+1));
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbuf1[i*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]);
    }
   //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      MPI_Recv(msgbuf1, NX*NY*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK-1), MPI_MSG_ZHI, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_ZHI from %d\n",PROCID,mpi_tile2procid(TI,TJ,TK-1));
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbuf1[i*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]);

    }
 
#endif
}

int
mpi_senddata()
{  
  int i,j,k,iv;
  int verbose=1;
  ldouble temp;
#ifdef MPI
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbuf1[(i+NG)*NY*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      if(verbose) printf("%d attempting to send MPI_MSG_XLO to %d\n",PROCID,mpi_tile2procid(TI-1,TJ,TK));
      MPI_Send(msgbuf1, NY*NZ*NV*NG, MPI_DOUBLE,
	       mpi_tile2procid(TI-1,TJ,TK), MPI_MSG_XLO, MPI_COMM_WORLD);
      if(verbose) printf("%d sent MPI_MSG_XLO to %d\n",PROCID,mpi_tile2procid(TI-1,TJ,TK));
      
    }

  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbuf1[(i-NX)*NY*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      if(verbose) printf("%d attempting to send MPI_MSG_XHI to %d\n",PROCID,mpi_tile2procid(TI+1,TJ,TK));
      MPI_Send(msgbuf1,NV*NY, MPI_DOUBLE,
	       mpi_tile2procid(TI+1,TJ,TK), MPI_MSG_XHI, MPI_COMM_WORLD);
      if(verbose) printf("%d sent MPI_MSG_XHI to %d\n",PROCID,mpi_tile2procid(TI+1,TJ,TK));

    }
return 0;
  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=-NG;j<0;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbuf1[i*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Send(msgbuf1, NX*NG*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ-1,TK), MPI_MSG_YLO, MPI_COMM_WORLD);
    }
  //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbuf1[i*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Send(msgbuf1, NX*NG*NZ*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ+1,TK), MPI_MSG_YHI, MPI_COMM_WORLD);
    }
  //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbuf1[i*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Send(msgbuf1, NX*NY*NG*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK-1), MPI_MSG_ZLO, MPI_COMM_WORLD);
    }
  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbuf1[i*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Send(msgbuf1, NX*NY*NG*NV, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK+1), MPI_MSG_ZHI, MPI_COMM_WORLD);
    }
#endif

  return 0;
}

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
