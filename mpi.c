#include "ko.h"

int
mpi_recvdata()
{
#ifdef MPI
  MPI_Status status;
  double temp;
  int verbose=1;
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      MPI_Recv(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI-1,TJ,TK), MPI_MSG_XHI, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_XHI from %d\n",PROCID,mpi_tile2procid(TI-1,TJ,TK));
    }
  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      MPI_Recv(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI+1,TJ,TK), MPI_MSG_XLO, MPI_COMM_WORLD, &status);
       if(verbose) printf("%d received MPI_MSG_XLO from %d\n",PROCID,mpi_tile2procid(TI+1,TJ,TK));
    }
  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      MPI_Recv(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ-1,TK), MPI_MSG_YHI, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_YHI from %d\n",PROCID,mpi_tile2procid(TI,TJ-1,TK));
    }
  //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      MPI_Recv(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ+1,TK), MPI_MSG_YLO, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_YLO from %d\n",PROCID,mpi_tile2procid(TI,TJ+1,TK));
    }
  //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      MPI_Recv(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK-1), MPI_MSG_ZHI, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_ZHI from %d\n",PROCID,mpi_tile2procid(TI,TJ,TK-1));
    }
  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      MPI_Recv(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK+1), MPI_MSG_ZLO, MPI_COMM_WORLD, &status);
      if(verbose) printf("%d received MPI_MSG_ZLO from %d\n",PROCID,mpi_tile2procid(TI,TJ,TK+1));
    }
#endif
}

int
mpi_senddata()
{
#ifdef MPI
  double temp;
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      MPI_Send(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI-1,TJ,TK), MPI_MSG_XLO, MPI_COMM_WORLD);
    }
  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      MPI_Send(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI+1,TJ,TK), MPI_MSG_XHI, MPI_COMM_WORLD);
    }
  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      MPI_Send(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ-1,TK), MPI_MSG_YLO, MPI_COMM_WORLD);
    }
  //upper x
  if(mpi_isitBC(YBCHI)==0)
    {
      MPI_Send(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ+1,TK), MPI_MSG_YHI, MPI_COMM_WORLD);
    }
  //lower x
  if(mpi_isitBC(ZBCLO)==0)
    {
      MPI_Send(&temp, 1, MPI_DOUBLE,
	       mpi_tile2procid(TI,TJ,TK-1), MPI_MSG_ZLO, MPI_COMM_WORLD);
    }
  //upper x
  if(mpi_isitBC(ZBCHI)==0)
    {
      MPI_Send(&temp, 1, MPI_DOUBLE,
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
