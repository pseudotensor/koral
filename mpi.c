#include "ko.h"

void
mpi_senddata()
{


}

void
mpi_recvdata()
{


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

void
mpi_tile2procid(int tilei, int tilej, int tilek, int *procid)
{
  *procid = tilek * NTX * NTY + tilej * NTX + tilei;
}
