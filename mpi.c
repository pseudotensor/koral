#include "ko.h"

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
