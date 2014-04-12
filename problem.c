//KORAL - problem.c
//problem-related but not problem-specific routines

 #include "ko.h"

/*********************************************/
/* returns locations of cell boundaries */
/*********************************************/
ldouble
calc_xb(int i,int idim)
{
  ldouble c0,c1,dc,xb;
  if(idim==0)
    {
      c0=MINX;
      c1=MAXX;
      dc=(c1-c0)/(ldouble)TNX;
      xb=c0+(ldouble)(i+TOI)*dc;
    }
  if(idim==1)
    {
      c0=MINY;
      c1=MAXY;
      dc=(c1-c0)/(ldouble)TNY;
      xb=c0+(ldouble)(i+TOJ)*dc;
    }
  if(idim==2)
    {
      c0=MINZ;
      c1=MAXZ;
      dc=(c1-c0)/(ldouble)TNZ;
      xb=c0+(ldouble)(i+TOK)*dc;
    }

  return xb;
}

/*********************************************/
/* returns problem specific BC defined in PROBLEMS/XXX/bc.c */
/*********************************************/
int
calc_bc(int ix,int iy,int iz,ldouble t,
	ldouble *uu,ldouble *pp,int ifinit,int BCtype)
{

#include PR_BC
  
  return 0;
}

/*********************************************/
/* sets the initial profile defined in PROBLEMS/XXX/init.c */
/*********************************************/
int
set_initial_profile()
{
  if(PROCID==0) {printf("Initializing problem... ");fflush(stdout);}
  int ix,iy,iz;


#pragma omp parallel for private(ix,iy,iz) schedule (dynamic)
  for(ix=0;ix<NX;ix++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {

#include PR_INIT

	    }
	}
    }

  #ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
  #endif

  if(PROCID==0) printf("done!\n");

  return 0;
}



/*********************************************/
/* the humanoidal finger */
/*********************************************/
int my_finger(ldouble t)
{

#ifdef PR_FINGER
#include PR_FINGER
#endif

  return 0;
}

/*********************************************/
/* suplementary routines */
/* uses PROBLEMS/XXX/tools.c */   
/*********************************************/
#include PR_TOOLS



