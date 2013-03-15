//KORAL - problem.c
//problem-related but not problem-specific routines

 #include "ko.h"

/*********************************************/
/* returns locations of cell boundaries */
/*********************************************/
ldouble
calc_xb(int i,int idim)
{
#ifdef LOGXGRID
  //logarithimic in x
  ldouble c0,c1,dc,xb;
  if(idim==0)
    {
      c0=1.;
      c1=logl(MAXX/MINX);
      dc=(c1)/(ldouble)NX; 
      xb=(ldouble)i*dc;
     xb=MINX*expl(xb);
    }
  if(idim==1)
    {
      c0=MINY;
      c1=MAXY;
      dc=(c1-c0)/(ldouble)NY; 
      xb=c0+(ldouble)i*dc;
    }
  if(idim==2)
    {
      c0=MINZ;
      c1=MAXZ;
      dc=(c1-c0)/(ldouble)NZ;
      xb=c0+(ldouble)i*dc;
    }

  return xb;
#elif defined (LOGXGRIDREF)
  //logartithmic refined to center
  //x=log((r-p1)/p2)
  ldouble c0,c1,dc,xb;
  if(idim==0)
    {
      c0=logl((MINX-LOGPAR1)/LOGPAR2);
      c1=logl((MAXX-LOGPAR1)/LOGPAR2);
      dc=(c1-c0)/(ldouble)NX; 
      xb=c0+(ldouble)(i)*dc;
      xb=LOGPAR1+2.*expl(xb);
    }
  if(idim==1)
    {
      c0=MINY;
      c1=MAXY;
      dc=(c1-c0)/(ldouble)NY; 
      xb=c0+(ldouble)i*dc;
    }
  if(idim==2)
    {
      c0=MINZ;
      c1=MAXZ;
      dc=(c1-c0)/(ldouble)NZ;
      xb=c0+(ldouble)i*dc;
    }

  return xb;
#else

  //default - uniform grid
  ldouble c0,c1,dc;
  if(idim==0)
    {
      c0=MINX;
      c1=MAXX;
      dc=(c1-c0)/(ldouble)NX;
    }
  if(idim==1)
    {
      c0=MINY;
      c1=MAXY;
      dc=(c1-c0)/(ldouble)NY;
    }
  if(idim==2)
    {
      c0=MINZ;
      c1=MAXZ;
      dc=(c1-c0)/(ldouble)NZ;
    }

  ldouble xb=c0+(ldouble)i*dc;

  return xb;
#endif

}

/*********************************************/
/* returns problem specific BC defined in PROBLEMS/XXX/bc.c */
/*********************************************/
int
calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp)
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
  int ix,iy,iz;
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {

#include PR_INIT

	    }
	}
    }

  return 0;
}

/*********************************************/
/* may be useful some day */
/*********************************************/
int
initialize_problem()
{
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
/* suplementary routines
/* uses PROBLEMS/XXX/tools.c */   
/*********************************************/
#include PR_TOOLS



