
//KORAL - mstep.c
//routines related to multi-time-stepping

#include "ko.h"


int mstep_cell_levels[NX][NY][NZ];
int mstep_current_counts[NUMMSTEPS];

/* initiate counters etc. */
int 
mstep_init(void)
{
  int i;
  int ix;int iy;int iz;
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	mstep_cell_levels[ix][iy][iz]=0;

  ldouble pow2=1;
  for(i=0;i<NUMMSTEPS;i++)
    {
      mstep_current_counts[i]=0;
      mstep_level_multiplier[i]=pow2;
      pow2*=2.;
    }

  
  return 0;
}

/* iterates the counters */
int 
mstep_iterate(void)
{
  int i;
 
  for(i=1;i<NUMMSTEPS;i++) //[0]=0 - always evolved
    {
      mstep_current_counts[i]++;
      if(mstep_current_counts[i]==mstep_level_multiplier[i]) 
	mstep_current_counts[i]=0;
    }
  
  return 0;
}

/* determine the time step level */
int
mstep_calc_level(ldouble dt, ldouble dtmin)
{
  int il,level;
  ldouble div;

  level=0;
  div=dt/dtmin;

  for(il=0;il<NUMMSTEPS;il++)
    {
      if(div < mstep_level_multiplier[il+1]) 
	{
	  level=il;
	  break;
	}
    }

  if(il==NUMMSTEPS)
    level=NUMMSTEPS-1;

  return level;
}

/* say if given cell is evolved in this time step */
int
mstep_is_cell_active(int ix, int iy, int iz)
{
  if(mstep_current_counts[mstep_cell_levels[ix][iy][iz]]==0)
    return 1;
  else
    return 0;
}

/* say if given level is evolved in this time step */
int
mstep_is_level_active(int level)
{
  if(mstep_current_counts[level]==0)
    return 1;
  else
    return 0;
}

/* tests */
int
mstep_test(int ix, int iy, int iz)
{
  mstep_init();

  mstep_cell_levels[0][0][0]=mstep_calc_level(1.5,1.);
  mstep_cell_levels[1][0][0]=mstep_calc_level(2.5,1.);
  mstep_cell_levels[2][0][0]=mstep_calc_level(4.5,1.);
  mstep_cell_levels[3][0][0]=mstep_calc_level(9.5,1.);

  int iter,i,j;

  for(i=0;i<NUMMSTEPS;i++)
    printf("%2d ",mstep_level_multiplier[i]);
  printf("|   ");
  for(i=0;i<4;i++)
    printf("%d ",mstep_cell_levels[i][0][0]);
  
  printf("\n\n");

  for(iter=0;iter<65;iter++)
    {
      for(i=0;i<NUMMSTEPS;i++)
	printf("%2d ",mstep_current_counts[i]);
      printf("   |   ");
      for(i=0;i<4;i++)
	printf("%d ",mstep_is_cell_active(i,0,0));
      printf("\n");
      mstep_iterate();
    }

  return 0;
}
