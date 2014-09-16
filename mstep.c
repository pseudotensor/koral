
//KORAL - mstep.c
//routines related to multi-time-stepping

#include "ko.h"


int mstep_cell_levels[NX][NY][NZ];
int mstep_current_counts[NUMMSTEPLEVELS];

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

  for(i=0;i<NUMMSTEPLEVELS;i++)
    {
      mstep_current_counts[i]=0;
    }

  ldouble pow2=1.;
  for(i=0;i<NUMMSTEPLEVELS+1;i++)
    {
      mstep_multiplier[i]=pow2;
      pow2*=2.;
    }

  
  return 0;
}

/* iterates the counters */
int 
mstep_iterate(void)
{
  int i;
 
  for(i=1;i<NUMMSTEPLEVELS;i++) //[0]=0 - always evolved
    {
      mstep_current_counts[i]++;
      if(mstep_current_counts[i]==(int)mstep_multiplier[i]) 
	mstep_current_counts[i]=0;
    }
  
  return 0;
}

/* prints_levels_to screen */
int
mstep_print_levels()
{
  int ix,iy,iz;
  
  for(ix=0;ix<NX;ix++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {
	      printf("%1d ",mstep_cell_levels[ix][iy][iz]);
	    }
	}
    }
  printf("\n");

  return 0;
}

/* update levels for a cell with new wavespeeds */
int
mstep_update_levels()
{
  ldouble dtmin = 1./tstepdenmax,dt;
  int levelnew,levelold;
  int ix,iy,iz,update,il;
  
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  if(mstep_is_cell_active(ix,iy,iz)==0) //only active cells can be updated
	    continue;

	  update=0;
	  dt=1./get_u_scalar(cell_tsteps,ix,iy,iz); //new time step
	  levelold=mstep_cell_levels[ix][iy][iz];
	  levelnew=mstep_calc_level(dt,dtmin);

	  //printf("%d %e %e > %d %d\n",ix,dt,dtmin,levelnew,levelold);
	  
	  if(levelnew>levelold) //can switch to longer time only when this one active
	    {
	      for(il=levelnew;il>levelold;il--) //tries from the longest
		{
		  if(mstep_is_level_active(il)==1) break;
		}
	      levelnew=il;
	    }

	  if(levelnew!=levelold) update=1; //and can decrease the timestep at any time
  
	  if(update)
	      mstep_cell_levels[ix][iy][iz]=levelnew;
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

  for(il=0;il<NUMMSTEPLEVELS;il++)
    {
      if(div < mstep_multiplier[il+1]) 
	{
	  level=il;
	  break;
	}
    }

  if(il==NUMMSTEPLEVELS)
    level=NUMMSTEPLEVELS-1;

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

/* say if given cell or its neighour in dim is evolved in this time step */
int
mstep_is_cell_or_neighbour_active(int ix, int iy, int iz,int idim)
{
  if(mstep_current_counts[mstep_cell_levels[ix][iy][iz]]==0)
    return 1;

  if(idim==0)
    {
      if(ix>0 && mstep_current_counts[mstep_cell_levels[ix-1][iy][iz]]==0)
	return 1;
      if(ix<NX-1 && mstep_current_counts[mstep_cell_levels[ix+1][iy][iz]]==0)
	return 1;
    }

  if(idim==1)
    {
      if(iy>0 && mstep_current_counts[mstep_cell_levels[ix][iy-1][iz]]==0)
	return 1;
      if(iy<NY-1 && mstep_current_counts[mstep_cell_levels[ix][iy+1][iz]]==0)
	return 1;
    }

  if(idim==2)
    {
      if(iz>0 && mstep_current_counts[mstep_cell_levels[ix][iy][iz-1]]==0)
	return 1;
      if(iz<NZ-1 && mstep_current_counts[mstep_cell_levels[ix][iy][iz+1]]==0)
	return 1;
    }

  //otherwise
  return 0;
}

/* say if given face in dim is needed in this time step */
int
mstep_is_face_active(int ix, int iy, int iz,int idim)
{
  if(idim==0)
    {
      if(ix>0 && mstep_current_counts[mstep_cell_levels[ix-1][iy][iz]]==0)
	return 1;
      if(ix<NX && mstep_current_counts[mstep_cell_levels[ix][iy][iz]]==0)
	return 1;
    }

  if(idim==1)
    {
      if(iy>0 && mstep_current_counts[mstep_cell_levels[ix][iy-1][iz]]==0)
	return 1;
      if(iy<NY && mstep_current_counts[mstep_cell_levels[ix][iy][iz]]==0)
	return 1;
    }

  if(idim==2)
    {
      if(iz>0 && mstep_current_counts[mstep_cell_levels[ix][iy][iz-1]]==0)
	return 1;
      if(iz<NZ && mstep_current_counts[mstep_cell_levels[ix][iy][iz]]==0)
	return 1;
    }

  //otherwise
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

/* gets number of cells at each level */
int
mstep_count_levels(int* nlev)
{
  int i;
  for(i=0;i<NUMMSTEPLEVELS;i++)
    nlev[i]=0;

  int ix,iy,iz;
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	nlev[mstep_cell_levels[ix][iy][iz]]++;

  return 0;
}

/* gets time step multiplier for a given cell */
ldouble
mstep_get_cell_multiplier(int ix,int iy,int iz)
{
  return mstep_multiplier[mstep_cell_levels[ix][iy][iz]];
}

/* gets time step multiplier for a given face in idim, i.e., shorter from left/right */
ldouble
mstep_get_face_multiplier(int ix,int iy,int iz,int dim)
{
  ldouble ml,mr,m;
  ml=mr=-1.;

  if(dim==0)
    {
      if(ix>0)
	ml=mstep_multiplier[mstep_cell_levels[ix-1][iy][iz]];
      if(ix<NX)
	mr=mstep_multiplier[mstep_cell_levels[ix][iy][iz]];
    }
  if(dim==1)
    {
      if(iy>0)
	ml=mstep_multiplier[mstep_cell_levels[ix][iy-1][iz]];
      if(iy<NY)
	mr=mstep_multiplier[mstep_cell_levels[ix][iy][iz]];
    }
  if(dim==2)
    {
      if(iz>0)
	ml=mstep_multiplier[mstep_cell_levels[ix][iy][iz-1]];
      if(iz<NZ)
	mr=mstep_multiplier[mstep_cell_levels[ix][iy][iz]];
    }

  if(ml<0. && mr<0.) 
    {
      printf("wrong multiplier in mstep_get_face_multiplier() for %d %d %d %d\n",ix,iy,iz,dim);
      m=1.;
    }
  else if(ml<0.)
    {
      m=mr;
    }
  else if(mr<0.)
    {
      m=ml;
    }
  else
    m=my_min(ml,mr);

  return m;
}

/* tests */
int
mstep_test()
{
  mstep_cell_levels[0][0][0]=mstep_calc_level(1.5,1.);
  mstep_cell_levels[1][0][0]=mstep_calc_level(2.5,1.);
  mstep_cell_levels[2][0][0]=mstep_calc_level(4.5,1.);
  mstep_cell_levels[3][0][0]=mstep_calc_level(9.5,1.);

  int iter,i,j;

  for(i=0;i<NUMMSTEPLEVELS;i++)
    printf("%2.0f ",mstep_multiplier[i]);
  printf("|   ");
  for(i=0;i<4;i++)
    printf("%d ",mstep_cell_levels[i][0][0]);
  
  printf("\n\n");

  for(iter=0;iter<65;iter++)
    {
      for(i=0;i<NUMMSTEPLEVELS;i++)
	printf("%2d ",mstep_current_counts[i]);
      printf("   |   ");
      for(i=0;i<4;i++)
	printf("%d ",mstep_is_cell_active(i,0,0));
      printf("\n");
      mstep_iterate();
    }

  return 0;
}


