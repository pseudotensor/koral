//KORAL - ana.c
//damp files postprocessing

#include "ko.h"

int 
main(int argc, char **argv)
{  
  #ifdef MPI
  printf("avg works on one core only, do not use MPI, please\n");
  exit(-1);
  #endif

  //which files to read
  int no1,no2,nostep;
  if(argc!=5)
    {
      printf("Not enough input arguments. Asks for ./avg no1 no2 nostep doingavg\n");
      return -1;
    }
  else
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
      doingavg=atoi(argv[4]);
    }

  int i;

  //currently gsl is not used
  gsl_set_error_handler_off();
  
  //random number gen. initialization
  srand ( time(NULL) );

  //preparing arrays
  initialize_arrays();

  //sets the grid
  set_grid(&min_dx,&min_dy,&min_dz,&max_dt);

  //precalculates metric etc.
  calc_metric();

  //folder to write in
  char folder[100],bufor[100];
  sprintf(folder,"analysis");

  //opens the scalar file
  sprintf(bufor,"%s/avgscalars.dat",folder);
  fout_scalars=fopen(bufor,"w");

  //arrays for averaging of primitives
  ldouble *pavgtot=(ldouble*)malloc((SX)*(SY)*(SZ)*(NV+NAVGVARS)*sizeof(ldouble));
  for(i=0;i<(SX)*(SY)*(SZ)*(NV+NAVGVARS);i++)
    pavg[i]=pavgtot[i]=0.;

  if(doingavg)
    {
      printf("working on avg files #%04d to #%04d with %d step \n",no1,no2,nostep);
    }
  else
    {
      printf("working on res files #%04d to #%04d with %d step \n",no1,no2,nostep);
    }
    

  int ifile,readret;
  ldouble t,ttot; ldouble scalars[NSCALARS];
  ttot=0.;
  t=global_time;

  int ix,iy,iz,iv;
  ldouble pp[NV],uu[NV];

  for(ifile=no1;ifile<=no2;ifile+=nostep)
    {
      if(doingavg)
	{
	  //reading avg file
	  readret=fread_avgfile(ifile,"dumps",pavg,&dt);
	}
      else
	{
	  //reading res file
	  readret=fread_restartfile(ifile,"dumps",&t);
	  dt=1.;
	  //rewrite primitives to pavg
	  for(iz=0;iz<NZ;iz++)
	    for(iy=0;iy<NY;iy++)
	      for(ix=0;ix<NX;ix++)
		for(iv=0;iv<NV;iv++)
		  set_uavg(pavg,iv,ix,iy,iz,get_u(p,iv,ix,iy,iz));
	}

      add_u_core(1.,pavgtot,dt,pavg,pavgtot,(SX)*(SY)*(SZ)*(NV+NAVGVARS));
	
      ttot+=dt;
    }

  //average primitives and averaged quantities
  copy_u_core(1./ttot,pavgtot,pavg,(SX)*(SY)*(SZ)*(NV+NAVGVARS));

  //rewrite primitives to p
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  for(iv=0;iv<NV;iv++)
	    set_u(p,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iz));
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	}

  //projects on ghost cells
  set_bc(t,0);

  char prefix[40];
  char suffix[10];
  if(doingavg)
    sprintf(suffix,"");
  else
    sprintf(suffix,"res");
  
  
  //dumps dumps to analysis analysis
#if(RADOUTPUT==1)
  sprintf(prefix,"radavg%s%04d-",suffix,no1);
  fprint_radprofiles(t,no2,"analysis",prefix);
#endif
  
#if(OUTOUTPUT==1)
  sprintf(prefix,"outavg%s%04d-",suffix,no1);
  fprint_outfile(t,no2,0,"analysis",prefix);
#endif
#if(SILOOUTPUT==1)
#ifndef NOSILO
  sprintf(prefix,"silavg%s%04d-",suffix,no1);
  fprint_silofile(t,no2,"analysis",prefix);
#endif
#endif
  
#if(SIMOUTPUT!=0)	  
  sprintf(prefix,"simavg%s%04d-",suffix,no1);
  fprint_simplecart(t,no2,"analysis",prefix);
#endif
  
  return 0;
}

