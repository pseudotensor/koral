//KORAL - ana.c
//damp files postprocessing

#include "ko.h"

int 
main(int argc, char **argv)
{  
  //which files to read
  int no1,no2,nostep;
  if(argc!=4)
    {
      printf("Not enough input arguments. Asks for ./avg no1 no2 nostep\n");
      return -1;
    }
  else
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
    }

  int i;
  doingavg=1;

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

  //sprintf(bufor,"rm %s/*",folder);
  //i=system(bufor);

  //opens the scalar file
  sprintf(bufor,"%s/scalars.dat",folder);
  fout_scalars=fopen(bufor,"w");

  //arrays for averaging of primitives
  ldouble *pavgtot=(ldouble*)malloc((SX)*(SY)*(SZ)*(NV+NAVGVARS)*sizeof(ldouble));

  for(i=0;i<(SX)*(SY)*(SZ)*(NV+NAVGVARS);i++)
    pavg[i]=pavgtot[i]=0.;

  int ifile,readret;
  ldouble t,ttot; ldouble scalars[NSCALARS];
  ttot=0.;
  t=global_time;

  printf("working on files #%04d to #%04d with %d step \n",no1,no2,nostep);

  int ix,iy,iz,iv;
  ldouble pp[NV],uu[NV];

  for(ifile=no1;ifile<=no2;ifile+=nostep)
    {
      //reading avg file
      readret=fread_avgfile(ifile,pavg,&dt);
      add_u_core(1.,pavgtot,dt,pavg,pavgtot,(SX)*(SY)*(SZ)*(NV+NAVGVARS));
      ttot+=dt;
    }

  //average primitives and averaged quantities
  copy_u_core(1./ttot,pavgtot,pavg,(SX)*(SY)*(SZ)*(NV+NAVGVARS));

  //rewrite primitives to p
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	for(iv=0;iv<NV;iv++)
	   set_u(p,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iz));

  //projects on ghost cells
  set_bc(t,0);

  char prefix[40];
  
  //dumps dumps to analysis analysis
#if(RADOUTPUT==1)
  sprintf(prefix,"radavg%04d-",no1);
  fprint_radprofiles(t,no2,"analysis",prefix);
#endif
#if(OUTOUTPUT==1)
  sprintf(prefix,"outavg%04d-",no1);
  fprint_outfile(t,no2,0,"analysis",prefix);
#endif
#if(SILOOUTPUT==1)
  sprintf(prefix,"silavg%04d-",no1);
  fprint_silofile(t,no2,"analysis",prefix);
#endif
#if(SIMOUTPUT==1)	  
  sprintf(prefix,"simavg%04d-",no1);
  fprint_simplecart(t,no2,"analysis",prefix);
#endif
  
  return 0;
}

