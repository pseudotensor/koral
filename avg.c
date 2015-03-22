//KORAL - ana.c
//damp files postprocessing

#include "ko.h"

int 
main(int argc, char **argv)
{  
  #ifdef MPI
  printf("avg works on shared memory only, do not use MPI, please\n");
  exit(-1);
  #endif

  #ifdef OMP
  omp_myinit();  
  #endif

  //which files to read
  int no1,no2,nostep,procotg,ifphiavg;
  if(argc!=4 && argc!=7)
    {
      printf("Not enough input arguments. Asks for ./avg no1 no2 nostep [doingavg=1 procotg=0 ifphiavg=0]\n");
      return -1;
    }
  else
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
      if(argc==7)
	{
	  doingavg=atoi(argv[4]);
	  procotg=atoi(argv[5]);      
	  ifphiavg=atoi(argv[6]);      
	}
      else
	{
	  doingavg=1;
	  procotg=0;
	  ifphiavg=0;
	}
    }

  doingpostproc=1;

  int i;

  //currently gsl is not used
  gsl_set_error_handler_off();
  
  //random number gen. initialization
  srand ( time(NULL) );

  //preparing arrays
  initialize_arrays();

  //sets the grid
  set_grid(&min_dx,&min_dy,&min_dz,&max_dt);
  alloc_loops(1,0.,0.);

  //precalculates metric etc.
  calc_metric();

  #ifdef COORDOUTPUT
  fprint_coordfile("analysis","coord");
  #endif

  //folder to write to
  char folder[100];
    //folder to write from
  char folderin[100],bufor[100];
  sprintf(folder,"analysis");

  if(ifphiavg==0)
    sprintf(folderin,"%s","dumps");
  else if(ifphiavg==1)
    sprintf(folderin,"%s","dumps_phiavg");
  else if(ifphiavg==2)
    my_err("using phisli(ced) data with ./avg makes no sense.\n");

  if(procotg)
    {
      //opens the scalar file
      sprintf(bufor,"%s/avgscalars.dat",folder);
      fout_scalars=fopen(bufor,"w");
#if(BOXOUTPUT==1)
      sprintf(bufor,"analysis/avgboxscalars.dat");
      fout_boxscalars=fopen(bufor,"w");
#endif
#if(VAROUTPUT==1)
      sprintf(bufor,"analysis/avgvarscalars.dat");
      fout_varscalars=fopen(bufor,"w");
#endif
    }


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
	  readret=fread_avgfile(ifile,folderin,pavg,&dt,&t);
	}
      else
	{
	  //reading res file
	  readret=fread_restartfile(ifile,folderin,&t);
	  dt=1.;
	  //rewrite primitives to pavg
	  for(iz=0;iz<NZ;iz++)
	    for(iy=0;iy<NY;iy++)
	      for(ix=0;ix<NX;ix++)
		for(iv=0;iv<NV;iv++)
		  set_uavg(pavg,iv,ix,iy,iz,get_u(p,iv,ix,iy,iz));
	  //printf("1 %e %e\n",get_u(p,EE,2*NX/3,NY/2,0),get_uavg(pavg,EE,2*NX/3,NY/2,0));
	}

      global_time=t;

      if(procotg)
	{
	  //calculates scaleheight etc.
	  calc_avgs_throughout();
	  
	  //sets bc
	  set_bc(t,0);
	  
	  //calculate scalars
	  calc_scalars(scalars,t);

	  fprint_scalars(t,scalars,NSCALARS);

	  
#if(BOXOUTPUT==1)
	  fprint_boxscalars(t);
#endif

#if(VAROUTPUT==1)
	  fprint_varscalars(t);
#endif
	}

      add_u_core(1.,pavgtot,dt,pavg,pavgtot,(SX)*(SY)*(SZ)*(NV+NAVGVARS));

	
      ttot+=dt;
    }
  printf("ttot: %f\n",ttot);

  //average primitives and averaged quantities
  copy_u_core(1./ttot,pavgtot,pavg,(SX)*(SY)*(SZ)*(NV+NAVGVARS));


  //avarage of average files
#ifdef AVGAVGOUTPUT
  avgtime=1./ttot;
  sprintf(bufor,"avgavg%04d-",no1);
  nfout2=no2;
  fprint_avgfile(0.,"analysis",bufor);
#endif

  //rewrite primitives to p
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  for(iv=0;iv<NV;iv++)
	    set_u(p,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iz));
	  struct geometry geom; //but what if pavg was written in BL?
	  fill_geometry(ix,iy,iz,&geom);
	  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	}



  //calculates scaleheight etc.
  calc_avgs_throughout();

  //projects on ghost cells
  set_bc(t,0);

  char prefix[40];
  char suffix[10];
  if(doingavg)
    sprintf(suffix,"");
  else
    sprintf(suffix,"res");

  if(ifphiavg)
    sprintf(suffix,"%sphiavg",suffix);
 

  
  //dumps dumps to analysis analysis
#if(RADOUTPUT==1)
  sprintf(prefix,"radavg%s%04d-",suffix,no1);
  fprint_radprofiles(t,no2,"analysis",prefix);
#endif
  
#if(THOUTPUT==1)
  sprintf(prefix,"thavg%s%04d-",suffix,no1);
  fprint_thprofiles(t,no2,"analysis",prefix);
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
  fprint_simplefile(t,no2,"analysis",prefix);
#endif
  
  if(procotg)
    {
#if(BOXOUTPUT==1)
  fclose(fout_boxscalars);
#endif

#if(VAROUTPUT==1)
  fclose(fout_varscalars);
#endif

  fclose(fout_scalars);
    }
  return 0;
}

