//KORAL - ana.c
//damp files postprocessing

#include "ko.h"

int 
main(int argc, char **argv)
{  
  doingavg=0;
  doingpostproc=1;

  #ifdef MPI
  printf("ana works on shared memory only, do not use MPI, please\n");
  exit(-1);
  #endif

  #ifdef OMP
  omp_myinit();  
  #endif

  //which files to read
  int no1,no2,nostep,ifphiavg;
  if(argc<4 || argc>5)
    {
      printf("Not enough input arguments. Asks for ./ana no1 no2 nostep [ifphiavg=0]\n");
      return -1;
    }
  else
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);

      if(argc==5)
	ifphiavg=atof(argv[4]);
      else
	ifphiavg=0;
    }

  char folder[100],bufer[100];
  if(!ifphiavg)
    sprintf(folder,"%s","dumps");
  else
    sprintf(folder,"%s","dumps_phiavg");

  int i;
  
  //no gsl error messages
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


  //precalculating problem related numbers
  /*
    #ifdef PR_PREPINIT
    #include PR_PREPINIT
    #endif
  */

  //getch();

  //opens the scalar file
  sprintf(bufer,"analysis/scalars.dat");
  fout_scalars=fopen(bufer,"w");

  int ifile,itot=0,readret;
  ldouble t,ttot; ldouble scalars[NSCALARS];
  ttot=0.;

  printf("working on files #%04d to #%04d with %d step \n",no1,no2,nostep);

  ldouble pp[NV],uu[NV];

  for(ifile=no1;ifile<=no2;ifile+=nostep)
    {
      itot++;

      //reading restart file
      readret=fread_restartfile(ifile,folder,&t);
      nfout1=ifile;

      //calculates scaleheight etc.
      calc_avgs_throughout();
      
      //sets bc
      #pragma omp parallel
      set_bc(t,0);
     
      //calculate scalars
      calc_scalars(scalars,t);

      //dumps dumps to analysis analysis
      
      char prefix[40];
      char suffix[10];
      sprintf(suffix,"");

      if(ifphiavg)
	sprintf(suffix,"%sphiavg",suffix);
 


#if(SCAOUTPUT==1)
      fprint_scalars(t,scalars,NSCALARS);
#endif
#if(RADOUTPUT==1)
      sprintf(prefix,"rad%s",suffix);  
      fprint_radprofiles(t,nfout1,"analysis",prefix);
#endif

#if(THOUTPUT==1)
      sprintf(prefix,"th%s",suffix);  
      fprint_thprofiles(t,nfout1,"analysis",prefix);
#endif

#if(OUTOUTPUT==1)
      sprintf(prefix,"out%s",suffix);  
      fprint_outfile(t,nfout1,0,"analysis",prefix);
#endif
#if(SILOOUTPUT==1)
#ifndef NOSILO
      sprintf(prefix,"sil%s",suffix);  
      fprint_silofile(t,nfout1,"analysis",prefix);
#endif
#endif
#if(SIMOUTPUT!=0)	  
      sprintf(prefix,"sim%s",suffix);  
      fprint_simplefile(t,nfout1,"analysis",prefix);
#endif
  

    }


  fclose(fout_scalars);
  
  return 0;
}

