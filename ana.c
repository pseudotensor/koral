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
  if(ifphiavg==0)
    sprintf(folder,"%s","dumps");
  else if(ifphiavg==1)
    sprintf(folder,"%s","dumps_phiavg");
  else if(ifphiavg==2)
    sprintf(folder,"%s","dumps_phisli");


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

  //opens the scalars file
  sprintf(bufer,"analysis/scalars.dat");
  fout_scalars=fopen(bufer,"w");

  if(ifphiavg<2)
    {
#if(BOXOUTPUT==1)
  sprintf(bufer,"analysis/boxscalars.dat");
  fout_boxscalars=fopen(bufer,"w");
  #endif
    }

  if(ifphiavg==2)
    {
 #if(VAROUTPUT==1)
  sprintf(bufer,"analysis/varscalars.dat");
  fout_varscalars=fopen(bufer,"w");
  #endif
    }


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
      set_bc(t,0);
     
      //calculate scalars
      calc_scalars(scalars,t);
     

      //dumps dumps to analysis analysis
      
      char prefix[40];
      char suffix[10];
      sprintf(suffix,"");

      if(ifphiavg==1)
	sprintf(suffix,"%sphiavg",suffix);
      if(ifphiavg==2)
	sprintf(suffix,"%sphisli",suffix);
 
      if(ifphiavg==2) //phisliced - only these below make sense for phi-slices
	{
#if(VAROUTPUT==1)
      fprint_varscalars(t);
#endif
#if(SILOOUTPUT==1)
#ifndef NOSILO
      sprintf(prefix,"sil%s",suffix);  
      fprint_silofile(t,nfout1,"analysis",prefix);
#endif
#endif
	}
      else
	{ //regular

#if(BOXOUTPUT==1)
      fprint_boxscalars(t);
#endif


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
  

    }

#if(BOXOUTPUT==1)
  fclose(fout_boxscalars);
#endif

#if(BOXOUTPUT==1)
  fclose(fout_varscalars);
#endif

  fclose(fout_scalars);
  
  return 0;
}

