//KORAL - ana.c
//damp files postprocessing

#include "ko.h"

int 
main(int argc, char **argv)
{  
  #ifdef MPI
  printf("ana works on one core only, do not use MPI, please\n");
  exit(-1);
  #endif

  //which files to read
  int no1,no2,nostep;
  if(argc!=4)
    {
      printf("Not enough input arguments. Asks for ./ana no1 no2 nostep\n");
      return -1;
    }
  else
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
    }

  char folder[100],bufer[100];
  sprintf(folder,"%s","dumps");

  int i;
  doingavg=0;
  doingpostproc=1;

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


  //precalculating problem related numbers
#ifdef PR_PREPINIT
#include PR_PREPINIT
#endif

  //opens the scalar file
  sprintf(bufer,"analysis/scalars.dat");
  fout_scalars=fopen(bufer,"w");

  //arrays for averaging of primitives

  //ldouble *pavg=(ldouble*)malloc((SX)*(SY)*(SZ)*(NV+NAVGVARS)*sizeof(ldouble));
  ldouble *pavgtot=(ldouble*)malloc((SX)*(SY)*(SZ)*(NV+NAVGVARS)*sizeof(ldouble));

  for(i=0;i<(SX)*(SY)*(SZ)*(NV+NAVGVARS);i++)
    pavg[i]=pavgtot[i]=0.;

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
#if(SCAOUTPUT==1)
      fprint_scalars(t,scalars,NSCALARS);
#endif
#if(RADOUTPUT==1)
      fprint_radprofiles(t,nfout1,"analysis","rad");
#endif
#if(OUTOUTPUT==1)
      fprint_outfile(t,nfout1,0,"analysis","out");
#endif
#if(SILOOUTPUT==1)
#ifndef NOSILO
      fprint_silofile(t,nfout1,"analysis","sil");
#endif
#endif
#if(SIMOUTPUT!=0)	  
      fprint_simplefile(t,nfout1,"analysis","sim");
#endif
  

    }


  fclose(fout_scalars);
  
  return 0;
}

