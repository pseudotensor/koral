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
      printf("Not enough input arguments. Asks for ./ana no1 no2 nostep\n");
      return -1;
    }
  else
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
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

  sprintf(bufor,"rm %s/*",folder);
  i=system(bufor);

  //opens the scalar file
  sprintf(bufor,"%s/scalars.dat",folder);
  fout_scalars=fopen(bufor,"w");

  //arrays for averaging of primitives

  ldouble *pavg=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
  ldouble *uavg=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
  for(i=0;i<(SX)*(SY)*(SZ)*NV;i++)
    pavg[i]=uavg[i]=0.;

  int ifile,itot=0;
  ldouble t; ldouble scalars[NSCALARS];
  for(ifile=no1;ifile<=no2;ifile+=nostep)
    {
      itot++;

      fread_restartfile(ifile,&t);

      //adding up to avg array
      add_u(1.,p,1.,pavg,pavg);

      nfout1--; //correcting index
  
      //sets bc
      set_bc(t,1);
     
      //calculate scalars
      calc_scalars(scalars,t);

      //dumps dumps
#if(SCAOUTPUT==1)
      fprint_scalars(t,scalars,NSCALARS,"dumps");
#endif
#if(RADOUTPUT==1)
      fprint_radprofiles(t,nfout1,"dumps","rad");
#endif
#if(OUTOUTPUT==1)
      fprint_outfile(t,nfout1,0,"dumps","out");
#endif
#if(SILOOUTPUT==1)
      fprint_silofile(t,nfout1,"dumps","sil");
#endif
#if(SIMOUTPUT==1)	  
      fprint_simplecart(t,folder);
#endif

    }

  fclose(fout_scalars);

  

  return 0;
}

