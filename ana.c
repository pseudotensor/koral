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

  ldouble tstart;
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


  int ifile;
  ldouble t; ldouble scalars[NSCALARS];
  for(ifile=no1;ifile<=no2;ifile+=nostep)
    {
      fread_dumpfile(ifile,&t);
  
      //sets bc
      set_bc(tstart);
     
      //calculate scalars
      calc_scalars(scalars,t);

      fprint_profiles(t,scalars,NSCALARS,0,folder);
    }

  fclose(fout_scalars);

  return 0;
}

