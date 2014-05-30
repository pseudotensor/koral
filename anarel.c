//KORAL - ana.c
//damp files postprocessing

#include "ko.h"

int 
main(int argc, char **argv)
{  
  #ifdef MPI
  printf("anarel works on one core only, do not use MPI, please\n");
  exit(-1);
  #endif

 int i,j;

   //which files to read
  int no1res,no2res,nostepres;
  int no1avg,no2avg,nostepavg;
  if(argc!=7)
    {
      printf("Not enough input arguments. Asks for ./anarel (no1 no2 nostep)_avg (no1 no2 nostep)_res\n");
      return -1;
    }
  else
    {      
      no1avg=atof(argv[1]);
      no2avg=atof(argv[2]);
      nostepavg=atof(argv[3]);
      no1res=atof(argv[4]);
      no2res=atof(argv[5]);
      nostepres=atof(argv[6]);
    }

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

  #ifdef COORDOUTPUT
  fprint_coordfile("analysis","coord");
  #endif

  //folder to write in
  char folder[100],bufor[100];
  sprintf(folder,"analysis");

  /*********************/
  //over avg files first
  /*********************/

  //arrays for averaging of primitives
  ldouble *pavgtot=(ldouble*)malloc((SX)*(SY)*(SZ)*(NV+NAVGVARS)*sizeof(ldouble));
  for(i=0;i<(SX)*(SY)*(SZ)*(NV+NAVGVARS);i++)
    pavg[i]=pavgtot[i]=0.;

  printf("using avg files #%04d to #%04d with %d step \n",no1avg,no2avg,nostepavg);
  
  int ifile,readret;
  ldouble t,ttot;
  ttot=0.;
  t=global_time;

  doingpostproc=1;
  int ix,iy,iz,iv;
  ldouble pp[NV],uu[NV];
  
  doingavg=1;
  for(ifile=no1avg;ifile<=no2avg;ifile+=nostepavg)
    {
      //reading avg file
      readret=fread_avgfile(ifile,"dumps",pavg,&dt,&t);

      add_u_core(1.,pavgtot,dt,pavg,pavgtot,(SX)*(SY)*(SZ)*(NV+NAVGVARS));
	
      ttot+=dt;
    }

  //average primitives and averaged quantities
  copy_u_core(1./ttot,pavgtot,pavg,(SX)*(SY)*(SZ)*(NV+NAVGVARS));
 
  /*********************/
  //at this point pavg holds averaged quantities
  //over snap files now
  /*********************/
  ldouble profiles[NANARELRADPROFILES][NX];      
  ldouble profilesavg[NANARELRADPROFILES][NX];   
  ldouble tprev=-1.;

  for(i=0;i<NANARELRADPROFILES;i++)
    for(j=0;j<NX;j++)
      profilesavg[i][j]=0.;

  doingavg=0;ttot=0.;
  for(ifile=no1res;ifile<=no2res;ifile+=nostepres)
    {
      //reading res file
      readret=fread_restartfile(ifile,"dumps",&t);
      nfout1=ifile;

      if(tprev>0.)
	{
	  dt=t-tprev;
	  
	  for(i=0;i<NANARELRADPROFILES;i++)
	    for(j=0;j<NX;j++)
	      profilesavg[i][j]+=profiles[i][j]*dt;

	  ttot+=dt;
	  
	}
      
      //sets bc
      set_bc(t,0); 
						
      //calculates scale-height
      calc_avgs_throughout();      
      
      //calculates and dumps 
#if(ANARELRADOUTPUT==1)
      calc_anarelradialprofiles(profiles);
      fprint_anarelradprofiles(t,nfout1,"analysis","relrad",profiles);
#endif

      tprev=t;
    }
  //the last one ignored

  for(i=0;i<NANARELRADPROFILES;i++)
    for(j=0;j<NX;j++)
      profilesavg[i][j]/=ttot;

  //averaged dumps
  char prefix[40];
  char suffix[10];
  sprintf(suffix,"");
#if(ANARELRADOUTPUT==1)
  sprintf(prefix,"relradavg%s%04d-",suffix,no1res);
  fprint_anarelradprofiles(t,no2res,"analysis",prefix,profilesavg);
#endif

  
  return 0;
}

