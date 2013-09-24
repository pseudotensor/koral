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

  int ifile,itot=0,readret;
  ldouble t; ldouble scalars[NSCALARS];

  printf("working on files #%04d to #%04d with %d step \n",no1,no2,nostep);

  for(ifile=no1;ifile<=no2;ifile+=nostep)
    {
      itot++;

      readret=fread_restartfile(ifile,&t);

      //adding up to avg array
      add_u(1.,p,1.,pavg,pavg);

      nfout1--; //correcting index
  
      //sets bc
      set_bc(t,0);
     
      //calculate scalars
      calc_scalars(scalars,t);

      //dumps dumps to analysis analysis
#if(SCAOUTPUT==1)
      fprint_scalars(t,scalars,NSCALARS,"analysis");
#endif
#if(RADOUTPUT==1)
      fprint_radprofiles(t,nfout1,"analysis","rad");
#endif
#if(OUTOUTPUT==1)
      fprint_outfile(t,nfout1,0,"analysis","out");
#endif
#if(SILOOUTPUT==1)
      fprint_silofile(t,nfout1,"analysis","sil");
#endif
#if(SIMOUTPUT==1)	  
      fprint_simplecart(t,nfout1,"analysis","sim");
#endif
  

    }

  fclose(fout_scalars);

  //preparing to dump the avg files
  copy_u(1./(ldouble)itot,pavg,p);
  //projects on ghost cells
  set_bc(t,0);
  //calculating conserved
  int ix,iy,iz,iv;
  ldouble pp[NV],uu[NV];
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      struct geometry geom;
	      fill_geometry(ix,iy,iz,&geom);
	      for(iv=0;iv<NV;iv++)
		pp[iv]=get_u(p,iv,ix,iy,iz);
	      p2u(pp,uu,&geom);
	      for(iv=0;iv<NV;iv++)
		set_u(u,iv,ix,iy,iz,uu[iv]);
	    }
	}
    }

  char prefix[40];
  
  //dumps dumps to analysis analysis
#if(RADOUTPUT==1)
  sprintf(prefix,"radavg%04d-",no1);
  fprint_radprofiles(t,nfout1,"analysis",prefix);
#endif
#if(OUTOUTPUT==1)
  sprintf(prefix,"outavg%04d-",no1);
  fprint_outfile(t,nfout1,0,"analysis",prefix);
#endif
#if(SILOOUTPUT==1)
  sprintf(prefix,"silavg%04d-",no1);
  fprint_silofile(t,nfout1,"analysis",prefix);
#endif
#if(SIMOUTPUT==1)	  
  sprintf(prefix,"simavg%04d-",no1);
  fprint_simplecart(t,nfout1,"analysis",prefix);
#endif


  return 0;
}

