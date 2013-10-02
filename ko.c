//KORAL - ko.c
//radiative hydrodynamical code

#include "ko.h"

int 
main(int argc, char **argv)
{  
#if(PROBLEM==60)
  //requires no rad. viscosity!
  test_jon_solve_implicit_lab();
  exit(0);
#endif

  ldouble tstart;
  int i;

  //print scalings GU->CGS
  //print_scalings();
  
  //gsl errors off
  gsl_set_error_handler_off();
  
  //random number gen. initialization
  srand ( time(NULL) );

  //preparing arrays
  initialize_arrays();

  if(argc!=NUM_INPUTARG+1)
    {
      my_err("Not enough input arguments.\n");
      return -1;
    }
  else
    {      
      for (i=0;i<NUM_INPUTARG;i++)
	{
	  inputarg[i]=atof(argv[i+1]);
	  printf("%d: %f\n",i,inputarg[i]);
	}
    }

  //sets the grid
  set_grid(&min_dx,&min_dy,&min_dz,&max_dt);
  //print_grid(min_dx,min_dy,min_dz);

#if(GRIDOUTPUT==1)
  fprint_gridfile("dumps");
#endif

  //precalculates metric etc.
  calc_metric();

  //precalculating problem related numbers
#ifdef PR_PREPINIT
#include PR_PREPINIT
#endif

  int ifinit=1;
#ifdef RESTART
  ifinit=fread_restartfile(RESTARTNUM,&tstart);
  set_bc(tstart,1);
#endif

  //no restart or no restart file
  if(ifinit==1)
    {
      //or initialize new problem
      set_initial_profile();
      tstart=0.;
      set_bc(tstart,1);
#ifdef VECPOTGIVEN
      calc_BfromA();
#endif

#ifdef PR_POSTINIT
#include PR_POSTINIT
#endif
    }

  /*******/
  //tests of implicit solver
  //test_solve_implicit_lab(); exit(0);
  /*******/

  //prepares files
  fprint_openfiles("dumps");
  
  //copies initial primitives to pinit
  copy_u(1.,p,pinit);

  //zeros the avg array
  copy_u(0.,p,pavg);
  avgtime=0.;

  //calculates initial scalars
  calc_scalars(scalars,tstart);

  //prints initial profiles to out0000.dat
  if(ifinit==1)
    {
      fprint_restartfile(tstart,"dumps");			
      //dumps dumps
#if(SCAOUTPUT==1)
      fprint_scalars(tstart,scalars,NSCALARS,"dumps");
#endif
#if(RADOUTPUT==1)
      fprint_radprofiles(tstart,nfout1,"dumps","rad");
#endif
#if(OUTOUTPUT==1)
      fprint_outfile(tstart,nfout1,0,"dumps","out");
#endif
#if(SILOOUTPUT==1)
      fprint_silofile(tstart,nfout1,"dumps","sil");
#endif

      nfout1++;
    }
  
  //evolves
  solve_the_problem(tstart);

  //free_arrays();
  fprint_closefiles();
  return 0;
}

/******************************************************/
/***************** time integration ********************/
/******************************************************/

int
solve_the_problem(ldouble tstart)
{
  ldouble t = tstart, t1 = TMAX, dtsource, taim;
  ldouble totalmass=0.;
  ldouble dtout = DTOUT1, lasttout = -1.,fprintf_time = 0.;
  int i1,i2,i3,lasttout_floor,ix,iy,iz,iv;
  struct timespec temp_clock;
  struct rad_parameters rp;
   
  i1=i2=0.;
  global_int_slot[GLOBALINTSLOT_NCRITFAILURES]=0; //counting number of critical failures

  lasttout=0.;lasttout_floor=floor(t/dtout); dt=-1.;
  max_ws[0]=max_ws[1]=max_ws[2]=1.;
  if(NZ>1)
    tstepdenmax=max_ws[0]/min_dx + max_ws[1]/min_dy + max_ws[2]/min_dz;
  else if(NY>1)
    tstepdenmax=max_ws[0]/min_dx + max_ws[1]/min_dy;
  else
    tstepdenmax=max_ws[0]/min_dx;            

  tstepdenmax*=10.;

  //copy primitives to hold previous time steps
  copy_u(1.,p,ptm1); ttm1=t;
  copy_u(1.,p,ptm2); ttm2=t;

  //main time loop
  int nstep=0;

  while (t < t1 && nfout1<=NOUTSTOP && i1<NSTEPSTOP)
    {    
      global_time=t;
      nstep++;

      //calculates the primitives to copy to previous time steps
      int ii;     

#pragma omp parallel for private(ii,ix,iy,iz,iv) schedule (dynamic)
      for(ii=0;ii<Nloop_0;ii++) //domain only
	{
	  ix=loop_0[ii][0];
	  iy=loop_0[ii][1];
	  iz=loop_0[ii][2]; 
      
	  calc_primitives(ix,iy,iz,0); 
	}

      //holds previous time steps
      copy_u(1.,ptm1,ptm2); ttm2=ttm1;
      copy_u(1.,p,ptm1); ttm1=t;             
      
      //initial time mark
      my_clock_gettime(&temp_clock);

      ldouble start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
      ldouble imp_time1=0.,imp_time2=0.,tstepden;

      dt=TSTEPLIM*1./tstepdenmax;

      if(t+dt>t1) {dt=t1-t;}

      //reseting wavespeeds
      max_ws[0]=-1.;
      max_ws[1]=-1.;
      max_ws[2]=-1.;
      tstepdenmax=-1.;

      //iteration counter
      i1++;

      //**********************************************************************
      //**********************************************************************
      //**********************************************************************
      
     if(TIMESTEPPING==RK2K2)
	{
	  //******************************* RK2 **********************************
	  //1st
	  f_timeder (t,.5*dt,ut0); //updates u
	  //2nd
	  f_timeder (t,dt,ut1); //in ut1 midpoint
	  add_u(1.,u,-1.,ut1,ut2); //k2 in ut2
	  //together     
	  t+=dt;    
	  add_u(1.,ut0,1.,ut2,u);
	  //************************** end of RK2 **********************************
	}
     else if(TIMESTEPPING==RK2K1K2)
	{
	  //******************************* RK2 **********************************
	  //1st
	  f_timeder (t,dt,ut0); //updates u
	  add_u(1.,u,-1.,ut0,ut1); //k1 in ut1
	  //midpoint
	  add_u(1.,ut0,.5,ut1,u);
	  //2nd
	  t+=.5*dt;
	  f_timeder (t,dt,ut2); 
	  add_u(1.,u,-1.,ut2,ut2); //k2 in ut2
	  //together     
	  t+=.5*dt;    
	  add_u(1.,ut0,.5,ut1,u);
	  add_u(1.,u,.5,ut2,u);
	  //************************** end of RK2 **********************************
	}
     else if(TIMESTEPPING==RK3)
	{
	  //******************************* RK3 **********************************
	  //1st
	  f_timeder (t,dt,ut0);  
	  copy_u(1.,u,ut1);
	  //2nd
	  f_timeder (t,dt,ut2); 
	  add_u(1.,u,-1.,ut2,ut2);   
	  add_u(.75,ut0,.25,ut1,u);
	  add_u(1.,u,.25,ut2,u);      
	  //3rd
	  f_timeder (t,dt,ut2); 
	  add_u(1.,u,-1.,ut2,ut3);   
	  //together     
	  t+=dt;    
	  add_u(1./3.,ut0,2./3.,ut2,u);
	  add_u(1.,u,2./3.,ut3,u);      
	  //************************** end of RK3 **********************************/
	}
     else 
       my_err("wrong time stepping specified\n");


      //**********************************************************************
      //************************* finger  ************************************
      //**********************************************************************

      my_finger(t);

      //**********************************************************************
      //************************* outputs ************************************
      //**********************************************************************
   

      //time mark
      my_clock_gettime(&temp_clock);    

      ldouble end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
     
      //performance
      ldouble znps=NX*NY*NZ/(end_time-start_time);

      //average number of iterations in the implicit solver
      ldouble avimpit=global_slot[0]/global_slot[1];

      //save to avg arrays
      save_avg(dt);

      //output to a file
      if(lasttout_floor!=floor(t/dtout) || ALLSTEPSOUTPUT || t>.9999999*t1)
	{
	  printf("otpt (no #%6d) at t=%10.3e with dt=%.3e  (%.3f) (real time: %10.4f) znps: %.0f avimpit: %.3f\n"
		 ,nfout1,t,dt,max_ws[0],end_time-start_time,znps,avimpit);
	  
	  //projects primitives onto ghost cells
	  set_bc(t,0);

	  //calculate scalars
	  calc_scalars(scalars,t);
	  
	  //dumps restart
	  fprint_restartfile(t,"dumps");

	  //dumps dumps

	  //avg goes first so that what is later can use it
#if(AVGOUTPUT==1) 
	  copy_u_core(1./avgtime,pavg,pavg,SX*SY*SZ*(NV+NAVGVARS));
	  fprint_avgfile(t,"dumps");
	  avgtime=0.;
#endif
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
	  
	  nfout1++;

	  lasttout=t;
	  lasttout_floor=floor(t/dtout);	 
	}
      //or performance to screen only every second
      else if(end_time-fprintf_time>1.) 
	{
	  printf("step (it #%6d) at t=%10.3e with dt=%.3e  (%.3f) (real time: %10.4f) znps: %.0f avimpit: %.3f\n"
		 ,nstep,t,dt,max_ws[0],end_time-start_time,znps,avimpit);
	  fprintf_time=end_time;
	  i2=i1;
	}
    }
  return 0;
}

int
test_maginv()
{
  ldouble uu[NV];
  ldouble pp[NV];

  //geometries
  struct geometry geom;
  fill_geometry(0,0,0,&geom);

  print_metric(geom.gg);

  pp[RHO]=1.;
  pp[UU]=0.001;
  pp[VX]=pp[VY]=pp[VZ]=0.;

  pp[VX]=0.0;
  pp[VY]=0.0;
  pp[VZ]=0.0;

#ifdef MAGNFIELD
  pp[B1]=pp[B2]=pp[B3]=0.;

  pp[B1]=0.e-2;
  pp[B2]=0.e-4;
  pp[B3]=1.e-0;
#endif

  //entropy
  pp[5]=calc_Sfromu(pp[RHO],pp[UU]);

  print_NVvector(pp);
  p2u(pp,uu,&geom);
  print_NVvector(uu);

  int aa[2],bb[2],ret;
  pp[UU]*=1.1;
  ret=u2p(uu,pp,&geom,aa,bb,0);
  printf("u2p ret: %d\n",ret);
  print_NVvector(pp);

  return 0;
}

int
print_scalings()
{
  
  printf("scalings (GU->CGS):\nrho: %.16e\nlen: %.16e\ntime: %.16e\nugas:"
	 "%.16e\nT(1,1): %.16e\nkbt: %.16e\nkappa: %.16e\n",
	 rhoGU2CGS(1.),
	 lenGU2CGS(1.),
	 timeGU2CGS(1.),
	 endenGU2CGS(1.),
	 calc_PEQ_Tfromurho(1.,1.),
	 K_BOLTZ/MU_GAS/M_PROTON,
	 kappaCGS2GU(1.)
	 );
  exit(0);
  
  return 0;
}
