//KORAL - ko.c
//radiative hydrodynamical code

#include "ko.h"

int 
main(int argc, char **argv)
{  
  ldouble tstart;
  int i;

  //currently gsl is not used
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

  //precalculates metric etc.
  calc_metric();

  //read restart file
  //TODO: restarting does not work
#ifdef RESTART
  fread_restartfile(&tstart);
#else
  //or initialize new problem
  initialize_problem();
  fprint_openfiles();
  tstart=0.;
#endif

  //sets initial profile of primitives
#ifndef RESTART
  set_initial_profile();
#endif
  //sets bc
  set_bc(tstart);

  //testing
  /*
    ldouble GG[4][5];
    pick_G(10,0,0,GG);
    ldouble gg[4][5];
    pick_g(10,0,0,gg);
	  
    ldouble xx=get_x(10,0);
    printf("x: %f\n",xx);
    print_metric(gg);
    print_metric(GG);
    double ucon[4],ucov[4];
    calc_normalobs_4vel(GG,ucon);
    print_4vector(ucon);
    indices_21(ucon,ucov,gg);
    print_4vector(ucov);
    printf("%e\n",dot(ucon,ucov));
    getchar();
  */


  //evolves
  solve_all_problems_5(tstart);

  //free_arrays();
  fprint_closefiles();
  return 0;
}

/******************************************************/
/***************** time integration ********************/
/******************************************************/

int
solve_all_problems_5(ldouble tstart)
{
  ldouble t = tstart, t1 = TMAX, dt, dtsource, taim;
  ldouble totalmass=0.;
  ldouble dtout = DTOUT1, lasttout = -1.,fprintf_time = 0.;
  int i1,i2,i3,lasttout_floor,ix,iy,iz,iv;
  struct timespec temp_clock;
  struct rad_parameters rp;
   
  i1=i2=0.;

  //prints initial profiles to out0000.dat
#ifndef RESTART
  fprint_profiles(t,totalmass);				
#endif

  lasttout=0.;lasttout_floor=floor(t/dtout); dt=-1.;
  max_ws[0]=max_ws[1]=max_ws[2]=1.;

  //main time loop
  while (t < t1 && nfout1<NOUTSTOP && i1<NSTEPSTOP)
    {    
     //initial time mark
#ifndef SKIP_CLOCK
      clock_gettime(CLOCK_REALTIME,&temp_clock);
#endif
      ldouble start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
      ldouble imp_time1,imp_time2,tstepden;

#ifndef RADIATION //pure hydro
      if(NZ>1)
	tstepden=max_ws[0]/min_dx + max_ws[1]/min_dy + max_ws[2]/min_dz;
      else if(NY>1)
	tstepden=max_ws[0]/min_dx + max_ws[1]/min_dy;
      else
	tstepden=max_ws[0]/min_dx;            
#else //radiation included
      //TODO: what is below assumes wavespeed=1 but for thick flows the real characterstic speed may be lower and time step larger
      if(NZ>1)
	tstepden=(1./min_dx + 1./min_dy + 1./min_dz);
      else if(NY>1)
	tstepden=(1./min_dx + 1./min_dy);
      else 
	tstepden=(1./min_dx);          
#endif
      
      dt=TSTEPLIM*1./tstepden;

      //reseting wavespeeds
      max_ws[0]=-1.;
      max_ws[1]=-1.;
      max_ws[2]=-1.;

      //iteration counter
      i1++;

      //**********************************************************************
      //**********************************************************************
      //**********************************************************************

#ifdef RK2STEPPING
      //******************************* RK2 **********************************
      //1st
      f_timeder (t,dt,1.,u,1,ut0);  
      copy_u(1.,u,ut1);
      //2nd
      f_timeder (t,dt,1.,u,1,ut2); 
      add_u(1.,u,-1.,ut2,ut2);     
      //together     
      t+=dt;    
      add_u(.5,ut0,.5,ut1,u);
      add_u(1.,u,.5,ut2,u);      
     //************************** end of RK2 **********************************
#endif

#ifdef RK3STEPPING
      //******************************* RK3 **********************************
      //TODO : clean up, think it over
      //1st
      copy_u(1.,u,ut0);
      f_timeder (t,dt,1.,u,0,ut0);  
      copy_u(1.,u,ut1);
      //2nd
      copy_u(1.,u,ut2);       
      f_timeder (t,dt,1.,u,0,ut2); 
      add_u(1.,u,-1.,ut2,ut2);   
      add_u(.75,ut0,.25,ut1,u);
      add_u(1.,u,.25,ut2,u);      
      //3rd
      copy_u(1.,u,ut2);
      f_timeder (t,dt,1.,u,0,ut2); 
      add_u(1.,u,-1.,ut2,ut3);   
      //together     
      t+=dt;    
      add_u(1./3.,ut0,2./3.,ut2,u);
      add_u(1.,u,2./3.,ut3,u);      
     //************************** end of RK3 **********************************
#endif

#ifdef RK4STEPPING
      //******************************* RK4 **********************************
      //1st
      f_timeder (t,dt,1.,u,1,ut0);  
      add_u(1.,u,-1.,ut0,ut1);
      add_u(1.,ut0,.5,ut1,u);
      //2nd
      t+=.5*dt;     
      f_timeder (t,dt,1.,u,1,ut2); 
      add_u(1.,u,-1.,ut2,ut2);     
      add_u(1.,ut0,.5,ut2,u);
      //3rd
      f_timeder (t,dt,1.,u,1,ut3); 
      add_u(1.,u,-1.,ut3,ut3);     
      add_u(1.,ut0,1.,ut3,u);
      //4th
      t+=.5*dt;     
      f_timeder (t,dt,1.,u,1,ut4);     
      add_u(1.,u,-1.,ut4,ut4);
      //together     
      add_u(1.,ut0,1./6.,ut1,u);
      add_u(1.,u,2./6.,ut2,u);
      add_u(1.,u,2./6.,ut3,u);
      add_u(1.,u,1./6.,ut4,u);

     //************************** end of RK4 **********************************

#endif


#ifndef SKIP_CLOCK
      //time mark
      clock_gettime(CLOCK_REALTIME,&temp_clock);    
#endif
      ldouble cons_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

      //**********************************************************************
      //**********************************************************************
      //**********************************************************************
   
      //integral of rho
      totalmass= calc_totalmass();

#ifndef SKIP_CLOCK
      //time mark
      clock_gettime(CLOCK_REALTIME,&temp_clock);    
#endif
      ldouble end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
     
      //performance
      ldouble znps=NX*NY*NZ/(end_time-start_time);


      //output to a file
      if(lasttout_floor!=floor(t/dtout) || ALLSTEPSOUTPUT)
	{
	  printf("otpt (no #%6d) at t=%10.3e with dt=%.3e  (%.3f) (real time: %10.4f|%10.4f|%10.4f|%10.4f) mass: %e znps: %f\n",nfout1,t,dt,max_ws[0],
		 cons_time-start_time-imp_time1-imp_time2,imp_time1+imp_time2,end_time-cons_time,end_time-start_time,totalmass,znps);

	  fprint_profiles(t,totalmass);
	  lasttout=t;
	  lasttout_floor=floor(t/dtout);	 
	}
      //or performance to screen only every second
      else if(end_time-fprintf_time>1.) 
	{
	  printf("step (it #%6d) at t=%10.3e with dt=%.3e  (%.3f) (real time: %10.4f|%10.4f|%10.4f|%10.4f) mass: %e znps: %f\n",i1,t,dt,max_ws[0],
		 cons_time-start_time-imp_time1-imp_time2,imp_time1+imp_time2,end_time-cons_time,end_time-start_time,totalmass,znps);
	  fprintf_time=end_time;
	  i2=i1;
	}
     }
  return 0;
}

