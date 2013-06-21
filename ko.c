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
  //print_grid(min_dx,min_dy,min_dz);

  //precalculates metric etc.
  calc_metric();


#ifdef RESTART
  fread_dumpfile(RESTARTNUM,&tstart);
#else
  //or initialize new problem
  set_initial_profile();
  my_finger(0.);
  tstart=0.;
#endif

  //prepares files
  fprint_openfiles("dumps");

  //sets bc
  set_bc(tstart);

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

  ldouble scalars[NSCALARS];
  calc_scalars(scalars,tstart);

  //prints initial profiles to out0000.dat
#ifndef RESTART
  fprint_profiles(t,scalars,NSCALARS,1,"dumps");			
#endif

  lasttout=0.;lasttout_floor=floor(t/dtout); dt=-1.;
  max_ws[0]=max_ws[1]=max_ws[2]=1.;

  //copy primitives to hold previous time steps
  //copy_u(1.,p,ptm1); ttm1=t;
  //copy_u(1.,p,ptm2); ttm2=t;

  //main time loop
  int nstep=0;
  while (t < t1 && nfout1<NOUTSTOP && i1<NSTEPSTOP)
    {    
      nstep++;
      //calculates the primitives to copy to previous time steps
      /*
      int ii;
#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
      for(ii=0;ii<Nloop_0;ii++) //domain only
	{
	  ix=loop_0[ii][0];
	  iy=loop_0[ii][1];
	  iz=loop_0[ii][2]; 
      
	  calc_primitives(ix,iy,iz);
	}
      //holds in previous time steps
      copy_u(1.,ptm1,ptm2); ttm2=ttm1;
      copy_u(1.,p,ptm1); ttm1=t;       
      */
      
      //initial time mark
#ifndef SKIP_CLOCK
      clock_gettime(CLOCK_REALTIME,&temp_clock);
#endif
      ldouble start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
      ldouble imp_time1=0.,imp_time2=0.,tstepden;

      //#ifndef RADIATION //pure hydro
      if(NZ>1)
	tstepden=max_ws[0]/min_dx + max_ws[1]/min_dy + max_ws[2]/min_dz;
      else if(NY>1)
	tstepden=max_ws[0]/min_dx + max_ws[1]/min_dy;
      else
	tstepden=max_ws[0]/min_dx;            
      /*
#else //radiation included
      //TODO: what is below assumes wavespeed=1 but for thick flows the real characterstic speed may be lower and time step larger
      if(NZ>1)
	tstepden=(1./min_dx + 1./min_dy + 1./min_dz);
      else if(NY>1)
	tstepden=(1./min_dx + 1./min_dy);
      else 
	tstepden=(1./min_dx);          
#endif
      */

      dt=TSTEPLIM*1./tstepden;

      if(t+dt>t1) {dt=t1-t;}

      //reseting wavespeeds
      max_ws[0]=-1.;
      max_ws[1]=-1.;
      max_ws[2]=-1.;

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
   

#ifndef SKIP_CLOCK
      //time mark
      clock_gettime(CLOCK_REALTIME,&temp_clock);    
#endif
      ldouble end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
     
      //performance
      ldouble znps=NX*NY*NZ/(end_time-start_time);


      //output to a file
      if(lasttout_floor!=floor(t/dtout) || ALLSTEPSOUTPUT || t>.9999999*t1)
	{
	  printf("otpt (no #%6d) at t=%10.3e with dt=%.3e  (%.3f) (real time: %10.4f) mass: %e znps: %f\n"
		 ,nfout1,t,dt,max_ws[0],end_time-start_time,totalmass,znps);
	  
	  //projects primitives onto ghost cells
	  set_bc(t);

	  //calculate scalars
	  calc_scalars(scalars,t);

	  fprint_profiles(t,scalars,NSCALARS,1,"dumps");
	  lasttout=t;
	  lasttout_floor=floor(t/dtout);	 
	}
      //or performance to screen only every second
      else if(end_time-fprintf_time>1.) 
	{
	  printf("step (it #%6d) at t=%10.3e with dt=%.3e  (%.3f) (real time: %10.4f) mass: %e znps: %f\n"
		  ,nstep,t,dt,max_ws[0],end_time-start_time,totalmass,znps);
	  fprintf_time=end_time;
	  i2=i1;
	}
     }
  return 0;
}

