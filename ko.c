//KORAL - ko.c
//radiative hydrodynamical code

#include "ko.h"

int 
main(int argc, char **argv)
{  
  mpi_myinit(argc,argv);

  //requires no rad. viscosity!
  //test_solve_implicit_lab();
  //exit(0);

  ldouble tstart;
  int i; char folder[100],bufer[100];
  sprintf(folder,"%s","dumps");
  #ifdef MPI
  #ifdef OUTPUTPERCORE
  sprintf(folder,"%s/%d",folder,PROCID);
  sprintf(bufer,"mkdir %s",folder);
  i=system(bufer);
  #endif
  #endif

  //this is not avg.c
  doingavg=0;
  //neither ana.c anarel.c ...
  doingpostproc=0;

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
  fprint_gridfile(folder);
#endif

  //precalculates metric etc.
  calc_metric();

  //prepare angular grid for radiative solver
#ifdef RADIATION
#if(RADCLOSURE==VETCLOSURE)
  zero_readangles();
#endif
#endif

  //**************
  //tests
  //**************
  /*
    ldouble pp[NV]={1.,1.,0.5,0.,0.,-1.,1.,.5,0.,0.};  
    struct geometry geom; ldouble Rij[4][4];
    fill_geometry(0,0,0,&geom);
    calc_Rij(pp,&geom,Rij);
    indices_2221(Rij,Rij,geom.gg);
    print_tensor(Rij);
    exit(-1);
  */
  
  //print scalings GU->CGS and quit
  //print_scalings(); exit(-1);

  //**************
  //end of tests
  //**************

  
  //precalculating problem related numbers
#ifdef PR_PREPINIT
#include PR_PREPINIT
#endif

  int ifinit=1;
#ifdef RESTART
  ifinit=fread_restartfile(RESTARTNUM,folder,&tstart);
  if(!ifinit)
    {
      //exchange initial state
      mpi_exchangedata();  
      calc_avgs_throughout();
      set_bc(tstart,1);
    }
#endif

  //no restart or no restart file
  if(ifinit==1)
    {
      //or initialize new problem
      set_initial_profile();

      tstart=0.;
      //exchange initial state
      if(PROCID==0) {printf("Sending initial data... ");fflush(stdout);}
      mpi_exchangedata();
      calc_avgs_throughout();
      set_bc(tstart,1);
      #ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      if(PROCID==0) {printf("done!\n");fflush(stdout);}
#ifdef VECPOTGIVEN
      if(PROCID==0) {printf("Calculating magn. field... ");fflush(stdout);}
      calc_BfromA(p,1);
      //exchange magn. field calculated in domain
      set_bc(tstart,1);
      mpi_exchangedata();
      calc_avgs_throughout();
      #ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      if(PROCID==0) {printf("done!\n");fflush(stdout);}
#endif

#ifdef PR_POSTINIT
#include PR_POSTINIT
#endif
    }

  //prepares files  
  fprint_openfiles(folder);
  
  //copies initial primitives to pinit
  copy_u(1.,p,pinit);

  //zeros the avg array
  copy_u(0.,p,pavg);
  avgtime=0.;

  //prints initial profiles to out0000.dat
  if(ifinit==1)
    {
      fprint_restartfile(tstart,folder);
			
      //dumps dumps only for shared memory
      #ifndef MPI
#if(SCAOUTPUT==1)
      fprint_scalars(tstart,scalars,NSCALARS);
#endif
#if(RADOUTPUT==1)
      fprint_radprofiles(tstart,nfout1,folder,"rad");
#endif
#if(OUTOUTPUT==1)
      fprint_outfile(tstart,nfout1,0,folder,"out");
#endif
#if(SILOOUTPUT==1)
#ifndef NOSILO
      fprint_silofile(tstart,nfout1,folder,"sil");
#endif
#endif
#if(SIMOUTPUT!=0)
      fprint_simplefile(tstart,nfout1,folder,"sim");

#endif
      #endif

      nfout1++;
    }
    
  //evolves
  solve_the_problem(tstart, folder);

  //free_arrays();
  fprint_closefiles();
  
  mpi_myfinalize();

  return 0;
}

/******************************************************/
/***************** time integration ********************/
/******************************************************/

int
solve_the_problem(ldouble tstart, char* folder)
{
  ldouble t = tstart, t1 = TMAX, dtsource, taim;
  ldouble totalmass=0.;
  ldouble dtout = DTOUT1;
  ldouble dtoutavg = DTOUT2;
  int lasttout_floor;
  int lasttoutavg_floor;
  int i,j;
  
  ldouble fprintf_time = 0.;
  int i1,i2,i3,ix,iy,iz,iv;
  struct timespec temp_clock;
  struct rad_parameters rp;
   
  i1=i2=0.;
  global_int_slot[GLOBALINTSLOT_NTOTALCRITFAILURES]=0; //counting number of critical failures
    
  lasttout_floor=floor(t/dtout); 
  lasttoutavg_floor=floor(t/dtoutavg);
 
  dt=-1.;
  max_ws[0]=max_ws[1]=max_ws[2]=10000.;
  if(NZ>1)
    tstepdenmax=max_ws[0]/min_dx + max_ws[1]/min_dy + max_ws[2]/min_dz;
  else if(NY>1)
    tstepdenmax=max_ws[0]/min_dx + max_ws[1]/min_dy;
  else
    tstepdenmax=max_ws[0]/min_dx;

  //chooses the smalles timestep etc.
  mpi_synchtiming(&t);

  //copy primitives to hold previous time steps
  //copy_u(1.,p,ptm1); ttm1=t;
  //copy_u(1.,p,ptm2); ttm2=t;

  //main time loop
  int nstep=0;

  while (t < t1 && nfout1<=NOUTSTOP && i1<NSTEPSTOP)
    {    
      global_time=t;
      nstep++;
      
      //initial time mark
      my_clock_gettime(&temp_clock);

      start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
      ldouble tstepden;
      
      //chooses the smalles timestep etc.
      mpi_synchtiming(&t);

      //dt based on the estimate from the last midpoint
      dt=TSTEPLIM*1./tstepdenmax;
 
      if(t+dt>t1) {dt=t1-t;}

      //reseting wavespeeds
      max_ws[0]=-1.;
      max_ws[1]=-1.;
      max_ws[2]=-1.;
      max_ws_ph=-1.;
      tstepdenmax=-1.;

      //iteration counter
      i1++;

      //**********************************************************************
      //**********************************************************************
      //**********************************************************************
    
      if(TIMESTEPPING==RK2IMEX)
	{
	  ldouble gamma=1.-1./sqrt(2.);
	  op_implicit (t,dt*gamma,ut0); //U(n) in *ut0;  U(1) in *u
	  add_u(1./(dt*gamma),u,-1./(dt*gamma),ut0,drt1); //R(U(1)) in *drt1;
	  op_explicit (t,dt,ut1); //U(1) in *ut1; 
	  add_u(1./dt,u,-1./dt,ut1,dut1); //F(U(1)) in *dut1;
	  add_u_3(1.,ut0,dt,dut1,dt*(1.-2.*gamma),drt1,u); //(U(n) + dt F(U(1)) + dt (1-2gamma) R(U(1))) in *u
	  op_implicit (t,gamma*dt,uforget); //U(2) in *u
	  add_u(1./(dt*gamma),u,-1./(dt*gamma),uforget,drt2); //R(U(2)) in *drt2;
	  op_explicit (t,dt,ut2); //U(2) in *ut2; 
	  add_u(1./dt,u,-1./dt,ut2,dut2); //F(U(2)) in *dut2;
	  add_u_3(1.,ut0,dt/2.,dut1,dt/2.,dut2,u); //U(n) + dt/2 (F(U(1)) + F(U(2))) in *u
	  add_u_3(1.,u,dt/2.,drt1,dt/2.,drt2,u); //u += dt/2 (R(U(1)) + R(U(2))) in *u
	  t+=dt;	 
	}
     else if(TIMESTEPPING==RK2)
       { 
	 //******************************* RK2 **********************************
	 //1st
	 op_explicit (t,.5*dt,ut0); 
	 op_implicit (t,.5*dt,ut3); 
	 //2nd
	 op_explicit (t,dt,ut1); 
	 op_implicit (t,dt,ut3); 
	 add_u(1.,u,-1.,ut1,ut2); //k2 in ut2
	 //together     
	 t+=dt;    
	 add_u(1.,ut0,1.,ut2,u);
	 //************************** end of RK2 **********************************
	}
     else if(TIMESTEPPING==RK2HEUN)
       { 
	 //******************************* RK2 **********************************
	 //1st
	 //todo:
	 op_explicit (t,1.*dt,ut0); 
	 op_implicit (t,1.*dt,uforget); 
	 add_u(1.,u,-1.,ut0,ut2); 
	 //2nd
	 op_explicit (t,dt,ut1); 
	 op_implicit (t,dt,uforget); 
	 add_u(1.,u,-1.,ut1,ut3); 
	 //together     
	 t+=dt;    
	 add_u_3(1.,u,1./2.,ut2,1./2.,ut3,u); //u += dt/2 (R(U(1)) + R(U(2))) in *u
	 //************************** end of RK2 **********************************
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
   

      #ifdef RADIATION
      //average number of iterations in the implicit solver
      ldouble avimpitloc[5],avimpit[5];
      int impnumsloc[7],impnums[7];
      
      avimpitloc[0]=global_int_slot[GLOBALINTSLOT_NIMPENERMHD]==0 ? 0. : (ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENERMHD]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERMHD];
      avimpitloc[1]=global_int_slot[GLOBALINTSLOT_NIMPENERRAD]==0 ? 0. : (ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENERRAD]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERRAD];
      avimpitloc[2]=global_int_slot[GLOBALINTSLOT_NIMPENTRMHD]==0 ? 0. : (ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENTRMHD]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENTRMHD];
      avimpitloc[3]=global_int_slot[GLOBALINTSLOT_NIMPENTRRAD]==0 ? 0. : (ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENTRRAD]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENTRRAD];
      avimpitloc[4]=global_int_slot[GLOBALINTSLOT_NIMPLTE]==0 ? 0. : (ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPLTE]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPLTE];
            
      impnumsloc[0]=global_int_slot[GLOBALINTSLOT_NIMPENERMHD];
      impnumsloc[1]=global_int_slot[GLOBALINTSLOT_NIMPENERRAD];
      impnumsloc[2]=global_int_slot[GLOBALINTSLOT_NIMPENTRMHD];
      impnumsloc[3]=global_int_slot[GLOBALINTSLOT_NIMPENTRRAD];
      impnumsloc[4]=global_int_slot[GLOBALINTSLOT_NIMPLTE];
      impnumsloc[5]=global_int_slot[GLOBALINTSLOT_NRADFIXUPS];
      impnumsloc[6]=global_int_slot[GLOBALINTSLOT_NCRITFAILURES]; 

      #ifdef MPI
      MPI_Allreduce(impnumsloc, impnums, 7, MPI_INT, MPI_SUM,
		    MPI_COMM_WORLD);  
      MPI_Allreduce(avimpitloc, avimpit, 5, MPI_LDOUBLE, MPI_MAX,
                   MPI_COMM_WORLD);  
      #else
      for(i=0;i<5;i++) avimpit[i]=avimpitloc[i];
      for(i=0;i<7;i++) impnums[i]=impnumsloc[i];
      #endif
      #endif
      
      //save to avg arrays
      save_avg(dt);

      //time mark
      my_clock_gettime(&temp_clock);    

      end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

      //performance
      ldouble znps=TNX*TNY*TNZ/(end_time-start_time);

      //avg files
#if(AVGOUTPUT==1) 
      if(lasttoutavg_floor!=floor(t/dtoutavg))
	{
	  if(PROCID==0)
	    printf("%d > avg file no #%6d dumped\n",PROCID,nfout2);
	  
	  //avg goes first so that what is later can use it
	  copy_u_core(1./avgtime,pavg,pavg,SX*SY*SZ*(NV+NAVGVARS));
	  fprint_avgfile(t,folder,"avg");
	  avgtime=0.;
  
	  nfout2++;

	  lasttoutavg_floor=floor(t/dtoutavg);	 
	}
#endif

      //snapshots
      if(lasttout_floor!=floor(t/dtout) || ALLSTEPSOUTPUT || t>.9999999*t1)
	{
	  if(PROCID==0)
	    printf("%d > snap file no #%6d dumped\n",PROCID,nfout1);
	  
	  //projects primitives onto ghost cells
	  set_bc(t,0);

	  //restart files
	  fprint_restartfile(t,folder);

	  //dumps dumps
	  #ifndef MPI //comment if you want .silo etc files per core for OUTPUTPERCORE
#if(SCAOUTPUT==1)
	  fprint_scalars(t,scalars,NSCALARS);
#endif
#if(RADOUTPUT==1)
	  fprint_radprofiles(t,nfout1,folder,"rad");
#endif
#if(OUTOUTPUT==1)
	  fprint_outfile(t,nfout1,0,folder,"out");
#endif
#if(SILOOUTPUT==1)
#ifndef NOSILO
	  fprint_silofile(t,nfout1,folder,"sil");
#endif
#endif
#if(SIMOUTPUT!=0)
	   fprint_simplefile(tstart,nfout1,folder,"sim");
#endif
	  #endif
	  
	  nfout1++;

	  lasttout_floor=floor(t/dtout);	 
	}


      //performance to screen only every second
      if(end_time-fprintf_time>1. && PROCID==0) 
	{
	  printf("(%d) step #%6d t=%10.3e dt=%.3e (real time: %7.2f | %7.6f) znps: %.0f "
		 ,PROCID,nstep,t,dt,end_time-start_time,2*maxmp_time,znps);
#ifdef RADIATION
	  printf("#:%d %d %d %d %d %d %d | %.1f %.1f %.1f %.1f %.1f\n",
		 impnums[0],impnums[1],impnums[2],impnums[3],impnums[4],impnums[5],impnums[6],
		 avimpit[0],avimpit[1],avimpit[2],avimpit[3],avimpit[4]);
#else
	  printf("\n");
#endif

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
  
  printf("BH mass: %.6f\nspin: %.6f\n\nscalings  (GU->CGS):\nrho: %.16e\nmdot: %.16e\nsigma: %.16e\nlen: %.16e\ntime: %.16e\nenden:"
	 "%.16e\nT(1,1): %.16e\nkbt: %.16e\nkappa: %.16e\n\nrhorizonBL: %.6f\nrISCOBL: %.6f\netaNT: %.6f\n\nmdotEdd: %.16e\n",
	 MASS,BHSPIN,
	 rhoGU2CGS(1.),
	 rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.),
	 rhoGU2CGS(1.)*lenGU2CGS(1.),
	 lenGU2CGS(1.),
	 timeGU2CGS(1.),
	 endenGU2CGS(1.),
	 calc_PEQ_Tfromurho(1.,1.),
	 K_BOLTZ/MU_GAS/M_PROTON,
	 kappaCGS2GU(1.),
	 rhorizonBL,
	 rISCOBL,
	 etaNT,
	 rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)/calc_mdotEdd()
	 );
  
  return 0;
}
