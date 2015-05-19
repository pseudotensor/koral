//KORAL - ko.c

//radiative hydrodynamical code

#include "ko.h"

int 
main(int argc, char **argv)
{  
  #ifdef MPI
  mpi_myinit(argc,argv);
  #else
  omp_myinit();  
  #endif
  mstep_init();

  ldouble tstart;
  int i,j,k; char folder[100],bufer[100];

  //dumping to
  sprintf(folder,"%s","dumps");
  //this is not avg.c
  doingavg=0;
  //neither ana.c anarel.c ...
  doingpostproc=0;
  
  global_time=0.;

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

  alloc_loops(1,0.,0.);

#if(GRIDOUTPUT==1)
  fprint_gridfile(folder);
#endif

  //precalculates metric etc.
  calc_metric();

#ifdef RADIATION
  //prepare angular grid for radiative solver
#if(RADCLOSURE==VETCLOSURE)
  zero_init();						
  //ZEROtest_oldmain();
#endif

  //prepare arrays for accelerating radiative viscosity
#if(RADVISCOSITY==SHEARVISCOSITY)
  reset_radviscaccel();
#endif
#endif

  //**************
  //tests
  //**************

  /*
  ldouble vec[4]={0.,1.,1./10.,1./10.};
  ldouble xx[4]={0.,10.,3*M_PI/4.,0.}; 
  printf("> %f %f %f\n",xx[1],xx[2],xx[3]);
  print_4vector(vec);
  trans2_coco(xx,vec,vec,SPHCOORDS,CYLCOORDS);
  print_4vector(vec);
  coco_N(xx,xx,SPHCOORDS,CYLCOORDS);
  printf("> %f %f %f\n",xx[1],xx[2],xx[3]);
  trans2_coco(xx,vec,vec,CYLCOORDS,SPHCOORDS);
  print_4vector(vec);
  coco_N(xx,xx,CYLCOORDS,SPHCOORDS);
  printf("> %f %f %f\n",xx[1],xx[2],xx[3]);
 
  exit(1);
  */

  //test_metric(); exit(1);

  //test_inversion_nonrel(); exit(1);

  //test_solve_implicit_lab(); exit(1);
  //test_Giff();  exit(-1);

  //print scalings GU->CGS and quit
  if(PROCID==0) print_scalings();

  //test_opacities(); exit(1);

  //**************
  //end of tests
  //**************

 

  int ifinit=1;
#ifdef RESTART
  ifinit=fread_restartfile(RESTARTNUM,folder,&tstart);
  if(!ifinit) global_time=tstart;
  //todo: read intensities from file!
#if (RADCLOSURE==VETCLOSURE)
#ifdef RADSTARTWITHM1INTENSITIES
  calc_M1intensities();
#endif
#endif

 
  //precalculating problem related numbers
#ifdef PR_PREPINIT
#include PR_PREPINIT
#endif

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
#ifdef MAGNFIELD
#ifdef VECPOTGIVEN
      if(PROCID==0) {printf("Calculating magn. field... ");fflush(stdout);}
      calc_BfromA(p,1);
      //exchange magn. field calculated in domain
      mpi_exchangedata();
      calc_avgs_throughout();
      set_bc(tstart,1);
       #ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      if(PROCID==0) {printf("done!\n");fflush(stdout);}
#endif
#endif
#ifdef RADSTARTWITHM1INTENSITIES
      calc_M1intensities();
      mpi_exchangedata();
      calc_avgs_throughout();
      set_bc(tstart,1);
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
  copy_u_core(0.,pavg,pavg,SX*SY*SZ*(NV+NAVGVARS));	  
  avgtime=0.;
  #ifdef SELFTIMESTEP
  int iz,iy,ix;
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	set_u_scalar(avgselftime,ix,iy,iz,0.);
  #endif
   

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
  ldouble dtoutbox = DTOUT3;
  ldouble dtoutvar = DTOUT4;
  ldouble lasttout_floor;
  ldouble lasttoutavg_floor;
  ldouble lasttoutbox_floor;
  ldouble lasttoutvar_floor;
  int i,j,ii;
  int loopsallociter;
  int spitoutput,lastzone;
  int nentr[4],nentr2[4];
  
  ldouble fprintf_time = 0.;
  int fprintf_nstep=0;
  int i1,i2,i3,ix,iy,iz,iv;
  struct timespec temp_clock;
  struct rad_parameters rp;
   
 
  i1=i2=0.;
  global_int_slot[GLOBALINTSLOT_NTOTALCRITFAILURES]=0; //counting number of critical failures
  global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]=0; //and fixups requests
  global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]=0; //and fixups requests


  lasttout_floor=floor(t/dtout); 
  lasttoutavg_floor=floor(t/dtoutavg);
  lasttoutbox_floor=floor(t/dtoutbox);
  lasttoutvar_floor=floor(t/dtoutvar);
 
  dt=-1.;
  max_ws[0]=max_ws[1]=max_ws[2]=10000.;
  if(NZ>1)
    tstepdenmax=max_ws[0]/min_dx + max_ws[1]/min_dy + max_ws[2]/min_dz;
  else if(NY>1)
    tstepdenmax=max_ws[0]/min_dx + max_ws[1]/min_dy;
  else
    tstepdenmax=max_ws[0]/min_dx;

  tstepdenmax/=TSTEPLIM;
  tstepdenmin=tstepdenmax;

  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      int ix,iy,iz;
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      set_u_scalar(cell_tstepden,ix,iy,iz,tstepdenmax);
      set_u_scalar(cell_dt,ix,iy,iz,1./tstepdenmax);
    }

  //chooses the smalles timestep etc.
  mpi_synchtiming(&t);

  //copy primitives to hold previous time steps
  //copy_u(1.,p,ptm1); ttm1=t;
  //copy_u(1.,p,ptm2); ttm2=t;

  //main time loop
  nstep=0;



  #ifdef SUBZONES
  currentzone=-1;
  currentzone=alloc_loops(1,t,dt);
  loopsallociter=0;
  currentzonetime=t;

  set_bc(t,0);
  int ii,jj;
  for(ii=0;ii<Nloop_02;ii++) //domain + gc
    {
      ix=loop_02[ii][0];
      iy=loop_02[ii][1];
      iz=loop_02[ii][2]; 
      PLOOP(jj)
      {
	set_u(u_bak_subzone,jj,ix,iy,iz,get_u(u,jj,ix,iy,iz));
	set_u(p_bak_subzone,jj,ix,iy,iz,get_u(p,jj,ix,iy,iz));
      }
    }
  
#endif
 
  while (t < t1 && nfout1<=NOUTSTOP && nstep<NSTEPSTOP)
    {   
      #ifdef METRICTIMEDEPENDENT
      calc_metric();
      #endif
 
      spitoutput=0;
      global_time=t;
      nstep++;

      //verify if broken
      if(global_int_slot[GLOBALINTSLOT_NTOTALCRITFAILURES]>1.e3 ||
	 global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]>1.e3 ||
	 global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]>1.e3)
	{
	  /*
	  printf("exceeded # of failures (%d %d %d) - exiting.\n",
		 global_int_slot[GLOBALINTSLOT_NTOTALCRITFAILURES],global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS],global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]);
		 exit(-1);*/
	}

      
      //chooses the smalles timestep etc.
      mpi_synchtiming(&t);

      //initial time mark
      my_clock_gettime(&temp_clock);

      start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
      ldouble tstepden;
      
      //reallocate loops to allow for sub-zones, but don't do it every step
      #ifdef SUBZONES

      loopsallociter++;
      if(loopsallociter>SUBZONES_NSTEPSTEP)
	{
	  lastzone=currentzone;
	  currentzone=alloc_loops(0,t,dt);
	  if(lastzone!=currentzone) //zone switched
	    {
	      //recomputes the timestep though wavespeeds
	      for(ii=0;ii<Nloop_1;ii++) //domain plus some ghost cells
		{
		  ix=loop_1[ii][0]; iy=loop_1[ii][1]; iz=loop_1[ii][2]; ldouble aaa[12],max_lws[3];
		  calc_wavespeeds_lr(ix,iy,iz,aaa); save_wavespeeds(ix,iy,iz,aaa,max_lws);
		}
              #ifdef OUTPUTAFTERSUBZONES
	      if(lastzone!=currentzone) spitoutput=1;
	      #endif	  
	    }
	  //if(lastzone==1 && currentzone!=1) spitoutput=1; //finished with the innermost one, output
	  loopsallociter=0;
	}
      #endif



      //dt based on the estimate from the last timestep
      dt=1./tstepdenmax;
#ifdef SELFTIMESTEP
      dt=1./tstepdenmin;
#endif
      global_dt=dt;

#ifdef MSTEP
      //adjusts the time levels if possible
      if(nstep>1)
	{
	  mstep_iterate();
	  mstep_update_levels();
	  //reallocate loops to adjust to the set of cells evolved at this time step
	  alloc_loops(0,t,dt);
 	  //mstep_print_state();
	}
#endif

      if(t+dt>t1) {dt=t1-t;}

   
      //reseting wavespeeds
      max_ws[0]=-1.;
      max_ws[1]=-1.;
      max_ws[2]=-1.;
      max_ws_ph=-1.;
      tstepdenmax=-1.;
      tstepdenmin=BIG;


      #if (RADVISCOSITY==SHEARVISCOSITY)
      calc_Rij_visc_total();
      #endif


      //**********************************************************************
      //**********************************************************************
      //**********************************************************************

      if(TIMESTEPPING==-100) //skip evolution completely
	{
	  save_timesteps(); 
	  t+=dt;
	}
      else
	if(TIMESTEPPING==RK2IMEX)
	  {
	    ldouble gamma=1.-1./sqrt(2.),dtcell;
	    save_timesteps(); 

	    dtcell=dt;

	    /******* 1st implicit **********/
	    copyi_u(1.,u,ut0);
	    count_entropy(&nentr[0],&nentr2[0]); copy_entropycount(); 
	    op_implicit (t,dt*gamma); //U(n) in *ut0;  U(1) in *u	  
	    
	    for(ii=0;ii<Nloop_0;ii++) { ix=loop_0[ii][0];      iy=loop_0[ii][1];      iz=loop_0[ii][2];
#ifdef SELFTIMESTEP
	      dtcell=get_u_scalar(cell_dt,ix,iy,iz);
#endif	   
	      PLOOP(iv) set_u(drt1,iv,ix,iy,iz,(1./(dtcell*gamma))*get_u(u,iv,ix,iy,iz)+(-1./(dtcell*gamma))*get_u(ut0,iv,ix,iy,iz)); 
	    } 
	    //addi_u(1./(dt*gamma),u,-1./(dt*gamma),ut0,drt1); //R(U(1)) in *drt1;

	    /******* 1st explicit **********/
	    copyi_u(1.,u,ut1);
      
	    calc_u2p();
#pragma omp barrier
	    do_correct();
#pragma omp barrier

	    count_entropy(&nentr[1],&nentr2[1]); 

	    op_explicit (t,dt); //U(1) in *ut1; 

	    for(ii=0;ii<Nloop_0;ii++)  { ix=loop_0[ii][0];      iy=loop_0[ii][1];      iz=loop_0[ii][2];
#ifdef SELFTIMESTEP
	      dtcell=get_u_scalar(cell_dt,ix,iy,iz);
#endif
	      PLOOP(iv) set_u(dut1,iv,ix,iy,iz,(1./(dtcell))*get_u(u,iv,ix,iy,iz)+(-1./(dtcell))*get_u(ut1,iv,ix,iy,iz)); 
	    }
	    //addi_u(1./dt,u,-1./dt,ut1,dut1); //F(U(1)) in *dut1;


	    /******* 1st together **********/
	    for(ii=0;ii<Nloop_0;ii++)  { ix=loop_0[ii][0];      iy=loop_0[ii][1];      iz=loop_0[ii][2];
#ifdef SELFTIMESTEP
	      dtcell=get_u_scalar(cell_dt,ix,iy,iz);
#endif
	      PLOOP(iv) set_u(u,iv,ix,iy,iz,get_u(ut0,iv,ix,iy,iz)+(dtcell)*get_u(dut1,iv,ix,iy,iz)+(dtcell*(1.-2.*gamma))*get_u(drt1,iv,ix,iy,iz)); 
	    }
	    //addi_u_3(1.,ut0,dt,dut1,dt*(1.-2.*gamma),drt1,u); //(U(n) + dt F(U(1)) + dt (1-2gamma) R(U(1))) in *u
	    
	    /******* 2nd implicit **********/
	    copyi_u(1.,u,uforget);
	    calc_u2p();
#pragma omp barrier
	    do_correct();
#pragma omp barrier

	    count_entropy(&nentr[2],&nentr2[2]); 
	    op_implicit (t,gamma*dt); //U(2) in *u

	    	    #if(AVGOUTPUT==1) 
	    //save to avg arrays	    
	    save_avg(dt);
	    #endif


	    for(ii=0;ii<Nloop_0;ii++) { ix=loop_0[ii][0];      iy=loop_0[ii][1];      iz=loop_0[ii][2];
#ifdef SELFTIMESTEP
	      dtcell=get_u_scalar(cell_dt,ix,iy,iz);
#endif
	      PLOOP(iv) set_u(drt2,iv,ix,iy,iz,(1./(dtcell*gamma))*get_u(u,iv,ix,iy,iz)+(-1./(dtcell*gamma))*get_u(uforget,iv,ix,iy,iz)); 
	    }
	    //addi_u(1./(dt*gamma),u,-1./(dt*gamma),uforget,drt2); //R(U(2)) in *drt2;
	    
	    /******* 2nd explicit **********/
	    copyi_u(1.,u,ut2);
	    calc_u2p();
#pragma omp barrier
	    do_correct();
#pragma omp barrier

	    count_entropy(&nentr[3],&nentr2[3]);
	    op_explicit (t,dt); //U(2) in *ut2; 

	    for(ii=0;ii<Nloop_0;ii++)  { ix=loop_0[ii][0];      iy=loop_0[ii][1];      iz=loop_0[ii][2];
#ifdef SELFTIMESTEP
	      dtcell=get_u_scalar(cell_dt,ix,iy,iz);
#endif
	      PLOOP(iv) set_u(dut2,iv,ix,iy,iz,(1./(dtcell))*get_u(u,iv,ix,iy,iz)+(-1./(dtcell))*get_u(ut2,iv,ix,iy,iz)); 
	    }
	    //addi_u(1./dt,u,-1./dt,ut2,dut2); //F(U(2)) in *dut2;
	    
            /******* explicit together **********/
	    for(ii=0;ii<Nloop_0;ii++)  { ix=loop_0[ii][0];      iy=loop_0[ii][1];      iz=loop_0[ii][2];
#ifdef SELFTIMESTEP
	      dtcell=get_u_scalar(cell_dt,ix,iy,iz);
#endif	      
	      PLOOP(iv) set_u(u,iv,ix,iy,iz,get_u(ut0,iv,ix,iy,iz)+(dtcell/2.)*get_u(dut1,iv,ix,iy,iz)+(dtcell/2.)*get_u(dut2,iv,ix,iy,iz)); 
	    }
	    //addi_u_3(1.,ut0,dt/2.,dut1,dt/2.,dut2,u); //U(n) + dt/2 (F(U(1)) + F(U(2))) in *u
	    

	    /******* implicit together ***********/
	    for(ii=0;ii<Nloop_0;ii++)  { ix=loop_0[ii][0];      iy=loop_0[ii][1];      iz=loop_0[ii][2];
#ifdef SELFTIMESTEP
	      dtcell=get_u_scalar(cell_dt,ix,iy,iz);
#endif
	      PLOOP(iv) set_u(u,iv,ix,iy,iz,get_u(u,iv,ix,iy,iz)+(dtcell/2.)*get_u(drt1,iv,ix,iy,iz)+(dtcell/2.)*get_u(drt2,iv,ix,iy,iz));
	    }
	    //addi_u_3(1.,u,dt/2.,drt1,dt/2.,drt2,u); //u += dt/2 (R(U(1)) + R(U(2))) in *u
	    //	    calc_u2p();
	


#ifdef EVOLVEINTENSITIES //only 1st order!
	    update_intensities(t, dt);
#endif
	    //printf("nstep: %d\n",nstep);

	    calc_u2p();
#pragma omp barrier
	    do_correct();
#pragma omp barrier
	
	    t+=dt;	 

	  }
	else if(TIMESTEPPING==RK2)
	  { 
	    /******************************* RK2 **********************************/

	    save_timesteps();

	    //1st
	    copyi_u(1.,u,ut0);
	    count_entropy(&nentr[0],&nentr2[0]); copy_entropycount(); do_correct();
	    op_explicit (t,0.5*dt); 

	    calc_u2p();
	    
#pragma omp barrier
	    count_entropy(&nentr[1],&nentr2[1]);
	    op_implicit (t,0.5*dt); 

	    copyi_u(1.,u,ut1);

#ifdef EVOLVEINTENSITIES
	    copyi_intensities(1.,Ibeam,Ibeam0);
	    update_intensities(t, 0.5*dt);
	    copyi_intensities(1.,Ibeam,Ibeam1);
#endif

	    //2nd
	    //calc_u2p();
#pragma omp barrier
	    count_entropy(&nentr[2],&nentr2[2]); 
	    op_explicit (t,dt); 

	    calc_u2p();
#pragma omp barrier
	    count_entropy(&nentr[3],&nentr2[3]); 
	    op_implicit (t,dt);

	    #if(AVGOUTPUT==1) 
	    //save to avg arrays	    
	    save_avg(dt);
	    #endif
 

	    //together
	    addi_u(1.,u,-1.,ut1,ut2); //k2 in ut2
	    addi_u(1.,ut0,1.,ut2,u);

#ifdef EVOLVEINTENSITIES
	    update_intensities(t, dt);
	    addi_intensities(1.,Ibeam,-1.,Ibeam1,Ibeambak);
	    addi_intensities(1.,Ibeam0,1.,Ibeambak,Ibeam);
#endif
	  
	    t+=dt;
	  
	    /************************** end of RK2 **********************************/
	  }
	else if(TIMESTEPPING==RK2HEUN)
	  { 
	    /******************************* RK2 Heun **********************************/
	    save_timesteps();

	    //1st	 
	    copyi_u(1.,u,ut0);
	    calc_u2p();
	    //#pragma omp barrier
	    count_entropy(&nentr[0],&nentr2[0]); copy_entropycount(); 
	    do_correct();
	    //#pragma omp barrier
	    op_explicit (t,1.*dt); 
#ifdef RADIATION
	    calc_u2p();
	    //#pragma omp barrier
	    do_correct();
#endif
	    //#pragma omp barrier
	    count_entropy(&nentr[1],&nentr2[1]); 
	    op_implicit (t,1.*dt); 
	    addi_u(1.,u,-1.,ut0,ut2); 

	    //2nd
	    copyi_u(1.,u,ut1);
	    calc_u2p();
	    //#pragma omp barrier
	    do_correct();
	    //#pragma omp barrier
	    count_entropy(&nentr[2],&nentr2[2]); 
	    //#pragma omp barrier
	    op_explicit (t,dt); 
#ifdef RADIATION
	    calc_u2p();
	    //#pragma omp barrier
	    do_correct();
#endif
	    //#pragma omp barrier
	    count_entropy(&nentr[3],&nentr2[3]); 
	    op_implicit (t,dt);

	    #if(AVGOUTPUT==1) 
	    //save to avg arrays	    
	    save_avg(dt);
	    #endif
 
	    addi_u(1.,u,-1.,ut1,ut3); 

	    //together     
	    addi_u_3(1.,ut0,1./2.,ut2,1./2.,ut3,u); //u += dt/2 (R(U(1)) + R(U(2))) in *u

	    //	  getch();
	    t+=dt;
	    /************************** end of RK2 **********************************/
	  }
	else 
	my_err("wrong time stepping specified\n");
      
 
      
      //**********************************************************************
      //************************* updating intensities *************************
      //**********************************************************************
 
      #ifdef EVOLVEINTENSITIES1STORDER
      update_intensities(t, dt);
#endif


      //**********************************************************************
      //************************* outputs ************************************
      //**********************************************************************
   
 
      
      


      #ifdef RADIATION
      //average number of iterations in the implicit solver
      ldouble avimpitloc[5],avimpit[5];
      int impnumsloc[7],impnums[7];
      
      avimpitloc[0]=global_int_slot[GLOBALINTSLOT_NIMPENERMHD]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENERMHD]
	/((ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERMHD]+(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERMHDFF]);
      avimpitloc[1]=global_int_slot[GLOBALINTSLOT_NIMPENERRAD]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENERRAD]
	/((ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERRAD]+(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERRADFF]);
      avimpitloc[2]=global_int_slot[GLOBALINTSLOT_NIMPENTRMHD]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENTRMHD]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENTRMHD];
      avimpitloc[3]=global_int_slot[GLOBALINTSLOT_NIMPENTRRAD]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENTRRAD]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENTRRAD];
      avimpitloc[4]=global_int_slot[GLOBALINTSLOT_NIMPLTE]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPLTE]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPLTE];
            
      impnumsloc[0]=global_int_slot[GLOBALINTSLOT_NIMPENERMHD];
      impnumsloc[1]=global_int_slot[GLOBALINTSLOT_NIMPENERRAD];
      impnumsloc[2]=global_int_slot[GLOBALINTSLOT_NIMPENERMHDFF];
      impnumsloc[3]=global_int_slot[GLOBALINTSLOT_NIMPENERRADFF];
      
      impnumsloc[4]=global_int_slot[GLOBALINTSLOT_NIMPENTRMHD];
      impnumsloc[5]=global_int_slot[GLOBALINTSLOT_NIMPENTRRAD];
      impnumsloc[6]=global_int_slot[GLOBALINTSLOT_NIMPLTE]+global_int_slot[GLOBALINTSLOT_NRADFIXUPS]+global_int_slot[GLOBALINTSLOT_NCRITFAILURES];

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
      
 
      //time mark
      my_clock_gettime(&temp_clock);    

      end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

      //performance
      ldouble znps=TNX*TNY*TNZ/(end_time-start_time);

      //avg files
#if(AVGOUTPUT==1) 
      //save to avg arrays	    
      //save_avg(dt);
      //dump avg file?
      if(lasttoutavg_floor!=floor(t/dtoutavg))
	{
	  if(PROCID==0)
	    printf("%d > avg file no #%6d dumped\n",PROCID,nfout2);
	  
	  //avg goes first so that what is later can use it
	  #ifndef SELFTIMESTEP
	  copy_u_core(1./avgtime,pavg,pavg,SX*SY*SZ*(NV+NAVGVARS));
	  #else
	  for(iz=0;iz<NZ;iz++)
	    for(iy=0;iy<NY;iy++)
	      for(ix=0;ix<NX;ix++)
		{
		  for(iv=0;iv<(NV+NAVGVARS);iv++)
		    set_uavg(pavg,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iz) / get_u_scalar(avgselftime,ix,iy,iz));
		  set_u_scalar(avgselftime,ix,iy,iz,0.);
		}              
	  #endif
	  fprint_avgfile(t,folder,"avg");
	  //zeros avg
	  copy_u_core(0.,pavg,pavg,SX*SY*SZ*(NV+NAVGVARS));
	  avgtime=0.;
  
	  nfout2++;

	  lasttoutavg_floor=floor(t/dtoutavg);	 
	}
#endif

     //boxscalars.dat file
#if(BOXOUTPUT==1) 
      if(lasttoutbox_floor!=floor(t/dtoutbox))
	{
	  fprint_boxscalars(t);
	  lasttoutbox_floor=floor(t/dtoutbox);	 
	}
#endif

 //varscalars.dat file
#if(VAROUTPUT==1) 
      if(lasttoutvar_floor!=floor(t/dtoutvar))
	{
	  fprint_varscalars(t);
	  lasttoutvar_floor=floor(t/dtoutvar);	 
	}
#endif

      //snapshots
      if(lasttout_floor!=floor(t/dtout) || ALLSTEPSOUTPUT || t>.9999999*t1 || spitoutput==1)
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
#if(THOUTPUT==1)
	  fprint_thprofiles(t,nfout1,folder,"th");
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
	  //znps = TNX*TNY*TNZ*(nstep-fprintf_nstep);
	  
	  printf("(%d) step #%6d t=%10.3e dt=%.3e (tot.time: %.2e (%.2e|%3d > %.2e|%3d) mp: %7.6f) znps: %.0f "
		 ,PROCID,nstep,t,dt,end_time-start_time,2.*max_u2ptime,max_u2ptime_loc,2.*min_u2ptime,min_u2ptime_loc,2*maxmp_time,znps);

#ifdef BALANCEENTROPYWITHRADIATION
	  printf("#: %d|%d %d|%d %d|%d %d|%d ",nentr[0],nentr2[0],nentr[1],nentr2[1],nentr[2],nentr2[2],nentr[3],nentr2[3]);
#else
	  printf("#: %d %d %d %d ",nentr[0],nentr[1],nentr[2],nentr[3]);
#endif

#ifdef RADIATION
#ifndef SKIPRADSOURCE
	  printf("#:%d %d %d %d %d %d %d | %.1f %.1f %.1f %.1f %.1f ",
		 impnums[0],impnums[1],impnums[2],impnums[3],impnums[4],impnums[5],impnums[6],
		 avimpit[0],avimpit[1],avimpit[2],avimpit[3],avimpit[4]);
#endif
#endif

#ifdef MSTEP
	  int nlev[NUMMSTEPLEVELS];
	  mstep_count_levels(nlev);
	  printf("| ");
	  int ilev,ilevshow=NUMMSTEPLEVELS-1;
	  for(ilev=0;ilev<ilevshow;ilev++) 
	    printf("%d ",nlev[ilev]);
	  for(ilev=ilevshow+1;ilev<NUMMSTEPLEVELS;ilev++) 
	    nlev[ilevshow]+=nlev[ilev];
	  printf(" +%d ",nlev[ilevshow]);
	  //getch();
#endif

	  printf("\n");

	  fflush(stdout);

	  fprintf_time=end_time;
	  fprintf_nstep = nstep;
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

  int aa[3],bb[2],ret;
  pp[UU]*=1.1;
  ret=u2p(uu,pp,&geom,aa,bb,0);
  printf("u2p ret: %d\n",ret);
  print_NVvector(pp);

  return 0;
}

int
print_scalings()
{
  printf("\n ***************************************\n\n");
  printf("BH mass: %.6f\nspin: %.6f\n\nscalings  (GU->CGS):\nrho: %.6e\nmdot: %.6e\nsigma: %.6e\nlen: %.6e\ntime: %.6e\nenden:"
	 "%.6e\nflux: %.6e\nT(1,1): %.6e\nkbt: %.6e\nkb/me: %.6e\nsigma_rad: %.6e\nkappa: %.6e\nGt: %.6e\n\n"
	 "rhorizonBL: %.6f\nrISCOBL: %.6f\netaNT: %.6f\n\n->mdotEdd: %.6e\n->lumEdd: %.6e\n\nGMc2: %.6e\nGMc3: %.6e\n",
	 MASS,BHSPIN,
	 rhoGU2CGS(1.),
	 rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.),
	 rhoGU2CGS(1.)*lenGU2CGS(1.),
	 lenGU2CGS(1.),
	 timeGU2CGS(1.),
	 endenGU2CGS(1.),
	 fluxGU2CGS(1.),
	 calc_PEQ_Tfromurho(1.,1.),
	 K_BOLTZ/MU_GAS/M_PROTON,
	 K_BOLTZ/M_ELECTR,
	 SIGMA_RAD,
	 kappaCGS2GU(1.),
	 kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.)*CCC,
	 rhorizonBL,
	 rISCOBL,
	 etaNT,
	 rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)/calc_mdotEdd(),
	 rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)*CCC0*CCC0/calc_lumEdd(),
	 GMC2,
	 GMC3
	 );
    printf("\n ***************************************\n\n");
  return 0;
}
