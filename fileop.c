//KORAL - fileop.c
//file operations

#include "ko.h"

/*********************************************/
/* opens files etc. */
/*********************************************/
int 
fprint_openfiles(char* folder)
{
  char bufor[100];

#ifndef RESTART
  sprintf(bufor,"rm %s/*",folder);
  int i=system(bufor);
  nfout1=0;
#endif

  sprintf(bufor,"%s/scalars.dat",folder);
  fout_scalars=fopen(bufor,"a");

  return 0;
}


/*********************************************/
/* closes file handles */
/*********************************************/
int 
fprint_closefiles()
{
  fclose(fout_scalars);
  return 0;
}

/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_gridfile(char* folder)
{
  FILE* out;
  char bufor[50];
  sprintf(bufor,"%s/grid.dat",folder);
  out=fopen(bufor,"w");

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

	      ldouble xxcar[4],xxsph[4];

	      coco_N(geom.xxvec,xxcar,MYCOORDS,MINKCOORDS); 
	      coco_N(geom.xxvec,xxsph,MYCOORDS,KERRCOORDS); 

	      fprintf(out,"%d %d %d %f %f %f %f %f %f %f %f %f\n",
		      ix,iy,iz, geom.xxvec[1],geom.xxvec[2],geom.xxvec[3], xxcar[1],xxcar[2],xxcar[3],xxsph[1],xxsph[2],xxsph[3]);

	    }
	}
    }

  fclose(out);

  return 0;
}
/*********************************************/
/*********************************************/
/*********************************************/
/* prints dumps to files .dat and calls gnuplot */
/* codeprim == 1 - prints out code primitives, only coordinates converted to OUTCOORDS */
/* codeprim == 0 - prints ZAMO frame etc primitives - post processing, called by ana.c */
/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_profiles(ldouble t, ldouble *scalars, int nscalars, int codeprim, char* folder)
{
  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/rad%04d.dat",folder,nfout1);
  fout_radprofiles=fopen(bufor,"w");

  sprintf(bufor,"%s/out%04d.dat",folder,nfout1);
  fout1=fopen(bufor,"w");
  
  //header
  //## nout time problem NX NY NZ
  fprintf(fout1,"## %d %e %d %d %d %d\n",nfout1,t,PROBLEM,NX,NY,NZ);

  sprintf(bufor2,"%s/out%04d.%s",folder,nfout1,IMAGETYPE);  
  int ix,iy,iz,iv;
  int gclx,gcrx,gcly,gcry,gclz,gcrz;

  //whether print ghost cells or not - default values
  gclx=gcly=gclz=0;
  gcrx=gcry=gcrz=0;
#ifdef PRINTGC_LEFT
  gclx=1;//gcly=1;
#endif
#ifdef PRINTGC_RIGHT
  gcrx=1;
#endif
#ifdef PRINTXGC_LEFT
  gclx=1;
#endif
#ifdef PRINTXGC_RIGHT
  gcrx=1;
#endif
#ifdef PRINTYGC_LEFT
  gcly=1;
#endif
#ifdef PRINTYGC_RIGHT
  gcry=1;
#endif
#ifdef PRINTZGC_LEFT
  gclz=1;
#endif
#ifdef PRINTZGC_RIGHT
  gcrz=1;
#endif
#ifdef PRINTZONEMORE
  gcrz=1;
#endif
  
  
  //printing scalars
  fprintf(fout_scalars,"%e ",t);
  for(iv=0;iv<nscalars;iv++)
    fprintf(fout_scalars,"%e ",scalars[iv]);
  fprintf(fout_scalars,"\n");
  fflush(fout_scalars);
  
  if(MYCOORDS == BLCOORDS || MYCOORDS == KSCOORDS || MYCOORDS == MKS1COORDS)
    {
      //calculating radial profiles
      ldouble profiles[NRADPROFILES][NX];
      calc_radialprofiles(profiles);
      //printing radial profiles  
      for(ix=0;ix<NX;ix++)
	{
	  ldouble xx[4],xxout[4];
	  get_xx(ix,0,0,xx);
	  coco_N(xx,xxout,MYCOORDS,BLCOORDS); 
	  if(xxout[1]<r_horizon_BL(BHSPIN) && codeprim==0) continue;
	  fprintf(fout_radprofiles,"%e ",xxout[1]);
	  for(iv=0;iv<NRADPROFILES;iv++)
	    fprintf(fout_radprofiles,"%e ",profiles[iv][ix]);
	  fprintf(fout_radprofiles,"\n");
	}
      fflush(fout_radprofiles);
    }
  fclose(fout_radprofiles);
  
  /**************************/  
  /** writing order *********/  
  /**************************/  
 
#ifdef YZXDUMP
  for(iy=0;iy<NY;iy++)
    {
      for(iz=-gclz*NG;iz<NZ+gcrz*NG;iz++)
	{
	  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
	    {
#elif defined(ZXYDUMP)
	      for(iz=0;iz<NZ;iz++)
		{
		  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
		    {
		      for(iy=-gcly*NG;iy<NY+gcry*NG;iy++)
			{
#elif defined(YSLICE)
			  for(iy=YSLICE;iy<YSLICE+1;iy++)
			    {
			      for(iz=0;iz<NZ+gcrz*NG;iz++)
				{
				  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
				    {
#elif defined(ZSLICE)
				      for(iz=ZSLICE;iz<ZSLICE+1;iz++)
					{
					  for(iy=0;iy<NY;iy++)
					    {
					      for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
						{
#elif defined(YZSLICE)
						  for(iy=NY/2;iy<NY/2+1;iy++)
						    {
						      for(iz=NZ/2;iz<NZ/2+1;iz++)
							{
							  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
							    {
#else
							      for(iz=0;iz<NZ;iz++)
								{
								  for(iy=-gcly*NG;iy<NY+gcry*NG;iy++)
								    {
								      for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
									{
#endif

									  /**************************/  
									  /**************************/  
									  /**************************/  


									  //within domain:
									  if(if_indomain(ix,iy,iz)==0 && if_outsidegc(ix,iy,iz)==1) continue;
						  
									  struct geometry geom;
									  fill_geometry(ix,iy,iz,&geom);

									  ldouble mx,my,mz,E,e,xx,yy,zz,phipot,xxx[4],dx[3],vv[10],a0,a1,a2,v1,v2,dphidx,v3,Tgas,Trad,v4,v5,v6,v7,Fx,Fy,Fz;
									  ldouble gg[4][5],GG[4][5];
									  ldouble pp[NV],uu[NV];
									  int i,j;

									  v1=v2=v3=v4=v5=v6=v7=0.;

									  ldouble xxvec[4],xxvecout[4];

									  //transforming code coordinates to output coordinates
									  get_xx(ix,iy,iz,xxvec);
						      
									  coco_N(xxvec,xxvecout,MYCOORDS,OUTCOORDS);
									  //						  coco_N(xxvec,xxvecout,MYCOORDS,BLCOORDS);

									  xx=xxvecout[1];
									  yy=xxvecout[2];
									  zz=xxvecout[3];

									  if(codeprim==0)
									    {
#ifndef PRINTINSIDEBH						  
									      if((OUTCOORDS==KERRCOORDS || OUTCOORDS==BLCOORDS) && xx<1.*r_horizon_BL(BHSPIN)) continue;
#endif
									    }

									  xxx[0]=t;
									  xxx[1]=xx;
									  xxx[2]=yy;
									  xxx[3]=zz;

									  pick_g(ix,iy,iz,gg);
									  pick_G(ix,iy,iz,GG);
									  ldouble gdet=gg[3][4];

									  dx[0]=get_size_x(ix,0)*sqrt(gg[1][1]);
									  dx[1]=get_size_x(iy,1)*sqrt(gg[2][2]);
									  dx[2]=get_size_x(iz,2)*sqrt(gg[3][3]);   

									  //						  calc_primitives(ix,iy,iz);

									  ldouble pporg[NV];
									  for(iv=0;iv<NV;iv++)
									    {
									      uu[iv]=get_u(u,iv,ix,iy,iz);
									      pp[iv]=get_u(p,iv,ix,iy,iz);
									      pporg[iv]=get_u(p,iv,ix,iy,iz);
									    }	 

						  
									  ldouble tup[4][4],tlo[4][4];
									  pick_T(tmuup,ix,iy,iz,tup);
									  pick_T(tmulo,ix,iy,iz,tlo);	    
									  ldouble eup[4][4],elo[4][4];
									  pick_T(emuup,ix,iy,iz,eup);
									  pick_T(emulo,ix,iy,iz,elo);
									  //to transform primitives between coordinates if necessary
									  ldouble ggout[4][5],GGout[4][5];
									  struct geometry geomout;
									  calc_g_arb(xxvecout,ggout,OUTCOORDS);
									  calc_G_arb(xxvecout,GGout,OUTCOORDS);
									  fill_geometry_arb(ix,iy,iz,&geomout,OUTCOORDS);

						  

									  if(MYCOORDS!=OUTCOORDS && codeprim==0)
									    {
#ifdef RADIATION
									      trans_prad_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);
#endif
									      trans_pmhd_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);

									      //from now on geom,gg, GG, tup, etc. defined in OUTCOORDS!
									      fill_geometry_arb(ix,iy,iz,&geom,OUTCOORDS);
									      for(i=0;i<4;i++)
										for(j=0;j<5;j++)
										  { gg[i][j]=ggout[i][j]; GG[i][j]=GGout[i][j]; }

									      calc_tetrades(gg,tup,tlo,OUTCOORDS);
									      calc_ZAMOes(gg,eup,elo,OUTCOORDS);
									    }

									  ldouble rho=pp[0];
									  ldouble uint=pp[1];
									  ldouble S=pp[5];
									  ldouble pre=(GAMMA-1.)*uint;
									  gdet=gg[3][4];
									  ldouble ut=uu[0]/gdet/rho;
									  Tgas=pre*MU_GAS*M_PROTON/K_BOLTZ/rho;

									  ldouble vx=pp[2];
									  ldouble vy=pp[3];
									  ldouble vz=pp[4];
									  ldouble vrel[4]={0,vx,vy,vz};	

									  if(codeprim==0)
									    {
									      conv_vels(vrel,vrel,VELPRIM,VEL4,gg,GG);						  
									      trans2_cc2on(vrel,vrel,tup);
									      //outvel - ortonormal VEL4
									      vx=vrel[1];
									      vy=vrel[2];
									      vz=vrel[3];
									    }
						  
						  

									  /******************/
									  /* extra lines to calculate v1...v4 from PROBLEMS/XXX/dump.c */

#include PR_DUMP
									  /******************/

#ifdef RADIATION

									  if(codeprim==0)
									    {
#ifndef LABRADFLUXES
#ifndef EDDINGTON_APR						

#ifdef RADOUTPUTINFF
									  prad_lab2ff(pp,pp,&geom);
#endif

#ifdef RADOUTPUTINZAMO //to print  radiation primitives in ZAMO
									  prad_lab2on(pp,pp,&geom);
									  //prad_lab2ff(pp,pp,&geom);
									  //prad_ff2zamo(pp,pp,geomout.gg,geomout.GG,geomout.eup);
#endif

#ifdef RADOUTPUTVELS
									  ldouble vrelrad[4]={0,pp[7],pp[8],pp[9]};
									  conv_vels(vrelrad,vrelrad,VELPRIMRAD,VEL4,gg,GG);						  
									  trans2_cc2on(vrelrad,vrelrad,tup);
									  //rad outvel - ortonormal VEL4
									  pp[7]=vrelrad[1];
									  pp[8]=vrelrad[2];
									  pp[9]=vrelrad[3];	  
#endif

#endif
#endif
									    }
#endif
						  
									  //**********************************************************************
									  //**********************************************************************
									  //**********************************************************************

									  //summing up multifluids
									  int irf;
									  E=Fx=Fy=Fz=0.;
									  for(irf=0;irf<NRF;irf++)
									    {

									      //						      irf=0;
									      E+=pp[EE(irf)];
									      Fx+=pp[FX(irf)];
									      Fy+=pp[FY(irf)];
									      Fz+=pp[FZ(irf)];
									      //						      break;
									    }
#ifdef LABRADFLUXES
									  E*=-1.;
#endif

									  Trad=calc_LTE_TfromE(fabs(E));
									  if(E<EEFLOOR || isnan(E)) E=EEFLOOR;

									  fprintf(fout1,"%.4e %.4e %.4e "
										  "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e "
										  "%.3e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e "
										  "%.6e %.6e ",
										  xx,     //1
										  yy,     //2
										  zz,     //3		      
										  uu[0],  //4
										  uu[1],  //5
										  uu[2],  //6
										  uu[3],  //7
										  uu[4],  //8
										  uu[5],  //9
										  uu[6],  //10
										  uu[7],  //11
										  uu[8],  //12
										  uu[9],  //13

										  rho,    //14
										  uint, 
										  vx,     //16
										  vy,     //17
										  vz,     //18
										  S,      //19
#ifdef MAGNFIELD
										  pp[B1],
										  pp[B2],
										  pp[B3],
										  0.
#else
										  E,
										  Fx,
										  Fy,
										  Fz
#endif
										  );


									  fprintf(fout1,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
										  v1,     //24
										  v2,     //25
										  v3,     //26 
										  v4,     //27
										  v5,     //28
										  v6,     //29 
										  v7      //30 
										  );

									}
								      fprintf(fout1,"\n");
								    }
								  fprintf(fout1,"\n\n");
								}
							      fflush(fout1);
							      fclose(fout1);

							      //calling gnuplot to produce gifs
#ifdef YZSLICE
							      convert_out2gif_1d(bufor,bufor2,nfout1,t);
#else
							      if(NY>3 || NZ>3)
								convert_out2gif_2d(bufor,bufor2,nfout1,t);
							      else
								convert_out2gif_1d(bufor,bufor2,nfout1,t);
#endif
  
							      
							      return 0;

							    }
/*********************************************/
/*********************************************/
/*********************************************/
/* prints restart files */
/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_restartfile(ldouble t, char* folder)
{
  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);
  fout1=fopen(bufor,"w");
  
  //header
  //## nout time problem NX NY NZ
  fprintf(fout1,"## %d %e %d %d %d %d\n",nfout1,t,PROBLEM,NX,NY,NZ);

  /***********************************/  
  /** writing order is fixed  ********/  
  /***********************************/  
 
  int ix,iy,iz,iv;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      fprintf(fout1,"%d %d %d ",ix,iy,iz);
	      for(iv=0;iv<NV;iv++)
		{
		  pp[iv]=get_u(p,iv,ix,iy,iz);
		}	 
	     
	      for(iv=0;iv<NV;iv++)
		fprintf(fout1,"%.10e ",pp[iv]);
	      fprintf(fout1,"\n");
	    }
	}
    }

  fflush(fout1);
  fclose(fout1);

  return 0;
}
							  
/*********************************************/
/*********************************************/
/*********************************************/
/* reads dump file */
/* puts conserved into the memory */
/* converts them to primitives */
/*********************************************/
/*********************************************/
/*********************************************/
int 
fread_restartfile(int nout1, ldouble *t)
{
  //opening dump file
  int i,ret;
  char fname[40];
  sprintf(fname,"dumps/res%04d.dat",nout1);
  nfout1=nout1+1; //global file no.
  FILE *fdump=fopen(fname,"r");

  //reading parameters, mostly time
  int intpar[5];
  ret=fscanf(fdump,"## %d %lf %d %d %d %d\n",&intpar[0],t,&intpar[1],&intpar[2],&intpar[3],&intpar[4]);
  printf("dump file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	 fname,intpar[0],*t,intpar[1],intpar[2],intpar[3],intpar[4]); 


  int ix,iy,iz,iv;
  //reading conserved
  int gclx,gcrx,gcly,gcry,gclz,gcrz;
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      ret=fscanf(fdump,"%*f %*f %*f ");
									  
	      struct geometry geom;
	      fill_geometry(ix,iy,iz,&geom);

	      ldouble xxvec[4],xxvecout[4];

	      /**************************/  
	      /**************************/  
	      /**************************/  
	      ldouble uu[NV],pp[NV],ftemp;
	      char c;
	      
	      //reading primitives from file
	      for(i=0;i<NV;i++)
		{
		  ret=fscanf(fdump,"%lf ",&pp[i]);
		}
	      
	      p2u(pp,uu,&geom);

	      //saving primitives
	      for(iv=0;iv<NV;iv++)    
		{
		  set_u(u,iv,ix,iy,iz,uu[iv]);
		  set_u(p,iv,ix,iy,iz,pp[iv]);
		}
	    								  
	    }
	}
    }

    return 0;
}

/*********************************************/
/*********************************************/
/*********************************************/
/* reads dump file */
/* puts conserved into the memory */
/* converts them to primitives */
/*********************************************/
/*********************************************/
/*********************************************/
int 
fread_dumpfile(int nout1, ldouble *t)
{
  //opening dump file
  int i,ret;
  char fname[40];
  sprintf(fname,"dumps/out%04d.dat",nout1);
  nfout1=nout1+1; //global file no.
  FILE *fdump=fopen(fname,"r");

  //reading parameters, mostly time
  int intpar[5];
  ret=fscanf(fdump,"## %d %lf %d %d %d %d\n",&intpar[0],t,&intpar[1],&intpar[2],&intpar[3],&intpar[4]);
  printf("dump file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	 fname,intpar[0],*t,intpar[1],intpar[2],intpar[3],intpar[4]); 


  int ix,iy,iz,iv;
  //reading conserved
  int gclx,gcrx,gcly,gcry,gclz,gcrz;

  /**************************************************/  
  /** reading size  *********************************/  
  /** must be the same as in fprint_profile *********/  
  /**************************************************/  
  gclx=gcly=gclz=0;
  gcrx=gcry=gcrz=0;
#ifdef PRINTGC_LEFT
  gclx=1;//gcly=1;
#endif
#ifdef PRINTGC_RIGHT
  gcrx=1;
#endif
#ifdef PRINTXGC_LEFT
  gclx=1;
#endif
#ifdef PRINTXGC_RIGHT
  gcrx=1;
#endif
#ifdef PRINTYGC_LEFT
  gcly=1;
#endif
#ifdef PRINTYGC_RIGHT
  gcry=1;
#endif
#ifdef PRINTZGC_LEFT
  gclz=1;
#endif
#ifdef PRINTZGC_RIGHT
  gcrz=1;
#endif
#ifdef PRINTZONEMORE
  gcrz=1;
#endif

  /**************************************************/  
  /** reading order *********************************/  
  /** must be the same as in fprint_profile *********/  
  /**************************************************/  


#ifdef YZXDUMP
  for(iy=0;iy<NY;iy++)
    {
      for(iz=-gclz*NG;iz<NZ+gcrz*NG;iz++)
	{
	  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
	    {
#elif defined(ZXYDUMP)
	      for(iz=0;iz<NZ;iz++)
		{
		  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
		    {
		      for(iy=-gcly*NG;iy<NY+gcry*NG;iy++)
			{
#elif defined(YSLICE)
			  for(iy=YSLICE;iy<YSLICE+1;iy++)
			    {
			      for(iz=0;iz<NZ;iz++)
				{
				  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
				    {
#elif defined(ZSLICE)
				      for(iz=ZSLICE;iz<ZSLICE+1;iz++)
					{
					  for(iy=0;iy<NY;iy++)
					    {
					      for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
						{
#elif defined(YZSLICE)
						  for(iy=NY/2;iy<NY/2+1;iy++)
						    {
						      for(iz=NZ/2;iz<NZ/2+1;iz++)
							{
							  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
							    {
#else
							      for(iz=0;iz<NZ;iz++)
								{
								  for(iy=-gcly*NG;iy<NY+gcry*NG;iy++)
								    {
								      for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
									{
#endif 
									  struct geometry geom;
									  fill_geometry(ix,iy,iz,&geom);

									  //*****************
									  //verify consistency here with fprint_profiles!
									  //*****************
									  //within domain:
									  if(if_indomain(ix,iy,iz)==0) continue;
									  ldouble xxvec[4],xxvecout[4];

									  //transforming code coordinates to output coordinates
									  get_xx(ix,iy,iz,xxvec);
									  coco_N(xxvec,xxvecout,MYCOORDS,OUTCOORDS);
									  ldouble xx = xxvecout[1];									  

									  /**************************/  
									  /**************************/  
									  /**************************/  
									  ldouble uu[NV],pp[NV],ftemp;
									  char c;
	      
									  //reading conserved and primitives from file
									  ret=fscanf(fdump,"%*f %*f %*f ");
									  
									  //10 hardcoded here
									  for(i=0;i<10;i++)
									    {
									      ret=fscanf(fdump,"%lf ",&ftemp);
									      if(i<NV) uu[i]=ftemp;
									    }
									  for(i=0;i<10;i++)
									    {
									      ret=fscanf(fdump,"%lf ",&ftemp);
									      if(i<NV) pp[i]=ftemp;
									    }

									  //rest of line
									  do
									    {
									      if(feof(fdump)) return 0;
									      c = fgetc(fdump);
									    }
									  while(c != '\n');

									  //saving primitives
									  for(iv=0;iv<NV;iv++)    
									    {
									      set_u(u,iv,ix,iy,iz,uu[iv]);
									      set_u(p,iv,ix,iy,iz,pp[iv]);
									    }

									  //sanity check using pp as initial guess
									  calc_primitives(ix,iy,iz);									  
									}
								    }
								}

    return 0;
}

/*********************************************/
/*********************************************/
/*********************************************/
/* prints in ASCII indices, cart coordinates,*/
/* primitives, velocities in cartesian       */
/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_simplecart(ldouble t, char* folder)
{
  if(PROBLEM==55) //G2STATIC
    {
      char bufor[50];
      sprintf(bufor,"%s/sim%04d.dat",folder,nfout1);
      fout1=fopen(bufor,"w");
  
      //header
      //## nout time problem NX NY NZ
      fprintf(fout1,"## %d %e %d %d %d %d\n",nfout1,t,PROBLEM,NX,NY,NZ);

      /***********************************/  
      /** writing order is fixed  ********/  
      /***********************************/  
 
      int ix,iy,iz,iv;
      ldouble pp[NV];
      for(iz=0;iz<NZ;iz++)
	{
	  for(iy=0;iy<NY;iy++)
	    {
	      for(ix=0;ix<NX;ix++)
		{
		  for(iv=0;iv<NV;iv++)
		    {
		      pp[iv]=get_u(p,iv,ix,iy,iz);
		    }
		  struct geometry geom,geomcart,geomout,geomsph;
		  fill_geometry(ix,iy,iz,&geom);
		  fill_geometry_arb(ix,iy,iz,&geomcart,MINKCOORDS);
		  fill_geometry_arb(ix,iy,iz,&geomout,OUTCOORDS);
		  fill_geometry_arb(ix,iy,iz,&geomsph,SPHCOORDS);

		  ldouble dx[3];
		  dx[0]=get_size_x(ix,0);
		  dx[1]=get_size_x(iy,1);
		  dx[2]=get_size_x(iz,2);
		  ldouble gdet=geom.gdet;
		  ldouble volume=dx[0]*dx[1]*dx[2]*gdet;
		  trans_pmhd_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomout);
		  ldouble rho=rhoGU2CGS(pp[RHO]);
		  ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
		  ldouble tracer=pp[TRA];
		  ldouble vel[4]={0,pp[VX],pp[VY],pp[VZ]};	
		  ldouble vx,vy,vz;
		  conv_vels(vel,vel,VELPRIM,VEL4,geomout.gg,geomout.GG);						  
		  trans2_cc2on(vel,vel,geomout.tup);
		  //transform to cartesian
		  if (MYCOORDS==SCHWCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS)
		    {
		      ldouble r=geomsph.xx;
		      ldouble th=geomsph.yy;
		      ldouble ph=geomsph.zz;

		      vx = sin(th)*cos(ph)*vel[1] 
			+ cos(th)*cos(ph)*vel[2]
			- sin(ph)*vel[3];

		      vy = sin(th)*sin(ph)*vel[1] 
			+ cos(th)*sin(ph)*vel[2]
			+ cos(ph)*vel[3];

		      vz = cos(th)*vel[1] 
			- sin(th)*vel[2];
		    }
	     
		  fprintf(fout1,"%d %d %d ",ix,iy,iz);

		  fprintf(fout1,"%.5e %.5e %.5e ",geomcart.xx,geomcart.yy,geomcart.zz);

		  fprintf(fout1,"%.5e %.5e %.5e ",rho,temp,tracer);

		  fprintf(fout1,"%.5e %.5e %.5e ",vx,vy,vz);

		  fprintf(fout1,"%.5e \n",volume);
		}
	    }
	}

      fflush(fout1);
      fclose(fout1);
    }

  return 0;
}
