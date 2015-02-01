//KORAL - fileop.c
//file operations

#include "ko.h"

/*********************************************/
/*  adds up current quantities to the pavg array */
/*********************************************/
int
save_avg(ldouble dt)
{
  int ix,iy,iz,iv,ii;

  //commented out to work with SUBZONES too
  //  for(ii=0;ii<Nloop_0;ii++) //domain 
  //ix=loop_0[ii][0];
  //iy=loop_0[ii][1];
  //iz=loop_0[ii][2]; 

#pragma omp parallel private(ix,iy,iz,iv) 
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++) 
	    {
     
      
	      ldouble avg[NV+NAVGVARS];
	      p2avg(ix,iy,iz,avg);

	      for(iv=0;iv<NV+NAVGVARS;iv++)
		{
		  set_uavg(pavg,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iz)+avg[iv]*dt);
		}
	    }
	}
    }

  avgtime+=dt;
  
  return 0;
}


/*********************************************/
/* opens files etc. */
/*********************************************/
int 
fprint_openfiles(char* folder)
{
  char bufor[100];

#ifndef RESTART
  if(PROCID==0)
    {
      sprintf(bufor,"rm %s/*",folder);
      int i=system(bufor);
    }
  nfout1=0;
  nfout2=0;
#endif

#ifndef MPI
  sprintf(bufor,"%s/scalars.dat",folder);
  fout_scalars=fopen(bufor,"a");

  sprintf(bufor,"%s/failures.dat",folder);
  fout_fail=fopen(bufor,"a");
#endif

  return 0;
}


/*********************************************/
/* closes file handles */
/*********************************************/
int 
fprint_closefiles()
{
#ifndef MPI
  fclose(fout_scalars);
  fclose(fout_fail);  
#endif
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
/* prints scalar quantities to scalars.dat */
/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_scalars(ldouble t, ldouble *scalars, int nscalars)
{
  #ifndef MPI
  calc_scalars(scalars,t);

  int iv;
  //printing scalars
  fprintf(fout_scalars,"%e ",t);
  for(iv=0;iv<nscalars;iv++)
    fprintf(fout_scalars,"%e ",scalars[iv]);
  fprintf(fout_scalars,"\n");
  fflush(fout_scalars);
  #endif

  return 0;
}

/*********************************************/
/*********************************************/
/*********************************************/
/* prints radial profiles to radNNNN.dat */
/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_radprofiles(ldouble t, int nfile, char* folder, char* prefix)
{
  //#ifdef BHDISK_PROBLEMTYPE 
      char bufor[50],bufor2[50];
      sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

      fout_radprofiles=fopen(bufor,"w");

      ldouble mdotscale = (rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.))/calc_mdotEdd();
      ldouble lumscale = (fluxGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.))/calc_lumEdd();

      fprintf(fout_radprofiles,"# mdotGU2Edd: %e lumGU2Edd: %e\n",mdotscale,lumscale);

      int ix,iv;
      //calculating radial profiles
      ldouble profiles[NRADPROFILES][NX];
      calc_radialprofiles(profiles);
      //printing radial profiles  
      for(ix=0;ix<NX;ix++)
	{
	  ldouble xx[4],xxout[4];
	  get_xx(ix,0,0,xx);
	  coco_N(xx,xxout,MYCOORDS,BLCOORDS); 
	  if(xxout[1]<rhorizonBL) continue;
	  fprintf(fout_radprofiles,"%e ",xxout[1]);
	  for(iv=0;iv<NRADPROFILES;iv++)
	    fprintf(fout_radprofiles,"%e ",profiles[iv][ix]);
	  fprintf(fout_radprofiles,"\n");
	}
      fclose(fout_radprofiles);
      //#endif
  
  return 0;
}
 
/*********************************************/
/*********************************************/
/*********************************************/
/* prints theta profiles to thNNNN.dat */
/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_thprofiles(ldouble t, int nfile, char* folder, char* prefix)
{
  //#ifdef BHDISK_PROBLEMTYPE 
  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

  FILE *fout_thprofiles=fopen(bufor,"w");

  int ix,iy,iv;

  //search for appropriate radial index
  ldouble xx[4],xxBL[4];
  ldouble radius=1.e3;
  #ifdef THPROFRADIUS
  radius=THPROFRADIUS;
  #endif
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  //calculating theta profiles
  ldouble profiles[NTHPROFILES][NY];
  calc_thetaprofiles(profiles);
  //printing th profiles  
  for(iy=0;iy<NY;iy++)
    {
      get_xx(ix,iy,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS); 
      
      fprintf(fout_thprofiles,"%e ",xxBL[2]);
      for(iv=0;iv<NTHPROFILES;iv++)
	fprintf(fout_thprofiles,"%e ",profiles[iv][iy]);
      fprintf(fout_thprofiles,"\n");
    }
  fclose(fout_thprofiles);
  
  return 0;
}
 

/*********************************************/
/*********************************************/
/*********************************************/
/* prints radial profiles to anarelradNNNN.dat   */
/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_anarelradprofiles(ldouble t, int nfile, char* folder, char* prefix, ldouble profiles[][NX])
{
#ifdef BHDISK_PROBLEMTYPE 
      char bufor[50],bufor2[50];
      sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

      fout_radprofiles=fopen(bufor,"w");

      int ix,iv;
      //printing radial profiles  
      for(ix=0;ix<NX;ix++)
	{
	  ldouble xx[4],xxout[4];
	  get_xx(ix,0,0,xx);
	  coco_N(xx,xxout,MYCOORDS,BLCOORDS); 
	  if(xxout[1]<rhorizonBL) continue;
	  fprintf(fout_radprofiles,"%e ",xxout[1]);
	  for(iv=0;iv<NANARELRADPROFILES;iv++)
	    fprintf(fout_radprofiles,"%e ",profiles[iv][ix]);
	  fprintf(fout_radprofiles,"\n");
	}
      fclose(fout_radprofiles);
#endif
  
  return 0;
}
 

/*********************************************/
/*********************************************/
/*********************************************/
/* prints dumps to files outNNNN.dat and calls gnuplot */
/* codeprim == 1 - prints out code primitives, only coordinates converted to OUTCOORDS */
/* codeprim == 0 - prints ZAMO frame etc primitives - post processing, called by ana.c */
/*********************************************/
/*********************************************/
/*********************************************/
int
fprint_outfile(ldouble t, int nfile, int codeprim, char* folder, char *prefix)
{
  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);
  fout1=fopen(bufor,"w");
  
 
  //header
  //## nout time problem NX NY NZ
  fprintf(fout1,"## %d %e %d %d %d %d\n",nfout1,t,PROBLEM,NX,NY,NZ);

  sprintf(bufor2,"%s/%s%04d.%s",folder,prefix,nfile,IMAGETYPE);  

  int ix,iy,iz,iv;
  int gclx,gcrx,gcly,gcry,gclz,gcrz;

  //whether print ghost cells or not - default values
  gclx=gcly=gclz=0;
  gcrx=gcry=gcrz=0;
#ifdef PRINTGC_LEFT
  gclx=1;
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
									  //if(if_indomain(ix,iy,iz)==0 && if_outsidegc(ix,iy,iz)==1) continue;
						  

									  


									  struct geometry geom;
									  fill_geometry(ix,iy,iz,&geom);

									  ldouble mx,my,mz,E,e,xx,yy,zz,phipot,xxx[4],dx[3],vv[10],a0,a1,a2,v1,v2,dphidx,v3,Tgas,Trad,v4,v5,v6,v7,v8,v9,v10,v11,v12,Fx,Fy,Fz;
									  ldouble gg[4][5],GG[4][5];
									  ldouble pp[NV],uu[NV];
									  int i,j;

									  v1=v2=v3=v4=v5=v6=v7=v8=v9=v10=v11=v12=0.;

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
									      if((OUTCOORDS==KERRCOORDS || OUTCOORDS==BLCOORDS) && xx<1.*rhorizonBL) continue;
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


									  ldouble pporg[NV];
									  for(iv=0;iv<NV;iv++)
									    {
									      uu[iv]=get_u(u,iv,ix,iy,iz);
									      pp[iv]=get_u(p,iv,ix,iy,iz);
									      pporg[iv]=get_u(p,iv,ix,iy,iz);
									    }	 



						  
									  ldouble tup[4][4],tlo[4][4];    
									  ldouble eup[4][4],elo[4][4];

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
									      //TODO
									      //tetrads sie zesraly
									      //trans2_cc2on(vrel,vrel,tup);
									      #ifdef BHDISK_PROBLEMTYPE
									      vrel[2]*=xx;
									      vrel[3]*=xx*sin(yy);
									      #endif

									      //outvel - ortonormal VEL4
									      vx=vrel[1];
									      vy=vrel[2];
									      vz=vrel[3];
									    }
						  
						  


#ifdef RADIATION

									  if(codeprim==0)
									    {
#ifdef RADOUTPUTINFF
									      prad_lab2ff(pp,pp,&geom);
#endif

#ifdef RADOUTPUTINZAMO //to print  radiation primitives in ZAMO, sometimes may not tetrad properly
									      prad_lab2on(pp,pp,&geom);
#endif

#ifdef RADOUTPUTVELS
									      ldouble vrelrad[4]={0,pp[7],pp[8],pp[9]};
									      conv_vels(vrelrad,vrelrad,VELPRIMRAD,VEL4,gg,GG);						  
									      //trans2_cc2on(vrelrad,vrelrad,tup);
									      //rad outvel - ortonormal VEL4
									      pp[7]=vrelrad[1];
									      pp[8]=vrelrad[2];
									      pp[9]=vrelrad[3];	  
#endif

									    }
#endif
						  


									  //**********************************************************************
									  //**********************************************************************
									  //**********************************************************************

									  //summing up multifluids
									  int irf;
									  E=Fx=Fy=Fz=0.;
									  //for(irf=0;irf<NRF;irf++)
									    {

									      //						      irf=0;
									      E+=pp[EE];
									      Fx+=pp[FX];
									      Fy+=pp[FY];
									      Fz+=pp[FZ];
									      //						      break;
									    }

									  #ifndef NCOMPTONIZATION
									  Trad=calc_LTE_TfromE(E);
									  #else
									  Trad=calc_ncompt_Thatrad_full(pp,&geomout);
									  #endif
									  if(E<EEFLOOR || isnan(E)) E=EEFLOOR;
									  

									  /******************/
									  /* extra lines to calculate v1...v4 from PROBLEMS/XXX/dump.c */

#include PR_DUMP
									  /******************/


									  fprintf(fout1,"%.4e %.4e %.4e "
										  "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e "
										  "%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e "
										  "%.9e %.9e ",
										  xx,     //1
										  yy,     //2
										  zz,     //3		      
										  uu[0],  //4
										  uu[1],  //5
										  uu[2],  //6
										  uu[3],  //7
										  uu[4],  //8
										  uu[5],  //9
										  uu[EE],  //10
										  uu[FX],  //11
										  uu[FY],  //12
										  uu[FZ],  //13

										  rho,    //14
										  uint, 
										  vx,     //16
										  vy,     //17
										  vz,     //18
										  S,      //19
#ifdef RADIATION
										  E,
										  Fx,
										  Fy,
										  Fz
#elif defined(MAGNFIELD)
										  pp[B1],
										  pp[B2],
										  pp[B3],
										  0.
#else
										  0.,
										  0.,
										  0.,
										  0.
#endif
										  );
									  
								

									  fprintf(fout1,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
										  v1,     //24
										  v2,     //25
										  v3,     //26 
										  v4,     //27
										  v5,     //28
										  v6,     //29 
										  v7,     //30 
										  v8,     //31 
										  v9,     //32 
										  v10,     //33
										  v11,     //34 
										  v12     //35
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
  #ifdef MPI

  fprint_restartfile_mpi(t,folder);

  #else 

  fprint_restartfile_bin(t,folder); 

  #endif
  
  return 0;
}

/*********************************************/
/*********************************************/

int //parallel output to a single file
fprint_restartfile_mpi(ldouble t, char* folder)
{
  #ifdef MPI
  char bufor[250];

  //header
  if(PROCID==0)
    {
       sprintf(bufor,"%s/res%04d.head",folder,nfout1);
       fout1=fopen(bufor,"w"); 
       //## nout time problem NX NY NZ
       sprintf(bufor,"## %5d %5d %10.6e %5d %5d %5d %5d\n",nfout1,nfout2,t,PROBLEM,TNX,TNY,TNZ);
       fprintf(fout1,"%s",bufor);
       fclose(fout1);
    }
  
  //maybe not needed
  MPI_Barrier(MPI_COMM_WORLD);

  //body
  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

  int rc = MPI_File_open( MPI_COMM_WORLD, bufor, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", bufor );fflush(stdout); exit(-1);
    }

  /***** first write all the indices ******/

  int indices[NX*NY*NZ*3];

  int ix,iy,iz,iv;
  int gix,giy,giz;

  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  indices[ix*NY*NZ*3+iy*NZ*3+iz*3+0]=gix;
	  indices[ix*NY*NZ*3+iy*NZ*3+iz*3+1]=giy;
	  indices[ix*NY*NZ*3+iy*NZ*3+iz*3+2]=giz;
	}

  //set the initial location at each process for indices
  MPI_Offset pos;
  pos=PROCID*NX*NY*NZ*(3*sizeof(int));  
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 

  //write all indices
  MPI_File_write( cFile, indices, NX*NY*NZ*3, MPI_INT, &status );
  
  /***** then primitives in the same order ******/

  //now let's try manually
  pos=TNX*TNY*TNZ*(3*sizeof(int)) + PROCID*NX*NY*NZ*(NV*sizeof(ldouble)); 
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 


  ldouble pout[NX*NY*NZ*NV];
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	for(iv=0;iv<NV;iv++)
	  pout[ix*NY*NZ*NV+iy*NZ*NV+iz*NV+iv]=get_u(p,iv,ix,iy,iz);

  MPI_File_write( cFile, pout, NX*NY*NZ*NV, MPI_LDOUBLE, &status );
  
  MPI_File_close( &cFile );


  //maybe not needed
  MPI_Barrier(MPI_COMM_WORLD);

  if(PROCID==0)
    {
      //sprintf(bufor,"cp %s/res%04d.dat %s/reslast.dat",folder,nfout1,folder);
      sprintf(bufor,"rm %s/reslast.dat",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.dat %s/reslast.dat",nfout1,folder);
      iv=system(bufor);
      //sprintf(bufor,"cp %s/res%04d.head %s/reslast.head",folder,nfout1,folder);
      sprintf(bufor,"rm %s/reslast.head",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.head %s/reslast.head",nfout1,folder);
      iv=system(bufor);
    }

#endif
  return 0;
}

/*********************************************/
/*********************************************/

int //serial binary output
fprint_restartfile_bin(ldouble t, char* folder)
{
  char bufor[250];
  
  //header
  if(PROCID==0)
    {
       sprintf(bufor,"%s/res%04d.head",folder,nfout1);
       fout1=fopen(bufor,"w"); 
       //## nout time problem NX NY NZ
       sprintf(bufor,"## %5d %5d %10.6e %5d %5d %5d %5d\n",nfout1,nfout2,t,PROBLEM,TNX,TNY,TNZ);
       fprintf(fout1,"%s",bufor);
       fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);
  fout1=fopen(bufor,"wb"); 

  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  //indices first
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  fwrite(&gix,sizeof(int),1,fout1);
	  fwrite(&giy,sizeof(int),1,fout1);
	  fwrite(&giz,sizeof(int),1,fout1);
	}

  //then, in the same order, primitives
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  fwrite(&get_u(p,0,ix,iy,iz),sizeof(ldouble),NV,fout1);
	}

  fclose(fout1);

  if(PROCID==0)
    {
      //sprintf(bufor,"cp %s/res%04d.dat %s/reslast.dat",folder,nfout1,folder)
      sprintf(bufor,"rm %s/reslast.dat",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.dat %s/reslast.dat",nfout1,folder);
      iv=system(bufor);
      //sprintf(bufor,"cp %s/res%04d.head %s/reslast.head",folder,nfout1,folder);
      sprintf(bufor,"rm %s/reslast.head",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.head %s/reslast.head",nfout1,folder);
      iv=system(bufor);    
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
fread_restartfile(int nout1, char* folder,ldouble *t)
{
  int ret;
  char bufor[250];
  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);

  #ifdef MPI

  ret=fread_restartfile_mpi(nout1,folder,t);

  #else //no MPI 

  ret=fread_restartfile_bin(nout1,folder,t);
  
  #endif
  
  return ret;
}

/*********************************************/
/*********************************************/

int 
fread_restartfile_bin(int nout1, char *folder, ldouble *t)
{
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz;
  char fname[40],fnamehead[40];
  if(nout1>=0)
    {
      sprintf(fname,"%s/res%04d.dat",folder,nout1);
      #ifdef MPI
      sprintf(fnamehead,"%s/../0/res%04d.head",folder,nout1);
      #else
      sprintf(fnamehead,"%s/res%04d.head",folder,nout1);
      #endif
    }
  else
    {
      sprintf(fname,"%s/reslast.dat",folder);
      #ifdef MPI
      sprintf(fnamehead,"%s/../0/reslast.head",folder);
      #else
      sprintf(fnamehead,"%s/reslast.head",folder);
      #endif
    }

  FILE *fdump;

  /***********/
  //header file
  fdump=fopen(fnamehead,"r");
  if(fdump==NULL) return 1; //request start from scratch

  //reading parameters, mostly time
  int intpar[6];
  ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],t,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
  printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	 fname,intpar[0],*t,intpar[2],intpar[3],intpar[4],intpar[5]); 
  nfout1=intpar[0]+1; //global file no.
  nfout2=intpar[1]; //global file no. for avg
  fclose(fdump);

  /***********/
  //body file
  fdump=fopen(fname,"rb");
 
  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;
  int indices[NX*NY*NZ][3];

  //first indices
  for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ret=fread(&gix,sizeof(int),1,fdump);
      ret=fread(&giy,sizeof(int),1,fdump);
      ret=fread(&giz,sizeof(int),1,fdump);

      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);

      indices[ic][0]=ix;
      indices[ic][1]=iy;
      indices[ic][2]=iz;
    }

  //then primitives
   for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ret=fread(pp,sizeof(ldouble),NV,fdump);

      ix=indices[ic][0];
      iy=indices[ic][1];
      iz=indices[ic][2];

      fill_geometry(ix,iy,iz,&geom);
      p2u(pp,uu,&geom);

      //saving primitives
      for(iv=0;iv<NV;iv++)    
	{
	  set_u(u,iv,ix,iy,iz,uu[iv]);
	  set_u(p,iv,ix,iy,iz,pp[iv]);
	}
    }

  fclose(fdump);

  return 0;
}

/*********************************************/
/*********************************************/

int 
fread_restartfile_mpi(int nout1, char *folder, ldouble *t)
{
  #ifdef MPI
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz,tix,tiy,tiz;
  char fname[40],fnamehead[40];
  if(nout1>=0)
    {
      sprintf(fname,"%s/res%04d.dat",folder,nout1);
      sprintf(fnamehead,"%s/res%04d.head",folder,nout1);
    }
  else
    {
      sprintf(fname,"%s/reslast.dat",folder);
      sprintf(fnamehead,"%s/reslast.head",folder);
    }

  FILE *fdump;

  /***********/
  //header file
  fdump=fopen(fnamehead,"r");
  if(fdump==NULL) return 1; //request start from scratch

  //reading parameters, mostly time
  int intpar[6];
  ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],t,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
  printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	 fname,intpar[0],*t,intpar[2],intpar[3],intpar[4],intpar[5]); 
  nfout1=intpar[0]+1; //global file no.
  nfout2=intpar[1]; //global file no. for avg
  fclose(fdump);


  //maybe not needed
  MPI_Barrier(MPI_COMM_WORLD);

  /***********/
  //body file
  struct geometry geom;
  ldouble uu[NV],pp[NV],ftemp;

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

  int rc = MPI_File_open( MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", fname );fflush(stdout); exit(-1);
    }

  /***** first read all the indices ******/

  
  //first read the indices
  #ifdef RESTARTGENERALINDICES
  int indices[TNX*TNY*TNZ*3];
  int len=TNX*TNY*TNZ;
  #else
  int indices[NX*NY*NZ*3];
  int len=NX*NY*NZ;
  #endif

  //set the initial location
  MPI_Offset pos;
  #ifdef RESTARTGENERALINDICES
  pos=0;
  #else
  pos=PROCID*NX*NY*NZ*(3*sizeof(int));  
  #endif
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  //read them
  MPI_File_read( cFile, indices, 3*len, MPI_INT, &status );

  //convert to local
  for(ic=0;ic<len;ic++)
    {
      gix=indices[ic*3+0];
      giy=indices[ic*3+1];
      giz=indices[ic*3+2];
      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);
      indices[ic*3+0]=ix;
      indices[ic*3+1]=iy;
      indices[ic*3+2]=iz;
    }

  /***** then primitives in the same order ******/

  //new location in the second block
  #ifdef RESTARTGENERALINDICES
  pos=TNX*TNY*TNZ*(3*sizeof(int));
  ldouble pout[TNX*TNY*TNZ*NV];
  #else
  pos=TNX*TNY*TNZ*(3*sizeof(int)) + PROCID*NX*NY*NZ*(NV*sizeof(ldouble)); 
  ldouble pout[NX*NY*NZ*NV];
  #endif
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  //so far manually

  MPI_File_read( cFile, pout, len*NV, MPI_LDOUBLE, &status );
 
  //rewriting to p
  int ppos;
  for(ic=0;ic<len;ic++)
    {
      ix=indices[ic*3+0];
      iy=indices[ic*3+1];
      iz=indices[ic*3+2];

      ppos=ic*NV;

      if(if_indomain(ix,iy,iz))
	{
	  fill_geometry(ix,iy,iz,&geom);

	  PLOOP(iv)
	    set_u(p,iv,ix,iy,iz,pout[ppos+iv]);

	  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	}
    }

  MPI_File_close( &cFile );
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  return 0;
}
							    
/*********************************************/
/*********************************************/
/*********************************************/
/* prints avg files */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_avgfile(ldouble t, char* folder,char* prefix)
{
  #ifdef MPI

  fprint_avgfile_mpi(t,folder,prefix);

  #else

  fprint_avgfile_bin(t,folder,prefix); 

  #endif
  
  return 0;
}

/*********************************************/
/*********************************************/

int //parallel output to a single file
fprint_avgfile_mpi(ldouble t, char* folder, char* prefix)
{
  #ifdef MPI
  char bufor[250];
  //header
  if(PROCID==0)
    {
      sprintf(bufor,"%s/%s%04d.head",folder,prefix,nfout2);
      fout1=fopen(bufor,"w"); 
      sprintf(bufor,"## %5d %10.6e %10.6e %10.6e\n",nfout2,t-avgtime,t,avgtime);
      fprintf(fout1,"%s",bufor);
      fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfout2);

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

  int rc = MPI_File_open( MPI_COMM_WORLD, bufor, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", bufor );fflush(stdout); exit(-1);
    }

  //set the initial location
  MPI_Offset pos;
  if(PROCID==0) pos=0;
  else
    pos=PROCID*NX*NY*NZ*(3*sizeof(int)+(NV+NAVGVARS)*sizeof(ldouble));
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  int ix,iy,iz,iv,idx[3];
  int gix,giy,giz;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  
	  idx[0]=gix;idx[1]=giy;idx[2]=giz;
	  MPI_File_write( cFile, idx, 3, MPI_INT, &status );
	  MPI_File_write( cFile, &get_uavg(pavg,0,ix,iy,iz), NV+NAVGVARS, MPI_LDOUBLE, &status );
	}

  MPI_File_close( &cFile );

  #endif
  return 0;
}

/*********************************************/
/*********************************************/

int //serial binary output
fprint_avgfile_bin(ldouble t, char* folder,char *prefix)
{
  char bufor[250];
  
  //header
  if(PROCID==0)
    {
      sprintf(bufor,"%s/%s%04d.head",folder,prefix,nfout2);
      fout1=fopen(bufor,"w"); 
      sprintf(bufor,"## %5d %10.6e %10.6e %10.6e\n",nfout2,t-avgtime,t,avgtime);
      fprintf(fout1,"%s",bufor);
      fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfout2);
  fout1=fopen(bufor,"wb"); 

  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  fwrite(&gix,sizeof(int),1,fout1);
	  fwrite(&giy,sizeof(int),1,fout1);
	  fwrite(&giz,sizeof(int),1,fout1);
	  fwrite(&get_uavg(pavg,0,ix,iy,iz),sizeof(ldouble),NV+NAVGVARS,fout1);
	}

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
fread_avgfile(int nout1, char* folder,ldouble *pavg, ldouble *dt,ldouble *t)
{
  char bufor[250];

  #ifdef MPI
  
  fread_avgfile_mpi(nout1,folder,pavg,dt,t);
  
  #else //no MPI

  fread_avgfile_bin(nout1,folder,pavg,dt,t);

  #endif
  
  return 0;
}

/*********************************************/
/*********************************************/


/*********************************************/
/*********************************************/

int 
fread_avgfile_bin(int nout1, char *folder,ldouble *pavg, ldouble *dt,ldouble *t)
{
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz;
  char fname[40],fnamehead[40];

  sprintf(fname,"%s/avg%04d.dat",folder,nout1);
#ifdef MPI
  sprintf(fnamehead,"%s/../0/avg%04d.head",folder,nout1);
#else
  sprintf(fnamehead,"%s/avg%04d.head",folder,nout1);
#endif
  

  FILE *fdump;

  /***********/
  //header file
  fdump=fopen(fnamehead,"r");

  //reading parameters, mostly time
  int intpar[5];
  ldouble ldpar[5];
  ret=fscanf(fdump,"## %d %lf %lf %lf\n",&intpar[0],&ldpar[0],&ldpar[1],&ldpar[2]);
  if(PROCID==0) printf("avg file (%s) read no. %d at times: %.6e to %.6e (dt=%.6e)\n",
	 fname,intpar[0],ldpar[0],ldpar[1],ldpar[2]); 
  
  *t=.5*(ldpar[0]+ldpar[1]);
  *dt=ldpar[2];
  fclose(fdump);
 
  /***********/
  //body file

  fdump=fopen(fname,"rb");
 
  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;
  for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ret=fread(&gix,sizeof(int),1,fdump);
      ret=fread(&giy,sizeof(int),1,fdump);
      ret=fread(&giz,sizeof(int),1,fdump);
      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);
      ret=fread(&get_uavg(pavg,0,ix,iy,iz),sizeof(ldouble),NV+NAVGVARS,fdump);
    }

  fclose(fdump);

  return 0;
}

/*********************************************/
/*********************************************/

int 
fread_avgfile_mpi(int nout1, char *folder,ldouble *pavg, ldouble *dt,ldouble *t)
{
   #ifdef MPI
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz;
  char fname[40],fnamehead[40];

  sprintf(fname,"%s/avg%04d.dat",folder,nout1);
  sprintf(fnamehead,"%s/avg%04d.head",folder,nout1);

  FILE *fdump;

  /***********/
  //header file
  fdump=fopen(fnamehead,"r");

  //reading parameters, mostly time
  int intpar[5];
  ldouble ldpar[5];
  ret=fscanf(fdump,"## %d %lf %lf %lf\n",&intpar[0],&ldpar[0],&ldpar[1],&ldpar[2]);
  if(PROCID==0) printf("avg file (%s) read no. %d at times: %.6e to %.6e (dt=%.6e)\n",
	 fname,intpar[0],ldpar[0],ldpar[1],ldpar[2]); 
  *t=.5*(ldpar[0]+ldpar[1]);
  *dt=ldpar[2];
  fclose(fdump);

  /***********/
  //body file
  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

  int rc = MPI_File_open( MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", fname );fflush(stdout); exit(-1);
    }

  //set the initial location
  MPI_Offset pos;
  if(PROCID==0) pos=0;
  else
    pos=PROCID*NX*NY*NZ*(3*sizeof(int)+(NV+NAVGVARS)*sizeof(ldouble));
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  
	  MPI_File_read( cFile, &gix, 1, MPI_INT, &status );
	  MPI_File_read( cFile, &giy, 1, MPI_INT, &status );
	  MPI_File_read( cFile, &giz, 1, MPI_INT, &status );

	  mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);

	  MPI_File_read( cFile, &get_uavg(pavg,0,ix,iy,iz), NV+NAVGVARS, MPI_LDOUBLE, &status );

	}

  MPI_File_close( &cFile );
#endif

  return 0;
}

/*********************************************/
/*********************************************/
/*********************************************/
/* wrapper for coordinate output */
/*********************************************/
/*********************************************/
/*********************************************/
int fprint_coordfile(char* folder,char* prefix)
{
#if (COORDOUTPUT==1)
  fprint_coordBL(folder,prefix);
#endif
#if (COORDOUTPUT==2)
  fprint_coordBL_shell(folder,prefix);
#endif

  return 0;
}


/*********************************************/
/*********************************************/
/*********************************************/
/* prints BL coordinates,  */
/*********************************************/
/*********************************************/
/*********************************************/
int fprint_coordBL(char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%sBL.dat",folder,prefix);
   FILE* fout1=fopen(bufor,"w");

   int ix,iy,iz,iv;
   ldouble pp[NV];
   for(iz=0;iz<NZ;iz++)
     {
       for(iy=0;iy<NY;iy++)
	 {
	   for(ix=0;ix<NX;ix++)
	     {
	       struct geometry geom,geomBL;
	       fill_geometry(ix,iy,iz,&geom);
	       fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	       ldouble r=geomBL.xx;
	       ldouble th=geomBL.yy;
	       ldouble ph=geomBL.zz;
	     
	       fprintf(fout1,"%d %d %d ",ix,iy,iz);

	       fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph);

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
/* prints BL polar/azimuthal coordinates on a shell  */
/*********************************************/
/*********************************************/
/*********************************************/
int fprint_coordBL_shell(char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%sBL.dat",folder,prefix);
   FILE* fout1=fopen(bufor,"w");

   int ix,iy,iz,iv;
   ldouble pp[NV];

   ix=NX-1;

   for(iz=0;iz<NZ;iz++)
     {
       for(iy=0;iy<NY;iy++)
	 {
	   struct geometry geom,geomBL;
	   fill_geometry(ix,iy,iz,&geom);
	   fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	   ldouble r=geomBL.xx;
	   ldouble th=geomBL.yy;
	   ldouble ph=geomBL.zz;
	     
	   fprintf(fout1,"%d %d %d ",ix,iy,iz);

	   fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph);

	   fprintf(fout1,"\n");
	 }
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }

/*********************************************/
/*********************************************/
/*********************************************/
/* wrapper for simple output */
int fprint_simplefile(ldouble t, int nfile, char* folder,char* prefix)
{
#if (SIMOUTPUT==1)
  fprint_simplecart(t,nfile,folder,prefix);
#endif

#if (SIMOUTPUT==2)
  fprint_simplesph(t,nfile,folder,prefix);
#endif

#if (SIMOUTPUT==3)
  fprint_simplebondi(t,nfile,folder,prefix);
#endif

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
int fprint_simplecart(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);
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
	       trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomout);
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

	       fprintf(fout1,"%.5e %.5e ",rho,temp);

	       fprintf(fout1,"%.5e %.5e %.5e ",vx,vy,vz);

	       fprintf(fout1,"%.5e ",volume);

	       #ifdef RADIATION
	       ldouble Rtt,ehat,Rij[4][4];
	       ldouble ugas[4],Fx,Fy,Fz;
	       if(doingavg==0)
		{
		  calc_ff_Rtt(pp,&Rtt,ugas,&geomout);
		  ehat=-Rtt;  
		  calc_Rij(pp,&geomout,Rij); //calculates R^munu in OUTCOORDS
		  indices_2221(Rij,Rij,geomout.gg);	      							  
		}
	      else
		{
		  ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  int i,j;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		}

	       //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
		{
		  ldouble r=geomsph.xx;
		  ldouble th=geomsph.yy;
		  ldouble ph=geomsph.zz;

		  Rij[2][0]*=r;
		  Rij[3][0]*=r*sin(th);

		  Fx = sin(th)*cos(ph)*Rij[1][0] 
		    + cos(th)*cos(ph)*Rij[2][0]
		    - sin(ph)*Rij[3][0];

		  Fy = sin(th)*sin(ph)*Rij[1][0] 
		    + cos(th)*sin(ph)*Rij[2][0]
		    + cos(ph)*Rij[3][0];

		  Fz = cos(th)*Rij[1][0] 
		    - sin(th)*Rij[2][0];
		}
	       
	      fprintf(fout1,"%.5e %.5e %.5e %.5e",endenGU2CGS(ehat),fluxGU2CGS(Fx),fluxGU2CGS(Fy),fluxGU2CGS(Fz));
#endif

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
/* prints in ASCII & BL coordinates,  */
/*********************************************/
/*********************************************/
/*********************************************/
int fprint_simplesph(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

   fout1=fopen(bufor,"w");
  
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
	       struct geometry geom,geomBL;
	       fill_geometry(ix,iy,iz,&geom);
	       fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	       ldouble dx[3];
	       dx[0]=get_size_x(ix,0);
	       dx[1]=get_size_x(iy,1);
	       dx[2]=get_size_x(iz,2);
	       ldouble gdet=geom.gdet;
	       ldouble volume=dx[0]*dx[1]*dx[2]*gdet;
	       trans_pall_coco(pp, pp, MYCOORDS,BLCOORDS, geom.xxvec,&geom,&geomBL);

	       ldouble rho=rhoGU2CGS(pp[RHO]);
	       ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
	       ldouble vel[4]={0,pp[VX],pp[VY],pp[VZ]};	
	       conv_vels(vel,vel,VELPRIM,VEL4,geomBL.gg,geomBL.GG);						  
	       ldouble r=geomBL.xx;
	       ldouble th=geomBL.yy;
	       ldouble ph=geomBL.zz;
	     
	       fprintf(fout1,"%d %d %d ",ix,iy,iz);

	       fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph);

	       fprintf(fout1,"%.5e %.5e ",rho,temp);

	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",vel[0],vel[1],vel[2],vel[3]);

	       fprintf(fout1,"%.5e ",volume);// (13)

	       #ifdef RADIATION
	       ldouble Rtt,ehat,Rij[4][4];
	       ldouble ugas[4],Fx,Fy,Fz;
	       ldouble Gi[4],Giff[4]={0.,0.,0.,0.};
	       ldouble Gic[4],Gicff[4]={0.,0.,0.,0.};
	       if(doingavg==0)
		{
		  calc_ff_Rtt(pp,&Rtt,ugas,&geomBL);
		  ehat=-Rtt;  
		  calc_Rij(pp,&geomBL,Rij); //calculates R^munu in OUTCOORDS
		  indices_2221(Rij,Rij,geomBL.gg);

		  //four fource
		  calc_Gi(pp,&geomBL,Gi,1); 
		  boost2_lab2ff(Gi,Giff,pp,geomBL.gg,geomBL.GG);
                  #ifdef COMPTONIZATION
		  ldouble kappaes=calc_kappaes(pp,&geomBL);
		  calc_Compt_Gi(pp,&geomBL,Gic,ehat,temp,kappaes,vel);
		  boost2_lab2ff(Gic,Gicff,pp,geomBL.gg,geomBL.GG);
                  #endif 
		  
		  //test
		  //calc_Gi(pp,&geomBL,Gi,0); 
		}
	      else
		{
		  ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  int i,j;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 
		  for(j=0;j<4;j++)
		    Giff[j]=get_uavg(pavg,AVGGHAT(j),ix,iy,iz);
                  #ifdef COMPTONIZATION
		  for(j=0;j<4;j++)
		    Gicff[j]=get_uavg(pavg,AVGGHATCOMPT(j),ix,iy,iz);
		  #endif		  
		}
	       
	       //flux
	       Fx=Rij[1][0];
	       Fy=Rij[2][0];
	       Fz=Rij[3][0];
	       	       
	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",endenGU2CGS(ehat),fluxGU2CGS(Fx),fluxGU2CGS(Fy),fluxGU2CGS(Fz)); //(14) - (17)
	       
	       ldouble conv=kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.)*CCC; //because (cE-4piB) in non-geom
	       fprintf(fout1,"%.5e %.5e ",Giff[0]*conv,Gicff[0]*conv);
#endif

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
/* prints in ASCII & BL coordinates,  */
/*********************************************/
/*********************************************/
/*********************************************/
int fprint_simplebondi(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

   fout1=fopen(bufor,"w");
  
   /***********************************/  
   /** writing order is fixed  ********/  
   /***********************************/  
 
   int ix,iy,iz,iv;
   ldouble pp[NV];
   struct struct_of_state state;
   for(iz=0;iz<NZ;iz++)
     {
       for(iy=0;iy<NY;iy++)
	 {
	   for(ix=0;ix<NX;ix++)
	     {
	       for(iv=0;iv<NV;iv++)
		 pp[iv]=get_u(p,iv,ix,iy,iz);
	       struct geometry geom,geomBL;
	       fill_geometry(ix,iy,iz,&geom);
	       fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	       trans_pall_coco(pp, pp, MYCOORDS,BLCOORDS, geom.xxvec,&geom,&geomBL);

	       get_state(pp,&geomBL,&state);

	       ldouble r=geomBL.xx;
	       ldouble th=geomBL.yy;
	       ldouble ph=geomBL.zz;
	       ldouble radlum,totlum;
	       calc_lum(r,3,&radlum,&totlum);
	       ldouble luminosity = 1.e-40+radlum/calc_lumEdd()*(rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)*velGU2CGS(1.)*velGU2CGS(1.));
	       ldouble mdotscale = rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.);
	       ldouble mdot = calc_mdot(r,0) * mdotscale / calc_mdotEdd();
	       ldouble Fr=state.Rij[1][0];
	       ldouble convGt=kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.)*CCC; //because (cE-4piB) in non-geom
	       ldouble tausca = kappaGU2CGS(state.kappaes/state.rho) *rhoGU2CGS(state.rho) * lenGU2CGS(r);
	       ldouble tauabs = kappaGU2CGS(state.kappa/state.rho) * rhoGU2CGS(state.rho) *lenGU2CGS(r);
	       //ldouble tausca = (state.kappaes) * (r);
	       //ldouble tauabs = (state.kappa) * (r);


	       fprintf(fout1,"%.5e ",r);

	       fprintf(fout1,"%.5e ",rhoGU2CGS(state.rho));

	       fprintf(fout1,"%.5e ",endenGU2CGS(state.uint));

	       fprintf(fout1,"%.5e ",tempGU2CGS(state.Tgas)); //4

	       fprintf(fout1,"%.5e ",(state.entr));

	       fprintf(fout1,"%.5e ",state.ucon[1]);

	       fprintf(fout1,"%.5e ",mdot); //7
	       	       
	       fprintf(fout1,"%.5e ",endenGU2CGS(state.Ehat)); //8

	       fprintf(fout1,"%.5e ",tempGU2CGS(state.Trad)); //9

	       fprintf(fout1,"%.5e ",(state.radentr)); //10
	       	       
	       fprintf(fout1,"%.5e ",luminosity); //11

	       fprintf(fout1,"%.5e ",state.Giff[0]*convGt); //12

	       fprintf(fout1,"%.5e ",state.Gicff[0]*convGt); //13

	       fprintf(fout1,"%.5e ",tausca); //14

	       fprintf(fout1,"%.5e ",tauabs); //15

	       fprintf(fout1,"\n");


	     }
	 }
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }
