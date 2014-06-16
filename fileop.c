//KORAL - fileop.c
//file operations

#include "ko.h"

/*********************************************/
/* adds up current quantities to the pavg array
/*********************************************/
int
save_avg(ldouble dt)
{
  int ix,iy,iz,iv,ii;

#pragma omp parallel for private(ix,iy,iz,iv) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 
      
      ldouble avg[NV+NAVGVARS];
      p2avg(ix,iy,iz,avg);

      for(iv=0;iv<NV+NAVGVARS;iv++)
	{
	  set_uavg(pavg,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iz)+avg[iv]*dt);
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
/* prints scalar quantities to scalars.dat
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
/* prints radial profiles to radNNNN.dat
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
/* prints radial profiles to anarelradNNNN.dat
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
									  for(irf=0;irf<NRF;irf++)
									    {

									      //						      irf=0;
									      E+=pp[EE(irf)];
									      Fx+=pp[FX(irf)];
									      Fy+=pp[FY(irf)];
									      Fz+=pp[FZ(irf)];
									      //						      break;
									    }

									  #ifndef NCOMPTONIZATION
									  Trad=calc_LTE_TfromE(E);
									  #else
									  Trad=calc_ncompt_Thatrad(pp,&geomout,E);
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
  #ifdef RESOUTPUT_ASCII //ASCII output

  fprint_restartfile_ascii(t,folder);

  #else //binary output

  #ifdef MPI

  #ifdef OUTPUTPERCORE

  fprint_restartfile_bin(t,folder); 

  #else

  fprint_restartfile_mpi(t,folder);

  #endif

  #else 

  fprint_restartfile_bin(t,folder); 

  #endif

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

  //body
  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);

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
    pos=PROCID*NX*NY*NZ*(3*sizeof(int)+NV*sizeof(ldouble));
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  //what is below should be combined into less MPI_File_writes

  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  
	  MPI_File_write( cFile, &gix, 1, MPI_INT, &status );
	  MPI_File_write( cFile, &giy, 1, MPI_INT, &status );
	  MPI_File_write( cFile, &giz, 1, MPI_INT, &status );
	  MPI_File_write( cFile, &get_u(p,0,ix,iy,iz), NV, MPI_LDOUBLE, &status );
	  /*
	  //somehow hangs up
	  MPI_File_iwrite( cFile, &gix, 1, MPI_INT, &req );
	  MPI_File_iwrite( cFile, &giy, 1, MPI_INT, &req );
	  MPI_File_iwrite( cFile, &giz, 1, MPI_INT, &req );
	  MPI_File_iwrite( cFile, &get_u(p,0,ix,iy,iz), NV, MPI_LDOUBLE, &req );
	  */
	}

  MPI_File_close( &cFile );

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
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  fwrite(&gix,sizeof(int),1,fout1);
	  fwrite(&giy,sizeof(int),1,fout1);
	  fwrite(&giz,sizeof(int),1,fout1);
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
//serial writing in ascii per core 
int
fprint_restartfile_ascii(ldouble t, char* folder)
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

  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);
  fout1=fopen(bufor,"w");

  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  fprintf(fout1,"%d %d %d %.2e %.2e %.2e ",gix,giy,giz,get_x(ix,0),get_x(iy,0),get_x(iz,0));
	  for(iv=0;iv<NV;iv++)
	    {
	      pp[iv]=get_u(p,iv,ix,iy,iz);
	      fprintf(fout1,"%.12e ",pp[iv]);
	    }
	  fprintf(fout1,"\n");
	}

  fclose(fout1);

  //sprintf(bufor,"cp %s/res%04d.dat %s/reslast.dat",folder,nfout1,folder);
  sprintf(bufor,"rm %s/reslast.dat",folder);
  iv=system(bufor);
  sprintf(bufor,"ln -s res%04d.dat %s/reslast.dat",nfout1,folder);
  iv=system(bufor);

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

  #ifdef RESOUTPUT_ASCII

  ret=fread_restartfile_ascii(nout1,folder,t);

  #else //binary output

  #ifdef MPI

  #ifdef OUTPUTPERCORE //each process dumps independent files
  
  ret=fread_restartfile_bin(nout1,folder,t);

  #else 

  ret=fread_restartfile_mpi(nout1,folder,t);

  #endif

  #else //no MPI 

  ret=fread_restartfile_bin(nout1,folder,t);
  
  #endif
  #endif
  
  return ret;
}

/*********************************************/
/*********************************************/

int 
fread_restartfile_ascii(int nout1, char *folder, ldouble *t)
{
  //opening dump file
  int ret;
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
  
  //header
  FILE *fdump=fopen(fnamehead,"r");

  if(fdump==NULL) return 1; //request start from scratch

  //reading parameters, mostly time
  int intpar[6];
  ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],t,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
  printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	 fname,intpar[0],*t,intpar[2],intpar[3],intpar[4],intpar[5]); 

  nfout1=intpar[0]+1; //global file no.
  nfout2=intpar[1]; //global file no. for avg
  fclose(fdump);

  //body
  fdump=fopen(fname,"r");

  int ix,iy,iz,iv,i,ic,gix,giy,giz;

  //reading primitives from file
  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;
  for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ret=fscanf(fdump,"%d %d %d %*f %*f %*f ",&gix,&giy,&giz);
      for(i=0;i<NV;i++)
      ret=fscanf(fdump,"%lf ",&pp[i]);
      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);
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
  for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ret=fread(&gix,sizeof(int),1,fdump);
      ret=fread(&giy,sizeof(int),1,fdump);
      ret=fread(&giz,sizeof(int),1,fdump);
      ret=fread(pp,sizeof(ldouble),NV,fdump);

      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);

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

#ifndef RESTARTGENERALINDICES
  if(PROCID==0) pos=0;
  else
    pos=PROCID*NX*NY*NZ*(3*sizeof(int)+NV*sizeof(ldouble));
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  for(tiz=0;tiz<NZ;tiz++)
    for(tiy=0;tiy<NY;tiy++)
      for(tix=0;tix<NX;tix++)
#else //going through the whole file, useful when changing number of cores
  pos=0;
  MPI_File_seek( cFile, pos, MPI_SEEK_SET );
  for(tiz=0;tiz<TNZ;tiz++)
    for(tiy=0;tiy<TNY;tiy++)
      for(tix=0;tix<TNX;tix++)
#endif
	{
	  MPI_File_read( cFile, &gix, 1, MPI_INT, &status );
	  MPI_File_read( cFile, &giy, 1, MPI_INT, &status );
	  MPI_File_read( cFile, &giz, 1, MPI_INT, &status );

	  mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);

	  MPI_File_read( cFile, pp, NV, MPI_LDOUBLE, &status );

	  if(if_indomain(ix,iy,iz))
	    {
	      fill_geometry(ix,iy,iz,&geom);
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,pp[iv]);
	      p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
	}

  MPI_File_close( &cFile );
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
  #ifdef RESOUTPUT_ASCII

  fprint_avgfile_ascii(t,folder,prefix);

  #else //binary output

  #ifdef MPI
  
  #ifdef OUTPUTPERCORE //each process dumps independent files
  
  fprint_avgfile_bin(t,folder,prefix); 

  #else 

  fprint_avgfile_mpi(t,folder,prefix);

  #endif

  #else

  fprint_avgfile_bin(t,folder,prefix); 

  #endif
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
  
  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  
	    MPI_File_write( cFile, &gix, 1, MPI_INT, &status );
	    MPI_File_write( cFile, &giy, 1, MPI_INT, &status );
	    MPI_File_write( cFile, &giz, 1, MPI_INT, &status );
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

//serial writing in ascii per core 
int
fprint_avgfile_ascii(ldouble t, char* folder,char *prefix)
{
  char bufor[250];

  //header
  if(PROCID==0)
    {
      sprintf(bufor,"%s/%s%04d.head",prefix,folder,nfout2);
      fout1=fopen(bufor,"w");
      sprintf(bufor,"## %5d %10.6e %10.6e %10.6e\n",nfout2,t-avgtime,t,avgtime);
      fprintf(fout1,"%s",bufor);
      fclose(fout1);
    }

  sprintf(bufor,"%s/%s%04d.dat",prefix,folder,nfout2);
  fout1=fopen(bufor,"w");

  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  fprintf(fout1,"%d %d %d ",gix,giy,giz);
	  for(iv=0;iv<NV+NAVGVARS;iv++)
	    {
	      fprintf(fout1,"%.12e ",get_uavg(pavg,iv,ix,iy,iz));
	    }
	  fprintf(fout1,"\n");
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

  #ifdef RESOUTPUT_ASCII

  fread_avgfile_ascii(nout1,folder,pavg,dt,t);

  #else //binary output

  #ifdef MPI
  
  #ifdef OUTPUTPERCORE //each process dumps independent files
  
  fread_avgfile_bin(nout1,folder,pavg,dt,t);

  #else //MPI-IO, each process writes in parallel to the same file

  fread_avgfile_mpi(nout1,folder,pavg,dt,t);
  
  #endif

  #else //no MPI

  fread_avgfile_bin(nout1,folder,pavg,dt,t);

  #endif
  #endif
  
  return 0;
}

/*********************************************/
/*********************************************/

int 
fread_avgfile_ascii(int nout1, char *folder,ldouble *pavg, ldouble *dt,ldouble *t)
{
  //opening dump file
  int ret;
  char fname[40],fnamehead[40];


  sprintf(fname,"%s/avg%04d.dat",folder,nout1);
#ifdef MPI
  sprintf(fnamehead,"%s/../0/avg%04d.head",folder,nout1);
#else
  sprintf(fnamehead,"%s/avg%04d.head",folder,nout1);
#endif
   
  
  //header
  FILE *fdump=fopen(fnamehead,"r");

  if(fdump==NULL) return 1; //request start from scratch

  //reading parameters, mostly time
  int intpar[5];
  ldouble ldpar[5];
  ret=fscanf(fdump,"## %d %lf %lf %lf\n",&intpar[0],&ldpar[0],&ldpar[1],&ldpar[2]);
  if(PROCID==0) printf("avg file (%s) read no. %d at times: %.6e to %.6e (dt=%.6e)\n",
	 fname,intpar[0],ldpar[0],ldpar[1],ldpar[2]); 
  *t=.5*(ldpar[0]+ldpar[1]);
 *dt=ldpar[2];
 
  fclose(fdump);

  //body
  fdump=fopen(fname,"r");

  int ix,iy,iz,iv,i,ic,gix,giy,giz;

  //reading primitives from file
  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;
  for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ret=fscanf(fdump,"%d %d %d ",&gix,&giy,&giz);
      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);
      for(i=0;i<NV+NAVGVARS;i++)
	ret=fscanf(fdump,"%lf ",&get_uavg(pavg,i,ix,iy,iz));
    }

  fclose(fdump);

  return 0;
}

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
/* wrapper for simple output */
int fprint_simplefile(ldouble t, int nfile, char* folder,char* prefix)
{
#if (SIMOUTPUT==1)
  fprint_simplecart(t,nfile,folder,prefix);
#endif

#if (SIMOUTPUT==2)
  fprint_simplesph(t,nfile,folder,prefix);
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
	       trans_pmhd_coco(pp, pp, MYCOORDS,BLCOORDS, geom.xxvec,&geom,&geomBL);

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

	       fprintf(fout1,"%.5e %.5e %.5e ",vel[1],vel[2],vel[3]);

	       fprintf(fout1,"%.5e ",volume);

	       #ifdef RADIATION
	       ldouble Rtt,ehat,Rij[4][4];
	       ldouble ugas[4],Fx,Fy,Fz;
	       if(doingavg==0)
		{
		  calc_ff_Rtt(pp,&Rtt,ugas,&geomBL);
		  ehat=-Rtt;  
		  calc_Rij(pp,&geomBL,Rij); //calculates R^munu in OUTCOORDS
		  indices_2221(Rij,Rij,geomBL.gg);	      							  
		}
	      else
		{
		  ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  int i,j;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		}
	       
	       //flux
	       Fx=Rij[1][0];
	       Fy=Rij[2][0];
	       Fz=Rij[3][0];

	       //four fource
	       ldouble Gi[4],Giff[4]={0.,0.,0.,0.};
	       ldouble Gic[4],Gicff[4]={0.,0.,0.,0.};
	       calc_Gi(pp,&geomBL,Gi); 
	       boost2_lab2ff(Gi,Giff,pp,geomBL.gg,geomBL.GG);
	       #ifdef COMPTONIZATION
	       ldouble kappaes=calc_kappaes(rho,temp,-1.,-1.,-1.);  
	       calc_Compt_Gi(pp,&geomBL,Gic,ehat,temp,kappaes,vel);
	       boost2_lab2ff(Gic,Gicff,pp,geomBL.gg,geomBL.GG);
               #endif 

	       /*
	       if(ix==NX/2 && iy==NY/2)
		 {
		   print_4vector(Gi);
		   print_4vector(Giff);
		   print_4vector(Gic);
		   print_4vector(Gicff);
		   getchar();
		 }
	       */
	       	       
	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",endenGU2CGS(ehat),fluxGU2CGS(Fx),fluxGU2CGS(Fy),fluxGU2CGS(Fz));
	       
	       ldouble conv=kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.);
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
