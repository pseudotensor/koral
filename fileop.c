//KORAL - fileop.c
//file operations

#include "ko.h"

/*********************************************/
/* opens files etc. */
/*********************************************/
int 
fprint_openfiles()
{
  fout_totmass=fopen("0log.dat","w");
  char bufor[100];
  sprintf(bufor,"rm dumps/*dat gifs/*");
  int i=system(bufor);

  nfout1=0;
  return 0;
}

/*********************************************/
/* reads restart file */
/* TODO: currently uncompatible! */
/*********************************************/
int 
fread_restartfile(ldouble *t)
{
#ifdef RESTART
  //reading time
  int i,ret;
  ldouble time;
  fout_totmass=fopen("0log.dat","r+");
  for(i=0;i<RESTART_NUM-1;i++)
    {
      ret=fscanf(fout_totmass,"%*f %*f\n");
    }
  ret=fscanf(fout_totmass,"%Lf %*f\n",&time);
  *t=time;
  printf("restart no. %d at time: %Lf\n",RESTART_NUM,time); 
  //reading conserved
  char fname[40];
  sprintf(fname,"dumps/out%04d.dat",RESTART_NUM);
  nfout1=RESTART_NUM+1;
  FILE *frestart=fopen(fname,"r");
  int ix,iy,iz,iv;

  for(iz=-0;iz<NZ+0;iz++)
    {
      for(iy=-0;iy<NY+0;iy++)
	{
	  for(ix=-0;ix<NX+0;ix++)
	    {
	      ldouble uu[10],pp[10];
	      
	      //reading conserved from file
	      ret=fscanf(frestart,"%*f %*f %*f "
			 "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf "
			 "%Lf %Lf %*f %*f %*f %*f %*f %*f %*f "
			 "%*f %*f %*f %*f\n",
			 &uu[0],
			 &uu[1],
			 &uu[2],
			 &uu[3],
			 &uu[4],
			 &uu[5],
			 &uu[6],
			 &uu[7],
			 &uu[8],
			 &uu[9],
			 &pp[0],
			 &pp[1]);

	      //setting conserved
	      for(iv=0;iv<NV;iv++)    
		set_u(u,iv,ix,iy,iz,uu[iv]);

	      //calculating primitives
	      ldouble gg[4][5];
	      pick_g(ix,iy,iz,gg);    

	      //saving primitives
	      for(iv=0;iv<NV;iv++)    
		set_u(p,iv,ix,iy,iz,pp[iv]);	      
	    }
	}
    }
#endif


  return 0;
}

/*********************************************/
/* closes file handles */
/*********************************************/
int 
fprint_closefiles()
{
  fclose(fout_totmass);
  return 0;
}

/*********************************************/
/* prints dumps to files and calls gnuplot */
/*********************************************/
int
fprint_profiles(ldouble t, ldouble totmass)
{
  char bufor[50],bufor2[50];
  sprintf(bufor,"dumps/out%04d.dat",nfout1);
  sprintf(bufor2,"gifs/out%04d.%s",nfout1,IMAGETYPE);  
  fout1=fopen(bufor,"w");

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

  if(totmass!=0)
    fprintf(fout_totmass,"%e %e\n",t,totmass);
  fflush(fout_totmass);


  /**************************/  

#ifdef YZXDUMP
  for(iy=0;iy<NY;iy++)
    {
      for(iz=-gclz*NG;iz<NZ+gcrz*NG;iz++)
	{
	  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
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
						  //within domain:
						  if(if_indomain(ix,iy,iz)==0 && if_outsidegc(ix,iy,iz)==1) continue;

						  calc_primitives(ix,iy,iz);

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

						  xx=xxvecout[1];
						  yy=xxvecout[2];
						  zz=xxvecout[3];						  

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
						  dx[0]=get_size_x(ix,0);

						  calc_primitives(ix,iy,iz);

						  for(iv=0;iv<NV;iv++)
						    {
						      uu[iv]=get_u(u,iv,ix,iy,iz);
						      pp[iv]=get_u(p,iv,ix,iy,iz);
						    }	 

						  ldouble tup[4][4],tlo[4][4];
						  pick_T(tmuup,ix,iy,iz,tup);
						  pick_T(tmulo,ix,iy,iz,tlo);	    
						  ldouble eup[4][4],elo[4][4];
						  pick_T(emuup,ix,iy,iz,eup);
						  pick_T(emulo,ix,iy,iz,elo);

						  //to transform primitives between coordinates if necessary
						  if(MYCOORDS!=OUTCOORDS)
						    {
						      ldouble ggout[4][5],GGout[4][5];
						      struct geometry geomout;
						      calc_g_arb(xxvecout,ggout,OUTCOORDS);
						      calc_G_arb(xxvecout,GGout,OUTCOORDS);
						      fill_geometry_arb(ix,iy,iz,&geomout,OUTCOORDS);

#ifdef RADIATION
						      trans_prad_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,gg,GG,ggout,GGout);
#endif
						      trans_phd_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,gg,GG,ggout,GGout);

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
						  ldouble vx=pp[2];
						  ldouble vy=pp[3];
						  ldouble vz=pp[4];
						  ldouble vrel[4]={0,vx,vy,vz};
						  
						  conv_vels(vrel,vrel,VELPRIM,OUTVEL,gg,GG);
						  vx=vrel[1];
						  vy=vrel[2];
						  vz=vrel[3];
						  ldouble S=pp[5];
						  ldouble p=(GAMMA-1.)*uint;
						  ldouble ut=uu[0]/rho;
						  Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;

#ifdef RADIATION						
#ifdef RADOUTPUTINFF
						  prad_lab2ff(pp,pp,&geom);
#elif defined(RADOUTPUTINZAMO) //to print  radiation primitives in ZAMO
						  prad_lab2ff(pp,pp,&geom);
						  prad_ff2zamo(pp,pp,gg,GG,eup); 
#endif

 						  E=pp[6];
						  Fx=pp[7];
						  Fy=pp[8];
						  Fz=pp[9];
						  Trad=calc_LTE_TfromE(fabs(E));
#endif

						  /******************/
						  /* extra lines to calculate v1...v4 from PROBLEMS/XXX/dump.c */

                                                  #include PR_DUMP
						  /******************/
						  
						  //**********************************************************************
						  //**********************************************************************
						  //**********************************************************************

						  fprintf(fout1,"%.4e %.4e %.4e "
							  "%.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e "
							  "%.7e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e "
							  "%.10e %.10e ",
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
#ifdef CGSOUTPUT
							  rhoGU2CGS(rho),    //14
							  endenGU2CGS(uint),   //15
							  (vx),     //16
							  (vy),     //17
							  (vz),     //18
							  S,      //19
							  endenGU2CGS(E),     //20
							  fluxGU2CGS(Fx),     //21
							  fluxGU2CGS(Fy),     //22
							  fluxGU2CGS(Fz),     //23
#else		    
							  rho,    //14
							  uint, 
							  vx,     //16
							  vy,     //17
							  vz,     //18
							  S,      //19
							  E,      //20
							  Fx,     //21
							  Fy,     //22
							  Fz     //23
#endif
							  );

#ifdef MULTIRADFLUID
						  int irf;
						  for(irf=1;irf<NRF;irf++)
						    {
						      fprintf(fout1,"%.10e %.10e %.10e %.10e ",
#ifdef CGSOUTPUT
							      endenGU2CGS(pp[EE(irf)]),     //20
							      fluxGU2CGS(pp[FX(irf)]),     //21
							      fluxGU2CGS(pp[FY(irf)]),     //22
							      fluxGU2CGS(pp[FZ(irf)]),     //23						      
#else
							      pp[EE(irf)],      //20
							      pp[FX(irf)],     //21
							      pp[FY(irf)],     //22
							      pp[FZ(irf)]     //23
							      );
						    }
#endif
						  
#endif

						  fprintf(fout1,"%.10e %.10e %.10e %.10e\n",
							  v1,     //24
							  v2,     //25
							  v3,     //26 
							  v4      //27 
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
  
				      nfout1++;
				      return 0;

				    }
	  
