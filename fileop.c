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
  int i=system("rm dumps/*dat gifs/*gif");

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
	      //u2p(uu,pp,gg);

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
  sprintf(bufor2,"gifs/out%04d.gif",nfout1);  
  fout1=fopen(bufor,"w");

  int ix,iy,iz,iv;
  int gcl,gcr;

  //whether print ghost cells or not - default values
  gcl=gcr=0;

#ifdef PRINTGC_LEFT
  gcl=1;
#endif
#ifdef PRINTGC_RIGHT
  gcr=1;
#endif

  if(totmass!=0)
    fprintf(fout_totmass,"%Le %Le\n",t,totmass);
  fflush(fout_totmass);


  /**************************/  

#ifdef YZXDUMP
  for(iy=0;iy<NY;iy++)
    {
      for(iz=-gcl*NG;iz<NZ+gcr*NG;iz++)
	{
	  for(ix=-gcl*NG;ix<NX+gcr*NG;ix++)
	    {
#elif defined(YSLICE)
	      for(iy=YSLICE;iy<YSLICE+1;iy++)
		{
		  for(iz=0;iz<NZ;iz++)
		    {
		      for(ix=-gcl*NG;ix<NX+gcr*NG;ix++)
			{
#elif defined(YZSLICE)
			  for(iy=NY/2;iy<NY/2+1;iy++)
			    {
			      for(iz=NZ/2;iz<NZ/2+1;iz++)
				{
				  for(ix=-gcl*NG;ix<NX+gcr*NG;ix++)
				    {
#else
				      for(iz=0;iz<NZ;iz++)
					{
					  for(iy=-0*NG;iy<NY+0*NG;iy++)
					    {
					      for(ix=-gcl*NG;ix<NX+gcr*NG;ix++)
						{
#endif
						  //within domain:
						  if(if_indomain(ix,iy,iz)==0 && if_outsidegc(ix,iy,iz)==1) continue;

						  ldouble mx,my,mz,E,e,xx,yy,zz,phipot,xxx[4],dx[3],vv[10],a0,a1,a2,v1,v2,dphidx,v3,Tgas,Trad,v4,v5,v6,v7,Fx,Fy,Fz;
						  ldouble gg[4][5],GG[4][5];
						  ldouble pp[NV],uu[NV];

						  v1=v2=v3=v4=v5=v6=v7=0.;

						  xx=get_x(ix,0);
						  yy=get_x(iy,1);
						  zz=get_x(iz,2);

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

						  for(iv=0;iv<NV;iv++)
						    {
						      uu[iv]=get_u(u,iv,ix,iy,iz);
						      pp[iv]=get_u(p,iv,ix,iy,iz);
						    }	      

						  ldouble rho=pp[0];
						  ldouble uint=pp[1];
						  ldouble vx=pp[2];
						  ldouble vy=pp[3];
						  ldouble vz=pp[4];
						  ldouble S=pp[5];
						  ldouble p=(GAMMA-1.)*uint;
						  ldouble ut=uu[0]/rho;
						  Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;

#ifdef RADIATION
						  ldouble tup[4][4],tlo[4][4];
						  pick_T(tmuup,ix,iy,iz,tup);
						  pick_T(tmulo,ix,iy,iz,tlo);	    
						  ldouble eup[4][4],elo[4][4];
						  pick_T(emuup,ix,iy,iz,eup);
						  pick_T(emulo,ix,iy,iz,elo);	    

						  prad_lab2ff(pp,pp,gg,GG,tup);
#ifdef RADOUTPUTINZAMO
						  prad_ff2zamo(pp,pp,gg,eup); //to print out radiation primitives in ZAMO
#endif
						  E=pp[6];
						  Fx=pp[7];
						  Fy=pp[8];
						  Fz=pp[9];
						  Trad=calc_LTE_TfromE(E);
#endif

						  /******************/
						  /* extra lines to calculate v1...v4 analytically from PROBLEMS/XXX/anasol.c */
						  
                                                  #include PR_ANASOL

						  v1=vv[0];
						  v2=vv[1];
						  v3=vv[2];
						  v4=vv[3];
						  /******************/

						  /******************/
						  /* extra lines to calculate v1...v4 from PROBLEMS/XXX/dump.c */

                                                  #include PR_DUMP
						  /******************/


						  //**********************************************************************
						  //**********************************************************************
						  //**********************************************************************

						  fprintf(fout1,"%.4Le %.4Le %.4Le "
							  "%.7Le %.7Le %.7Le %.7Le %.7Le %.7Le %.7Le %.7Le %.7Le "
							  "%.7Le %.10Le %.10Le %.10Le %.10Le %.10Le %.10Le %.10Le %.10Le "
							  "%.10Le %.10Le %.10Le %.10Le %.10Le %.10Le\n",
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
							  velGU2CGS(vx*sqrtl(gg[1][1])),     //16
							  velGU2CGS(vy*sqrtl(gg[2][2])),     //17
							  velGU2CGS(vz*sqrtl(gg[3][3])),     //18
							  S,      //19
							  endenGU2CGS(E),     //20
							  fluxGU2CGS(Fx),     //21
							  fluxGU2CGS(Fy),     //22
							  fluxGU2CGS(Fz),     //23
#else		    
							  rho,    //14
							  uint, 
							  vx*sqrtl(gg[1][1]),     //16
							  vy*sqrtl(gg[2][2]),     //17
							  vz*sqrtl(gg[3][3]),     //18
							  S,      //19
							  E,      //20
							  Fx,     //21
							  Fy,     //22
							  Fz,     //23
#endif


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
	  
