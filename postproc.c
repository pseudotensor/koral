//KORAL - misc.c
//routines for postprocessing 
//used both on the go and separately

#include "ko.h"

/*********************************************/
/* calculates radial profiles - L(r) etc. */
/* uses mostly primitives, but may use avg quantities as well */
/* however, this part would work only when postprocessing with */
/* ./avg, otherwise pavg hold non-normalized numbers */

//surface density (2) (column)
//rest mass flux (3)
//rho-weighted minus radial velocity (4)
//<u_phi>* (5)
//Keplerian u_phi (6)
//abs optical depth (7)
//tot optical depth (8)
//net accretion rate at given radius (9)
//inflow accretion rate at given radius (10)
//outflow accretion rate at given radius (11)
//luminosity at given radius opt. thin (12)
//location of the photosphere (13)
//total mhd energy flux (14)
//outflowin mhd energy flux (15)
//jet mhd energy flux (16)
//total rad energy flux (17)
//outflowin rad energy flux (18)
//jet rad energy flux (19)
//outflowin mass flux (20)
//jet mass flux (21)
//luminosity at given radius everywhere (22)
//surface density in the inflow (23)
//rho-weighted minus radial velocity in the inflow (24)
//opt thin mhd energy flux (25)
//opt. thin rad energy flux (26)
//opt. thin mass flux (27)
//rho-weighted q-theta (28)
//rho-weighted temperature (29)
//rho-weighted magn.field angle <sqrt(grr gphph)b^r b^ph> / <bsq> (30)
//scale-height (31)
//rho-weighted beta (32)
//rho-wighted prad/pgas (33)
//alpha (34)
//rad. viscosity energy flux (35)
//rho-weighted minus radial velocity in the outflow (36)
//conserved flux rho ur transformed to BLCOORDS (37)
//conserved flux rho ur in MYCOORDS (38)
//conserved flux rhout+Trt in MYCOORDS (39)
//conserved flux for Rrt int MYCOORDS(40)
//surface density of energy = int (Ttt+rhout+Rtt) dz (41)
//rho-weighted radial velocity in the jet (42)
//magnetic flux in the jet (43)
//kinetic flux in the jet (44)
//radial velocity close to the axis (45)                                                                                                                                                                              //Bernoulli close to the axis (46)                                                                                                                                                                                                                                    


/*********************************************/
int calc_radialprofiles(ldouble profiles[][NX])
{
  //adjust NRADPROFILES in problem.h

  //calculates scale-height etc.
  calc_avgs_throughout();     
  int ix;

  //search for appropriate radial index
#pragma omp parallel for
  for(ix=0;ix<NX;ix++)
    {
      int iy,iz,iv,i,j;
      ldouble x0[3],x0l[3],x0r[3],xm1[3],xp1[3];
      ldouble dx0, dxm2, dxm1, dxp1, dxp2;  
      ldouble xx[4],xxBL[4],dx[3],dxph[3],dxcgs[3],mdot,rho,rhouconrcons,uint,temp,ucon[4],utcon[4],ucon3[4];
      ldouble rhoucont,enden,rhouconr,Tij[4][4],Tij22[4][4],Rij[4][4],Rviscij[4][4],Trt,Fluxx[NV],Rrt,Rviscrt,bsq,bcon[4],bcov[4],Qtheta;
      ldouble Trtmagn,Trtkin;
      ldouble ucov[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5],Ehat;
      ldouble tautot,tautotloc,tauabs,tauabsloc;
      ldouble avgsums[NV+NAVGVARS][NX];
      ldouble Bangle1,Bangle2,brbphi;
      ldouble Sigmagdet;
      ldouble diffflux[NV];
      ldouble fd_u0[NV],fd_up1[NV],fd_up2[NV],fd_um1[NV],fd_um2[NV];
      ldouble fd_p0[NV],fd_pp1[NV],fd_pp2[NV],fd_pm1[NV],fd_pm2[NV],fd_pm3[NV],fd_pp3[NV];
      ldouble fd_pl[NV],fd_pr[NV],fd_plm1[NV],fd_prm1[NV],fd_plp1[NV],fd_prp1[NV];
      ldouble fd_ul[NV],fd_ur[NV],fd_ulm1[NV],fd_urm1[NV],fd_ulp1[NV],fd_urp1[NV];
      ldouble du[NV],dul[NV],dur[NV],aaa[12],ahd,arad;
      int injet;

      //vertically integrated/averaged profiles

      for(iv=0;iv<NRADPROFILES;iv++)
	profiles[iv][ix]=0.;

      //outside horizon?
      struct geometry geomBLtemp;
      fill_geometry_arb(ix,0,0,&geomBLtemp,BLCOORDS);
      if(geomBLtemp.xx<=1.1*rhorizonBL) continue; //to avoid working inside horizon
      
      Bangle1=Bangle2=0.;
      Sigmagdet=0.;
      ldouble jetsigma=0.;

      for(iv=0;iv<NAVGVARS;iv++)
	avgsums[iv][ix]=0.;

      tautot=tauabs=0.;

      // #ifdef BHDISK_PROBLEMTYPE
      if(NZ==1) //phi-symmetry
	{
	  iz=0;
	  for(iy=0;iy<NY;iy++)
	    {
	      //metric
	      pick_g(ix,iy,iz,gg);
	      pick_G(ix,iy,iz,GG);

	      struct geometry geom;
	      fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	      struct geometry geomm1;
	      fill_geometry_arb(ix-1,iy,iz,&geomm1,MYCOORDS);

	      struct geometry geomp1;
	      fill_geometry_arb(ix+1,iy,iz,&geomp1,MYCOORDS);

	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	      struct geometry geomBLm1;
	      fill_geometry_arb(ix-1,iy,iz,&geomBLm1,BLCOORDS);

	      struct geometry geomBLp1;
	      fill_geometry_arb(ix+1,iy,iz,&geomBLp1,BLCOORDS);

	      struct geometry geoml;
	      fill_geometry_face(ix,iy,iz,0,&geoml);

	      struct geometry geomr;
	      fill_geometry_face(ix+1,iy,iz,0,&geomr);	
	      
	      struct geometry geomBLl;
	      fill_geometry_face_arb(ix,iy,iz,0,&geomBLl,BLCOORDS);

	      struct geometry geomBLr;
	      fill_geometry_face_arb(ix+1,iy,iz,0,&geomBLr,BLCOORDS);	      

	      ldouble gdetuBL=geomBL.gdet;

	      //	      if(PROCID==0) {printf("%d %d %d > %e %e %e\n",ix,iy,iz,geom.gg[1][1],geom.gdet,geomBL.gdet);getch();}
	      //exit(1);

#if (GDETIN==0) //gdet out of derivatives
	      gdetuBL=1.;
#endif

	      //coordinates
	      get_xx(ix,iy,iz,xx);	      
	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	      ldouble dxph[3];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,0);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;
	      dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);

	      //primitives at the cell - either averaged or original, in BL or MYCOORDS
	      for(iv=0;iv<NV;iv++)
		{
		  pp[iv]=get_u(p,iv,ix,iy,iz);
		}
	      
	      
	      //primitives at radial neighbours
	      for(i=0;i<NV;i++)
		{
		  fd_p0[i]=pp[i];
		  fd_pp1[i]=get_u(p,i,ix+1,iy,iz);
		  fd_pm1[i]=get_u(p,i,ix-1,iy,iz);
		  fd_pm2[i]=get_u(p,i,ix-2,iy,iz);
		  fd_pp2[i]=get_u(p,i,ix+2,iy,iz);
		}

	      //internal coordinates
	      x0[0]=get_x(ix,0);

	      x0l[0]=get_xb(ix,0);
	      xm1[0]=get_x(ix-1,0);
	      x0l[1]=xm1[1]=get_x(iy,1); 
	      x0l[2]=xm1[2]=get_x(iz,2);

	      x0r[0]=get_xb(ix+1,0);
	      xp1[0]=get_x(ix+1,0);
	      x0r[1]=xp1[1]=get_x(iy,1);
	      x0r[2]=xp1[2]=get_x(iz,2);

	      dx0=get_size_x(ix,0);    
	      dxm1=get_size_x(ix-1,0);    
	      dxp1=get_size_x(ix+1,0);  

	      //interpolation to get left/right biased interpolated primtives
	      avg2point(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2,0); 
	      //INTORDER==1!
	      avg2point(fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pp2,fd_plp1,fd_prp1,dxm1,dx0,dxp1,dxp2,dxp2,0);   
	      avg2point(fd_pm2,fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_plm1,fd_prm1,dxm2,dxm2,dxm1,dx0,dxp1,0);   
		  
	      //to BL, res-files and primitives in avg in MYCOORDS
	      trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);

	      //transforming interpolated primitives to BL
	      trans_pall_coco(fd_pm1,fd_pm1,MYCOORDS,BLCOORDS,xx,&geomm1,&geomBLm1);
	      trans_pall_coco(fd_pp1,fd_pp1,MYCOORDS,BLCOORDS,xx,&geomp1,&geomBLp1);
	      trans_pall_coco(fd_pl,fd_pl,MYCOORDS,BLCOORDS,xx,&geoml,&geomBLl);
	      trans_pall_coco(fd_pr,fd_pr,MYCOORDS,BLCOORDS,xx,&geomr,&geomBLr);
	      trans_pall_coco(fd_plp1,fd_plp1,MYCOORDS,BLCOORDS,xx,&geoml,&geomBLl);
	      trans_pall_coco(fd_prm1,fd_prm1,MYCOORDS,BLCOORDS,xx,&geomr,&geomBLr);

	      if(doingavg)
		{
		  rho=get_uavg(pavg,RHO,ix,iy,iz);
		  uint=get_uavg(pavg,UU,ix,iy,iz);
		  bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		  bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
		  bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
		  bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
		  bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
		  utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);

		  rhouconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
		  rhoucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz);

		  rhouconrcons=get_uavg(pavg,AVGRHOURDIFF,ix,iy,iz);
		  
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)

		  Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(i,j),ix,iy,iz)
		    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(i,j),ix,iy,iz)
		    + get_uavg(pavg,AVGBSQUCONUCOV(i,j),ix,iy,iz)
		    - get_uavg(pavg,AVGBCONBCOV(i,j),ix,iy,iz); 

		  for(i=0;i<NV;i++)
		    Fluxx[i]=get_uavg(pavg,AVGFLUXXL(i),ix,iy,iz);

		  Trt=Tij[1][0];
		  enden = Tij[0][0]+rhoucont;

		  Trtmagn = get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		    - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

		  Trtkin = get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz);

#ifdef RADIATION  
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 

		  enden += Rij[0][0];
		  Rrt = Rij[1][0];
		  Ehat = get_uavg(pavg,AVGEHAT,ix,iy,iz);

		  int derdir[3]={0,0,0};
		  calc_Rij_visc(pp,&geomBL,Rviscij,derdir);
      
		  Rviscrt = Rviscij[1][0];
#endif
		  
		  //no need of transforming interpolated primitives to BL, already there
		 
		}
	      else //on the go from the primitives
		{ 
		  rho=pp[0];
		  uint=pp[1];
		  utcon[1]=pp[2];
		  utcon[2]=pp[3];
		  utcon[3]=pp[4];
		  
#ifdef MAGNFIELD
		  calc_bcon_prim(pp,bcon,&geomBL);
		  indices_21(bcon,bcov,geomBL.gg); 
		  bsq = dot(bcon,bcov); 
#endif

		  conv_vels_both(utcon,utcon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
		  rhouconr=rhouconrcons=rho*utcon[1];

		  calc_Tij(pp,&geomBL,Tij22);
		  indices_2221(Tij22,Tij,geomBL.gg);

		  Trt = Tij[1][0];

		  Trtmagn = bsq*utcon[1]*ucov[0] - bcon[1]*bcov[0];
		  Trtkin =  rho*utcon[1]*ucov[0];
		  enden = Tij[0][0] + rho*utcon[0];

#ifdef RADIATION
		  calc_Rij(pp,&geomBL,Rij);
		  indices_2221(Rij,Rij,geomBL.gg);

		  Rrt = Rij[1][0];

		  ldouble Rtt,uconr[4];
		  calc_ff_Rtt(&get_u(p,0,ix,iy,iz),&Rtt,uconr,&geomBL);
		  Ehat=-Rtt; 	
		  enden+=Rij[0][0];

		  int derdir[3]={0,0,0}; 
		  calc_Rij_visc(pp,&geomBL,Rviscij,derdir);
      
		  Rviscrt = Rviscij[1][0];
#endif

		}

	      ldouble muBe;
	      muBe=-(Trt+rhouconr)/rhouconr;
	      #ifdef RADIATION
	      muBe+=-Rrt/rhouconr;
	      #endif
	      
	      int isjet;
	      if(muBe>0.05 && (xxBL[2]<M_PI/4. || xxBL[2]>3.*M_PI/4.))
		isjet=1;
	      else isjet=0;

	      //estimating the diffusive flux
	      //to conserved
	      p2u(fd_pl,fd_ul,&geomBLl);
	      p2u(fd_pr,fd_ur,&geomBLr);
	      p2u(fd_plp1,fd_ulp1,&geomBLr);
	      p2u(fd_prm1,fd_urm1,&geomBLl);

	      //gradient of conserved
	      PLOOP(iv)
	      {
		//getting rid of gdetu in conserved - integrated with gdet lateron
		dul[iv]=fd_ul[iv]-fd_urm1[iv];
		dur[iv]=fd_ulp1[iv]-fd_ur[iv];
		du[iv]=.5*(dul[iv]+dur[iv]);
		du[iv]/=gdetuBL;

		//test - right face only
		/*
		du[iv]=dur[iv];
		du[iv]/=gdetuBL; //de facto substracting ahd (getd_r (rho ut_rR - rho ut_rL))
		*/
	      }
	      
	      //test - right face only
	      /*
	      double ff1[NV],ff2[NV];
	      f_flux_prime(fd_pr,0,ix+1,iy,iz,ff1,0); 
	      f_flux_prime(fd_plp1,0,ix+1,iy,iz,ff2,0); 
	      rhouconr=.5*(ff1[0]+ff2[0])/gdetuBL; //de facto plotting gdet_r * rhour_r
	      Trt=.5*(ff1[1]-ff1[0]+ff2[1]-ff2[0])/gdetuBL; 
	      //Trt=.5*(ff1[1]+ff2[1])/gdetuBL; 
	      */

	      //wavespeeds
	      calc_wavespeeds_lr_pure(pp,&geomBL,aaa);
	      ahd=my_max(fabs(aaa[0]),fabs(aaa[1]));
	      arad=my_max(fabs(aaa[6]),fabs(aaa[7]));

	      //test
	      calc_wavespeeds_lr_pure(fd_pm1,&geomBLm1,aaa);
	      ahd=my_max(ahd,my_max(fabs(aaa[0]),fabs(aaa[1])));
	      arad=my_max(arad,my_max(fabs(aaa[6]),fabs(aaa[7])));

	      calc_wavespeeds_lr_pure(fd_pp1,&geomBLp1,aaa);
	      ahd=my_max(ahd,my_max(fabs(aaa[0]),fabs(aaa[1])));
	      arad=my_max(arad,my_max(fabs(aaa[6]),fabs(aaa[7])));

	      //diffusive flux
	      PLOOP(iv)
	      {
		if(iv<NVMHD) diffflux[iv]=-0.5*ahd*du[iv];
		else diffflux[iv]=-0.5*arad*du[iv];
	      }

	      ldouble pregas = GAMMAM1*uint;
	      ldouble ptot = pregas;
	      #ifdef RADIATION
	      ldouble prerad = Ehat/3.;
	      ptot+=prerad;
	      #endif
	      
	      //alpha 
	      boost22_lab2ff(Tij22,Tij22,pp,geomBL.gg,geomBL.GG);
	      ldouble alpha=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij[1][3]/ptot;

	      //temperature
	      ldouble temp=calc_PEQ_Tfromurho(uint,rho);

	      //angular velocity
	      ldouble Omega=utcon[3]/utcon[0];

	      //MRI resolution parameter
	      ldouble qtheta=calc_Qtheta(ix,iy,iz);//2.*M_PI/Omega/dx[1]*fabs(bcon[2])/sqrt(rho);

	      //to calculate magn. field angle
	      //bsq and brbphi taken from avg if neeeded
	      ldouble bconfake[4],bcovfake[4],bsqfake;
	      calc_angle_brbphibsq(ix,iy,iz,&brbphi,&bsqfake,bconfake,bcovfake);
	      Bangle1+=rho*brbphi*dxph[1];
	      Bangle2+=rho*bsq*dxph[1];

	      //optical depths
	      ldouble k1,k2,k3,k4;
	      ldouble tauabsloc = utcon[0]*calc_kappa(pp,&geomBL,&k1,&k2,&k3,&k4);
	      ldouble tautotloc = utcon[0]*calc_kappaes(pp,&geomBL);
	      tautot+=tautotloc*dxph[1];
	      tauabs+=tauabsloc*dxph[1];	

	      //alpha (34) (column)
	      profiles[32][ix]+=alpha*rho*dxph[1];

	      //surface density (2) (column)
	      profiles[0][ix]+=rho*dxph[1];

	      //surface energy density (41)
	      profiles[39][ix]+=enden*dxph[1];

	      //numerator of scale height (31) (column)
	      //profiles[29][ix]+=rho*dxph[1]*pow(tan(fabs(M_PI/2.-xxBL[2])),2.);

	      //surface density in the inflow (23)
	      if(utcon[1]<0.)
		profiles[21][ix]+=rho*dxph[1];

	      //rho-weighted q-theta (28)
	      profiles[26][ix]+=rho*qtheta*dxph[1];
	       
	      //rho-weighted temperature (29)
	      profiles[27][ix]+=rho*temp*dxph[1];
	      
	      //rho-weighted beta (32)
	      ldouble prermhd = GAMMAM1*uint;
	      #ifdef RADIATION
	      prermhd+=Ehat/3.;
	      #endif
	      ldouble ibeta=bsq/2./(prermhd+bsq/2.);
	      profiles[30][ix]+=rho*ibeta*dxph[1];

	      //rho-weighted prad/pgas (33)
	      #ifdef RADIATION
	      profiles[31][ix]+=rho*prerad/pregas*dxph[1];
	      #else
	      profiles[31][ix]+=0.;
	      #endif
              
    
	      //rest mass flux (3)
	      profiles[1][ix]+=-rhouconr*dx[1]*dx[2]*geomBL.gdet;


	      //conserved flux (rhour) transformed to BLCOORDS (may be imprecise) (37)
	      //profiles[35][ix]+=-rhouconrcons*dx[1]*dx[2];

	      
	      //test
	      profiles[35][ix]+=get_uavg(pavg,AVGFLUXYL(EE0),ix,iy,iz)*get_size_x(ix,0)*dx[2]+
		get_uavg(pavg,AVGFLUXYL(UU),ix,iy,iz)*get_size_x(ix,0)*dx[2];
	      
	      /*
	      if(iy==0) {printf("%d %d > %e %e\n", ix,iy,
			       get_uavg(pavg,AVGFLUXYL(EE0),ix,iy,iz),
				get_uavg(pavg,AVGFLUXYL(UU),ix,iy,iz));getchar();}
	      */

	      //conserved flux (gdet rhour) in MYCOORDS (38)
	      profiles[36][ix]+=-Fluxx[RHO]*get_size_x(iy,1)*dx[2];

	      //conserved flux for Trt (39) in MYCOORDS 
	      profiles[37][ix]+=-Fluxx[UU]*get_size_x(iy,1)*dx[2];
	      
	      //temporary surface density to normalize what is above
	      Sigmagdet+=rho*dx[1]*dx[2]*geomBL.gdet;

	      //total mhd energy flux (14)
	      profiles[12][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;

	      //magnetic mhd energy flux in jet (43) 
	      if(isjet==1)
		profiles[41][ix]+=(-Trtmagn)*dx[1]*dx[2]*geomBL.gdet;

	      //kinetic mhd energy flux in jet (43)  
	      if(isjet==1)
		profiles[42][ix]+=(-Trtkin)*dx[1]*dx[2]*geomBL.gdet;



	      
	      //opt thin mhd energy flux (25)
	      if(tautot<1.)
		profiles[23][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;

	      //outflowin mhd energy flux (15)
	      if(utcon[1]>0.)
		profiles[13][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;

	      //jet mhd energy flux (16)
	      if(isjet==1)
		profiles[14][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;

	      //jet gas velocity (42)
	      if(isjet==1)
		{
		  profiles[40][ix]+=rhouconr*dx[1]*dx[2]*geomBL.gdet;
		  jetsigma+=rho*dx[1]*dx[2]*geomBL.gdet;
		}

	      //gas velocity near the axis (45)
	      if(iy==NCCORRECTPOLAR+1 || iy==(NY-NCCORRECTPOLAR-2))
		profiles[43][ix]+=0.5*utcon[1];

	      //Bernoulli near the axis (46)
	      if(iy==NCCORRECTPOLAR+1 || iy==(NY-NCCORRECTPOLAR-2))
		profiles[44][ix]+=0.5*muBe;

	      //total rad energy flux (17)
	      profiles[15][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;

	      //conserved flux for Rrt in MYCOORDS (40)
	      profiles[38][ix]+=-Fluxx[EE0]*get_size_x(iy,1)*dx[2];

	      //rad viscosity energy flux (35)
	      profiles[33][ix]+=(-Rviscrt)*dx[1]*dx[2]*geomBL.gdet;

	      //opt. thin rad energy flux (26)
	      if(tautot<1.)
		profiles[24][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;

	      //outflowin rad energy flux (18)
	      if(utcon[1]>0.)
		profiles[16][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;

	      //jet rad energy flux (19)
	      if(isjet==1)
		profiles[17][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;
	      
	      //outflowin mass flux (20)
	      if(utcon[1]>0.)
		profiles[18][ix]+=(-rhouconr)*dx[1]*dx[2]*geomBL.gdet;

	      //opt. thin mass flux (27)
	      if(tautot<1.)
		profiles[25][ix]+=(-rhouconr)*dx[1]*dx[2]*geomBL.gdet;

	      //jet mass flux (21)
	      if(isjet==1)
		profiles[19][ix]+=(-rhouconr)*dx[1]*dx[2]*geomBL.gdet;

	      //rho-weighted minus radial velocity (4)
	      //profiles[2][ix]+=-rhouconr*dxph[1];

	      //rho-weighted minus radial velocity in the inflow (24)
	      //if(utcon[1]<0.)
	      //profiles[22][ix]+=-rhouconr*dxph[1];

	      //abs optical depth (7)
	      profiles[5][ix]=tauabs;	

	      //tot optical depth (8)
	      profiles[6][ix]=tautot;

	      for(iv=0;iv<NV+NAVGVARS;iv++)
		avgsums[iv][ix]+=get_uavg(pavg,iv,ix,iy,iz)*dxph[1];

	      //<(rho+u+bsq/2)u^r><u_phi> (5)
	      //profiles[3][ix]+=get_uavg(pavg,AVGWUCON(1),ix,iy,iz)*get_uavg(pavg,AVGUCOV(3),ix,iy,iz)*dxph[1];

	      if(doingavg)
		profiles[3][ix]+=get_uavg(pavg,AVGRHOUCOV(3),ix,iy,iz)*dxph[1];
	      else
		profiles[3][ix]+=rho*ucov[3]*dxph[1];
	    }


	 
	  //normalizing by sigma
	  profiles[22][ix]/=profiles[21][ix];
	  profiles[26][ix]/=profiles[0][ix];
	  profiles[27][ix]/=profiles[0][ix];
	  //profiles[29][ix]/=profiles[0][ix];
	  profiles[30][ix]/=profiles[0][ix];
	  profiles[31][ix]/=profiles[0][ix];
	  profiles[32][ix]/=profiles[0][ix];
	  //profiles[29][ix]=sqrt(profiles[29][ix]); //scale height
	  profiles[29][ix]=scaleth_otg[ix]; //scale height
	  profiles[40][ix]/=jetsigma;
	  
	  Bangle1/=profiles[0][ix];
	  Bangle2/=profiles[0][ix];

	  //rho-weighted magn.field angle -<sqrt(grr gphph)b^r b^ph> / <bsq> (30)
	  profiles[28][ix]=-Bangle1/Bangle2;
	     
	  //normalizing by <(rho+u+bsq/2)u^r>
	  //profiles[3][ix]/=avgsums[AVGWUCON(1)][ix];
	  profiles[3][ix]/=profiles[0][ix];
 
	  //Keplerian u_phi (6)
	  ldouble r=xxBL[1];
	  profiles[4][ix]=((r*r-2.*BHSPIN*sqrt(r)+BHSPIN*BHSPIN)/(sqrt(r*(r*r-3.*r+2.*BHSPIN*sqrt(r)))));  

	  //net accretion rate at given radius (9)
	  profiles[7][ix]=fabs(calc_mdot(xxBL[1],0));
	  //inflow accretion rate at given radius (10)
	  profiles[8][ix]=fabs(calc_mdot(xxBL[1],1));
	  //outflow accretion rate at given radius (11)
	  profiles[9][ix]=fabs(calc_mdot(xxBL[1],2));

	  //to get velocities
	  profiles[2][ix]=profiles[1][ix]/Sigmagdet;
	  profiles[22][ix]=profiles[8][ix]/Sigmagdet;
	  profiles[34][ix]=profiles[9][ix]/Sigmagdet;


	  //luminosity at given radius (12)
	  ldouble radlum,totallum;
	  calc_lum(xxBL[1],0,&radlum,&totallum);
	  profiles[10][ix]=radlum;
	  //luminosity at given radius (22)
	  calc_lum(xxBL[1],1,&radlum,&totallum);
	  profiles[20][ix]=radlum;
	  //location of the photosphere (13)
	  profiles[11][ix]=calc_photloc(ix);

	}

     //#endif
    }

  return 0;
}

/*********************************************/
/* calculates theta profiles  */

//total energy flux (2) (column)

/*********************************************/
int calc_thetaprofiles(ldouble profiles[][NY])
{
  //adjust NTHPROFILES in problem.h
  ldouble rho,uint,bsq,bcon[4],bcov[4],utcon[4],ucov[4],rhouconr,rhoucont;
  ldouble Tij[4][4],Tij22[4][4],Rij[4][4],Trt,Rrt,Ehat,Rviscij[4][4],Rviscrt;
  ldouble pp[NV];

  //choose radius where to extract from
  int ix,i,j,iv,iz;

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

  int iy;
#pragma omp parallel for
  for(iy=0;iy<NY;iy++)
    {
      for(iv=0;iv<NTHPROFILES;iv++)
	profiles[iv][iy]=0.;

      if(NZ==1) //phi-symmetry
	{
	  iz=0;
	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	  struct geometry geom;
	  fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);
	      
	  //primitives at the cell - either averaged or original, in BL or MYCOORDS
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  //to BL, res-files and primitives in avg in MYCOORDS
	  trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);

	  if(doingavg)
	    {
	      rho=get_uavg(pavg,RHO,ix,iy,iz);
	      uint=get_uavg(pavg,UU,ix,iy,iz);
	      bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
	      bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
	      bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
	      bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
	      bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
	      utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      rhouconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	      rhoucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz);
		  
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)

		  Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(i,j),ix,iy,iz)
		    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(i,j),ix,iy,iz)
		    + get_uavg(pavg,AVGBSQUCONUCOV(i,j),ix,iy,iz)
		    - get_uavg(pavg,AVGBCONBCOV(i,j),ix,iy,iz); 

	      Trt=Tij[1][0];

#ifdef RADIATION  
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 

	      Rrt = Rij[1][0];
	      Ehat = get_uavg(pavg,AVGEHAT,ix,iy,iz);

	      int derdir[3]={0,0,0};
	      calc_Rij_visc(pp,&geomBL,Rviscij,derdir);      
	      Rviscrt = Rviscij[1][0];
#endif
		  
	      //no need of transforming interpolated primitives to BL, already there
		 
	    }
	  else
	    { 
	      rho=pp[0];
	      uint=pp[1];
	      utcon[1]=pp[2];
	      utcon[2]=pp[3];
	      utcon[3]=pp[4];
		  
#ifdef MAGNFIELD
	      calc_bcon_prim(pp,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsq = dot(bcon,bcov); 
#endif

	      conv_vels_both(utcon,utcon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	      rhouconr=rho*utcon[1];

	      calc_Tij(pp,&geomBL,Tij22);
	      indices_2221(Tij22,Tij,geomBL.gg);

	      Trt = Tij[1][0];

#ifdef RADIATION
	      calc_Rij(pp,&geomBL,Rij);
	      indices_2221(Rij,Rij,geomBL.gg);

	      Rrt = Rij[1][0];

	      ldouble Rtt,uconr[4];
	      calc_ff_Rtt(&get_u(p,0,ix,iy,iz),&Rtt,uconr,&geomBL);
	      Ehat=-Rtt; 	

	      int derdir[3]={0,0,0}; 
	      calc_Rij_visc(pp,&geomBL,Rviscij,derdir);
      
	      Rviscrt = Rviscij[1][0];
#endif
	      
	    }
	  


	  int isconverged;
	  if(fabs(utcon[1])>xxBL[1]/(global_time/2.)) 
	    isconverged=1;
	  else
	    isconverged=0;
	  
	  ldouble fluxconv=fluxGU2CGS(1.); 

	  //total energy luminosity in cgs for converged region (2)
	  profiles[0][iy]=-fluxconv*(Rrt+rhouconr+Trt);
	      
	  //radiative luminosity in cgs for converged region (3)
	  profiles[1][iy]=-fluxconv*(Rrt);

	  //gas velocity (4)
	  profiles[2][iy]=utcon[1];

	  //converged? (5)
	  profiles[3][iy]=(ldouble)isconverged;

	  //optical depths
	  int iix;
	  ldouble tau1,tau2,Rphot;
	  ldouble k1,k2,k3,k4;
	  tau1=tau2=0.;
	  Rphot=-1.;


	  for(iix=NX-1;iix>=0;iix--)
	    {
	      struct geometry geomBL2;
	      fill_geometry_arb(iix,iy,iz,&geomBL2,BLCOORDS);
	      ldouble grr=geomBL2.gg[1][1];
	      //coordinates
	      ldouble dxph[3],dx[0];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,0);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      dxph[0]=dx[0]*sqrt(geomBL2.gg[1][1]);
	      ldouble rho2=get_u(p,RHO,iix,iy,iz);
              ldouble uint2=get_u(p,UU,iix,iy,iz);
              ldouble Tgas2=calc_PEQ_Tfromurho(uint2,rho2);
	      ldouble kabsloc = calc_kappa(&get_u(p,0,iix,iy,iz),&geomBL2,&k1,&k2,&k3,&k4);
	      ldouble kscaloc = calc_kappaes(&get_u(p,0,iix,iy,iz),&geomBL2);
	      if(doingavg)
                {
                  utcon[0]=get_uavg(pavg,AVGRHOUCON(0),iix,iy,iz)/get_uavg(pavg,RHO,iix,iy,iz);
                  utcon[1]=get_uavg(pavg,AVGRHOUCON(1),iix,iy,iz)/get_uavg(pavg,RHO,iix,iy,iz);
                  utcon[2]=get_uavg(pavg,AVGRHOUCON(2),iix,iy,iz)/get_uavg(pavg,RHO,iix,iy,iz);
                  utcon[3]=get_uavg(pavg,AVGRHOUCON(3),iix,iy,iz)/get_uavg(pavg,RHO,iix,iy,iz);
		  conv_vels_both(utcon,utcon,ucov,VELPRIM,VEL4,geomBL2.gg,geomBL2.GG);
		}
	      else
		{
		  //temporary
		  ucov[0]=-1.; ucov[1]=0.;
		}

	      tau1+=-(kabsloc+kscaloc)*(ucov[0]+ucov[1])*sqrt(grr)*dxph[0];
	      tau2=-(kabsloc+kscaloc)*(ucov[0]+ucov[1])*geomBL2.xx*sqrt(grr);
	      if(Rphot<0. && my_max(tau1,tau2)>2./3.)
		Rphot=geomBL2.xx;
	    }

	  //location of photosphere (6)
	  profiles[4][iy]=Rphot;

	  //optical depth along radius between boundaries (7)
	  profiles[5][iv]=tau1;


	}

    }

  return 0;
}


/*********************************************/
/* calculates scalar s - total mass, accretion rate etc. */
/*********************************************/
int calc_scalars(ldouble *scalars,ldouble t)
{
  //adjust NSCALARS in problem.h

  /*********************************************/
  //base for BHDISK problems
  /*********************************************/

  //total mass inside the domain (2nd column)
  scalars[0]=calc_totalmass();

  //accretion rate through horizon (3)
  ldouble mdotscale = rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.);

 ldouble rmdot=rhorizonBL;
#if(PROBLEM==7) //BONDI
  rmdot = RMIN;
#endif
  ldouble mdot=calc_mdot(rmdot,0);
  scalars[1]=-mdot;

  //accretion rate through horizon in Edd. units (6)
  scalars[4]=-mdot*mdotscale/calc_mdotEdd();

  ldouble xx[4],xxBL[4];
  get_xx(NX-1,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

  //luminosities 
  ldouble rlum=15.;
#if(PROBLEM==69) //INJDISK
  rlum=2./3.*DISKRCIR;
#endif
#if(PROBLEM==7) //BONDI
  rlum=RMAX;
#endif
  ldouble radlum,totallum;
  calc_lum(rlum,3,&radlum,&totallum);

  //radiative luminosity everywhere (4)
  scalars[2]=radlum*mdotscale*CCC0*CCC0/calc_lumEdd();

  //total energy at infinity (rho ur + Trt + Rrt) (12)
  calc_lum(rlum,1,&radlum,&totallum);
  scalars[10]=totallum;
  

  //mri resolution parameter Q_theta (7) at rmri
  ldouble rmri=xxBL[1]/2.;
#if(PROBLEM==69) //INJDISK
  rmri=20.;
#endif
  scalars[5]=calc_resmri(rmri);

  //rho-weighted temperature at rtemp (8)
  ldouble rtemp=15.;
#if(PROBLEM==69) //INJDISK
  rtemp=20.;
#endif
#if(PROBLEM==7) //BONDI
  rtemp=RMAX;
#endif
  scalars[6]=calc_meantemp(rtemp);

  //magnetic flux through horizon parameter (5)
  ldouble Bfluxquad;
  ldouble Bflux;
  calc_Bflux(rhorizonBL,0.,&Bflux,&Bfluxquad);
  scalars[3]=Bflux;

  //MAD parameter (9)
  scalars[7]=Bflux/sqrt(fabs(mdot))*sqrt(4.*M_PI)/2.;

  //MAD-quad parameter (11)
  scalars[9]=(Bfluxquad/sqrt(fabs(mdot))*sqrt(4.*M_PI)/2.);

  //scaleheight at rtemp (10)
  ldouble rscale=15.;
  scalars[8]=calc_scaleheight(rscale);
 

  /*********************************************/
  //Tgas Trad Egas Erad for testing Comptonization
  /*********************************************/

#ifdef TESTCOMPTINSCALARS4FLAT
  ldouble pp[NV];
  int iv;
  PLOOP(iv) pp[iv]=get_u(p,iv,0,0,0);
  struct geometry geom;
  fill_geometry(0,0,0,&geom);
  ldouble ugas[4],Tgas,Trad1,Trad2,Rtt,Ehatrad;
  int dominates;
  calc_ff_Rtt(pp,&Rtt,ugas,&geom);
  Ehatrad=-Rtt;
  Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
  #ifndef NCOMPTONIZATION
  Trad1=Trad2=calc_LTE_TfromE(Ehatrad);
  #else
  Trad1=calc_LTE_TfromE(Ehatrad);
  Trad2=calc_ncompt_Thatrad(pp,&geom,Ehatrad);
  #endif

  scalars[0]=Tgas;
  scalars[1]=Trad1;
  scalars[2]=Trad2;
  scalars[3]=pp[UU];
  scalars[4]=Ehatrad;
  scalars[5]=pp[NF0];
#endif

  /*********************************************/
  //L1 ERRRORS for some problems
  /*********************************************/

#ifdef CALCL1_RMHDWAVE
  //temporarily here: L1 error for RMHDWAVE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      calc_primitives(i,0,0,0,0);
      ldouble xx=get_x(i,0);
      ldouble dx=get_size_x(i,0);
      ldouble myrho=RHOZERO+DRRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DRIM/DRRE*sin(OMRE*t-KK*xx));
      //L1 in rho:
      L1+=fabs(get_u(p,RHO,i,0,0)-myrho)*dx;
    }
  scalars[0]=L1;///(ldouble)NX;
#endif

#ifdef CALCL1_HDWAVE
  //temporarily here: L1 error for HDWAVE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      calc_primitives(i,0,0,0,0);
      ldouble xx=get_x(i,0);
      ldouble om=1./CC*2.*Pi;
      ldouble myrho=RHOZERO*(1.+AAA*cos(KK*xx-om*t));
      ldouble myuint=UINT*(1.+GAMMA*AAA*cos(KK*xx-om*t));
      ldouble mycs=1./CC;
      ldouble myvx=AAA*cos(KK*xx-om*t)*mycs;
      //L1 in rho:
      L1+=fabs(get_u(p,0,i,0,0)-myrho);
    }
  scalars[0]=L1/(ldouble)NX;
#endif

#ifdef CALCL1_HUBBLE
  //temporarily here: L1 error for HUBBLE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      ldouble xx=get_x(i,0);      
      ldouble myrho=RHO0 / (1.+VPRIME*t);
      ldouble myuint=UINT0 / pow(1.+VPRIME*t,GAMMA);
      ldouble myvx=VPRIME*xx / (1.+VPRIME*t);
      //L1 in rho:
      L1+=fabs(get_u(p,0,i,0,0)-myrho);
    }
  scalars[0]=L1/(ldouble)NX;
#endif

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//integrates mass in the domain
ldouble
calc_totalmass()
{
  int ix,iy,iz;
  ldouble xx[4],dx[3],mass,rho,gdet;
  
  mass=0.;
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      get_xx(ix,iy,iz,xx);
	      dx[0]=get_size_x(ix,0);
	      dx[1]=get_size_x(iy,1);
	      dx[2]=get_size_x(iz,2);
	      gdet=calc_gdet(xx);
	      rho=get_u(p,0,ix,iy,iz);
	      mass+=rho*dx[0]*dx[1]*dx[2]*gdet;
	    }
	}
    }
  return mass;
}
	  
//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates the Eddington mass accretion rate
ldouble
calc_mdotEdd()
{
#if (PROBLEM==7) //spherical bondi
  ldouble mcgs=2.23/16.*1e18*MASS; //g/s
#else
  ldouble mcgs=1.09649*2.23e18*MASS*(0.057/etaNT); //g/s \propto 1/etaNT(a)
#endif

  return mcgs;
}

//**********************************************************************
//**********************************************************************
//*********************************************************************
//calculates the Eddington luminosity
ldouble
calc_lumEdd()
{
  ldouble Lcgs=1.25e38*MASS; //erg/s

  return Lcgs;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates luminosity by integrating positive flux from the axis up to tau=1 surface
//normalized to total sphere, taken at radius radius
int
calc_lum(ldouble radius,int type,ldouble *radlum, ldouble *totallum)
{


    int ix,iy,iz,iv,i,j;
    ldouble xx[4],xxBL[4],dx[3],pp[NV],Rrt,rhour,Tij[4][4],Trt;
    ldouble Rij[4][4],Rtt,ehat,ucongas[4];
    ldouble tautot[3],tau=0.;
    ldouble gdet;

 
  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }
  if(ix==NX) ix=NX-1;

  ldouble lum=0.,jet=0.;

  if(NY==1 && NZ==1) //spherical symmetry
    {
      iz=0; 
      iy=0;

      for(iv=0;iv<NV;iv++)
	pp[iv]=get_u(p,iv,ix,iy,iz);

      struct geometry geomBL;
      fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);
  
      if(doingavg)
	{
	  PLOOP(iv)
	    pp[iv]=get_uavg(pavg,iv,ix,iy,iz);

	  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	  ldouble ucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	  ldouble uconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);		  
	  rhour=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	  Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
		    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
		    + get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		    - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

#ifdef RADIATION

	  if(type==0) //R^r_t outside photosphere
	    {
	      Rrt=0.;
	    }
	  else if(type==1) //R^r_t everywhere
	    {
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
	      Rrt=Rij[1][0];// + ehat*uconr);
	      if(Rrt<0.) Rrt=0.;
	    }
	  else if(type==2) //R^r_t everywhere in outflow
	    {
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
	      Rrt=Rij[1][0];// + ehat*uconr);
	      if(uconr<0. || Rrt<0.) Rrt=0.;
	    }
	  else
	    Rrt=0.;
#else
	  Rrt=0.;
#endif

	  lum=-geomBL.gdet*Rrt*4.*M_PI;
	  jet=geomBL.gdet*(Trt+rhour+Rrt)*4.*M_PI;
	}
      else
	{
	  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	  
	  ucongas[1]=pp[2];
	  ucongas[2]=pp[3];
	  ucongas[3]=pp[4];	      
	  conv_vels(ucongas,ucongas,VELPRIM,VEL4,geom.gg,geom.GG);

	  rhour = pp[RHO]*ucongas[1];
	  
	  calc_Tij(pp,&geom,Tij);
	  indices_2221(Tij,Tij,geom.gg);
	  Trt=Tij[1][0];


#ifdef RADIATION	      
	  if(type==0) //R^r_t outside photosphere
	    {
	      Rrt=0.;
	    }
	  else if(type==1) //sum of positive R^r_t everywhere
	    {
	      //calc_ff_Rtt(pp,&Rtt,ucongas,&geom);
	      //ehat=-Rtt;
	      calc_Rij(pp,&geom,Rij); 
	      indices_2221(Rij,Rij,geom.gg);
	      Rrt=Rij[1][0];// + ehat*ucongas[1];
	      if(Rrt<0.)
	      Rrt=0.;
	    }
	  else if(type==2) //R^r_t in the outflow region
	    {
	      //calc_ff_Rtt(pp,&Rtt,ucongas,&geom);
	      //ehat=-Rtt;
	      calc_Rij(pp,&geom,Rij); 
	      indices_2221(Rij,Rij,geom.gg);
	      Rrt=Rij[1][0];// + ehat*ucongas[1];
	      if(Rrt<0. || ucongas[1]<0.)
		Rrt=0.;
	    }
	  else if(type==3) //sum of R^r_t everywhere
	    {
	      calc_Rij(pp,&geom,Rij); 
	      indices_2221(Rij,Rij,geom.gg);
	      Rrt=Rij[1][0];// + ehat*ucongas[1];	      
	    }
	  else
	    Rrt=0.;
#else
	  Rrt=0.;
#endif

	  lum=-geom.gdet*Rrt*4.*M_PI;
	  jet=geom.gdet*(rhour+Trt+Rrt)*4.*M_PI;

	  //printf("%e %e %e %e\n",xxBL[1],Rrt,geom.gdet,lum);
	}

      *radlum=lum;
      *totallum=jet;
      return 0.;
    }
  else if(NZ==1) //phi-symmetry only
    {
      iz=0; 
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=2.*M_PI;
	  gdet=geom.gdet;
	  ldouble dxph[3],dxBL[3];
	  ldouble xx1[4],xx2[4];
	  xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	  xx2[0]=0.;xx2[1]=get_xb(ix+1,1);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	  coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	  coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	  dxBL[0]=fabs(xx2[1]-xx1[1]);
	  xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	  xx2[0]=0.;xx2[1]=get_xb(ix,1);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	  coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	  coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	  dxBL[1]=fabs(xx2[2]-xx1[2]);
	  dxBL[2]=2.*M_PI;
	  dxph[0]=dxBL[0]*sqrt(geomBL.gg[1][1]);
	  dxph[1]=dxBL[1]*sqrt(geomBL.gg[2][2]);
	  dxph[2]=dxBL[2]*sqrt(geomBL.gg[3][3]);
	  
	  if(doingavg)
	    {
	      PLOOP(iv)
		pp[iv]=get_uavg(pavg,iv,ix,iy,iz);

	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	      calc_tautot(pp,&geomBL,dxph,tautot);

	      ldouble ucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ldouble uconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);		  
	      rhour=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	      
		  Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
		    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
		    + get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		    - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

	      tau+=ucont*tautot[1];

	      #ifdef RADIATION
	      if(type==0) //R^r_t outside photosphere
		{
		  if(tau>1.) break;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  Rrt=Rij[1][0];
		  if(Rrt<0.) Rrt=0.;
		}
	      else if(type==1) //positive R^r_t everywhere
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		  Rrt=Rij[1][0];// + ehat*uconr);
		  if(Rrt<0.) Rrt=0.;

		  
		}
	      else if(type==2) //R^r_t everywhere in outflow
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		  Rrt=Rij[1][0];// + ehat*uconr);
		  if(uconr<0. || Rrt<0.) Rrt=0.;
		}
	      else if(type==3) //any R^r_t everywhere
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		  Rrt=Rij[1][0];// + ehat*uconr);
		  
		}
	      else
		Rrt=0.;
	      #else
	      Rrt=0.;
	      #endif

	      lum+=-geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
	      jet+=geomBL.gdet*(rhour+Rrt+Rrt)*dxBL[1]*dxBL[2];
	    }
	  else
	    {
	      

	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	      calc_tautot(pp,&geomBL,dxph,tautot);

	      ucongas[1]=pp[2];
	      ucongas[2]=pp[3];
	      ucongas[3]=pp[4];	      
	      conv_vels(ucongas,ucongas,VELPRIM,VEL4,geom.gg,geom.GG);

	      rhour = pp[RHO]*ucongas[1];
	  
	      calc_Tij(pp,&geom,Tij);
	      indices_2221(Tij,Tij,geom.gg);
	      Trt=Tij[1][0];

	      tau+=ucongas[0]*tautot[1];
	      
#ifdef RADIATION
	      if(type==0) //R^r_t outside photosphere
		{
		  if(tau>1.) break;	  
		  //trans_prad_coco(pp,pp,MYCOORDS,KERRCOORDS,xx,&geom,&geomBL);
		  //prad_lab2on(pp,pp,&geomBL);
		  calc_Rij(pp,&geom,Rij); 
		  indices_2221(Rij,Rij,geom.gg);
		  Rrt=Rij[1][0];
		  if(Rrt<0.) Rrt=0.;
		}
	      else if(type==1) //R^r_t everywhere
		{
		  //calc_ff_Rtt(pp,&Rtt,ucongas,&geom);
		  //ehat=-Rtt;
		  calc_Rij(pp,&geom,Rij); 
		  indices_2221(Rij,Rij,geom.gg);
		  Rrt=Rij[1][0];// + ehat*ucongas[1];
		  if(Rrt<0.) Rrt=0.;
		}
	      else if(type==2) //R^r_t in the outflow region
		{
		  //calc_ff_Rtt(pp,&Rtt,ucongas,&geom);
		  //ehat=-Rtt;
		  calc_Rij(pp,&geom,Rij); 
		  indices_2221(Rij,Rij,geom.gg);
		  Rrt=Rij[1][0];// + ehat*ucongas[1];
		  if(Rrt<0. || ucongas[1]<0.)
		    Rrt=0.;
		}
	      else
		Rrt=0.;
	      #else
	      Rrt=0.;
	      #endif

	      lum+=-gdet*Rrt*dx[1]*dx[2];
	      jet+=gdet*(rhour+Trt+Rrt)*dx[1]*dx[2];
	    }

	  //#ifdef CGSOUTPUT
	  //never!
	  //Rrt=fluxGU2CGS(Rrt);
	  //dx[1]=lenGU2CGS(dx[1]);
	  //dx[2]=lenGU2CGS(dx[2]);
	  //#endif		  


	  //printf("%e %e %e -> %e\n",Rrt,gdet,dx[1],dx[2],lum);getch();
	}
      
      *radlum=lum;
      *totallum=jet;
      return 0.;
    }
  else


    return -1;
}
//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates MRI resolution parameter Q_theta ar rmri
ldouble
calc_resmri(ldouble radius)
{
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no MRI
#endif

#ifdef MAGNFIELD

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3];
  ldouble qtheta=0.,sigma=0.,rho;
 
  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  if(NZ==1)
    {
      iz=0;
      for(iy=0;iy<NY;iy++)
	{
	  dx[1]=get_size_x(iy,1);
	  rho=get_u(p,RHO,ix,iy,iz);

	  sigma+=rho*dx[1];
	  qtheta+=rho*calc_Qtheta(ix,iy,iz)*dx[1];
	}

      return qtheta/sigma;
    }
  else
    return -1.;

#endif
    return -1.;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates mean temperature at rmri
ldouble
calc_meantemp(ldouble radius)
{
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no disk no cry
#endif


  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3];
  ldouble mtemp=0.,sigma=0.,rho,ugas;
 
  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

   if(NZ==1)
    {
      iz=0;
      for(iy=5;iy<NY-5;iy++)
	{
	  dx[1]=get_size_x(iy,1);
	  rho=get_u(p,RHO,ix,iy,iz);
	  ugas=get_u(p,UU,ix,iy,iz);
	  ldouble temp=calc_PEQ_Tfromurho(ugas,rho);
	  sigma+=rho*dx[1];
	  mtemp+=rho*temp*dx[1];
	}

      return mtemp/sigma;
    }
  else
    return -1.;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates mean temperature at rmri
ldouble
calc_scaleheight(ldouble radius)
{
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no disk no cry
#endif


  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3];
  ldouble mtemp=0.,sigma=0.,rho,ugas;
 
  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  return scaleth_otg[ix];
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates theta corresponding to integrated tau from the axis
ldouble
calc_photloc(int ix)
{
  if(MYCOORDS != BLCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS && MYCOORDS != MKS2COORDS)
    return -1.; //no BH

  ldouble tau=0.,pp[NV],xx[4],xxBL[4],dx[3];

  int iz=0; int iy,iv; 

  if(NZ==1)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  ldouble tautot[3];

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);

	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);
	  dx[0]=dx[0]*sqrt(geom.gg[1][1]);
	  dx[1]=dx[1]*sqrt(geom.gg[2][2]);
	  dx[2]=2.*M_PI*sqrt(geom.gg[3][3]);

	  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	  calc_tautot(pp,&geom,dx,tautot);
	  tau+=tautot[1];
	  if(tau>1.) break;
	}
      return xxBL[2];
    }
  else
    return -1;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates rest mass flux through r=radius within range of thetas
//normalized to 2pi in phi
//type == 0 (net)
//type == 1 (inflow only)
//type == 2 (outflow only)
ldouble
calc_mdot(ldouble radius,int type)
{
  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3],mdot,gdet,rho,rhouconr,ucon[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];

  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  mdot=0.;

  if(NZ==1) //phi-symmetry
    {
      iz=0;
      for(iy=0;iy<NY;iy++)
	{
	  struct geometry geom;
	  fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);
	  
	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	  if(doingavg)
	    {
	      rho=get_uavg(pavg,RHO,ix,iy,iz);
	      ucon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ucon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ucon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      rhouconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);	      	      
	      gdet=geomBL.gdet;
	      ldouble dxph[3];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,1);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,1);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;
	    }
	  else
	    {
	      for(iv=0;iv<NVMHD;iv++)
		pp[iv]=get_u(p,iv,ix,iy,iz);

	      get_xx(ix,iy,iz,xx);
	      dx[0]=get_size_x(ix,0);
	      dx[1]=get_size_x(iy,1);
	      dx[2]=get_size_x(iz,2);
	      pick_g(ix,iy,iz,gg);
	      pick_G(ix,iy,iz,GG);

	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	      

	      rho=pp[0];
	      ucon[1]=pp[2];
	      ucon[2]=pp[3];
	      ucon[3]=pp[4];
	      
	      conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	      rhouconr=rho*ucon[1];
	      gdet=geom.gdet;	
	      dx[1]=dx[1];	     
	      dx[2]=2.*M_PI;
	      /*
	      trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);
	      ldouble dxph[3];
	      conv_vels(ucon,ucon,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	      rhouconr=rho*ucon[1];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,1);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,1);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;
	      gdet=geomBL.gdet;
	      */
	      
	    }

	  if(NY==1)
	    {
	      dx[1]=2.;
	      dx[2]=2.*M_PI;
	    }

	  if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
	    mdot+=gdet*rhouconr*dx[1]*dx[2];	     
	}
    }
  else
    return -1;

  return mdot;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates magnetic flux through r=radius within range of thetas
//normalized to 2pi in phi
//type == 0 (default)
int
calc_Bflux(ldouble radius,int type,ldouble *Bflux, ldouble* Bfluxquad)
{
  *Bflux=*Bfluxquad=0.;
  
  if(MYCOORDS != BLCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS && MYCOORDS != MKS2COORDS)
    {
      return -1.; //no BH
    }

  #ifdef MAGNFIELD

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3],Psi,Psiquad,rho,ucon[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];

  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  Psi=0.;
  Psiquad=0.;

  if(NZ==1) //phi-symmetry
    {
      iz=0;
      for(iy=0;iy<NY;iy++)
	{
	  struct geometry geom;
	  fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);
	  
	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	  

	  if(doingavg)
	    {
	      ldouble bcon[4]={get_uavg(pavg,AVGBCON(0),ix,iy,iz),
			       get_uavg(pavg,AVGBCON(1),ix,iy,iz),
			       get_uavg(pavg,AVGBCON(2),ix,iy,iz),
			       get_uavg(pavg,AVGBCON(3),ix,iy,iz)};
	      
	      ldouble ucon[4];

	      ucon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ucon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ucon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      conv_vels(ucon,ucon,VEL4,VEL4,geomBL.gg,geomBL.GG);

	      //Bcon[1]
	      ldouble Br = bcon[1]*ucon[0] - bcon[0]*ucon[1];
	      
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,1);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,1);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;

	      if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
		Psi+=geomBL.gdet*fabs(Br)*dx[1]*dx[2];
	  
	      if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
		Psiquad+=geomBL.gdet*Br*my_sign(geomBL.yy-M_PI/2.)*dx[1]*dx[2];

	    }
	  else
	    {
	      for(iv=0;iv<NVMHD;iv++)
		pp[iv]=get_u(p,iv,ix,iy,iz);

	      dx[0]=get_size_x(ix,0);
	      dx[1]=get_size_x(iy,1);
	      dx[2]=get_size_x(iz,2);

	      get_xx(ix,iy,iz,xx);
	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	      struct geometry geom;
	      fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	      ldouble Br=pp[B1];

	      dx[1]=dx[1];
	      dx[2]=2.*M_PI;

	  
	      if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
		Psi+=geom.gdet*fabs(Br)*dx[1]*dx[2];
	  
	      if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
		Psiquad+=geom.gdet*Br*my_sign(geomBL.yy-M_PI/2.)*dx[1]*dx[2];
	    }
	  

	}
    }
  else
    return -1;

  *Bflux = Psi;

  *Bfluxquad = Psiquad;
  return 0;
#else
  return -1;
#endif
}

/*********************************************/
//calculates quantites with respect to the averaged solution
//and projects them on the equatorial plane
/*********************************************/
int calc_anarelradialprofiles(ldouble profiles[][NX])
{
  //adjust NRADPROFILES in problem.h

  int ix,iy,iz,iv,i,j;
  ldouble rhoavg,uintavg,tempavg,uconavg[4],utconavg[4],bconavg[4],ppavg[NV],bsqavg;
  ldouble rho,uint,temp,ucon[4],ucov[4],bcon[4],bcov[4],bsq;
  ldouble Tij[4][4],alpha,pp[NV],Rij[4][4],Ehat,urcon[4];

  //loop over radius
  for(ix=0;ix<NX;ix++)
    {
      //vertically integrated/averaged profiles

      for(iv=0;iv<NANARELRADPROFILES;iv++)
	profiles[iv][ix]=0.;
      
 #ifdef BHDISK_PROBLEMTYPE
     if(NZ==1) //phi-symmetry
	{
	  iz=0;
	  for(iy=0;iy<NY;iy++)
	    {
	      struct geometry geom;
	      fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	      
	      ldouble dxph[3],dx[3];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,0);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;
	      dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);
	      
	      /*******************************/
	      //averaged quantities (always VEL4 BL!)
	      rhoavg=get_uavg(pavg,RHO,ix,iy,iz);
	      uintavg=get_uavg(pavg,UU,ix,iy,iz);
	      bsqavg=get_uavg(pavg,AVGBSQ,ix,iy,iz);
	      bconavg[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
	      bconavg[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
	      bconavg[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
	      bconavg[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
	      uconavg[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      uconavg[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      uconavg[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      uconavg[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      conv_vels(uconavg,utconavg,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
	      ppavg[RHO]=rhoavg;
	      ppavg[UU]=uintavg;
	      ppavg[VX]=utconavg[1];
	      ppavg[VY]=utconavg[2];
	      ppavg[VZ]=utconavg[3];
	      //magnetic field etc. empty!

	      /*******************************/
	      //snapshot quantities (always VEL4 BL!)
	      //use old ones instead
	      for(iv=0;iv<NV;iv++)
		  pp[iv]=get_u(p,iv,ix,iy,iz);
	      //to BL     
	      trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,geom.xxvec,&geom,&geomBL);
	      rho=pp[0];
	      uint=pp[1];
	      ucon[1]=pp[2];
	      ucon[2]=pp[3];
	      ucon[3]=pp[4];
	      conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
		  
              #ifdef MAGNFIELD
	      calc_bcon_prim(pp,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsq = dot(bcon,bcov); 
              #endif

              #ifdef RADIATION
	      calc_Rij(pp,&geomBL,Rij);
	      indices_2221(Rij,Rij,geomBL.gg);
	      calc_ff_Rtt(pp,&Ehat,urcon,&geomBL);
	      Ehat*=-1.;
              #endif


	      /*******************************/
	      //calculating alpha
	      
	      //stress energy tensor in lab frame
	      calc_Tij(pp,&geomBL,Tij);

	      //boosting it from lab frame to the average comoving frame
	      //what if there are secular trends? 
	      //one should take avg from around snapshot file
	      //or do as Penna+12 did, i.e., take averaged phi velocity=VZ
	      boost22_lab2ff(Tij,Tij,pp,geomBL.gg,geomBL.GG);

	      //pressure
	      ldouble ptot = GAMMAM1*uint + 1./3.*Ehat;
	      alpha=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij[1][3]/ptot;
	      
	      /*******************************/
	      //vertical integrals
	      //within scale-height
	      if(fabs(geomBL.yy-M_PI/2.) < scaleth_otg[ix])
		{
		  //surface density (2nd column)
		  profiles[2-2][ix]+=rho*dxph[1];

		  //alpha numerator (3) (sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij[1][3]) 
		  profiles[3-2][ix]+=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij[1][3]*rho*dxph[1];
		  
		  //alpha denominator (4) (ptot)
		  profiles[4-2][ix]+=ptot*rho*dxph[1];
		  
		  //direct alpha (5) 
		  profiles[5-2][ix]+=alpha*rho*dxph[1];
		}

	    }

	  //normalizing by sigma
	  profiles[1][ix]/=profiles[0][ix];
	  profiles[2][ix]/=profiles[0][ix];
	  profiles[3][ix]/=profiles[0][ix];
	 
	}

#endif
    }

  return 0;
}
