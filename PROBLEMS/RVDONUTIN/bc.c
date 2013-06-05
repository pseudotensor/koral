//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

/**********************/
//geometries
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecBL[4],xx,yy,zz;

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);

//working in BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,KERRCOORDS);
calc_G_arb(xxvecBL,GGBL,KERRCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);
/**********************/

//radius
if(ix>=NX) //analytical solution within the torus and atmosphere outside
  {
    ldouble podpierd=-(GGBL[0][0]-2.*ELL*GGBL[0][3]+ELL*ELL*GGBL[3][3]);
    ldouble ut=-1./sqrt(podpierd);
    ut/=UTPOT;
    ldouble uint,Vphi,rho,Vr;
    ldouble xx=get_x(ix,0);
    ldouble D,E,W,eps,uT,uphi,uPhi;
    if(ut<-1 || podpierd<0.|| NODONUT ) //outside torus
      {
	iix=NX-1;
	iiy=iy;
	iiz=iz;

	//ambient
	set_hdatmosphere(pp,xxvec,gg,GG,4);


#ifdef RADIATION
	ldouble ppatm[NV];
	ldouble ucon[4];
	
	//pure atmosphere
	set_radatmosphere(ppatm,xxvec,gg,GG,0);

	pp[6]=get_u(p,6,iix,iiy,iiz);
	pp[7]=get_u(p,7,iix,iiy,iiz);
	pp[8]=get_u(p,8,iix,iiy,iiz);
	pp[9]=get_u(p,9,iix,iiy,iiz);
		
	pp[6]=ppatm[6];
	pp[7]=ppatm[7];
	pp[8]=ppatm[8];
	pp[9]=ppatm[9];
	
	
#endif

	/*
	//BL free-fall velocity
	ldouble ucon[4];
	ldouble r=xxvecBL[1];
	ucon[0]=0.;
	ucon[1]=-sqrtl(2./r)*(1.-2./r);
	ucon[2]=ucon[3]=0.;
	conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);
	trans2_coco(xxvecBL,ucon,ucon,BLCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];
	*/

      }
    else
      {
	ldouble h=-1./ut;
	ldouble eps=(h-1.)/GAMMA;
	rho=pow(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.));
	uint=rho*eps;
	uphi=-ELL*ut;
	uT=GGBL[0][0]*ut+GGBL[0][3]*uphi;
	uPhi=GGBL[3][3]*uphi+GGBL[0][3]*ut;
	Vphi=uPhi/uT;
	Vr=0.;

	Vr=-URIN;

	//3-velocity in BL 
	ldouble ucon[4]={0.,Vr,0.,Vphi};
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggBL,GGBL);
   
	pp[2]=ucon[1]; 
	pp[3]=ucon[2];
	pp[4]=ucon[3];
	pp[0]=rho; pp[1]=uint; 


    ldouble P,aaa,bbb;
    P=GAMMAM1*uint;
    //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
    aaa=4.*SIGMA_RAD/3.;
    bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
    ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));

    //    ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;
    ldouble T4=-Sqrt((-4*cbrt(0.666666666666)*P)/naw1 + naw1/(cbrt(2.)*cbrt(9.)*aaa))/2. + Sqrt((4*cbrt(0.666666666666)*P)/naw1 - naw1/(cbrt(2.)*cbrt(9.)*aaa) + (2*bbb)/(aaa*Sqrt((-4*cbrt(0.666666666666)*P)/naw1 + naw1/(cbrt(2.)*cbrt(9.)*aaa))))/2.;

    E=calc_LTE_EfromT(T4);
    Fx=Fy=Fz=0.;
    uint=calc_PEQ_ufromTrho(T4,rho);

#ifdef HDDONUTASWITHRAD
    pp[1]=uint;
#endif

#ifdef RADIATION
    pp[1]=uint;
    pp[6]=E;
    pp[7]=Fx;
    pp[8]=Fy;
    pp[9]=Fz;


   //estimating flux: F = -1/chi E,i
    ldouble kappa,kappaes,chi;
    chi=calc_kappa(pp[0],calc_PEQ_Tfromurho(pp[1],pp[0]),xxvec[1],xxvec[2],xxvec[3])
      +
      calc_kappaes(pp[0],calc_PEQ_Tfromurho(pp[1],pp[0]),xxvec[1],xxvec[2],xxvec[3]);
    
    ldouble xxvectemp[4]={xxvec[0],xxvec[1],xxvec[2],xxvec[3]};
    ldouble pptemp[NV],E1,E2,ggt[4][5],GGt[4][5];
    int anret,anretmin=0;

    //r dimension
    xxvectemp[1]=1.01*xxvecBL[1];
    xxvectemp[2]=1.0*xxvecBL[2];
    xxvectemp[3]=1.0*xxvecBL[3];
    calc_g_arb(xxvectemp,ggt,KERRCOORDS);
    calc_G_arb(xxvectemp,GGt,KERRCOORDS);

    anret=donut_analytical_solution(pptemp,xxvectemp,ggt,GGt);
    if(anret<0) anretmin=-1;
    E1=pptemp[6];

    xxvectemp[1]=.99*xxvecBL[1];
    xxvectemp[2]=1.0*xxvecBL[2];
    xxvectemp[3]=1.0*xxvecBL[3];
    calc_g_arb(xxvectemp,ggt,KERRCOORDS);
    calc_G_arb(xxvectemp,GGt,KERRCOORDS);

    anret=donut_analytical_solution(pptemp,xxvectemp,ggt,GGt);
    if(anret<0) anretmin=-1;
    E2=pptemp[6];

    Fx=(E2-E1)/(.02*xxvecBL[1]*sqrt(ggBL[1][1]))/chi/3.;

    //th dimension
    xxvectemp[1]=1.0*xxvecBL[1];
    xxvectemp[2]=1.01*xxvecBL[2];
    xxvectemp[3]=1.0*xxvecBL[3];
    calc_g_arb(xxvectemp,ggt,KERRCOORDS);
    calc_G_arb(xxvectemp,GGt,KERRCOORDS);

    anret=donut_analytical_solution(pptemp,xxvectemp,ggt,GGt);
    if(anret<0) anretmin=-1;
    E1=pptemp[6];

    xxvectemp[1]=1.0*xxvecBL[1];
    xxvectemp[2]=0.99*xxvecBL[2];
    xxvectemp[3]=1.0*xxvecBL[3];
    calc_g_arb(xxvectemp,ggt,KERRCOORDS);
    calc_G_arb(xxvectemp,GGt,KERRCOORDS);

    anret=donut_analytical_solution(pptemp,xxvectemp,ggt,GGt);
    if(anret<0) anretmin=-1;
    E2=pptemp[6];

    Fy=(E2-E1)/(.02*xxvecBL[2]*sqrt(ggBL[2][2]))/chi/3.;

    //ph dimension - symmetry
    Fz=0.;

    if(anretmin<0) //one of the points outside the donut
      Fx=Fy=Fz=0.;
    else
      {
	ldouble Fl=sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
	if(Fl>.99*E)
	  {
	    Fx=Fx/Fl*0.99*E;
	    Fy=Fy/Fl*0.99*E;
	    Fz=Fz/Fl*0.99*E;
	  }
      }

     //saving ff values to pp[]
     pp[7]=Fx;
     pp[8]=Fy;
     pp[9]=Fz;

#ifdef NOINITFLUX
     pp[7]=0.;
     pp[8]=0.;
     pp[9]=0.;
#endif
     
     //transforming from BL lab radiative primitives to code non-ortonormal primitives
     prad_ff2lab(pp,pp,&geomBL);
 
     //if(ix==NX && iy==NY-1){print_Nvector(pp,NV);getchar();}

#endif
      //transforming all primitives from BL to MYCOORDS
     trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,xxvecBL,ggBL,GGBL,gg,GG);
      }     
    
   
    pp[5]=calc_Sfromu(pp[0],pp[1]);
    
    //testing if primitives make sense
    //    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,&geom);

    return 0.;
  }
 else if(ix<0) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     //extrapolating along MYCOORDS

     //gc
     ldouble r=xxvec[1];

     //iix=0
     ldouble xxout[4]={0.,get_x(iix,0),get_x(iiy,1),get_x(iiz,2)};
     //     coco_N(xxout,xxout,MYCOORDS,BLCOORDS);
     ldouble r0=xxout[1];   
     
     //iix=1
     iix=1;
     xxout[1]=get_x(iix,0);
     //     coco_N(xxout,xxout,MYCOORDS,BLCOORDS);
     ldouble r1=xxout[1];   
     
     
     /*
     //copying YLCOORDS quantities
     for(iv=0;iv<NV;iv++)
       { 
	 //unchanged primitives
	 pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }

     if(MYCOORDS==SPHCOORDS)
       {
	 pp[0]=get_u(p,0,iix,iiy,iiz);
	 pp[1]=get_u(p,1,iix,iiy,iiz);
	 pp[2]=-10.;
       }
     else
       {
	 pp[0]=get_u(p,0,iix,iiy,iiz)*pow(r/r0,-1.5);
	 pp[1]=get_u(p,1,iix,iiy,iiz)*pow(r/r0,-2.5);
       }
     */

     //linear extrapolation
      for(iv=0;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,0,iiy,iiz)+(get_u(p,iv,1,iiy,iiz)-get_u(p,iv,0,iiy,iiz))*(r-r0)/(r1-r0);
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       }
 
     //atmosphere
      //set_radatmosphere(pp,xxvec,gg,GG,0);

#ifdef RADIATION
     //imposing inflowing velocity of the normal observer
     ldouble ucon[4];
     calc_normalobs_4vel(GG,ucon);
     pp[7]=ucon[1];
     pp[8]=ucon[2];
     pp[9]=ucon[3];

     if(MYCOORDS==KERRCOORDS)
       pp[7]=-100.;

     //pure copy
  
   iix=0;
   pp[6]=get_u(p,6,iix,iiy,iiz);

     //copying with scalings
  

     //pp[6]=get_u(p,6,iix,iiy,iiz)*pow(r/r0,-2.5);
     //pp[7]=get_u(p,7,iix,iiy,iiz)*pow(r/r0, 1.);
    
     //this works only for Kerr
     //if(pp[7]>0.) pp[7]=0.;
     //if(pp[2]>0.) pp[2]=0.;

#endif
     
     //testing if interpolated primitives make sense
     //     check_floors_hd(pp,VELPRIM,gg,GG);
     //end of floor section

     p2u(pp,uu,&geom);
     return 0;
   }

//reflections/outflow in theta 
if(iy<0.) //spin axis 
  {      
    
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	//v_theta
#ifndef PUREAXISOUTFLOW
	if(iv==3 || iv==8)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
#endif
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
     
    
    //testing if interpolated primitives make sense
    //    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,&geom);
    return 0;
  }
if(iy>=NY) //equatorial plane
  {
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
  	  
    for(iv=0;iv<NV;iv++)
      {
	if(iv==3 || iv==8)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	
      }

 

    //testing if interpolated primitives make sense
    //    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,&geom); 
    return 0; 
  }
   
//periodic in phi:
iiz=iz;
iiy=iy;
iix=ix;
if(iz<0) iiz=iz+NZ;
if(iz>NZ-1) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

