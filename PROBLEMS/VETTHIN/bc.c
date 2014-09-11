//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecSPH[4],xx,yy,zz;

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecSPH,MYCOORDS,SPHCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);


struct geometry geomSPH;
fill_geometry_arb(ix,iy,iz,&geomSPH,SPHCOORDS);

gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
//working in SPH
ldouble ggSPH[4][5],GGSPH[4][5];
calc_g_arb(xxvecSPH,ggSPH,SPHCOORDS);
calc_G_arb(xxvecSPH,GGSPH,SPHCOORDS);
ldouble eupSPH[4][4],eloSPH[4][4];
ldouble tupSPH[4][4],tloSPH[4][4];
calc_tetrades(ggSPH,tupSPH,tloSPH,SPHCOORDS);
calc_ZAMOes(ggSPH,eupSPH,eloSPH,SPHCOORDS);
/**********************/


//radius
if(ix>=NX) //analytical solution at rout only
  {
    /*
    ldouble podpierd=-(GGSPH[0][0]-2.*ELL*GGSPH[0][3]+ELL*ELL*GGSPH[3][3]);
    ldouble ut=-1./sqrt(podpierd);
    ut/=UTPOT;
    ldouble uint,Vphi,rho,Vr;
    ldouble xx=get_x(ix,0);
    ldouble D,E,W,eps,uT,uphi,uPhi;

    
    if(1)
      {

	iix=NX-1;
	iiy=iy;
	iiz=iz;

	//ambient
	set_hdatmosphere(pp,xxvec,gg,GG,0);

#ifdef RADIATION       
	set_radatmosphere(pp,xxvec,gg,GG,0);
#endif
      }      
    */
    
    iix=NX-1;
    iiy=iy;
    iiz=iz;

    //copying MYCOORDS quantities
    for(iv=0;iv<NV;iv++)
      { 
	//unchanged primitives
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }
    
    if(pp[VX]<0.) pp[VX]=0.;
    if(pp[FX0]<0.) pp[FX0]=0.;
   
    p2u(pp,uu,&geom); 

    /*
    //calculate the intensities
    double RijM1[4][4];double M1[5];
    calc_Rij_M1(pp,&geom,RijM1);
    trans22_coco(geom.xxvec,RijM1,RijM1,MYCOORDS,SPHCOORDS);
    //input
    M1[0]=RijM1[0][0];
    M1[1]=RijM1[0][1];
    M1[2]=RijM1[0][2];
    M1[3]=RijM1[0][3];
    M1[4]=pp[EE0];
      
    ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
    */
    int il;
     for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
    return 0.;
  }
 else if(ix<0) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     ldouble r=xxvecSPH[1];
     ldouble xxout[4]={0.,get_x(iix,0),get_x(iiy,1),get_x(iiz,2)};
     coco_N(xxout,xxout,MYCOORDS,SPHCOORDS);
     ldouble r0=xxout[1];      
     

     //copying MYCOORDS quantities
     for(iv=0;iv<NV;iv++)
       { 
	 //unchanged primitives
	 pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }

   
     p2u(pp,uu,&geom); 

     int il;
     for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
     return 0;
   }

//reflections in theta 
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
	if(iv==VY || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
   
    p2u(pp,uu,&geom);
    
    //should reflect the intensities!!!!

/*
    int il;
    for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
	}
*/
#if(RADCLOSURE==VETCLOSURE)
    double reflect_direction[3] = {1.,0.,0.};
    reflectI(reflect_direction, &Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][0], &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
#endif   
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

    /*
    int il;
    for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
    */

#if(RADCLOSURE==VETCLOSURE)
    double reflect_direction[3] = {0.,0.,1.};

    reflectI(reflect_direction, &Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][0], &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
#endif


    //hot disk:
    ldouble rSPH=xxvecSPH[1];
    ldouble rin=6.;
#ifndef HOURGLASS
    if(rSPH>rin) //hot boundary
      {

	pp[EE0]=calc_LTE_EfromT(1.e11)*(1.-sqrt(rin/rSPH))/pow(rSPH,3.);
	pp[FX0]=pp[FZ0]=0.;
	pp[FY0]=-0.5*pp[EE0];

	//Keplerian gas
	ldouble Om=1./pow(rSPH,1.5)*OMSCALE;
	
	ldouble ucon[4]={0.,0.,0.,Om};
	
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggSPH,GGSPH);
		
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];	
	
	prad_ff2lab(pp,pp,&geomSPH);
	
	trans_pall_coco(pp, pp, SPHCOORDS, MYCOORDS,xxvecSPH,&geomSPH,&geom);


#if(RADCLOSURE==VETCLOSURE)
	//calculate the intensities
	double RijM1[4][4];double M1[5];
	calc_Rij_M1(pp,&geom,RijM1);
	//converting to RADCLOSURECOORDS
	trans22_coco(geom.xxvec, RijM1, RijM1, MYCOORDS, RADCLOSURECOORDS);
	//to ortonormal
	trans22_cc2on(RijM1,RijM1,geomSPH.tup);
	

	//input
	M1[0]=RijM1[0][0];
	M1[1]=RijM1[0][1];
	M1[2]=RijM1[0][2];
	M1[3]=RijM1[0][3];
	M1[4]=pp[EE0];
      
	ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
#endif
	


	

      }

#else
   if(rSPH>10. && rSPH<15.) //hot boundary
      {

	pp[EE0]=calc_LTE_EfromT(1.e11)*(1.-sqrt(rin/rSPH))/pow(rSPH,3.);
	pp[FZ0]=0.;
	pp[FY0]=-1.*pp[EE0];
	pp[FX0]=-1.*pp[EE0];
	  
	//Keplerian gas
	ldouble Om=1./pow(rSPH,1.5)*OMSCALE;	
	ldouble ucon[4]={0.,0.,0.,Om};	
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggSPH,GGSPH);		
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];	
	
	prad_ff2lab(pp,pp,&geomSPH);
	
	trans_pall_coco(pp, pp, SPHCOORDS, MYCOORDS,xxvecSPH,&geomSPH,&geom);

#if(RADCLOSURE==VETCLOSURE)
	//calculate the intensities
	double RijM1[4][4];double M1[5];
	calc_Rij_M1(pp,&geom,RijM1);
	//converting to RADCLOSURECOORDS
	trans22_coco(geom.xxvec, RijM1, RijM1, MYCOORDS, RADCLOSURECOORDS);
	//to ortonormal
	trans22_cc2on(RijM1,RijM1,geomSPH.tup);
	

	//input
	M1[0]=RijM1[0][0];
	M1[1]=RijM1[0][1];
	M1[2]=RijM1[0][2];
	M1[3]=RijM1[0][3];
	M1[4]=pp[EE0];
      
	ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
#endif
	


	

      }

#endif 
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

