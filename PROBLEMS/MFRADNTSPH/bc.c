//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
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

gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
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
if(ix>=NX) //analytical solution at rout only
  {

    int irf;    
    iix=NX-1;
    iiz=iz;
    iiy=iy;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);

	  for(irf=0;irf<NRF;irf++)
	    if(iv==FX(irf))
	      if(pp[iv]<0.) pp[iv]=0.;
      }  
    
    check_floors_hd(pp,VELPRIM,gg,GG);

    p2u(pp,uu,gg,GG);

    return 0.;
  }
 else if(ix<0) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     ldouble r=xxvecBL[1];
     ldouble xxout[4]={0.,get_x(iix,0),get_x(iiy,1),get_x(iiz,2)};
     coco_N(xxout,xxout,MYCOORDS,BLCOORDS);
     ldouble r0=xxout[1];      
     

     //copying MYCOORDS quantities
     for(iv=0;iv<NV;iv++)
       { 
	 //unchanged primitives
	 pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }

 
     //atmosphere
     //set_radatmosphere(pp,xxvec,gg,GG,0);

     //imposing inflowing velocity of the normal observer
     ldouble ucon[4];
     int irf;

     if(MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS)
       {
	 for(irf=0;irf<NRF;irf++)
	   if(pp[FX(irf)]>0.) pp[FX(irf)]=0.;
       }
     else
       {
	 pp[0]=get_u(p,0,iix,iiy,iiz)*pow(r/r0,-1.5);
	 pp[1]=get_u(p,1,iix,iiy,iiz)*pow(r/r0,-2.5);
	 calc_normalobs_4vel(GG,ucon);
	 for(irf=0;irf<NRF;irf++)
	   {
	     pp[FX(irf)]=ucon[1];
	     pp[EE(irf)]=get_u(p,EE(irf),iix,iiy,iiz)*pow(r/r0,-2.5);
	   }
       }
       
     //testing if interpolated primitives make sense
     check_floors_hd(pp,VELPRIM,gg,GG);
     //end of floor section
     p2u(pp,uu,gg,GG);
     return 0;
   }

//reflections in theta 
if(iy<0.) //spin axis
  {      
    int irf;
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);

	//reflection
	//v_theta
	if(iv==3)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);

	for(irf=0;irf<NRF;irf++)
	  if(iv==FY(irf))
	    //radial component
	    pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	
        }
   
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section
   

    p2u(pp,uu,gg,GG);

#ifdef SKIP_MULTIRADFLUID
    redistribute_radfluids(pp,uu,&geom);
    u2p_rad(uu,pp,&geom,&irf);
#ifdef MFCORRECTPHI
    mf_correct_in_azimuth(pp,uu,&geom,-1.);
#endif
#endif


    return 0;
  }
if(iy>=NY) //equatorial plane
  {
    int irf;
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  

    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
	if(iv==3)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	for(irf=0;irf<NRF;irf++)
	  if(iv==FY(irf))
	    //radial component
	    pp[iv]=-get_u(p,iv,iix,iiy,iiz);
      }

    ldouble rBL=xxvecBL[1];
    ldouble rin=6.;
    if(rBL>rin) //hot boundary
      {

	pp[6]=calc_LTE_EfromT(1.e11)*(1.-sqrt(rin/rBL))/pow(rBL,3.);
	pp[7]=pp[8]=pp[9]=0.;
	pp[8]=-0.*pp[6];

	//Keplerian gas
	ldouble Om=1./pow(rBL,1.5)*OMSCALE;
	
	ldouble ucon[4]={0.,0.,0.,Om};
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggBL,GGBL);
		
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];	
	
	int irf;
#ifdef MULTIRADFLUID

	for(irf=1;irf<NRF;irf++)
	  {
	    pp[FX(irf)]=pp[FY(irf)]=pp[FZ(irf)]=0.;
	    pp[EE(irf)]=EEFLOOR;
	  }

#endif	

	prad_ff2lab(pp,pp,&geomBL);
	
	trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,xxvecBL,ggBL,GGBL,gg,GG);
      }
    

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section
    
    p2u(pp,uu,gg,GG); 
#ifdef MULTIRADFLUID
    redistribute_radfluids(pp,uu,&geom);
    u2p_rad(uu,pp,&geom,&irf);
#ifdef MFCORRECTPHI
    mf_correct_in_azimuth(pp,uu,&geom,-1.);
#endif
#endif

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
check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

