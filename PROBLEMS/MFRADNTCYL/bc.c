//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecCYL[4],xx,yy,zz;

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecCYL,MYCOORDS,CYLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);


struct geometry geomCYL;
fill_geometry_arb(ix,iy,iz,&geomCYL,CYLCOORDS);

gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
//working in CYL
ldouble ggCYL[4][5],GGCYL[4][5];
calc_g_arb(xxvecCYL,ggCYL,CYLCOORDS);
calc_G_arb(xxvecCYL,GGCYL,CYLCOORDS);
ldouble eupCYL[4][4],eloCYL[4][4];
ldouble tupCYL[4][4],tloCYL[4][4];
calc_tetrades(ggCYL,tupCYL,tloCYL,CYLCOORDS);
calc_ZAMOes(ggCYL,eupCYL,eloCYL,CYLCOORDS);
/**********************/

if(iy>NY-1) 
  {      
	  int irf;    
    iiy=NY-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);

	  for(irf=0;irf<NRF;irf++)
	    if(iv==FY(irf))
	      if(pp[iv]<0.) pp[iv]=0.;
      }  
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);

#ifdef SKMULTIRADFLUID    
    redistribute_radfluids(pp,uu,&geom);
    u2p_rad(uu,pp,&geom,&irf);
#ifdef MFCORRECTPHI
    mf_correct_in_azimuth(pp,uu,&geom,-1.);
#endif
#endif
    return 0;
  }

if(ix>NX-1)
  {      
    	int irf;

    iiy=iy;
    iiz=iz;
    iix=NX-1;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
     
	for(irf=0;irf<NRF;irf++)
	  if(iv==FX(irf))
	    if(pp[iv]<0.)
	      pp[iv]=0.;
      }

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);

#ifdef SKMULTIRADFLUID   
    redistribute_radfluids(pp,uu,&geom);
    u2p_rad(uu,pp,&geom,&irf);
#ifdef MFCORRECTPHI
    mf_correct_in_azimuth(pp,uu,&geom,-1.);
#endif
#endif
    return 0;
  }
if(ix<0.) //spin axis
  {      
    
    iiy=iy;
    iiz=iz;
    iix=0;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);

	int irf;
	for(irf=0;irf<NRF;irf++)
	  if(iv==FX(irf))
	    //radial component
	    pp[iv]=-get_u(p,iv,iix,iiy,iiz);
      }
   
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);

#ifdef MULTIRADFLUID
    int irf;
    //print_Nvector(uu,NV);
    //    redistribute_radfluids_along_axes(pp,uu,&geom);
    redistribute_radfluids(pp,uu,&geom);
    u2p_rad(uu,pp,&geom,&irf);
    //    print_Nvector(uu,NV);
#ifdef MFCORRECTPHI
    mf_correct_in_azimuth(pp,uu,&geom,-1.);
#endif


    //    print_Nvector(uu,NV);getchar();
    //    u2p_rad(uu,pp,&geom,&irf);
#endif
    return 0;
  }
if(iy<0) //equatorial plane
  {
    
    iiy=0;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  

    for(iv=0;iv<NV;iv++)
      {
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    ldouble rCYL=xxvecCYL[1];
    ldouble rin=6.;
    if(rCYL>rin) //hot boundary
      {

	pp[6]=calc_LTE_EfromT(1.e11)*(1.-sqrt(rin/rCYL))/pow(rCYL,3.);
	pp[7]=pp[8]=pp[9]=0.;
	pp[8]=0.*pp[6];

	//Keplerian gas
	ldouble Om=1./pow(rCYL,1.5)*OMSCALE;
	
	ldouble ucon[4]={0.,0.,0.,Om};
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggCYL,GGCYL);
		
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
	prad_ff2lab(pp,pp,&geomCYL);
	
	trans_pall_coco(pp, pp, CYLCOORDS, MYCOORDS,xxvecCYL,ggCYL,GGCYL,gg,GG);
      }
    else
      {
	set_radatmosphere(pp,xxvec,gg,GG,0);
	iiy=-iy-1;
	iiz=iz;
	iix=ix;
	gdet_src=get_g(g,3,4,iix,iiy,iiz);  
	gdet_bc=get_g(g,3,4,ix,iy,iz);  

	for(iv=0;iv<NV;iv++)
	  {
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);

	    int irf;
	    for(irf=0;irf<NRF;irf++)
	      if(iv==FY(irf))
		//radial component
		pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	  }
      }

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section
    p2u(pp,uu,gg,GG); 
  
#ifdef MULTIRADFLUID
    int irf;
    //print_Nvector(uu,NV);
    //    print_Nvector(pp,NV);
    //redistribute_radfluids_along_axes(pp,uu,&geom);
    redistribute_radfluids(pp,uu,&geom);
    u2p_rad(uu,pp,&geom,&irf);
    //    if(uu[FZ(1)]!=0.){print_Nvector(uu,NV); getchar();}
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

