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
fill_geometry_arb(ix,iy,iz,&geomBL,CYLCOORDS);

gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
//working in BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,CYLCOORDS);
calc_G_arb(xxvecBL,GGBL,CYLCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,CYLCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,CYLCOORDS);
/**********************/

if(iy>NY-1) 
  {      
    
    iiy=NY-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
   
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);
    return 0;
  }

if(ix>NX-1)
  {      
    
    iiy=iy;
    iiz=iz;
    iix=NX-1;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
   
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);
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
	//v_r
	if(iv==2 || iv==7)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
   
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);
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

    ldouble rBL=xxvecBL[1];
    ldouble rin=6.;
    if(rBL>rin) //hot boundary
      {

	pp[6]=calc_LTE_EfromT(1.e11)*(1.-sqrt(rin/rBL))/pow(rBL,3.);
	pp[7]=pp[8]=pp[9]=0.;
	pp[8]=0.*pp[6];

	//Keplerian gas
	ldouble Om=1./pow(rBL,1.5)*OMSCALE;
	
	ldouble ucon[4]={0.,0.,0.,Om};
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggBL,GGBL);
		
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];	
	
	prad_ff2lab(pp,pp,&geomBL);
	
	trans_pall_coco(pp, pp, CYLCOORDS, MYCOORDS,xxvecBL,ggBL,GGBL,gg,GG);
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
	    //v_z
	    if(iv==3|| iv==8)
	      pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	    else
	      pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }
      }
    

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section
    
    p2u(pp,uu,gg,GG); 
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

