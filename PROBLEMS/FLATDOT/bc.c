//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecMINK[4],xx,yy,zz;

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecMINK,MYCOORDS,MINKCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomMINK;
fill_geometry_arb(ix,iy,iz,&geomMINK,MINKCOORDS);

gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
//working in MINK
ldouble ggMINK[4][5],GGMINK[4][5];
calc_g_arb(xxvecMINK,ggMINK,MINKCOORDS);
calc_G_arb(xxvecMINK,GGMINK,MINKCOORDS);
ldouble eupMINK[4][4],eloMINK[4][4];
ldouble tupMINK[4][4],tloMINK[4][4];
calc_tetrades(ggMINK,tupMINK,tloMINK,MINKCOORDS);
calc_ZAMOes(ggMINK,eupMINK,eloMINK,MINKCOORDS);
/**********************/

int src=0;

if(ix>=NX )
  {
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<NV;iv++)
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      

    p2u(pp,uu,gg,GG);

    return 0.;
}

if(ix<0)
  {
    iix=-ix-1;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);      
	if(iv==2 || iv==7) 
	  pp[iv]*=-1.;
      }

    p2u(pp,uu,gg,GG);

    return 0.;
}

if(iy>=NY)
  {
    iix=ix;
    iiy=NY-1;
    iiz=iz;
    
    for(iv=0;iv<NV;iv++)
      pp[iv]=get_u(p,iv,iix,iiy,iiz);     

    p2u(pp,uu,gg,GG);

    return 0.;
}

if(iy<0)
  {
#ifndef HOTBOUNDARY
    iix=ix;
    iiy=0;
    iiz=iz;
    
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);  
      }    

    pp[8]=0.;


    p2u(pp,uu,gg,GG);

    return 0.;
#else
    if(ix<IXDOTMIN || ix>IXDOTMAX)
      {
	iix=ix;
	iiy=0;
	iiz=iz;
    
	for(iv=0;iv<NV;iv++)
	  {
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);  
	  }    

	pp[8]=0.;


	p2u(pp,uu,gg,GG);

	return 0.;
    
      }
    else
      {

	pp[0]=1.;
	pp[1]=1.;
	pp[2]=0.;
	pp[3]=0.;
	pp[4]=0.;
	pp[5]=calc_Sfromu(pp[0],pp[1]);
	pp[6]=LTEFACTOR*calc_LTE_Efromurho(pp[0],pp[1]);
	pp[7]=0.;
	pp[8]=0.;
	pp[9]=0.;

	if(NZ==1)
	  pp[6]*=10000.;
	else
	  pp[6]*=100000000.;
	pp[7]=FXDOT*pp[6];
	pp[8]=FYDOT*pp[6];
	pp[9]=FZDOT*pp[6];
	prad_ff2lab(pp,pp,&geom);

	p2u(pp,uu,gg,GG);

	return 0;
      }

#endif
}

  
//periodic in z
iix=ix;
iiy=iy;

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

