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

if(ix>=NX || ix<0)
  {

    if(iy>.9*NY/2. && iy<1.1*NY/2.) src=1;
   
    ldouble uint,Vphi,rho,Vr;
    ldouble xx=get_x(ix,0);
    ldouble D,E,W,eps,uT,uphi,uPhi;

    //ambient
    pp[2]=pp[3]=pp[4]=0.;
    pp[0]=1.;
    pp[1]=0.1;

#ifdef RADIATION

    if(src==1 && 1)
      {	
	pp[6]=1.;
	pp[7]=pp[8]=pp[9]=0.;
	//pp[7]=-.5*pp[6]; //isotropic
	pp[7]=0.;

	//a'la Keplerian gas
	ldouble xx=xxvec[1];
	ldouble vy=1./pow(fabs(xx),1.5)*OMSCALE*xx;

	//	if(ix<0) vy*=-1.;

	ldouble ucon[4]={0.,0.,vy,0.};

	//	print_4vector(ucon);getchar();
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggMINK,GGMINK);
		
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];	
	
	prad_ff2lab(pp,pp,&geomMINK);

	pp[4]=pp[2]=pp[3]=0.;
	
	trans_pall_coco(pp, pp, MINKCOORDS, MYCOORDS,xxvecMINK,ggMINK,GGMINK,gg,GG);
      }
    else
      {
	pp[6]=0.001;
	pp[7]=pp[8]=pp[9]=0.;
     }
#endif
  
      
    pp[5]=calc_Sfromu(pp[0],pp[1]);
    
    check_floors_hd(pp,VELPRIM,gg,GG);

    p2u(pp,uu,gg,GG);

    return 0.;
}

if(iy>=NY || iy<0 )
  {

    if(ix>.9*NX/2. && ix<1.1*NX/2.) src=1;
   
    ldouble uint,Vphi,rho,Vr;
    ldouble xx=get_x(ix,0);
    ldouble D,E,W,eps,uT,uphi,uPhi;

    //ambient

    pp[2]=pp[3]=pp[4]=0.;
    pp[0]=1.;
    pp[1]=0.1;

#ifdef RADIATION

    if(src==1 && 1)
      {	
	pp[6]=1.;
	pp[7]=pp[8]=pp[9]=0.;
	//pp[7]=-.5*pp[6]; //isotropic
	pp[7]=0.;

	//a'la Keplerian gas
	ldouble yy=xxvec[2];

	ldouble vx=-1./pow(fabs(yy),1.5)*OMSCALE*yy;


	ldouble ucon[4]={0.,vx,0.,0.};
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggMINK,GGMINK);
		
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];	
	
	prad_ff2lab(pp,pp,&geomMINK);

	pp[4]=pp[2]=pp[3]=0.;
	
	trans_pall_coco(pp, pp, MINKCOORDS, MYCOORDS,xxvecMINK,ggMINK,GGMINK,gg,GG);
      }
    else
      {
	pp[6]=0.001;
	pp[7]=pp[8]=pp[9]=0.;
     }
#endif
  
      
    pp[5]=calc_Sfromu(pp[0],pp[1]);
    
    check_floors_hd(pp,VELPRIM,gg,GG);

    p2u(pp,uu,gg,GG);

    return 0.;
  }

   
//periodic in z:
iiz=iz;


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

