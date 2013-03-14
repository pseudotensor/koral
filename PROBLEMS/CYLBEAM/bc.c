//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecBL[4],xx,yy,zz;

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,MYCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);


struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,MYCOORDS);

gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
//working in BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,MYCOORDS);
calc_G_arb(xxvecBL,GGBL,MYCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,MYCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,MYCOORDS);
/**********************/


//radius
if(ix>=NX) //analytical solution at rout only
  {
   
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
	pp[2]=0.;
	pp[0]=1.;
	pp[1]=0.1;

#ifdef RADIATION
	pp[6]=calc_LTE_EfromT(1.e10);
	pp[6]=1.;
	pp[7]=pp[8]=pp[9]=0.;
	pp[7]=-.5*pp[6]; //isotropic

	//Keplerian gas
	ldouble rBL=xx;
	ldouble Om=1./pow(rBL,1.5)*1.;
	//Om=0.;

	ldouble ucon[4]={0.,0.,0.,Om};
	conv_vels(ucon,ucon,VEL3,VELPRIM,ggBL,GGBL);
		
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];	
	
	prad_ff2lab(pp,pp,&geomBL);

	pp[4]=0.;
	
	trans_pall_coco(pp, pp, MYCOORDS, MYCOORDS,xxvecBL,ggBL,GGBL,gg,GG);
	//	print_Nvector(pp,NV);getchar();
#endif
      }
      
    pp[5]=calc_Sfromu(pp[0],pp[1]);
    
    check_floors_hd(pp,VELPRIM,gg,GG);

    p2u(pp,uu,gg,GG);

    return 0.;
  }
 else if(ix<0) //cylindrical axis - reflection
   {
     iix=-ix-1;
     iiz=iz;
     iiy=iy;
     gdet_src=get_g(g,3,4,iix,iiy,iiz);  
     gdet_bc=get_g(g,3,4,ix,iy,iz);  
     for(iv=0;iv<NV;iv++)
       {
	 //radial component
	 //if(iv==9)
	 //pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	 //else
	 pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
   
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section
    
    //    printf("bc %d > %d",ix,iix);
    //    print_4vector(&get_u(p,6,iix,iiy,iiz));
    //print_Nvector(pp,NV); 

    p2u(pp,uu,gg,GG);

    return 0;
  }

   
//periodic in phi and z:
iiz=iz;
iiy=iy;
iix=ix;

if(iz<0) iiz=iz+NZ;
if(iz>NZ-1) iiz=iz-NZ;
if(iy<0) iiy=iy+NY;
if(iy>NY-1) iiy=iy-NY;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//testing if interpolated primitives make sense
check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

