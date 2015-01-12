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
ldouble rSPH=geomSPH.xx;

gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
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

    #ifdef myVET
    int il;
     for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
     #endif
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

     if(pp[VX]>0.) pp[VX]=0.;
     if(pp[FX0]>0.) pp[FX0]=0.;

   
     p2u(pp,uu,&geom); 

     #ifdef myVET
     int il;
     for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
     #endif
     return 0;
   }

//reflections in theta 
if(iy<0. || iy>=NY) //zero gradient
  {      
    iiy=0;
    iiz=iz;
    iix=ix;

    for(iv=0;iv<NV;iv++)
      {
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
   
    p2u(pp,uu,&geom);
    
   #ifdef myVET
    int il;
    for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
   #endif

    return 0;
  }
   

if(iz<0.) //inflow
  {

   set_hdatmosphere(pp,xxvec,gg,GG,0);
#ifdef RADIATION
   set_radatmosphere(pp,xxvec,gg,GG,0);
#endif
   pp[5]=calc_Sfromu(pp[0],pp[1]);


   if(rSPH>RBEAM1 && rSPH<RBEAM2) //hot boundary
      {
	pp[EE0]*=1.e3;
	pp[FX0]=0.;
	pp[FZ0]=2./rSPH;
	pp[FY0]=0.;
      }
 

#ifdef myVET
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

	return 0;
  }

if(iz>=NZ) //outflow
  {

    iix=ix;
    iiy=iy;
    iiz=NZ-1;

    //copying MYCOORDS quantities
    for(iv=0;iv<NV;iv++)
      { 
	//unchanged primitives
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }
    
    if(pp[VZ]<0.) pp[VZ]=0.;
    if(pp[FZ0]<0.) pp[FZ0]=0.;
   
    p2u(pp,uu,&geom); 

    #ifdef myVET
    int il;
     for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
     #endif
    return 0.;
  }


return 0;

