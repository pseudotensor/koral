//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

int iix,iiy,iiz,iv,il;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/

/***********************************************/
/***********************************************/
if(ix>=NX) //outflow
  {
   
    for (iv=0; iv < NV; iv++)
      {
	pp[iv]=get_u(p,iv,NX-1,iy,iz);
      }
     
    p2u(pp,uu,&geom);
 
    for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[NX-1+NGCX][iy+NGCY][iz+NGCZ][il];
      }

    return 0;
  }
 
if(ix<0) //bulk inflow
  {
    ldouble E,Fx;
    if ((geom.yy < 0.25) || (geom.yy > 0.35))
      {
	E=1.e-5;
	Fx=0.;
      }
    else
      {
	E=1000.;
	Fx=5.;
      }

    
    for (iv=0; iv < NV; iv++)
      {
	pp[iv]=get_u(p,iv,0,iy,iz);
      }
    pp[EE0]=E;
    pp[FX0]=Fx;
    pp[FY0]=pp[FZ0]=0.;
    //stress energy tensor
    double RijM1[4][4];double M1[5];
    calc_Rij_M1(pp,&geom,RijM1);
    //input
    M1[0]=RijM1[0][0];
    M1[1]=RijM1[0][1];
    M1[2]=RijM1[0][2];
    M1[3]=RijM1[0][3];
    M1[4]=pp[EE0];
      
    ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
      
      /*
    if ((geom.yy < 0.25) || (geom.yy > 0.35))
      {
	for (iv=0; iv < NV; iv++)
	  {
	    pp[iv]=get_u(p,iv,0,iy,iz);
	  }
	
	for(il=0;il<NUMANGLES;il++)
	  {
	    Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[0+NGCX][iy+NGCY][iz+NGCZ][il];
	  }

      }
    else
      {
	for (iv=0; iv < NV; iv++)
	  {
	    pp[iv]=get_u(p,iv,0,iy,iz);
	  }
	pp[EE0]=1000.;
	pp[FX0]=5.;
        pp[FY0]=pp[FZ0]=0.;
	//stress energy tensor
	double RijM1[4][4];double M1[5];
	calc_Rij_M1(pp,&geom,RijM1);
	//input
	M1[0]=RijM1[0][0];
	M1[1]=RijM1[0][1];
	M1[2]=RijM1[0][2];
	M1[3]=RijM1[0][3];
	M1[4]=pp[EE0];
      
	ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
	
      }
    */



    p2u(pp,uu,&geom);
    return 0;
  }
 
if(iy<0) //beam
  {
    ldouble E,Fy;
    if ((geom.xx < 0.25) || (geom.xx > 0.35))
      {
	E=1.e-5;
	Fy=0.;
      }
    else
      {
	E=1000.;
	Fy=5.;
      }

    
    for (iv=0; iv < NV; iv++)
      {
	pp[iv]=get_u(p,iv,ix,0,iz);
      }
    pp[EE0]=E;
    pp[FY0]=Fy;
    pp[FX0]=pp[FZ0]=0.;
	
    //stress energy tensor
    double RijM1[4][4];double M1[5];
    calc_Rij_M1(pp,&geom,RijM1);
    //input
    M1[0]=RijM1[0][0];
    M1[1]=RijM1[0][1];
    M1[2]=RijM1[0][2];
    M1[3]=RijM1[0][3];
    M1[4]=pp[EE0];
      
    ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
	

    p2u(pp,uu,&geom);
    return 0;
  }
 
if(iy>=NY) 
  {
    
   for (iv=0; iv < NV; iv++)
      {
	pp[iv]=get_u(p,iv,ix,NY-1,iz);
      }
     
     p2u(pp,uu,&geom);

     for(il=0;il<NUMANGLES;il++)
       {
	 Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[ix+NGCX][NY-1+NGCY][iz+NGCZ][il];
       }

    return 0;
  }

/***********************************************/
/***********************************************/
//periodic in z:
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

//and that is all
 
return 0;

