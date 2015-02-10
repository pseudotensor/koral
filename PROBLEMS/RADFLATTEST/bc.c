//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

iix=ix;
    iiy=iy;
    iiz=iz;

/**********************/

//radius
if(ix>=NX) //outflow or 2nd beam
  {    
    #ifndef TWOBEAMS
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    for(iv=0;iv<NV;iv++)
      { 
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      } 

    if(pp[VX]<0.) pp[VX]=0.;
    if(pp[FX0]<0.) pp[FX0]=0.;
    
#ifdef myVET
    int il;
    for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
#endif
#else //beam from the right as well
    iix=NX-1;
 iiy=iy;
    iiz=iz;
     for(iv=0;iv<NVMHD;iv++)
       { 
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       } 

     pp[0]=RHOAMB;
     pp[1]=UUAMB;
     pp[2]=0.;
     pp[3]=0.;
     pp[4]=0.;
     pp[5]=calc_Sfromu(pp[0],pp[1]);

#ifdef RADIATION
     pp[6]=.1;
     pp[7]=-VELRAD;
     pp[8]=0.;
     pp[9]=0.;

     #if(RADCLOSURE==VETCLOSURE)
	//calculate the intensities
	double RijM1[4][4];double M1[5];
	calc_Rij_M1(pp,&geom,RijM1);
	
	//input
	M1[0]=RijM1[0][0];
	M1[1]=RijM1[0][1];
	M1[2]=RijM1[0][2];
	M1[3]=RijM1[0][3];
	M1[4]=pp[EE0];
      
	ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
#endif
#endif	 
    #endif

    p2u(pp,uu,&geom);
 
    return 0.;
  }

 else if(ix<0) //outflow near BH
   {
     iix=0;
    iiy=iy;
    iiz=iz;
     for(iv=0;iv<NVMHD;iv++)
       { 
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       } 

     pp[0]=LEFTRHO;
     pp[1]=UUAMB;
     if(pp[VX]<0.) pp[VX]=0.;
     pp[3]=0.;
     pp[4]=0.;
     pp[5]=calc_Sfromu(pp[0],pp[1]);

#ifdef RADIATION
     pp[6]=.1;
     pp[7]=VELRAD;
     pp[8]=0.;
     pp[9]=0.;

     #if(RADCLOSURE==VETCLOSURE)
	//calculate the intensities
	double RijM1[4][4];double M1[5];
	calc_Rij_M1(pp,&geom,RijM1);
	
	//input
	M1[0]=RijM1[0][0];
	M1[1]=RijM1[0][1];
	M1[2]=RijM1[0][2];
	M1[3]=RijM1[0][3];
	M1[4]=pp[EE0];
      
	ZERO_decomposeM1(ix,iy,iz,M1, &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
#endif
#endif	 
    
     p2u(pp,uu,&geom);

     //     print_primitives(pp);
     //getch();
 
     return 0;
   }

//periodic
while(iiz<0)    iiz+=NZ;
while(iiz>=NZ)    iiz-=NZ; 
//periodic
while(iiy<0)    iiy+=NY;
while(iiy>=NY)    iiy-=NY; 


for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }
for(iv=0;iv<NUMANGLES;iv++)
    Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][iv]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][iv];

  
return 0;
  
