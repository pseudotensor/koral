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
 
	#ifdef myVET
    for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[NX-1+NGCX][iy+NGCY][iz+NGCZ][il];
      }
#endif

    return 0;
  }
 
if(ix<0) //bulk inflow
  {
    //if (geom.yy < 0.2)
    if (0)
      {
	for (iv=0; iv < NV; iv++)
	  {
	    pp[iv]=get_u(p,iv,0,iy,iz);
	  }
		#ifdef myVET
	for(il=0;il<NUMANGLES;il++)
	  {
	    Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[0+NGCX][iy+NGCY][iz+NGCZ][il];
	  }
	#endif

      }
        else
      {
	for (iv=0; iv < NV; iv++)
	  {
	    pp[iv]=get_u(p,iv,0,iy,iz);
	  }
	ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);

	pp[EE0]=1.0e4*calc_LTE_EfromT(temp);
//	pp[EE0]=0.01*pp[UU];

	pp[FX0]=1.;
        pp[FY0]=pp[FZ0]=0.;
	
	//stress energy tensor
	#ifdef myVET
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
      }

/*
	else
      {
	for (iv=0; iv < NV; iv++)
	  {
	    pp[iv]=get_u(p,iv,0,iy,iz);
	  }
	ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
	pp[EE0]=1.0e4*calc_LTE_EfromT(temp);


	pp[FX0]=1.;
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
    
    //if (geom.xx < 0.2)
    if (0)
      {
	for (iv=0; iv < NV; iv++)
	  {
	    pp[iv]=get_u(p,iv,ix,0,iz);
	  }
     

	#ifdef myVET
	for(il=0;il<NUMANGLES;il++)
	  {
	    Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[ix+NGCX][0+NGCY][iz+NGCZ][il];
	  }
	#endif
      }
    else
      {
	for (iv=0; iv < NV; iv++)
	  {
	    pp[iv]=get_u(p,iv,ix,0,iz);
	  }
	ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);

	pp[EE0]=1.0e4*calc_LTE_EfromT(temp);
//	pp[EE0]=0.01*pp[UU];

	pp[FY0]=1.;
        pp[FX0]=pp[FZ0]=0.;
	
	#ifdef myVET
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
#endif
      }

    p2u(pp,uu,&geom);
    return 0;
  }
 
if(iy>=NY) //bulk motion
  {
    
   for (iv=0; iv < NV; iv++)
      {
	pp[iv]=get_u(p,iv,ix,NY-1,iz);
      }
     
     p2u(pp,uu,&geom);

	#ifdef myVET
     for(il=0;il<NUMANGLES;il++)
       {
	 Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[ix+NGCX][NY-1+NGCY][iz+NGCZ][il];
       }
#endif
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

