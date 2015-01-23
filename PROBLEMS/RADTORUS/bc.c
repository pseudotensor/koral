//returns problem specific BC
//ix, iy, iz hold the ghost cell indices, e.g., (-1,0,0)
//BCtype gives the type of the boundary

int iix,iiy,iiz,iv;

/***********************************************/
//structure of geometry
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
/***********************************************/

//atmosphere
if(BCtype==XBCHI) 
  {

    //outflow for radiation
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }


#ifdef RADIATION
    if(pp[FX0]<0.)
      pp[FX0]=0.;
#endif

    //fixed hydro part
   pp[RHO]=RHOAMB; 
   pp[UU]=UUAMB; 
   pp[VZ]=0.;
   pp[VY]=0.;
   pp[VX]=0.;
   pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

  
 #ifdef myVET
    int il;
     for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
     #endif
     

   p2u(pp,uu,&geom);
   return 0;  
  }


//outflow:
if(BCtype==XBCLO) 
  {
    iix=0;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    //check for inflow
    if(pp[VX]>0.)
      {
	/*pp[RHO]=RHOAMB; 
	pp[UU]=UUAMB; 
	pp[VZ]=0.;
	pp[VY]=0.;
	*/
	pp[VX]=0.;
	//pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);
      }

#ifdef RADIATION
    if(pp[FX0]>0.)
      pp[FX0]=0.;
	//set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
	
#endif

     #ifdef myVET
    int il;
     for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
     #endif

    p2u(pp,uu,&geom);


    return 0;  
  }

//reflections in theta 
if(BCtype==YBCLO) 
  {      
    
    iiy=-iy-1;
    iiz=iz;
    iix=ix;

    for(iv=0;iv<NV;iv++)
      {
	//v_theta
	if(iv==VY || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
   
    p2u(pp,uu,&geom);
    
    


#ifdef myVET

 /*
 for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
 */   //in cart coords
       double reflect_direction[3] = {1.,0.,0.};
       reflectI(reflect_direction, &Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][0], &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
      
#endif   

    return 0;
  }

//reflections in theta 
if(BCtype==YBCHI) 
  {      
    
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;

    for(iv=0;iv<NV;iv++)
      {
	//v_theta
	if(iv==VY || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
   
    p2u(pp,uu,&geom);
    
    /*
    int il;
    for(il=0;il<NUMANGLES;il++)
      {
	Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il]=Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][il];
      }
    */

#ifdef myVET
    
    //in cart coords
       double reflect_direction[3] = {1.,0.,0.};

    reflectI(reflect_direction, &Ibeam[iix+NGCX][iiy+NGCY][iiz+NGCZ][0], &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
    /*
    double Itemp[NUMANGLES];int il;
       for(il=0;il<NUMANGLES;il++)
	 Itemp[il]=Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][il];
       reflect_direction[0] = 0.;
       reflect_direction[1] = 1.;
      

       reflectI(reflect_direction, &Itemp[0], &Ibeam[ix+NGCX][iy+NGCY][iz+NGCZ][0]);
    */
#endif   

    return 0;
  }
 
//periodic in phi, used under OMP
iiz=iz;
iiy=iy;
iix=ix;
if(BCtype==ZBCLO) iiz=iz+NZ;
if(BCtype==ZBCHI) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }


//and that is all

return 0;
