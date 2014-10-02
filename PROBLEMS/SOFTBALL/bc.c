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

    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

   pp[RHO]=RHOAMB; 
   pp[UU]=UUAMB; 
   pp[VZ]=0.;
   pp[VY]=0.;
   pp[VX]=0.;
   pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

   //#ifdef MAGNFIELD
   //pp[B1]=pp[B2]=pp[B3]=0.;
   //#endif

 //check for inflow
    if(pp[VX]<0.)
      {
	pp[RHO]=RHOAMB; 
	pp[UU]=UUAMB; 
	pp[VZ]=0.;
	pp[VY]=0.;
	pp[VX]=0.;
	pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);
      }

#ifdef RADIATION
    if(pp[FX0]<0.)
      pp[FX0]=0.;
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
    
    for(iv=0;iv<=NV;iv++)
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

    p2u(pp,uu,&geom);


    return 0;  
  }


//outflow:                                                                                                   
if(BCtype==YBCLO)
  {
    iix=ix;
    iiy=0;
    iiz=iz;

    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }    

    //check for inflow
    if(pp[VY]>0.)
      {
	pp[VY]=0.;
	/*	
	pp[RHO]=RHOAMB; 
	pp[UU]=UUAMB; 
	pp[VZ]=0.;
	pp[VX]=0.;
	pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);*/
      }

#ifdef RADIATION
    if(pp[FY0]>0.)
      pp[FY0]=0.;
#endif

    p2u(pp,uu,&geom);
    return 0;
  }
//outflow:                                                                                                                                                                                                    
if(BCtype==YBCHI)
  {
    iix=ix;
    iiy=NY-1;
    iiz=iz;

    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }
    
    //check for inflow
    if(pp[VY]<0.)
      {pp[VY]=0.;
	/*
	pp[RHO]=RHOAMB; 
	pp[UU]=UUAMB; 
	pp[VZ]=0.;
	pp[VY]=0.;
	pp[VX]=0.;
	pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);*/
      }

#ifdef RADIATION
    if(pp[FY0]<0.)
      pp[FY0]=0.;
#endif

    p2u(pp,uu,&geom);
    return 0;
  }

//and that is all

return 0;
