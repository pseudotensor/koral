//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

int iix,iiy,iiz,iv;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
/***********************************************/

/***********************************************/
/***********************************************/
if(ix>=NX) //Sgr A* atmosphere
  {
    if(ifinit==0)
      {
	for(iv=0;iv<NV;iv++)
	  {
	    pp[iv]=get_u(pinit,iv,ix,iy,iz);	
	  }
      }
    else
      {
	//flat atmosphere
	//set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,2);
	//Sgr A* atmosphere
#ifdef DONUT
	int anret=donut_analytical_solution(pp,geomBL.xxvec,geomBL.gg,geomBL.GG);
	if(anret<0) //atmosphere
	  {
	    //ambient
	    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
	  }
	else
	  {
	    //transforming primitives from BL to MYCOORDS
	    trans_phd_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,geomBL.gg,geomBL.GG,geom.gg,geom.GG);
     
	  }
#else
	set_sgradisk(pp,geom.xxvec,&geom,&geomBL);
#endif

#ifdef TRACER
	pp[TRA]=get_u(p,TRA,NX-1,iy,iz);
#endif

	//testing if interpolated primitives make sense
	check_floors_hd(pp,VELPRIM,&geom);
	//end of floor section
      }

    p2u(pp,uu,&geom);
 
    return 0;
  }
 
/***********************************************/
/***********************************************/
if(ix<0) //outflow at inner edge / fixed atmosphere
   {
     //Sgr A* atmosphere
     set_sgradisk(pp,geom.xxvec,&geom,&geomBL);

     
    //extrapolating along MYCOORDS
     //gc radialcoordinate
     ldouble r=geomBL.xx;
     //iix=0
     struct geometry geomBL0;
     fill_geometry_arb(0,iy,iz,&geomBL0,KERRCOORDS);     
     ldouble r0=geomBL0.xx;
     
     //for velocities or for everything
     for(iv=0;iv<NV;iv++)
       {
	 //linear extrapolation
	 //pp[iv]=get_u(p,iv,0,iy,iz)+(get_u(p,iv,1,iy,iz)-get_u(p,iv,0,iy,iz))*(r-r0)/(r1-r0);

	 /*
	 if(iv==0)
	   {
	     pp[0]=get_u(p,iv,0,iy,iz)*r0/r;
	   }
	 else if(iv==1)
	   {
	     pp[1]=get_u(p,iv,0,iy,iz)*r0*r0/r/r;
	   }
	 else if(iv==4)
	   {
	     pp[4]=get_u(p,iv,0,iy,iz)*sqrt(r0*r0*r0/r/r/r);
	   }
	 else
	 */
	 //constant
	 pp[iv]=get_u(p,iv,0,iy,iz);
       }

     //no inflow 
     ldouble ucon[4]={0.,pp[2],pp[3],pp[4]};
     conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
     trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,KERRCOORDS);
     if(ucon[1]>0.)
       {
	 ucon[1]=0.;
	 trans2_coco(geom.xxvec,ucon,ucon,KERRCOORDS,MYCOORDS);
	 conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	 pp[2]=ucon[1];
	 pp[3]=ucon[2];
	 pp[4]=ucon[3];	   
       }
 
     //testing if interpolated primitives make sense
     check_floors_hd(pp,VELPRIM,&geom);
     //end of floor section
    
     p2u(pp,uu,&geom);
     return 0;
   }

/***********************************************/
/***********************************************/
//transmissive near the axis
if(iy<0.) 
  {      
    iiy=-iy-1;
    iiz=iz+NZ/2;
    if(iiz>NZ-1) iiz-=NZ;
    iix=ix;

    for(iv=0;iv<NV;iv++)
      {
	//reflect theta, conserve phi - handles rotation properly
	if(iv==3)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
    
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section

    p2u(pp,uu,&geom);
    return 0;
  }

if(iy>NY-1) 
  {      
    iiy=NY-1-(iy-NY);
    iiz=iz+NZ/2;
    if(iiz>NZ-1) iiz-=NZ;
    iix=ix;

    for(iv=0;iv<NV;iv++)
      {
	//reflect theta, conserve phi - handles rotation properly
	if(iv==3)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
    
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section

    p2u(pp,uu,&geom);
    return 0;
  }

/***********************************************/
/***********************************************/
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

//and that is all
 
return 0;

