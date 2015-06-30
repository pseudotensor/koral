ldouble uu[NV]; //this will hold conserved quantities
ldouble pp[NV]; //this will hold primitive quantities

/***********************************************/
struct geometry geom; 
fill_geometry(ix,iy,iz,&geom); //this gives to knowledge about the given cell you are working on

/***********************************************/
double x=geom.xx;
double y=geom.yy;
double z=geom.zz;

//by default, let's set up ambient medium
pp[RHO]=AMBRHO;
pp[UU]=calc_PEQ_ufromTrho(AMBTEMP,AMBRHO);
pp[VX]=0.;
pp[VY]=0.;
pp[VZ]=0.;

//if we are in plane 1
if(x>PL1X1 && x<PL1X2 && y>PL1Y)
  {
    pp[RHO]=PL1RHO;
    pp[UU]=calc_PEQ_ufromTrho(PL1TEMP,PL1RHO);
    pp[VX]=PL1VX;
    pp[VY]=0.;
    pp[VZ]=0.;
  }


//if we are in plane 1
if(x>PL2X1 && x<PL2X2 && y<PL2Y)
  {
    pp[RHO]=PL2RHO;
    pp[UU]=calc_PEQ_ufromTrho(PL2TEMP,PL2RHO);
    pp[VX]=PL2VX;
    pp[VY]=0.;
    pp[VZ]=0.;
  }


#ifdef MAGNFIELD

//if VECPOTGIVEN is uncommented then what is below is actually vector potential
pp[B1]=pp[B2]=pp[B3]=0.;
//if we are in plane 1
if(x>PL1X1 && x<PL1X2 && y>PL1Y)
  {
    pp[B3]=exp(-(x-.5*(PL1X1+PL1X2))*(x-.5*(PL1X1+PL1X2))/(PL1X2-PL1X1)/(PL1X2-PL1X1)/2./2.);
  }
//if we are in plane 2
if(x>PL2X1 && x<PL2X2 && y<PL2Y)
  {
    pp[B3]=exp(-(x-.5*(PL2X1+PL2X2))*(x-.5*(PL2X1+PL2X2))/(PL2X2-PL2X1)/(PL2X2-PL2X1)/2./2.);
  }

#endif


/***********************************************/
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
p2u(pp,uu,&geom); //converst primitves to conserved
/***********************************************/

//print_primitives(pp); getch();

int iv;
//save to memory
for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }
