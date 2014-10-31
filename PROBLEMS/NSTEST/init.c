ldouble uu[NV];
ldouble pp[NV];

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

/***********************************************/

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;
ldouble ph=geomBL.zz;

ldouble rfac=pow(r/RMIN,-6.);

rfac=1.;
pp[RHO]=RHO_AMB*rfac;// + RHO_BLOB*exp(-((y-5.)*(y-5.) + (z-20.)*(z-20.) + (x*x))/10.);
pp[UU]=U_AMB*rfac;//*rfac;

pp[VX]=0.;
pp[VY]=0.;
pp[VZ]=0.;

#ifdef MAGNFIELD
//vector potential
pp[B1]=get_u(pproblem1,B1,ix,iy,iz);
pp[B2]=get_u(pproblem1,B2,ix,iy,iz);
pp[B3]=get_u(pproblem1,B3,ix,iy,iz);

//test
/*
pp[B1]=2.*cos(th)/r/r/r/sqrt(geom.gg[1][1])*BETANORMFACTOR;
pp[B2]=sin(th)/r/r/r/sqrt(geom.gg[2][2])*BETANORMFACTOR;
pp[B3]=0.;
*/
//test
//pp[B3]=pp[B1]=pp[B2]=0.;
#endif

/*
//velocity in cartesian
ldouble ucart[3]={0.,-.1,0.};

ldouble usph[4];
ldouble vx,vy,vz;
vx=vz=0.;
vy=-.1;

ldouble vr,vth,vph;
ldouble cosph,sinth,costh,sinph;
sinth=sin(th);
costh=cos(th);
sinph=sin(ph);
cosph=cos(ph);


vr =-((-(cosph*sinth*vx) - sinph*sinth*vy - Power(cosph,2)*costh*vz - 
          costh*Power(sinph,2)*vz)/
      ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
vth = -((-(cosph*costh*vx) - costh*sinph*vy + Power(cosph,2)*sinth*vz + 
          Power(sinph,2)*sinth*vz)/
        ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
vph = -((sinph*vx - cosph*vy)/(Power(cosph,2) + Power(sinph,2)));

vth /= r;
vph /= r*sinth;

usph[1]=vr;
usph[2]=vth;
usph[3]=vph;

if (pp[RHO]/RHO_AMB > 2.)
  {   
    pp[VX]=vr;
    pp[VY]=vth;
    pp[VZ]=vph;   
  }
*/

/***********************************************/
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
p2u(pp,uu,&geom);
/***********************************************/

//print_primitives(pp); getch();

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }
