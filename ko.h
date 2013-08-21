#define NDIM 4

//physics
#define GGG0 (6.674e-8)
#define CCC0 (2.998e10)
#define gTILDA (1.)
#define cTILDA (1.)
#define GGG (GGG0/gTILDA)
#define CCC (CCC0/cTILDA)

//GM/C2 for Msun in cm
#define MSUNCM 147700.

//conversions
#define tempCGS2GU(x)    (x)
#define tempGU2CGS(x)    (x)
#define lenCGS2GU(x)    (x/MASSCM)
#define lenGU2CGS(x)    (x*MASSCM)
#define timeCGS2GU(x)   (x/MASSCM*CCC)
#define timeGU2CGS(x)   (x*MASSCM/CCC)
#define velCGS2GU(x)    (x/CCC)
#define velGU2CGS(x)    (x*CCC)
#define rhoCGS2GU(x)    (x*GGG/CCC/CCC*MASSCM*MASSCM)
#define rhoGU2CGS(x)    (x/GGG*CCC*CCC/MASSCM/MASSCM)
#define massCGS2GU(x)    (x*GGG/CCC/CCC/MASSCM)
#define massGU2CGS(x)    (x/GGG*CCC*CCC*MASSCM)
#define kappaCGS2GU(x)  (x/GGG*CCC*CCC/MASSCM)
#define kappaGU2CGS(x)  (x*GGG/CCC/CCC*MASSCM)
#define endenCGS2GU(x) (x*GGG*MASSCM*MASSCM/CCC/CCC/CCC/CCC)
#define endenGU2CGS(x) (x/GGG/MASSCM/MASSCM*CCC*CCC*CCC*CCC)
#define fluxCGS2GU(x) (x*GGG*MASSCM*MASSCM/CCC/CCC/CCC/CCC/CCC)
#define fluxGU2CGS(x) (x/GGG/MASSCM/MASSCM*CCC*CCC*CCC*CCC*CCC)

//constants
#define K_BOLTZ (1.3806488e-16 * GGG / CCC / CCC / CCC / CCC / MASSCM)
#define M_PROTON massCGS2GU(1.67262158e-24)
#define SIGMA_RAD (5.67e-5 * GGG / CCC / CCC / CCC / CCC / CCC * MASSCM * MASSCM)
#define MU_GAS 1.
#define Z_RATIO (1.0)
#define Pi (3.141592654)     
#define KAPPA_ES_COEFF (kappaCGS2GU(0.4))

//other stuff
#include "problem.h"

#include "mdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>  
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multiroots.h>



#ifdef PR_DEFS
#include PR_DEFS
#endif

#define ldouble double
#define FTYPE ldouble
#define gSIZE 20 //size of metric arrays = 16 + 1 (gdet) + 3 (dlgdet)

//global variables
ldouble *u,*x,*xb,*du,*ut1,*ut2,*ut3,*ut4,*ut0,*u_bak,*p_bak,*u_step1,*u_step2,*u_step3,*u_step4,*ahdx,*ahdy,*ahdz,*aradx,*arady,*aradz,
  *ahdxl,*ahdyl,*ahdzl,*aradxl,*aradyl,*aradzl,  *ahdxr,*ahdyr,*ahdzr,*aradxr,*aradyr,*aradzr,*p,*pinit,*pproblem1,*pproblem2,*emf,
  *ptm1,*ptm2,*pt0,*px,*py,*pz,*s,*g,*gbx,*gby,*gbz,*Gbx,*Gby,*Gbz,
  *pbLx,*pbRx,*pbLy,*pbRy,*pbLz,*pbRz,*sbLx,*sbRx,*sbLy,*sbRy,*sbLz,*sbRz,*ubLx,*ubRx,*ubLy,*ubRy,*ubLz,*ubRz,
  *flbx,*flby,*flbz,*flLx,*flRx,*flLy,*flRy,*flLz,*flRz,*gKr,*gKrbx,*gKrby,*gKrbz,*G,
  *emuup,*emulo,*emuupbx,*emulobx,*emuupby,*emuloby,*emuupbz,*emulobz,
  *emuup2,*emulo2,*emuupbx2,*emulobx2,*emuupby2,*emuloby2,*emuupbz2,*emulobz2,
  *tmuup,*tmulo,*tmuupbx,*tmulobx,*tmuupby,*tmuloby,*tmuupbz,*tmulobz,
  *tmuup2,*tmulo2,*tmuupbx2,*tmulobx2,*tmuupby2,*tmuloby2,*tmuupbz2,*tmulobz2;
int *cellflag,**loop_0,**loop_1,**loop_2,**loop_3,**loop_4,**loop_02,Nloop_0,Nloop_1,Nloop_2,Nloop_02,Nloop_3,Nloop_4;

ldouble Kr_tmp[4][4][4],g_tmp[4][4];

ldouble inputarg[10];
int **gcidx;

ldouble max_ws[3],max_dt,ttm1,ttm2;
ldouble min_dx,min_dy,min_dz;
ldouble dt,tstepdenmax;
FILE *fout1,*fout_scalars,*fout_radprofiles;
int nfout1;

//some macros
#define my_max(x,y) (x>y?x:y)

//geometry structure
struct geometry
{
  int ix,iy,iz;
  ldouble xxvec[4];
  ldouble xx,yy,zz;
  ldouble gg[4][5];
  ldouble GG[4][5];
  ldouble gdet;
  ldouble alpha;
  ldouble tlo[4][4];
  ldouble tup[4][4];
  ldouble elo[4][4];
  ldouble eup[4][4];
};

//main.c
int solve_the_problem(ldouble);

gsl_odeiv2_step **odeiv2_step_1;
gsl_odeiv2_step **odeiv2_step_2;
struct evolve_fluxes_1_param
{
  ldouble ix,iy,iz,t,dt;
};
struct evolve_fluxes_2_param
{
  ldouble ix,iy,iz,t,dt;
  struct rad_parameters *rp;
};

//postproc.c
int calc_scalars(ldouble*,ldouble);
ldouble calc_mdotEdd();
ldouble calc_lumEdd();
int calc_radialprofiles(ldouble profiles[][NX]);
ldouble calc_totalmass();
ldouble calc_mdot(ldouble radius,int);
ldouble calc_lum(ldouble radius);
ldouble calc_photloc(int ix);

//misc.c
int my_clock_gettime(void* tsptr);
int print_NVvector(ldouble v[4]);
int print_tensor(ldouble T[][4]);
int print_metric(ldouble T[][5]);
int print_4vector(ldouble v[4]);
int print_Nvector(ldouble v[4],int);
ldouble calc_eigen_4x4(ldouble g[][4], ldouble *ev);
ldouble my_atan2(ldouble y, ldouble x);
void shuffle_loop(int **array, size_t n);
ldouble step_function(ldouble x,ldouble k);
int calc_stationary1d_solution()  ;
int initialize_arrays();
int free_arrays();
ldouble my_min(ldouble a, ldouble b);
//ldouble my_max(ldouble a, ldouble b);
ldouble my_sign(ldouble);
int find_eigenvalues3(ldouble[3][3],ldouble*);
int find_eigenvalues(ldouble*,int);
ldouble find_max_eigenvalue(ldouble *data,int N);
int find_max_eigenvalue_lr(ldouble *data,int N,ldouble*,ldouble*);
int my_err(char *);
int inverse_44matrix(ldouble a[][4], ldouble ia[][4]);
int convert_out2gif_1d(char *fname,char*,int niter,ldouble t);
int convert_out2gif_2d(char *fname,char*,int niter,ldouble t);
int getch(void);
int dosthelse(void);

//finite.c

int
solve_implicit_metric(int ix,int iy,int iz,ldouble dt,ldouble *ubase);
int cell_fixup_rad();
int cell_fixup_hd();
ldouble f_calc_fluxes_at_faces(int ix,int iy,int iz);
int f_timeder (ldouble t, ldouble dt,ldouble *);
int set_grid(ldouble*, ldouble*, ldouble*,ldouble*);
int print_grid(ldouble,ldouble,ldouble);
ldouble fd_flux_limiter(ldouble r);
ldouble minmod_fd_flux_limiter(ldouble ,ldouble,ldouble);
ldouble f_der_kurganovtadmor(int ix,int iy, int yz,ldouble*);
ldouble f_der_hlle_obsolete(int ix,int iy, int yz,ldouble*);
ldouble f_der_muscl(int ix,int iy, int yz,ldouble*);
int copy_u(ldouble,ldouble*,ldouble*);
int add_u(ldouble f1, ldouble* u1, ldouble f2, ldouble *u2, ldouble *u3);
ldouble f_timeder_source_term(ldouble t, const ldouble y[], ldouble f[],  void *params);

int set_bc(ldouble,int);
int if_indomain(int ix,int iy,int iz);
int if_outsidegc(int ix,int iy,int iz);
int if_outsidewave(int ix,int iy,int iz);
ldouble get_size_x(int ic, int idim);


//size of 3d arrays
#define SX (NX+2*NG)
#define NGCX NG
#define iX(ix) (ix)
#if(NY>1)
#define NGCY NG
#define SY (NY+2*NG)
#define iY(iy) (iy)
#else
#define NGCY 0
#define SY 1
#define iY(iy) (0)
#endif
#if(NZ>1)
#define NGCZ NG
#define SZ (NZ+2*NG)
#define iZ(iz) (iz)
#else
#define NGCZ 0
#define SZ 1
#define iZ(iz) (0)
#endif

//memory wrappers
//ldouble get_x(int,int);
int get_xx(int ix,int iy,int iz,ldouble *xx);
#define get_x(ic,idim) (idim==0 ? x[ic+NG] : (idim==1 ? x[ic+NG + NX+2*NG] : (idim==2 ? x[ic+NG + NX+2*NG + NY+2*NG ] : 0.)))
//ldouble get_xb(int,int);
#define get_xb(ic,idim) (idim==0 ? xb[ic+NG] : (idim==1 ? xb[ic+NG + NX+2*NG + 1] : (idim==2 ? xb[ic+NG + NX+2*NG +1 + NY+2*NG +1 ] : 0.)))
int set_x(int,int,ldouble);
int set_xb(int,int,ldouble);
//end of coordinates

#define get_emf(iv,ix,iy,iz) (emf[iv-1 + (ix)*3 + (iy)*(NX+1)*3 + (iz)*(NY+1)*(NX+1)*3])
#define set_emf(iv,ix,iy,iz,val) emf[iv-1 + (ix)*3 + (iy)*(NX+1)*3 + (iz)*(NY+1)*(NX+1)*3]=val
//end of emf

//ldouble get_u(ldouble*,int,int,int,int);
#define get_cflag(iflag,ix,iy,iz) (cellflag[iflag + (iX(ix)+NG)*NFLAGS + (iY(iy)+NGCY)*(SX)*NFLAGS + (iZ(iz)+NGCZ)*(SY)*(SX)*NFLAGS])
#define set_cflag(iflag,ix,iy,iz,val) cellflag[iflag + (iX(ix)+NG)*NFLAGS + (iY(iy)+NGCY)*(SX)*NFLAGS + (iZ(iz)+NGCZ)*(SY)*(SX)*NFLAGS]=val
#define get_u(uarr,iv,ix,iy,iz) (uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV])
//int set_u(ldouble*,int,int,int,int,ldouble);
#define set_u(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV]=val
#define get_u_scalar(uarr,ix,iy,iz) (uarr[(iX(ix)+NGCX) + (iY(iy)+NGCY)*(SX) + (iZ(iz)+NGCZ)*(SY)*(SX)])
//ldouble get_u_scalar(ldouble*,int,int,int);
//int set_u_scalar(ldouble*,int,int,int,ldouble);
#define set_u_scalar(uarr,ix,iy,iz,val) uarr[iX(ix)+NGCX + (iY(iy)+NGCY)*(SX) + (iZ(iz)+NGCZ)*(SY)*(SX)] = val
int set_ub(ldouble* uarr,int iv,int ix,int iy,int iz,ldouble value,int idim);
//#define set_ub(uarr,iv,ix,iy,iz,idim,val) (idim==0 ? uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX+1)*NV + (iZ(iz)+NGCZ)*(SY)*(SX+1)*NV]=val : (idim==1 ? uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY+1)*(SX)*NV]=val : (idim==2 ? uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV]=val : 0.)))
#define set_ubx(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX+1)*NV + (iZ(iz)+NGCZ)*(SY)*(SX+1)*NV]=val
#define set_uby(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY+1)*(SX)*NV]=val
#define set_ubz(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV]=val
//ldouble get_ub(ldouble* uarr,int iv,int ix,int iy,int iz,int idim);
#define get_ub(uarr,iv,ix,iy,iz,idim) (idim==0 ? uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX+1)*NV + (iZ(iz)+NGCZ)*(SY)*(SX+1)*NV] : (idim==1 ? uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY+1)*(SX)*NV] : (idim==2 ? uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV] : 0.)))
//ldouble get_g(ldouble* uarr, int i,int j, int iX(ix), int iy, int iz);
#define get_g(uarr,i,j,ix,iy,iz) uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX)*gSIZE + (iZ(iz)+NGCZ)*(SY)*(SX)*gSIZE]
int set_g(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value);
int set_T(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value);

#define get_T(uarr,i,j,ix,iy,iz) uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZ(iz)+NGCZ)*(SY)*(SX)*16]
#define get_Tb(uarr,i,j,ix,iy,iz,idim) (idim==0 ? uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX+1)*16 + (iZ(iz)+NGCZ)*(SY)*(SX+1)*16] : (idim==1 ? uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZ(iz)+NGCZ)*(SY+1)*(SX)*16] : (idim==2 ? uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZ(iz)+NGCZ)*(SY)*(SX)*16] : 0.)))
int set_Tb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim);
int set_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim);
//ldouble get_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,int idim);
#define get_gb(uarr,i,j,ix,iy,iz,idim) (idim==0 ? uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX+1)*gSIZE + (iZ(iz)+NGCZ)*(SY)*(SX+1)*gSIZE] : (idim==1 ? uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX)*gSIZE + (iZ(iz)+NGCZ)*(SY+1)*(SX)*gSIZE] : (idim==2 ? uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX)*gSIZE + (iZ(iz)+NGCZ)*(SY)*(SX)*gSIZE] : 0.)))
#define get_gKr(i,j,k,ix,iy,iz) gKr[i*4*4+j*4+k + (iX(ix)+NGCX)*64 + (iY(iy)+NGCY)*(SX)*64 + (iZ(iz)+NGCZ)*(SY)*(SX)*64]
#define get_gKrb(i,j,k,ix,iy,iz,idim) (idim==0 ? gKrbx[i*4*4+j*4+k + (iX(ix)+NGCX)*64 + (iY(iy)+NGCY)*(SX+1)*64 + (iZ(iz)+NGCZ)*(SY)*(SX+1)*64] : (idim==1 ? gKrby[i*4*4+j*4+k + (iX(ix)+NGCX)*64 + (iY(iy)+NGCY)*(SX)*64 + (iZ(iz)+NGCZ)*(SY+1)*(SX+1)*64] : (idim==2 ? gKrbz[i*4*4+j*4+k + (iX(ix)+NGCX)*64 + (iY(iy)+NGCY)*(SX)*64 + (iZ(iz)+NGCZ)*(SY)*(SX+1)*64] : 0.)))
#define set_gKr(i,j,k,ix,iy,iz,val) gKr[i*4*4+j*4+k + (iX(ix)+NGCX)*64 + (iY(iy)+NGCY)*(SX)*64 + (iZ(iz)+NGCZ)*(SY)*(SX)*64]=val
int set_Krb(int i,int j,int k,int ix,int iy,int iz,ldouble value,int idim);

//other wrappers
#define delta(i,j) (i==j ? 1 : 0)

//fileop.c
int fprint_restartfile(ldouble t, char* folder);
int fprint_simplecart(ldouble t, char* folder);
int fread_restartfile(int,ldouble*);
int fprint_gridfile(char* folder);
int fread_dumpfile(int,ldouble*);
int fread_dumpfile(int,ldouble*);
int fprint_openfiles(char *);
int fprint_closefiles();
int fprint_profiles(ldouble,ldouble*,int,int,char*);
int print_profiles();

//physics.c
struct rad_parameters
{
  ldouble f_Edd[3][3];
  ldouble F0[3];
  ldouble lambda;
  ldouble chi; 
  ldouble kappa;
  ldouble B;
  ldouble tau[3];
  ldouble tautot[3];
  ldouble T;
  ldouble dE0dx;
  ldouble R;
  ldouble f;
  ldouble x,y,z;
};
int f_metric_source_term_arb(ldouble *pp,void *ggg,ldouble *ss);
int calc_hd_shearviscosity(ldouble *pp,void* ggg,ldouble shear[][4],ldouble *nuret,ldouble *vdiff2ret);
int calc_shear_comoving(int ix,int iy,int iz,ldouble S[][4],int hdorrad);
int calc_shear_lab(int ix,int iy,int iz,ldouble S[][4],int hdorrad);

int calc_wavespeeds_lr_pure(ldouble *pp,void*,ldouble *aaa);
int calc_visc_Tij(ldouble *pp, void* ggg, ldouble T[][4]);
int calc_Tij( ldouble *p, void*, ldouble T[][4]);
ldouble max_eigen_Jac(ldouble *,ldouble*,int,void*);
int calc_wavespeeds(int,int,int,ldouble*,ldouble*,ldouble*,ldouble*,ldouble*,ldouble*);
int calc_wavespeeds_lr(int,int,int,ldouble*);
int calc_wavespeeds_lr_new(int,int,int,ldouble*);
int calc_wavespeeds_lr_faces( int,int,int,int, ldouble*,ldouble*);
int max_eigen_lr_Jac(ldouble *,ldouble*,int,void*,ldouble*,ldouble*);
int calc_Jac_num(ldouble *xx,ldouble *ujac, int idim,void *parameters,ldouble *fd_jac);
int set_initial_profile();
ldouble f_flux_prime(ldouble *uu, int,int,int,int,ldouble *ff);
ldouble f_diffusion_prime(ldouble *uu, ldouble *du, int iv,int,void*);
ldouble f_grav_potential(ldouble,ldouble,ldouble);
ldouble f_der_grav_potential(ldouble,ldouble,ldouble,int);
int f_source_term(int,int,int,ldouble *);
int f_fourforce_source_term(int,int,int,ldouble *);
int f_implicit_rhou(int ix, int iy, int iz, ldouble *rho, ldouble *uint, ldouble dt);
int f_other_source_term(int ix, int iy, int iz,ldouble *ss);
int f_metric_source_term(int ix, int iy, int iz,ldouble *ss);
int f_metric_source_term_face(int ix, int iy, int iz,int idim,int ifleft,ldouble *ss);
int f_source_term_face(int ix, int iy, int iz,int idim,int ifleft,ldouble *ss);
int initialize_problem();
int calc_maxgravspeed(int ix, int iy, int iz, ldouble dt);
int calc_maxwavespeed_obsolete(int ix, int iy, int iz,void*);
int calc_maxwavespeed_osbolete_cs(int ix, int iy, int iz,void*);
int rad_fld_factors(int ix,int iy, int iz,ldouble *,void *rad_param,int);
int rad_fld_factors_arb(int ix,int iy, int iz,ldouble*,void *rad_param);

int radx_flux_implicit(ldouble *uu,ldouble dtt);
int LTE_implicit(ldouble *uu,ldouble dtt);
ldouble max_eigen_radx(ldouble nx,ldouble ny,ldouble nz,int idim);
int max_eigen_lr_radx(ldouble nx,ldouble ny,ldouble nz,int idim,ldouble *al, ldouble *ar);
int gsl_poly_complex_solve_quartic (double a, double b, double c, double d,
                                gsl_complex * z0, gsl_complex * z1,
                                gsl_complex * z2, gsl_complex * z3);
ldouble f_der_hlle         (int ix,int iy,int iz, ldouble *fd_der);
ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z);
ldouble calc_kappaes(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z);
ldouble calc_ufromS(ldouble S,ldouble rho);
ldouble calc_Sfromu(ldouble S,ldouble u);
int
avg2point(ldouble *um2,ldouble *um1,ldouble *u0,ldouble *up1,ldouble *up2,ldouble*,ldouble*,ldouble dxm2,ldouble dxm1,ldouble dx0,ldouble dxp1,ldouble dxp2);

//problem.c
ldouble calc_xb(int i,int idim);
int calc_bc(int,int,int,ldouble,ldouble*, ldouble*,int,int);
int pr_tophat_inside(ldouble x,ldouble y,ldouble z);
int my_finger(ldouble);
int analytical_solution(ldouble t,int ix,int iy,int iz,ldouble *uu,ldouble *pp,ldouble *vv);

//metric.c

int fill_geometry(int ix,int iy,int iz,void *geom);
int fill_geometry_face(int ix,int iy,int iz,int,void *geom);
int fill_geometry_arb(int ix,int iy,int iz,void *geom,int COORDS);
int calc_metric();
int calc_tetrades(ldouble g[][5], ldouble tmuup[][4], ldouble tmulo[][4],int);
int calc_ZAMOes(ldouble g[][5], ldouble emuup[][4], ldouble emulo[][4],int);
int dxdx_KS2BL(ldouble *xx, ldouble dxdx[][4]);
int dxdx_BL2KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2MKS1(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MKS12KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MCYL12CYL(ldouble *xx, ldouble dxdx[][4]);
int dxdx_CYL2MCYL1(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MKER12KER(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KER2MKER1(ldouble *xx, ldouble dxdx[][4]);

int calc_g(ldouble*,ldouble[][5]);
int calc_G(ldouble*,ldouble[][5]);
int calc_Krzysie(ldouble*,ldouble[][4][4]);
int calc_Krzysie_at_center(int ix,int iy,int iz, ldouble Krzys[][4][4]);
int calc_g_arb(ldouble*,ldouble[][5],int);
int calc_G_arb(ldouble*,ldouble[][5],int);
int calc_Krzysie_arb(ldouble*,ldouble[][4][4],int);
int print_Krzysie(ldouble g[][4][4]);
int print_g(ldouble [][5]);
ldouble calc_gdet(ldouble *xx);
ldouble calc_dlgdet(ldouble *xx, int idim);
ldouble calc_gdet_arb(ldouble *xx,int);
ldouble calc_dlgdet_arb(ldouble *xx, int idim,int);
int coco_N(ldouble *x1, ldouble *x2,int CO1, int CO2);
int coco_BL2KS(ldouble *xBL, ldouble *xKS);
int coco_KS2BL(ldouble *xBL, ldouble *xKS);

//relele.c
int calc_normalobs_4vel(ldouble GG[][5], ldouble *ncon);
int set_hdatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype);
int calc_photonrad_4vel(ldouble gg[][5],ldouble GG[][5], ldouble *ucon);
int conv_velsinprims(ldouble *pp,int which1, int which2,ldouble gg[][5],ldouble GG[][5]);
#define dot(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3])
#define dot3(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define kron(i,j) (i == j ? 1. : 0.)
ldouble pick_gdet(int ix,int iy,int iz);
int pick_g(int ix,int iy,int iz,ldouble gg[][5]);
int pick_G(int ix,int iy,int iz,ldouble gg[][5]);
int pick_gb(int ix,int iy,int iz,int,ldouble gg[][5]);
int pick_Gb(int ix,int iy,int iz,int,ldouble gg[][5]);
int pick_T(ldouble *arr,int ix,int iy,int iz,ldouble T[][4]);
int pick_Tb(ldouble *arr,int ix,int iy,int iz,int,ldouble T[][4]);
int p2u_Sonly(ldouble *p, ldouble *u,ldouble[][5]);
int print_p(ldouble *p);
int print_u(ldouble *p);
int convert_uold2urel(ldouble *x,ldouble *u);
int calc_sourceterms(int,int,int);
ldouble r_horizon_BL(ldouble a);
ldouble r_mbound_BL(ldouble a);
ldouble r_photon_BL(ldouble a);
int update_entropy(int ix,int iy,int iz,int u2pflag);
int conv_vels(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);

//u2p.c
int calc_primitives_local(int ix,int iy,int iz,ldouble *pp);
int calc_primitives(int,int,int);
int check_floors_hd(ldouble *uu, int,void*);
int u2p(ldouble *uu, ldouble *pp, void*,int*,int*);
int u2p_hot(ldouble*,ldouble*,void*);
int u2p_entropy(ldouble*,ldouble*,void*);
int u2p_cold(ldouble*,ldouble*,void*);
int u2p_hotmax(ldouble*,ldouble*,void*);
int u2p_solver(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose);

//u2prad.c
int check_floors_rad(ldouble *uu, int,void*);
int u2p_rad(ldouble *uu, ldouble *pp, void*,int*);
int u2p_rad_onff(ldouble *uu, ldouble *pp, void*,int*);
int u2p_rad_urf(ldouble *uu, ldouble *pp,void* ggg, int *corrected);

//p2u.c
int calc_conserved(int ix,int iy,int iz);
int p2u(ldouble *p, ldouble *u,void*);
int pff2u(ldouble *p, ldouble *u,ldouble[][5],ldouble[][4],ldouble[][4]);
int p2u_rad(ldouble *p,ldouble *u,void*);

//frames.c
int trans_pmhd_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble *xxvec, void* ggg1,void* ggg2);
int boost2_lab2rf(ldouble A1[4],ldouble A2[4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5]);
int boost22_rf2lab(ldouble T1[][4],ldouble T2[][4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5]);
int boost22_ff2lab(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5]);
int boost22_lab2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5]);
int boost22_lab2rf(ldouble T1[][4],ldouble T2[][4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5]);
int boost2_lab2ff(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5]);
int boost2_ff2lab(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5]);
int boost22_ff2zamo(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4]);
int boost22_zamo2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4]);
int boost2_zamo2ff(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4]);
int boost2_ff2zamo(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4]);
int trans22_lab2zamo(ldouble T1[][4],ldouble T2[][4],ldouble tup[][4]);
int trans22_zamo2lab(ldouble T1[][4],ldouble T2[][4],ldouble tup[][4]);
int trans2_lab2zamo(ldouble *u1,ldouble *u2,ldouble e[][4]);
int trans2_zamo2lab(ldouble *u1,ldouble *u2,ldouble e[][4]);
int trans22_cc2on(ldouble T1[][4],ldouble T2[][4],ldouble tup[][4]);
int trans22_on2cc(ldouble T1[][4],ldouble T2[][4],ldouble tlo[][4]);
int trans2_cc2on(ldouble *u1,ldouble *u2,ldouble t[][4]);
int trans2_on2cc(ldouble *u1,ldouble *u2,ldouble t[][4]);
int indices_2221(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);
int indices_2122(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5]);
int indices_2211(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5]);
int indices_1122(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5]);
int indices_21(ldouble A1[4],ldouble A2[4],ldouble gg[][5]);
int indices_12(ldouble A1[4],ldouble A2[4],ldouble GG[][5]);
int prad_ff2lab(ldouble *pp1, ldouble *pp2, void* ggg);
int prad_lab2ff(ldouble *pp1, ldouble *pp2, void *ggg);
int prad_on2lab(ldouble *pp1, ldouble *pp2, void* ggg);
int prad_lab2on(ldouble *pp1, ldouble *pp2, void *ggg);
int prad_zamo2ff(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble GG[][5], ldouble elo[][4]);
int prad_ff2zamo(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble GG[][5], ldouble eup[][4]);
int multiply22(ldouble T1[][4],ldouble T2[][4],ldouble A[][4]);
int multiply2(ldouble *u1,ldouble *u2,ldouble A[][4]);
int trans22_coco(ldouble *xx,ldouble T1[][4],ldouble T2[][4],int CO1, int CO2);
int trans2_coco(ldouble *xx,ldouble *,ldouble *,int CO1, int CO2);
int trans_prad_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble*, void*,void*);
int trans_mhd_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble *xxvec, void*,void*);
int trans_pall_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble *xxvec, void* ggg1, void* ggg2);
int coco_3vector(ldouble A1[3],ldouble A2[3],int CO1,int CO2,void* ggg);


//rad.mf.c
int
calc_rad_wavespeeds_on_base_mf(ldouble *pp, ldouble*,ldouble *avaltop);
int redistribute_with_velocities(ldouble avals[6],ldouble A[NRF],ldouble,ldouble);
int calc_rad_wavespeeds_on(ldouble nx,ldouble ny,ldouble nz, ldouble *avals);
int mf_correct_in_azimuth(ldouble *pp, ldouble *uu, void* ggg,ldouble);
int set_radatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype);
int calc_Rij_mf(ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble Rij[][4][4]);
int calc_rad_wavespeeds_mf_total(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble tautot[3],ldouble *aval);
int calc_rad_wavespeeds_pure_mf_each(ldouble *pp,void*,ldouble aval[][6]);
int calc_Rij_ff_mf(ldouble *pp, ldouble  Rij[][4][4]);
int redistribute_radfluids(ldouble *pp, ldouble *uu0, void* ggg);
int redistribute_radfluids_m1(ldouble *pp, ldouble *uu0, void* ggg);
int redistribute_radfluids_m2(ldouble *pp, ldouble *uu0, void* ggg);
int redistribute_radfluids_m3(ldouble *pp, ldouble *uu0, void* ggg);
int redistribute_radfluids_along_axes(ldouble *pp, ldouble *uu0, void* ggg);
int redistribute_radfluids_at_cell(int ix,int iy,int iz);
int mf_correct_in_azimuth_at_cell(int ix,int iy,int iz,ldouble);



//rad.c
int solve_explicit_lab_core(ldouble *uu,ldouble *pp,void* ggg,ldouble dt,ldouble* deltas,int verbose);
int apply_rad_source_del4(int ix,int iy,int iz,ldouble *del4);
int test_if_rad_implicit(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5], ldouble *del4);
int implicit_lab_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5],ldouble tlo[][4], ldouble tup[][4],ldouble *pp);
ldouble calc_LTE_temp(ldouble *pp,void *ggg);
int test_solve_implicit_lab();
int test_jon_solve_implicit_lab();
int calc_LTE_state(ldouble *pp,ldouble *ppLTE,void *ggg);
int calc_ff_Rtt(ldouble *pp,ldouble *Rtt, ldouble* ucon,void* ggg);
int calc_rad_shearviscosity(ldouble *pp,void* ggg,ldouble shear[][4],ldouble *nuret,ldouble *vdiff2ret);
int calc_visc_Rij(ldouble *pp, void* ggg, ldouble T[][4], ldouble R[][4]);
int set_radatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype);
int calc_Rij(ldouble *pp, void*, ldouble Rij[][4]);
int calc_Rij_ff(ldouble *pp, ldouble  Rij[][4]);
int solve_explicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int);
int solve_implicit_ff(int ix,int iy,int iz,ldouble dt,ldouble* deltas);
ldouble calc_LTE_EfromT(ldouble);
ldouble calc_LTE_TfromE(ldouble);
ldouble calc_LTE_Efromurho(ldouble E,ldouble);
ldouble calc_PEQ_ufromTrho(ldouble,ldouble);
ldouble calc_PEQ_Tfromurho(ldouble,ldouble);
int calc_LTE_ff(ldouble,ldouble*,ldouble*,ldouble,int);
int solve_LTE_ff(int ix,int iy,int iz,ldouble dt);
int solve_LTE(int ix,int iy,int iz,ldouble dt);
int solve_radforce_ff(int ix,int iy,int iz,ldouble dt);
int solve_radforce(int ix,int iy,int iz,ldouble dt);
ldouble calc_chi(ldouble *pp, ldouble *xx);
int calc_tautot(ldouble *pp, ldouble *xx, ldouble *dl, ldouble *tautot);
int calc_tauabs(ldouble *pp, ldouble *xx, ldouble *dl, ldouble *tauabs);
int calc_Gi_ff(ldouble *pp, ldouble Gi[4]);
int calc_Gi(ldouble *pp, void*,ldouble Gi[4]);
int calc_rad_Jac_eval(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble *aval,int);
int calc_rad_wavespeeds(ldouble *pp,void*,ldouble tautot[3],ldouble *aval,int verbose);
int calc_rad_wavespeeds_pure(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble *aval);
int solve_implicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int verbose);

//magn.c
void calc_bcon_4vel(ldouble *pr, ldouble *ucon, ldouble *ucov, ldouble *bcon); 
int calc_BfromA();
ldouble calc_divB(int ix,int iy,int iz);
void calc_bcon_prim(double *pp, double *bcon, void* ggg);
void calc_Bcon_prim(double *pp, double *bcon,double *Bcon, void* ggg);


#include "mnemonics.h"

/*********************/
/*********************/
/*********************/
/*********************/
/*****  loops   ******/
/*********************/
/*********************/
/*********************/
/*********************/

/* loop over all primitives */
#define PLOOP(j) for(j=0;j<NV;j++)
/* loop over all Dimensions; second rank loop */
#define DLOOP(j,k) for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA(j) for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP(j,k) for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA(j) for(j=1;j<NDIM;j++)
