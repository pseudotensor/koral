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
#define surfdensCGS2GU(x)    (x*GGG/CCC/CCC*MASSCM)
#define surfdensGU2CGS(x)    (x/GGG*CCC*CCC/MASSCM)
#define massCGS2GU(x)    (x*GGG/CCC/CCC/MASSCM)
#define massGU2CGS(x)    (x/GGG*CCC*CCC*MASSCM)
#define kappaCGS2GU(x)  (x/GGG*CCC*CCC/MASSCM)
#define kappaGU2CGS(x)  (x*GGG/CCC/CCC*MASSCM)
//verify
#define endenCGS2GU(x) (x*GGG*MASSCM*MASSCM/CCC/CCC/CCC/CCC)
#define endenGU2CGS(x) (x/GGG/MASSCM/MASSCM*CCC*CCC*CCC*CCC)
#define fluxCGS2GU(x) (x*GGG*MASSCM*MASSCM/CCC/CCC/CCC/CCC/CCC)
#define fluxGU2CGS(x) (x/GGG/MASSCM/MASSCM*CCC*CCC*CCC*CCC*CCC)

//constants
#define K_BOLTZ (1.3806488e-16 * GGG / CCC / CCC / CCC / CCC / MASSCM)
#define M_PROTON massCGS2GU(1.67262158e-24)
#define M_ELECTR massCGS2GU(9.11e-28)
#define SIGMA_RAD (5.67e-5 * GGG / CCC / CCC / CCC / CCC / CCC * MASSCM * MASSCM)
#define A_RAD (4.*SIGMA_RAD)
#define MU_GAS 1.
#define Z_RATIO (1.0)
#define Pi (3.141592654)     
#define KAPPA_ES_COEFF (kappaCGS2GU(0.4))

//other stuff
#include "problem.h"
#include "mnemonics.h"

#include "mdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>

#ifdef MPI
#include <mpi.h>
#endif

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
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
#include <gsl/gsl_blas.h>
//#include <gsl/gsl_odeiv2.h>

#ifdef PR_DEFS
#include PR_DEFS
#endif

#define ldouble double

#define FTYPE ldouble
#define gSIZE 20 //size of metric arrays = 16 + 1 (gdet) + 3 (dlgdet)

//global variables, some of them to be distributed in some way over processes
ldouble global_time;
ldouble global_dt;
ldouble global_tstepdenmax;
ldouble start_time, end_time, mid1_time, mid2_time, maxmp_time;

ldouble avgtime,dt;
int nstep;
int global_int_slot[NGLOBALINTSLOT];
ldouble max_ws[3],max_dt,ttm1,ttm2,max_ws_ph;
ldouble tstepdenmax;
ldouble min_dx,min_dy,min_dz;

//for subzones
ldouble currentzonetime;
int currentzone;
int global_ix1,global_ix2;
int global_iy1,global_iy2;
int global_iz1,global_iz2;

//tile specific
int TI,TJ,TK; //tile order
int TOI,TOJ,TOK; //indices of the tile origin
int PROCID;
int NPROCS;

//arrays and stuff
ldouble **msgbufs;
ldouble *u,*x,*xb,*du,*ut1,*ut2,*ut3,*ut4,*ut0,*u_bak_fixup,*p_bak_fixup,
  *u_step1,*u_step2,*u_bak_subzone,*p_bak_subzone,
  *u_step3,*u_step4,*ahdx,*ahdy,*ahdz,*aradx,*arady,*aradz,*cell_tsteps,
  *dut0,*dut1,*dut2,*dut3,*uforget,*drt0,*drt1,*drt2,*drt3,
  *ahdxl,*ahdyl,*ahdzl,*aradxl,*aradyl,*aradzl,  *ahdxr,*ahdyr,*ahdzr,*aradxr,
  *aradyr,*aradzr,*p,*pinit,*pproblem1,*pproblem2,*emf,*ptemp1,*pvecpot,
  *ptm1,*ptm2,*pt0,*px,*py,*pz,*s,*g,*gbx,*gby,*gbz,*Gbx,*Gby,*Gbz,*pavg,
  *pbLx,*pbRx,*pbLy,*pbRy,*pbLz,*pbRz,*sbLx,*sbRx,*sbLy,*sbRy,*sbLz,*sbRz,*ubLx,*ubRx,*ubLy,*ubRy,*ubLz,*ubRz,
  *flbx,*flby,*flbz,*flLx,*flRx,*flLy,*flRy,*flLz,*flRz,*gKr,*gKrbx,*gKrby,*gKrbz,*G,
  *emuup,*emulo,*emuupbx,*emulobx,*emuupby,*emuloby,*emuupbz,*emulobz,
  *emuup2,*emulo2,*emuupbx2,*emulobx2,*emuupby2,*emuloby2,*emuupbz2,*emulobz2,
  *tmuup,*tmulo,*tmuupbx,*tmulobx,*tmuupby,*tmuloby,*tmuupbz,*tmulobz,
  *tmuup2,*tmulo2,*tmuupbx2,*tmulobx2,*tmuupby2,*tmuloby2,*tmuupbz2,*tmulobz2;
int *cellflag,**loop_0,**loop_02,**loop_1,**loop_2,**loop_3,**loop_4,**loop_5,**loop_6;
int Nloop_0,Nloop_1,Nloop_2,Nloop_02,Nloop_3,Nloop_4,Nloop_5,Nloop_6;

ldouble sigma_otg[TNX];
ldouble sigma_otg_temp[TNX];
ldouble scaleth_otg[TNX];
ldouble scaleth_otg_temp[TNX];


ldouble scalars[NSCALARS];
int doingavg;
int doingpostproc;
ldouble Kr_tmp[4][4][4],g_tmp[4][4];
ldouble inputarg[10];
int **gcidx;
FILE *fout1,*fout_scalars,*fout_radprofiles,*fout_fail;
int nfout1,nfout2;

//precalculated metric parameters
ldouble rhorizonBL,rISCOBL,rmboundBL,rphotonBL,etaNT;

//some macros
#define my_max(x,y) (x>y?x:y)
#define my_max3(x,y,z) (x>my_max(y,z)?x:my_max(y,z))

//geometry structure
struct geometry
{
  int coords;
  int ix,iy,iz;
  int ifacedim; //-1 - cell center, 0 - x face, 1 - y face, 2 - z face
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
  int par; //some parameter to be used by user
};

//main.c
int solve_the_problem(ldouble,char*);
int print_scalings(void);
 
//postproc.c
ldouble calc_resmri(ldouble);
ldouble calc_meantemp(ldouble);
ldouble calc_scaleheight(ldouble);
int calc_scalars(ldouble*,ldouble);
int calc_Bflux(ldouble radius,int type,ldouble*,ldouble*);
ldouble calc_mdotEdd();
ldouble calc_lumEdd();
int calc_radialprofiles(ldouble profiles[][NX]);
int calc_anarelradialprofiles(ldouble profiles[][NX]);
ldouble calc_totalmass();
ldouble calc_mdot(ldouble radius,int);
int calc_lum(ldouble radius,int,ldouble*,ldouble*);
ldouble calc_photloc(int ix);

//misc§<.c
int print_primitives(ldouble *u);
int print_conserved(ldouble *u);
int my_clock_gettime(void* tsptr);
int print_NVvector(ldouble v[4]);
int print_tensor(ldouble T[][4]);
int print_metric(ldouble T[][5]);
int print_4vector(ldouble v[4]);
int print_Nvector(ldouble v[4],int);
ldouble calc_eigen_4x4(ldouble g[][4], ldouble *ev);
ldouble calc_eigen_4x4symm(ldouble g[][4], ldouble *ev);
ldouble calc_eigen_3x3symm(ldouble g[][3], ldouble *ev);
int make_matrixposdef_3x3symm(ldouble g[][3]);
ldouble my_atan2(ldouble y, ldouble x);
void shuffle_loop(int **array, size_t n);
ldouble step_function(ldouble x,ldouble k);
int calc_stationary1d_solution()  ;
int initialize_arrays();
int free_arrays();
ldouble my_min(ldouble a, ldouble b);
ldouble my_min_N(ldouble* a, int);
ldouble my_max_N(ldouble* a, int );
//ldouble my_max(ldouble a, ldouble b);
ldouble my_sign(ldouble);
int find_eigenvalues3(ldouble[3][3],ldouble*);
int find_eigenvalues(ldouble*,int);
ldouble find_max_eigenvalue(ldouble *data,int N);
int find_max_eigenvalue_lr(ldouble *data,int N,ldouble*,ldouble*);
int my_err(char *);
int inverse_44matrix(ldouble a[][4], ldouble ia[][4]);
int inverse_matrix(ldouble *a, ldouble *ia, int N);
int convert_out2gif_1d(char *fname,char*,int niter,ldouble t);
int convert_out2gif_2d(char *fname,char*,int niter,ldouble t);
int getch(void);
int dosthelse(void);
ldouble opacity_BellLin(ldouble rhoc, ldouble Tc);

//finite.c

int calc_u2p();
int do_finger();
int correct_polaraxis();
int solve_implicit_metric(int ix,int iy,int iz,ldouble dt,ldouble *ubase);
int cell_fixup_rad();
int cell_fixup_hd();
ldouble f_calc_fluxes_at_faces(int ix,int iy,int iz);
int op_explicit(ldouble t, ldouble dt);
int op_implicit(ldouble t, ldouble dt);
int set_grid(ldouble*, ldouble*, ldouble*,ldouble*);
int alloc_loops(int,ldouble,ldouble);
int print_grid(ldouble,ldouble,ldouble);
ldouble fd_flux_limiter(ldouble r);
ldouble minmod_fd_flux_limiter(ldouble ,ldouble,ldouble);
ldouble f_der_kurganovtadmor(int ix,int iy, int yz,ldouble*);
ldouble f_der_hlle_obsolete(int ix,int iy, int yz,ldouble*);
ldouble f_der_muscl(int ix,int iy, int yz,ldouble*);
int copy_u(ldouble,ldouble*,ldouble*);
int copy_u_core(ldouble factor,ldouble *uu1,ldouble* uu2, int N );
int add_u(ldouble f1, ldouble* u1, ldouble f2, ldouble *u2, ldouble *u3);
int add_u_core(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3,int N);
int add_u_3(ldouble f1, ldouble* u1, ldouble f2, ldouble *u2, ldouble f3, ldouble *u3, ldouble *u4);
int add_u_core_3(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble f3, ldouble *uu3, ldouble *uu4,int N);
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

//ldouble get_x(int,int);
int get_xx(int ix,int iy,int iz,ldouble *xx);
int  get_xx_arb(int ix,int iy,int iz,ldouble *xx,int COORDSOUT);
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
#define get_cflag(iflag,ix,iy,iz) (cellflag[iflag + (iX(ix)+NGCX)*NFLAGS + (iY(iy)+NGCY)*(SX)*NFLAGS + (iZ(iz)+NGCZ)*(SY)*(SX)*NFLAGS])
#define set_cflag(iflag,ix,iy,iz,val) cellflag[iflag + (iX(ix)+NGCX)*NFLAGS + (iY(iy)+NGCY)*(SX)*NFLAGS + (iZ(iz)+NGCZ)*(SY)*(SX)*NFLAGS]=val
#define get_u(uarr,iv,ix,iy,iz) (uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV])
#define get_uavg(uarr,iv,ix,iy,iz) (uarr[iv + (iX(ix)+NGCX)*(NV+NAVGVARS) + (iY(iy)+NGCY)*(SX)*(NV+NAVGVARS) + (iZ(iz)+NGCZ)*(SY)*(SX)*(NV+NAVGVARS)])
//int set_u(ldouble*,int,int,int,int,ldouble);
#define set_u(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV]=val
#define set_uavg(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+NGCX)*(NV+NAVGVARS) + (iY(iy)+NGCY)*(SX)*(NV+NAVGVARS) + (iZ(iz)+NGCZ)*(SY)*(SX)*(NV+NAVGVARS)]=val
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
int save_avg(ldouble dt);

int fprint_coordfile(char* folder,char* prefix);
int fprint_coordBL(char* folder,char* prefix);
int fprint_restartfile(ldouble t, char* folder);
int fprint_restartfile_ascii(ldouble t, char* folder);
int fprint_restartfile_bin(ldouble t, char* folder);
int fprint_restartfile_mpi(ldouble t, char* folder);
int fread_restartfile(int,char *,ldouble*);
int fread_restartfile_bin(int,char *,ldouble*);
int fread_restartfile_ascii(int,char *,ldouble*);
int fread_restartfile_mpi(int,char *,ldouble*);

int fprint_avgfile(ldouble t, char* folder,char *);
int fprint_avgfile_ascii(ldouble t, char* folder,char *);
int fprint_avgfile_bin(ldouble t, char* folder,char *);
int fprint_avgfile_mpi(ldouble t, char* folder,char *);
int fread_avgfile(int,char *,ldouble*,ldouble *,ldouble *);
int fread_avgfile_bin(int,char *,ldouble*,ldouble *,ldouble *);
int fread_avgfile_ascii(int,char *,ldouble*,ldouble *,ldouble *);
int fread_avgfile_mpi(int,char *,ldouble*,ldouble *,ldouble *);

int fprint_simplefile(ldouble t, int nfile, char* folder, char* prefix);
int fprint_simplecart(ldouble t, int nfile, char* folder, char* prefix);
int fprint_simplesph(ldouble t, int nfile, char* folder, char* prefix);
int fprint_scalars(ldouble t, ldouble *scalars, int nscalars);
int fprint_radprofiles(ldouble t, int nfile, char* folder, char* prefix);
int fprint_anarelradprofiles(ldouble t, int nfile, char* folder, char* prefix, ldouble[NANARELRADPROFILES][NX]);
int fprint_gridfile(char* folder);
int fprint_openfiles(char *);
int fprint_closefiles();
int fprint_outfile(ldouble t, int nfile, int codeprim, char* folder, char *prefix);
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
ldouble calc_PEQ_rhofromTu(ldouble,ldouble);
ldouble calc_PEQ_ufromTrho(ldouble,ldouble);
ldouble calc_PEQ_csfromT(ldouble);
ldouble calc_PEQ_Tfromurho(ldouble,ldouble);
int f_metric_source_term_arb(ldouble *pp,void *ggg,ldouble *ss);
int calc_hd_shearviscosity(ldouble *pp,void* ggg,ldouble shear[][4],ldouble *nuret,ldouble *vdiff2ret);
int calc_shear_lab(ldouble *pp,void* ggg,ldouble S[][4],int hdorrad,int *);
int calc_shear_rad_lab(ldouble *pp,void* ggg,ldouble S[][4],int *);

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
int f_flux_prime(ldouble *uu, int,int,int,int,ldouble *ff,int lr);
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
ldouble calc_kappa(ldouble*,void*,ldouble*,ldouble*,ldouble*,ldouble*);
ldouble calc_kappaes(ldouble*,void*);
ldouble calc_ufromS(ldouble S,ldouble rho);
ldouble calc_Sfromu(ldouble S,ldouble u);
int
avg2point(ldouble *um2,ldouble *um1,ldouble *u0,ldouble *up1,ldouble *up2,ldouble*,ldouble*,ldouble dxm2,ldouble dxm1,ldouble dx0,ldouble dxp1,ldouble dxp2,int);

//problem.c
ldouble calc_xb(int i,int idim);
int calc_bc(int,int,int,ldouble,ldouble*, ldouble*,int,int);
int pr_tophat_inside(ldouble x,ldouble y,ldouble z);
int my_finger(ldouble);
int analytical_solution(ldouble t,int ix,int iy,int iz,ldouble *uu,ldouble *pp,ldouble *vv);

//metric.c

int fill_geometry(int ix,int iy,int iz,void *geom);
int fill_geometry_face(int ix,int iy,int iz,int,void *geom);
int fill_geometry_face_arb(int ix,int iy,int iz,int,void *geom,int);
int fill_geometry_arb(int ix,int iy,int iz,void *geom,int COORDS);
int calc_metric();
int calc_tetrades(ldouble g[][5], ldouble tmuup[][4], ldouble tmulo[][4],int);
int calc_ZAMOes(ldouble g[][5], ldouble emuup[][4], ldouble emulo[][4],int);
int dxdx_KS2BL(ldouble *xx, ldouble dxdx[][4]);
int dxdx_BL2KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2MKS1(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MKS12KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2MKS2(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MKS22KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MCYL12CYL(ldouble *xx, ldouble dxdx[][4]);
int dxdx_CYL2MCYL1(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MSPH12SPH(ldouble *xx, ldouble dxdx[][4]);
int dxdx_SPH2MSPH1(ldouble *xx, ldouble dxdx[][4]);
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
int coco_KS2MINK(ldouble *xBL, ldouble *xKS);

//relele.c
int calc_normalobs_4vel(ldouble GG[][5], ldouble *ncon);
int calc_normalobs_relvel(ldouble GG[][5], ldouble *ncon);
int set_hdatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype);
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
ldouble r_ISCO_BL(ldouble ac);
ldouble r_mbound_BL(ldouble a);
ldouble r_photon_BL(ldouble a);
int update_entropy(int ix,int iy,int iz,int u2pflag);
int conv_vels_both(ldouble *u1,ldouble *u2con,ldouble *u2cov,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
int conv_velscov(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
int conv_vels(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
int conv_vels_ut(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
int conv_vels_core(ldouble *u1,ldouble *u2,ldouble *u2cov,int which1,int which2,ldouble gg[][5],ldouble GG[][5],ldouble,int);

//u2p.c
int count_entropy(int *,int*);
int copy_entropycount();
int calc_primitives_local(int ix,int iy,int iz,ldouble *pp);
int calc_primitives(int,int,int,int,int);
int check_floors_mhd(ldouble *uu, int,void*);
int u2p(ldouble *uu, ldouble *pp, void*,int*,int*,int);
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
int p2avg(int,int,int,ldouble*);
int calc_conserved(int ix,int iy,int iz);
int p2u(ldouble *p, ldouble *u,void*);
int pff2u(ldouble *p, ldouble *u,ldouble[][5],ldouble[][4],ldouble[][4]);
int p2u_rad(ldouble *p,ldouble *u,void*);
int p2u_mhd(ldouble *p,ldouble *u,void*);

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
int indices_1121(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5]);
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

//rad.c
//rad-viscosity specific
ldouble Rijviscprev[SX][SY][SZ][4][4],radvisclasttime[SX][SY][SZ];

void reset_radviscaccel();
int update_intensities();
int calc_M1intensities(void);
ldouble calc_ncompt_nphlab(ldouble *pp, void* ggg);
ldouble calc_ncompt_Thatrad(ldouble *pp, void* ggg,ldouble);
ldouble calc_ncompt_Thatrad_full(ldouble *pp, void* ggg);
ldouble calc_ncompt_Thatrad_4vel(ldouble *pp, void* ggg,ldouble,ldouble *,ldouble *);
int calc_rad_visccoeff(ldouble *pp,void *ggg,ldouble *,ldouble *mfpret,ldouble *);
int f_flux_prime_rad( ldouble *pp, int idim, void *ggg,ldouble *ff);
int f_flux_prime_rad_total(ldouble *pp, void *ggg,ldouble Rij[][4],ldouble RijM1[][4], ldouble Rijvisc[][4]);
int solve_implicit_lab_4dprim(ldouble *uu00,ldouble *pp00,void *ggg,ldouble dt,ldouble* deltas,int verbose,int *params,ldouble *);
int solve_implicit_lab_4dcon(ldouble *uu00,ldouble *pp00,void *ggg,ldouble dt,ldouble* deltas,int verbose,ldouble *pp);
int solve_explicit_lab_core(ldouble *uu,ldouble *pp,void* ggg,ldouble dt,ldouble* deltas,int verbose);
int explicit_rad_source_term(int ix,int iy, int iz,ldouble dt);
int apply_rad_source_del4(int ix,int iy,int iz,ldouble *del4);
int test_if_rad_implicit(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5], ldouble *del4);
int implicit_lab_rad_source_term(int ix,int iy, int iz,ldouble dt);
ldouble calc_LTE_temp(ldouble *pp,void *ggg,int);
int test_solve_implicit_lab();
int test_jon_solve_implicit_lab();
int calc_LTE_state(ldouble *pp,ldouble *ppLTE,void *ggg);
int calc_LTE_state_temp(ldouble *pp,void *ggg);
int calc_ff_Rtt(ldouble *pp,ldouble *Rtt, ldouble* ucon,void* ggg);
int calc_normal_Rtt(ldouble *pp,ldouble *Rtt, ldouble* ucon,void* ggg);
int calc_rad_shearviscosity(ldouble *pp,void* ggg,ldouble shear[][4],ldouble *nuret,int *);
int calc_Rij_visc(ldouble *pp, void* ggg, ldouble T[][4],int *);
int calc_Rij_total(ldouble *pp, void* ggg, ldouble R[][4]);
int set_radatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype);
int calc_Rij(ldouble *pp, void*, ldouble Rij[][4]);
int calc_Rij_M1(ldouble *pp, void*, ldouble Rij[][4]);
int calc_Rij_M1_ff(ldouble *pp, ldouble  Rij[][4]);
int calc_Rij_Minerbo_ff(ldouble *pp, ldouble  Rij[][4]);
int radclosure_Edd(ldouble *pp, void *ggg, ldouble Rij[][4]);
int radclosure_Minerbo(ldouble *pp, void *ggg, ldouble Rij[][4]);
int radclosure_M1orto(ldouble *pp, void *ggg, ldouble Rij[][4]);
int radclosure_VET(ldouble *pp, void *ggg, ldouble Rij[][4]);
int solve_explicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int);
int solve_implicit_ff(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int);
int solve_implicit_ff_core(ldouble *uu0,ldouble *pp0,void* ggg,ldouble dt,ldouble* deltas,int verbose);
int solve_implicit_backup(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int);
int solve_implicit_backup_core(ldouble *uu0,ldouble *pp0,void* ggg,ldouble dt,ldouble* deltas,int verbose);
int solve_explicit_lab_core(ldouble *uu0,ldouble *pp0,void* ggg,ldouble dt,ldouble* deltas,int verbose);
ldouble calc_NFfromT(ldouble T);
ldouble calc_NFfromE(ldouble E);
ldouble calc_LTE_EfromT(ldouble);
ldouble calc_LTE_TfromE(ldouble);
ldouble calc_LTE_Efromurho(ldouble E,ldouble);
int calc_LTE_ff(ldouble,ldouble*,ldouble*,ldouble,ldouble,int);
int solve_LTE_ff(int ix,int iy,int iz,ldouble dt);
int solve_LTE(int ix,int iy,int iz,ldouble dt);
int solve_radforce_ff(int ix,int iy,int iz,ldouble dt);
int solve_radforce(int ix,int iy,int iz,ldouble dt);
ldouble calc_chi(ldouble *pp,void*);
int calc_tautot(ldouble *pp, void*, ldouble *dl, ldouble *tautot);
int calc_tauabs(ldouble *pp, void*, ldouble *dl, ldouble *tauabs);
int calc_Gi(ldouble *pp, void*,ldouble Gi[4],int);
int calc_Compt_Gi(ldouble *pp, void* ggg, ldouble *Gtc, ldouble Ehatrad, ldouble Tgas, ldouble kappaes, ldouble *ucon);
int calc_rad_Jac_eval(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble *aval,int);
int calc_rad_wavespeeds(ldouble *pp,void*,ldouble tautot[3],ldouble *aval,int verbose);
int calc_rad_wavespeeds_pure(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble *aval);
int solve_implicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int verbose);
ldouble calc_nsource(ldouble *pp, void* ggg);

//magn.c
int mimic_dynamo(ldouble);
int calc_curl(ldouble *ptu, ldouble idx, int ix, int iy, int iz, void* ggg, ldouble *curl);
int calc_angle_brbphibsq(int ix, int iy, int iz, ldouble *brbphi, ldouble *bsq,ldouble*, ldouble*);
int calc_angle_bpbphibsq(int ix, int iy, int iz, ldouble *brbphi, ldouble *bsq,ldouble*, ldouble*);
ldouble calc_BpBphi(int ix, int iy, int iz);
int adjust_fluxcttoth_emfs();
void calc_bcon_4vel(ldouble *pr, ldouble *ucon, ldouble *ucov, ldouble *bcon); 
int calc_BfromA(ldouble*,int);
int calc_BfromA_core();
ldouble calc_divB(int ix,int iy,int iz);
void calc_bcon_prim(double *pp, double *bcon, void* ggg);
void calc_Bcon_prim(double *pp, double *bcon,double *Bcon, void* ggg);
int flux_ct();
ldouble calc_Qtheta(int ix, int iy, int iz);


//silo.c
int fprint_silofile(ldouble time, int num, char* folder, char* prefix);

//mpi.c
int
mpi_isitBC(int BCtype);
void
mpi_global2localidx(int gix,int giy, int giz, int *lix, int *liy, int *liz);
void
mpi_local2globalidx(int lix,int liy, int liz, int *gix, int *giy, int *giz);
int calc_avgs_throughout();
int mpi_exchangedata();
void
mpi_myinit(int argc, char *argv[]);
void
mpi_myfinalize();
void
mpi_synchtiming(ldouble *time);


//zeroshort.c
//Define what is contained in each node
struct bsptree
{
	int angIndex;  //specify which angle index this node corresponds to
	int iter;	//specify which axis we do the median cutting

	struct bsptree *lower;
	struct bsptree *upper;
};

void reflectI(double reflect_direction[3], double I_start[NUMANGLES], double I_return[NUMANGLES]); 
int get_angDualIndex(double targetAng[3], double angGridCoords[NUMDUALANGLES][3]);
void bspGetNearestDualNeighbor(double targetAng[3], double angGridCoords[NUMDUALANGLES][3], struct bsptree *bspCurrentLoc, double *bestDistance, int *bestIndex);
int ZEROtest_oldmain();
int zero_readangles();
struct bsptree *angGridRoot;
struct bsptree *angDualGridRoot;
void setupInterpWeights_cart2D(int ix, int iy, int iz, double angGridCoords[NUMANGLES][3], int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4], double intersectGridWeights[SX][SY][SZ][NUMANGLES][4], double intersectDistances[SX][SY][SZ][NUMANGLES]);
void setupInterpWeights_cart3D(int ix, int iy, int iz, double angGridCoords[NUMANGLES][3], int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4], double intersectGridWeights[SX][SY][SZ][NUMANGLES][4], double intersectDistances[SX][SY][SZ][NUMANGLES]);
void setupInterpWeights_sph3D(int ix, int iy, int iz, double angGridCoords[NUMANGLES][3], int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4], double intersectGridWeights[SX][SY][SZ][NUMANGLES][4], double intersectDistances[SX][SY][SZ][NUMANGLES]);
void setupInterpWeights_sph2D(int ix, int iy, int iz, double angGridCoords[NUMANGLES][3], int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4], double intersectGridWeights[SX][SY][SZ][NUMANGLES][4], double intersectDistances[SX][SY][SZ][NUMANGLES], double intersectGridPhi[SX][SY][SZ][NUMANGLES]);


void ZERO_shortCharI(int,int,int,double delta_t, double I_Data[3][3][3][NUMANGLES], double source_Data[3][3][3][4],
		     double I_return[NUMANGLES],int verbose);

void initAngIndex(double angGridCoords[NUMANGLES][3], double angDualGridCoords[NUMDUALANGLES][3], int angGridIndexSort[NUMANGLES][3], int angDualGridIndexSort[NUMANGLES][3]);
void splitAngGrid(int numAvailAnglesInit, int angGridIndexInit[NUMANGLES][3], int iter, double angGridCoords[NUMANGLES][3], struct bsptree **node);
void splitDualAngGrid(int numAvailAnglesInit, int angGridIndexInit[NUMDUALANGLES][3], int iter, double angGridCoords[NUMDUALANGLES][3], struct bsptree **node);
void
ZERO_calcVET(int,int,int,double I_time[NUMANGLES], double eddingtonFactor[3][3], double angGridCoords[NUMANGLES][3]);

void transformI(int,int,int,double I_return[NUMANGLES], double*);
void transformI_stretch(int,int,int,double I_return[NUMANGLES], double*);
void transformI_quad(int,int,int,double I_return[NUMANGLES], double*);
int transformI_Lagrange(int,int,int,double I_return[NUMANGLES], double*);
void ZERO_decomposeM1(int,int,int,double M1_Data[5], double I_return[NUMANGLES]);
double
calc_stretchFactor(void *argsin);


double Ibeam[SX][SY][SZ][NUMANGLES];                //specific intensities at cell centers
double Ibeam2[SX][SY][SZ][NUMANGLES];                //specific intensities at cell centers _ auxiliary
double angGridCoords[NUMANGLES][3];  		//Store xyz locations of angle grid
double angDualGridCoords[NUMDUALANGLES][3]; 	//Store xyz locations of dual angle grid
int dualAdjacency[NUMDUALANGLES][3]; 		//Store index information for adjacent angles
int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4];	//indices for gridLoc of 4 points bounding intersection
double intersectGridWeights[SX][SY][SZ][NUMANGLES][4];	//weights corresponding to 4 points bounding intersection
double intersectDistances[SX][SY][SZ][NUMANGLES];		//distance to intersection point from center
double intersectGridPhi[SX][SY][SZ][NUMANGLES];		//phi value of intersection point used in axisymmetryc in spherical-like
double carttetrad[SX][SY][SZ][3][3]; //cartesian components of the local RADCLOSURECOORDS tetrad

#ifdef MPI
MPI_Group mpi_all_group;
#define MPI_LDOUBLE MPI_DOUBLE
void mpi_procid2tile(int procid, int* tilei, int* tilej, int* tilek);
int mpi_tile2procid(int tilei, int tilej, int tilek);
void mpi_tileorigin(int ti, int tj, int tk, int* toi, int* toj, int* tok);
void mpi_synchtiming(ldouble*);
int mpi_senddata(MPI_Request *reqs, int *nreqs);
int mpi_recvdata(MPI_Request *reqs, int *nreqs);
int mpi_savedata();
int mpi_isitBC(int BCtype);
void mpi_myinit(int argc, char *argv[]);
void mpi_myfinalize();
void mpi_global2localidx(int gix,int giy, int giz, int *lix, int *liy, int *liz);
void mpi_local2globalidx(int lix,int liy, int liz, int *gix, int *giy, int *giz);
/*
#ifdef CALCSIGMAONTHEGO
MPI_Group mpi_inttotal_group[NTX], mpi_intbelow_group[NTX]; 
MPI_Comm mpi_inttotal_comm[NTX], mpi_intbelow_comm[NTX]; 
#endif
*/
#endif

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

//multisteps
int mstep_cell_levels[NX][NY][NZ];
ldouble mstep_multiplier[NUMMSTEPLEVELS+1];
int mstep_current_counts[NUMMSTEPLEVELS];

//mstep.c
int mstep_calc_level(ldouble,ldouble);
int mstep_is_cell_active(int,int,int);
int mstep_is_level_active(int);
int mstep_init();
int mstep_iterate(void);
int mstep_test();
int mstep_update_levels();
int mstep_print_levels();
int mstep_is_cell_or_neighbour_active(int ix, int iy, int iz,int idim);
int mstep_is_face_active(int ix, int iy, int iz,int idim);
ldouble mstep_get_cell_multiplier(int ix,int iy,int iz);
ldouble mstep_get_face_multiplier(int ix,int iy,int iz,int dim);


