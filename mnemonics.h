//mnemonical definitions

//rad vs hydro
#define RAD 1
#define MHD 2 

//velocities
#define VEL4 1 //lab four-velocity u^i
#define VEL3 2 //lab three-velocity u^i/u^t
#define VELR 3 //relative velocity \tilde u^i

//primitive/conserved ordering
#define RHO 0
#define UU 1
#define VX 2
#define VY 3
#define VZ 4
#define ENTR 5
#define TRA (NVHD-1)

#define B1 (NVHD+0)
#define B2 (NVHD+1)
#define B3 (NVHD+2)

//multi rad-fluid macros
#define EE(nf) (NVMHD+nf*4)
#define FX(nf) (NVMHD+nf*4+1)
#define FY(nf) (NVMHD+nf*4+2)
#define FZ(nf) (NVMHD+nf*4+3)

//single fluid macros
#define EE0 EE(0)
#define FX0 FX(0)
#define FY0 FY(0)
#define FZ0 FZ(0)

//u2p inversion types
#define U2P_HOT 0
#define U2P_ENTROPY 1
#define U2P_HOTMAX 2
#define U2P_COLD 3

//coordinates/metric
#define BLCOORDS 1
#define SCHWCOORDS 1
#define KERRCOORDS 1
#define KSCOORDS 2
#define MKSCOORDS 3
#define MINKCOORDS 4
#define CYLCOORDS 5
#define SPHCOORDS 6
#define MKS1COORDS 7
#define MKS2COORDS 10
#define MCYL1COORDS 8
#define MKER1COORDS 9

//type of boundary
#define XBCLO 1
#define XBCHI 2
#define YBCLO 3
#define YBCHI 4
#define ZBCLO 5
#define ZBCHI 6

//cell flags
#define NFLAGS 5

//boolean
#define ENTROPYFLAG 0
#define RADSOURCEWORKEDFLAG 1
#define HDFIXUPFLAG 2
#define RADFIXUPFLAG 3

//values for RADFIXUP
#define RADSOURCETYPEFLAG 4
#define RADSOURCETYPEEXPLICIT 1
#define RADSOURCETYPEEXPLICITSUBSTEP 2
#define RADSOURCETYPEIMPLICITLAB 3
#define RADSOURCETYPEIMPLICITFF 10

//fluxes
#define LAXF_FLUX 0
#define HLL_FLUX 1

//timestepping
#define RK2IMEX 1
#define RK2 0
#define RK2HEUN 2

//types of hd/rad viscosity
#define NOVISCOSITY 0
#define SIMPLEVISCOSITY 1
#define SHEARVISCOSITY 2

//rad.implicit solver parameters
#define RADIMPLICIT_ENERGYEQ 0
#define RADIMPLICIT_ENTROPYEQ 1
#define RADIMPLICIT_LTEEQ 2

#define RADIMPLICIT_LAB 0
#define RADIMPLICIT_FF 1

//global integer slots
#define NGLOBALINTSLOT 13
#define GLOBALINTSLOT_NIMPENERRAD 0
#define GLOBALINTSLOT_ITERIMPENERRAD 1
#define GLOBALINTSLOT_NIMPENERMHD 2
#define GLOBALINTSLOT_ITERIMPENERMHD 3
#define GLOBALINTSLOT_NIMPENTRRAD 4
#define GLOBALINTSLOT_ITERIMPENTRRAD 5
#define GLOBALINTSLOT_NIMPENTRMHD 6
#define GLOBALINTSLOT_ITERIMPENTRMHD 7
#define GLOBALINTSLOT_NIMPLTE 8
#define GLOBALINTSLOT_ITERIMPLTE 9
#define GLOBALINTSLOT_NRADFIXUPS 10
#define GLOBALINTSLOT_NCRITFAILURES 11
#define GLOBALINTSLOT_NTOTALCRITFAILURES 12

//frames
#define ZAMOFRAME 0
#define FFFRAME 1

//avg quantities
#define AVGBSQ (NV+0)
#define AVGUCON(i) (NV+1+i)
#define AVGUCOV(i) (NV+5+i)
#define AVGBCON(i) (NV+9+i)
#define AVGBCOV(i) (NV+13+i)
#define AVGRHOUCON(i) (NV+17+i)
#define AVGRHOUCOV(i) (NV+21+i)
#define AVGUUUCON(i) (NV+25+i)
#define AVGUUCOV(i) (NV+29+i)
#define AVGBSQUCON(i) (NV+33+i)
#define AVGBSQUCOV(i) (NV+37+i)
#define AVGRHOUCONUCOV(i,j) (NV+41+i*4+j)
#define AVGUUUCONUCOV(i,j) (NV+57+i*4+j)
#define AVGBSQUCONUCOV(i,j) (NV+73+i*4+j)
#define AVGBCONBCOV(i,j) (NV+89+i*4+j)
#define AVGWUCON(i) (NV+105+i)
#define AVGEHAT (NV+109)
#define AVGRIJ(i,j) (NV+110+i*4+j)
#define AVGEHATUCON(i) (NV+126+i)
#define AVGEHATUCOV(i) (NV+130+i)
#define AVGURFCON(i) (NV+134+i)
#define AVGURFCOV(i) (NV+138+i)


//MPI mnemonics
#define MPI_MSG_TIME 100
