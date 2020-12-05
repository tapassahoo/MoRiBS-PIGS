#ifndef _MC_PIQMC_H
#define _MC_PIQMC_H 1
#include "rngstream.h"
#include "omprng.h"

void MCMolecularMove(int);
void MCMolecularMovePIGS(int);
void MCMolecularMoveNaive(int);
void MCTransLinStepPIGS(int,int,int,int,double,double,double,double,double &,double &);
void MCMolecularMoveNaiveNorm(int);
void MCTransLinStepPIGSNorm(int,int,int,int,double,double &,double &);
double GetTransDensityPIGS(int, int, int, int, int, double **);
double Gauss(double, double);	
#ifdef GAUSSIANMOVE
void MCMolecularMoveGauss(int);
#endif
void MCBisectionMove(int,int);
void MCBisectionMovePIGS(int,int);
void MCRotationsMove(int);
void MCRotationsMoveIndex(int);
void MCRotationsMoveCL(int);
void MCRotLinStep(int,int,int,int,double,double,double,double,double &,double &);
void MCRotLinStepPIMC(int,int,int,int,double,double,double,double,double &,double &);
void MCRotLinStepIndex(int,int,int,int,double,double,double,double &,double &);
void MCRotLinStepPIGS(int,int,int,int,double,double,double,double,double &,double &);
void MCRotLinStepCL(int,int,double,double,double,int,double,double &,double &, RngStream *);
int ClusterGrowth(int,double *,int,int,int,int, RngStream *, double **);
void MCRotLinStepSwap(int,int,int,int,double,double,double,double,double &,double &, string);
void MCRotLinStepSwapBroken(int,int,int,int,double,double,double,double,double &,double &);
void MCSwap(int, double, string &);
double PotRotEnergyPIMC(int, double *,int );   
double PotRotEnergyPIGS(int, double *,int , int );   
double PotRotEnergySwap(int,int,const double *,int , string);   
double PotRotE3DSwap(int,int,double *,int , string);   
double PotRotE3DBrokenPath(int,int,double *,int);   
double PotRotEnergySwapBroken(int, double *,int);   
//double PotRotEnergySwap(int,double **,int it, int );   
// Toby adds rotation move for nonlinear rotor
void MCRotations3D(int);
void MCRot3Dstep(int, int, int, int, double,double,double,double,double,int, int, double &, double &);
void MCRot3DstepPIGS(int, int, int, int, double,double,double,double,double,int, int, double &, double &);
void MCRot3DstepSwap(int, int, int, int, double,double,double,double,double,int, int, double &, double &, string);
void MCRot3DstepBrokenPath(int, int, int, int, double,double,double,double,double,int, int, double &, double &);
void Reflect_MF_XZ(void);
void Reflect_MF_YZ(void);
void Reflect_MF_XY(void);
void RotSymConfig(void);

void MCMolecularMoveExchange(int);
void MCBisectionMoveExchange(int,int);

double PotEnergy(int,double **);
double PotEnergy(int,double **,int);
double PotEnergyPIGS(int,double **);
double PotEnergyPIGS(int,double **,int);

double PotRotEnergy(int,double **,int it);   
double PotRotE3D(int,double *,int it);   
double PotRotE3DPIGS(int,double *,int it);   
double PotRotE3DHF(int,double **,int it);   

extern double  **MCTotal;  // MC counters (total number of moves)
extern double  **MCAccep;  // MC counters (number of accepted moves)
extern double  **MCTotalEndBeads; 
extern double  **MCAccepEndBeads; 
extern double MCTransChunkAcp;    // total number of translational moves for one chunk loop
extern double MCTransChunkTot;    // total accept number of translational moves for one chunk loop

// counters for parallel 3d rotation move.  They are supposed to be declared for all CPUs
extern double MCRotChunkAcp;    // total number of rotational moves for one chunk loop
extern double MCRotChunkTot;    // total accept number of rotational moves for one chunk loop
extern double MCRotTot; // the sum of all MCRotChunkTot from all CPUs
extern double MCRotAcp; // the sum of all MCRotChunkAcp from all CPUs


void ResetMCCounts(void);
void MemAllocMCCounts(void);
void MFreeMCCounts(void);

extern int PrintYrfl; // integer flag for printing reflected coordinates
extern int PrintXrfl; // integer flag for printing reflected coordinates
extern int PrintZrfl; // integer flag for printing reflected coordinates

//Last two lines added by Tapas Sahoo
double DotProduct(double *, double *); 
void CrossProduct(double *, double *, double *);
int myRand(double *, double );
int findCeil(double *, double);
double GetDensityENT(string, int, int, int, int, int, int, int, int, int, int, int, int, double, double, double **);
double GetDensity3DENT(string, int, int, int, int, int, int, int, int, int, double *, double *, double *);
double GetDensity3DPIGS(int, double *, double *, double *);
void CodeExit(int );
#endif  //mc_pimc.h
