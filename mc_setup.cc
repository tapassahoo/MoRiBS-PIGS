#include "mc_setup.h"
#include "mc_utils.h"
#include "mc_confg.h"
#include "mc_const.h"
#include "rngstream.h"
#include "omprng.h"
#include <cstdlib>
#include <omp.h>

#include <math.h>

//------------- MC FLAGS ---------------------------------

bool    WORM     = false;      // use the worm algorithm 

//------------- MC SYSTEM ---------------------------------

bool    IMPURITY = false;    // set true if there is a molecule in the system
bool    MINIMAGE = false;    // set true to apply minimum image convention 
bool    CRYSTAL = false;     // set true to use crystal geometry
bool    FCC = false;
bool    HCP = false;

int     IMTYPE   = -1;       // atom type for dopant molecule (rotational degrees of freedom)

int     ISPHER   = 0;        // whether to treat the asymmetric top dopant as spherical particle. 0: no; 1; yes
int     IREFLY   = 0;        // whether to reflect all particles wrt the xz plane of the dopant, 0: no; 1: yes
int     IREFLX   = 0;        // whether to reflect all particles wrt the yz plane of the dopant, 0: no; 1: yes
int     IREFLZ   = 0;        // whether to reflect all particles wrt the xy plane of the dopant, 0: no; 1: yes
int     IROTSYM  = 0;        // whether to rotate the dopants by their body-fixed axis, 0: no; 1: yes
int     NFOLD_ROT;           // foldness of rotational symmetry of the dopant

bool    ROTATION = false;    // set to 1 to account for the rotational degrees of freedom
bool    TRANSLATION = false;    // set to 1 to account for the translational degrees of freedom

bool    BOSONS   = false;    // true if there're bosons in the system
int     BSTYPE   = -1;       // atom type for bosons
 
bool    FERMIONS = false;    // true if there're fermions in the system
int     FERMTYPE = -1;       // atom type for fermions

int     NUMB_ATOMS = 0;      // total number of atoms  
int     NUMB_MOLCS = 0;      // total number of molecules  

int     NUMB_ATOMTYPES = 0;  // total number of atoms types 
int     NUMB_MOLCTYPES = 0;  // total number of molecules types
 
int     NDIM;
double  Temperature;

double Distance;
int    NDIMinit;
double DipoleMoment;
double DipoleMomentAU2;
double RR;
int RefAtom;
bool PIGS_SIM = false;
bool PIMC_SIM = false;
bool ENT_SIM = false;
string ENT_ENSMBL;
string ENT_ALGR;

#ifdef EWALDSUM
double prefSelf, prefBfun, alpha, alpha2, prefUk1, prefUk2, prefUk3, boxLength;
int    KMAX;
#endif

double  Density;
double *BoxSize;             // size of the entire box (in NDIM dimensions)
double  UnitCell[3];         // size of a single crystal unit cell (containing 4 atoms)
int     N1d[3];      

double  MCBeta;
double  MCTau;      // imaginary time step 

double  MCRotTau;   // imaginary time step for rotational degrees of freedom

int RotDenType = 0; // Rotational Density Type, exact (default 0) or rattle and shake(1)
int RotOdEvn = 0; // rotor symmetry info: -1 distinguishable, 0 pH2O type, 1 oH2O type
double RotEoff = 0.0; // offset for the rotational energy estimator taken from Noya's formula at relative Euler angles = 0 0 0
double X_Rot = 0.0; // rotational constant A in the unit of cm-1
double Y_Rot = 0.0; // rotational constant B in the unit of cm-1
double Z_Rot = 0.0; // rotational constant C in the unit of cm-1
int    RNratio = 1; // ratio between RS and Noya steps in hybrid rotational energy estimation

int    NumbRotLim = 100; // limit of number of one type of rotors

int     NumbAtoms;  // total number of atoms and molecules
int     NumbTypes;  // Number of particles' types

TParticle MCAtom[MAX_NUMBER_TYPES];  // size should be NumbTypes+1

// -------------- MC TABLES -------------------------------   

int  * MCType;     // convert atom number into atom type
int  * PIndex;     // permutation index
int  * RIndex;     // inverse permutation index

//------------- MC PARAMETERS -----------------------------

int NumbTimes;     // number of time slices (number of slices per atom)
int MaxnTimes;     // NumbTimes*NumbAtoms - total number of slices 

int NumbRotTimes;  // number of time slices for rotational degrees of freedom
int RotRatio;      // RotRatio = NumbTimes/NumbRotTimes;

long int NumberOfMCPasses; // number of steps within a block
long int NumberOfMCBlocks; // number of blocks
long int InitialBlock; // starting block number
long int NumberOfEQBlocks; // number of equilibration blocks 

// number of MC steps to skip ...

int MCSKIP_RATIO = 100000;    //  to save information regarding the accept ratio
int MCSKIP_TOTAL = 10000;     //  to save accumulated average
int MCSKIP_AVERG = 1;         //  to evaluate averages

//------------- MC STATUS ----------------------------------

long int MCStartBlock;

//---------------- MPI PARALL ---------------
int NProcs; // the number of processors as a global variable
int chunksize;  // the size of a chunk of rotational time slices treated by MPI
int tagrunning; // tag for the running MPI tag, which runs with the loop but alway within 0-32767
long int SEED; // random seed that depends on CPU id

//---------------- OpenMP PARLL --------------
int NThreads; // the number of threads as a global variable

//------------ MC DATA STORAGE -----------------------------

double ** MCCoords;   // translational degrees of freedom
double ** MCCosine;   // orientational cosines
double ** MCCosinex;   // orientational cosines
double ** MCCosiney;   // orientational cosines
double ** MCAngles;
double ** DipoleCoords;   // translational degrees of freedom

#ifdef PROPOSED
double dcost;
double dphi;
int NPHI = 361;
int NCOST = 101;
double *costProposed;
double *phiProposed;
#endif

//------------ Initial MCCoords and MCAngles;
double * MCCooInit;   // store the read in MCCoords
double * MCAngInit;   // store the read in MCAngles

//------------ test the data structure of MCCoords ---------
//double ** TZMAT; // just a matrix

double ** newcoords;  // buffer for new coordinates
#ifdef PROPOSED
double ** tempcoords;  // buffer for new coordinates
#endif
int    *  atom_list;  // buffer for atom labels

double *  rhoprp;     // rotatinal propagator for non-linear rotor
double *  erotpr;     // rotational energy estimator for non-linear rotor
double *  erotsq;     // rotational energy square estimator for non-linear rotor

int       InitMCCoords; // integer flag for read in MCCoords;

#ifdef GAUSSIANMOVE
double ** gausscoords;
int size_Gauss;
double * AMat;
double * Eigen;
#endif

//----------------------------------------------------------
#ifdef CHAINCONFIG
void initChain_config(double **);
#endif
#ifdef LATTICECONFIG
void initLattice_config(double **);
#endif
#ifdef WATERCLUSTER
void initCluster_config(double **, double **);
#endif
void replInitial_config(double **);

void MCMemAlloc(void)  // allocate  memmory 
{
	BoxSize   = new double [NDIM];

  MCCoords  = doubleMatrix (NDIM,NumbAtoms*NumbTimes);  
  newcoords = doubleMatrix (NDIM,NumbAtoms*NumbTimes); 
  //DipoleCoords  = doubleMatrix (NDIM,NumbAtoms);  

  // TZMAT = doubleMatrix (NDIM,NDIM);
 
  atom_list = new int [NumbAtoms];

  MCCosine  = doubleMatrix (NDIM,NumbAtoms*NumbTimes); 
  //MCCosinex  = doubleMatrix (NDIM,NumbAtoms*NumbTimes); 
  //MCCosiney  = doubleMatrix (NDIM,NumbAtoms*NumbTimes); 
  MCAngles  = doubleMatrix (NDIM,NumbAtoms*NumbTimes); 

  MCCooInit = new double [NDIM*NumbAtoms*NumbTimes];
  MCAngInit = new double [NDIM*NumbAtoms*NumbTimes];

  // non-linear rotor
  rhoprp    = new double [SizeRotDen];
  erotpr    = new double [SizeRotDen];
  erotsq    = new double [SizeRotDen];

  // TABLES
  
  MCType    = new int [NumbAtoms];
  PIndex    = new int [NumbAtoms];
  RIndex    = new int [NumbAtoms];

/*
#ifdef PROPOSED
  	tempcoords = doubleMatrix (NDIM,NumbAtoms*NumbTimes); 
	costProposed = new double [NCOST*NPHI];
	phiProposed  = new double [NCOST*NPHI];
#endif
*/
#ifdef GAUSSIANMOVE
  	gausscoords = doubleMatrix (NDIM,NumbAtoms*NumbTimes); 
	AMat    = new double [NumbTimes*NumbTimes];
	Eigen   = new double [NumbTimes];
#endif
}

void MCMemFree(void)  //  free memory
{
	delete [] BoxSize;
  free_doubleMatrix(MCCoords);  
  free_doubleMatrix(newcoords); 
  //free_doubleMatrix(DipoleCoords);  

  free_doubleMatrix(MCCosine); 
  //free_doubleMatrix(MCCosinex); 
  //free_doubleMatrix(MCCosiney); 
  free_doubleMatrix(MCAngles); 

  delete [] atom_list;

  delete [] rhoprp;
  delete [] erotpr;
  delete [] erotsq;


  delete [] MCType;
  delete [] PIndex;
  delete [] RIndex;

/*
#ifdef PROPOSED
  	free_doubleMatrix(tempcoords); 
	delete [] costProposed;
	delete [] phiProposed;
#endif
*/
#ifdef GAUSSIANMOVE
  	free_doubleMatrix(gausscoords); 
	delete [] AMat;
	delete [] Eigen;
#endif
}

//------------ MC SYSTEM OF UNITS --------------------------

TSystemOfUnits Units;

//-------------------------
int   MPIsize;    // MPI
int   MPIrank;    // MPI
//-----------------------------

void MCSetUnits(void)
{
	Units.temperature = 1.0;   // Kelvin
	Units.energy      = 1.0;   // Kelvin
	Units.length      = 1.0;   // Angstrom
	Units.mass        = 1.0;   // amu
	Units.mass        = 1.0;   // amu

	Units.senergy     = "Kelvin";
	Units.slength     = "Angstrom";

	Temperature      /= Units.temperature;
	Density          *= (Units.length*Units.length*Units.length);

	double lambda     = 100.0*(HBAR*HBAR)/(AMU*K_B);   // \AA^2 K

	for (int atype=0;atype<NumbTypes;atype++)
    {
		MCAtom[atype].mcstep /= Units.length;
		MCAtom[atype].mass   /= Units.mass;
		MCAtom[atype].brot   /= Units.energy;

		MCAtom[atype].lambda  = 0.5*lambda/MCAtom[atype].mass;  // (hbar^2/2m)
    } 
}

void MCSetUnits_HO_TEST(void)
{
  Units.temperature = 1.0;   // Kelvin
  Units.energy      = 1.0;   // Kelvin
  Units.length      = 1.0;   // Angstrom
  Units.mass        = 1.0;   // amu

  Units.senergy     = "arb units";
  Units.slength     = "arb units";

  Temperature      /= Units.temperature;
  Density          *= (Units.length*Units.length*Units.length);

  double lambda     = 1.0;   // \AA^2 K

  for (int atype=0;atype<NumbTypes;atype++)
    {
      MCAtom[atype].mcstep /= Units.length;
      MCAtom[atype].mass   /= Units.mass;
      MCAtom[atype].brot   /= Units.energy;

      MCAtom[atype].lambda  = 0.5*lambda/MCAtom[atype].mass;  // (hbar^2/2m)
    } 
}

void MCInitParams(void)
{
	const char *_proc_=__func__;    // "MCInitParams()";

	double mass;
	double brot;

	for (int atype=0;atype<NumbTypes;atype++)
	{
    	string stype = MCAtom[atype].type;

      	if  (stype == HE4)
		{
          	mass = MASS_HE4; 
          	brot = 0.0; 
		}
      	else
		if (stype ==H2)
	  	{
	    	mass=MASS_H2;
	    	brot = 0.0; 
	  	}
		else
	  	if  (stype == OCS)
	    {
	      	mass =(MASS_O16 + MASS_C12 + MASS_S32);
	      	brot = B_OCS;
	    } 
	  	else
	    if  (stype == N2O)
	    {
			mass = (2.0*MASS_N14 + MASS_O16);
			brot = B_N2O;
	    } 
	    else 
	    if  (stype == CO2)
		{
		  	mass =(MASS_C12 + 2.0*MASS_O16);
		  	brot = B_CO2;
		} 
	    else
		if  (stype == CO)
		{
		    mass =(MASS_C12 + MASS_O16);
		    brot = B_CO;
		}
		else
		if  (stype == HCN)
		{
		    mass =(MASS_H1 + MASS_C12 + MASS_N14);
		    brot = B_HCN;
		} 
		else 
		if  (stype == HCCCN)
		{
			mass =(MASS_H1 + 3.0*MASS_C12 + MASS_N14);
			brot = B_HCCCN;
		} 
		else
		if  (stype == H2O)
		{
			mass =(2.0*MASS_H1 + MASS_O16);
		}
		else
		if (stype == CH3F)
		{	
		    mass=(MASS_C12+3.0*MASS_H1+MASS_F19);
		}	
		else 
		if  (stype == SO2)
			mass =(2.0*MASS_O16 + MASS_S32);
		else
		if  (stype == HCOOCH3)
		    mass =(4.0*MASS_H1 + 2.0*MASS_O16 + 2.0*MASS_C12);
		else 
		if (stype == HF)
		    mass =(MASS_H1 + MASS_F19);
		else 
		nrerror(_proc_,"Unknown atom/molecule type");

      	MCAtom[atype].mass = mass;
      	MCAtom[atype].brot = brot;
    }
}

void MCInit(void)  // only undimensional parameters in this function 
{
	const char *_proc_=__func__;    // "MCInit()";

	//  INITIALIZE MC TABLES ---------------------------------------

	NDIMinit = 0;
	int natom = 0;    // map atom number into atom type
	for (int type = 0; type < NumbTypes; type++)
	for (int atom = 0; atom < MCAtom[type].numb; atom++)
	{
		MCType[natom] = type;
		natom++;
	}

	for (int atom = 0; atom < NumbAtoms; atom++)
	{ 
		PIndex[atom] = atom;
		RIndex[atom] = atom;
	}
	// ------------------------------------------------------------

	// BoxSize  =  pow((double)NumbAtoms/Density,1.0/(double)NDIM); 
	// define a box size based on number of atoms only (molecules excluded)
  
   double volume = (double)(NUMB_ATOMS+NUMB_MOLCS)/Density;
   if(CRYSTAL)
   {
      if(FCC)
      {
         UnitCell[0] = pow(volume/N1d[0]/N1d[1]/N1d[2], 1.0/NDIM);
         UnitCell[1] = UnitCell[0];
         UnitCell[2] = UnitCell[0];
      }
      else if(HCP)
      {
         UnitCell[0] = pow(volume/N1d[0]/N1d[1]/N1d[2]/sqrt(8.0), 1.0/NDIM);
         UnitCell[1] = sqrt(3.0)*UnitCell[0];
         UnitCell[2] = sqrt(8.0/3.0)*UnitCell[0];
      }
      cout<<"UnitCell: "<<UnitCell[0]<<" "<<UnitCell[1]<<" "<<UnitCell[2]<<endl;

      for(int id=0;id<NDIM;id++)
      BoxSize[id] = UnitCell[id]*N1d[id];
   }
   else
   {
      for (int id=0;id<NDIM;id++)
      BoxSize[id] = pow(volume, 1.0/NDIM);
   }
   cout<<"BoxSize: "<<BoxSize[0];
   for (int id=1;id<NDIM;id++)
   cout<<" "<<BoxSize[id];
   cout<<endl;


	MCBeta   =  1.0/Temperature;

	if (PIMC_SIM) MCTau = MCBeta/(double)NumbTimes;
	else	      MCTau = MCBeta/((double)NumbTimes-1.0);

  	if (ROTATION)
	{
    	if (PIMC_SIM) MCRotTau = MCBeta/(double)NumbRotTimes;
    	else          MCRotTau = MCBeta/((double)NumbRotTimes-1.0);
	}

	RotRatio  = 1;  // div_t quot - it's important for the area estimator
	// even without rotations

	if (ROTATION)
	{
		int rt;
		if (PIMC_SIM) {
			RotRatio = NumbTimes/NumbRotTimes;  // div_t quot
			rt   = NumbTimes%NumbRotTimes;  // div_t rem
		} 
		else {
			RotRatio = (NumbTimes-1)/(NumbRotTimes-1);  // div_t quot
			rt   = (NumbTimes-1)%(NumbRotTimes-1);  // div_t rem
		}
#ifndef ROTS_TEST
		if (rt)
		nrerror (_proc_,"NumbTimes is not proportional to NumbRotTimes");
#endif
	}

	if (TRANSLATION)
	{
		for (int type=0;type<NumbTypes;type++)  
		{
			MCAtom[type].twave2 = 4.0*MCAtom[type].lambda * MCTau;   // thermal wavelength squared
			MCAtom[type].mlsegm = (int)pow(2.0,MCAtom[type].levels); // segmen size for multilevel
      
			if (MCAtom[type].mlsegm >= NumbTimes)
			nrerror (_proc_,"Segment size is larger then a number of time slices");
		}
	} 

	int bcount = 0;  // number of bosons' types
	int fcount = 0;  // number of fermions' types
	int icount = 0;  // number of molecules (impurities)

	BOSONS   = false;
	FERMIONS = false;

	for (int type=0;type<NumbTypes;type++)
	{
		if (MCAtom[type].stat == BOSE)
		{
			BOSONS = true;
			BSTYPE = type;
			bcount ++;  
		} 

		if (MCAtom[type].stat == FERMI)
		{
			FERMIONS = true;
			FERMTYPE = type;
			fcount ++;  
		} 

		if ((MCAtom[type].molecule == 1)||(MCAtom[type].molecule == 2)||(MCAtom[type].molecule == 3)||(MCAtom[type].molecule == 4))//Modified by Tapas Sahoo
		{
			IMTYPE = type;
			icount ++;  
		} 
    }

  if (bcount>1)
    nrerror (_proc_,"Too many boson atoms' types");

  if (fcount>1)
    nrerror (_proc_,"Too many fermion atoms' types");

  if (icount>1)
    nrerror (_proc_,"Too many dopant molecule' types");

  if ((icount == 0) && (IMPURITY))
    nrerror (_proc_,"Impurity type not defined");

  if ((BOSONS && FERMIONS) && (BSTYPE == FERMTYPE))
    nrerror (_proc_,"Wrong particle statistics");

  //   if (!WORM && BOSONS)
  //   nrerror (_proc_,"BE statistics: worm algorithm only");

/*
#ifndef HOSC_TEST 
	if (BOSONS && (IMTYPE < 0))  // need to define the reference axes for the area
	nrerror (_proc_,"Define the impurity type for the area estimator");
#endif 
*/
}

void MCConfigInit(void)
{
	const char *_proc_=__func__;    // "MCConfigInit()";
#ifndef HOSC_TEST
#ifdef CHAINCONFIG
	initChain_config(MCCoords); 
#endif
#ifdef LATTICECONFIG
	initLattice_config(MCCoords); 
#endif
#ifdef WATERCLUSTER
	initCluster_config(MCCoords,MCAngles); 
#endif
	replInitial_config(MCCoords);

	if (NDIM != 3) nrerror(_proc_,"Only 3D for rotational coordinates");
#else
	for (int id=0;id<NDIM;id++)	
	    for (int atom=0;atom<NumbAtoms;atom++)
    		for (int it=0;it<NumbTimes;it++)
				MCCoords[id][atom*NumbTimes+it] = 0.0;	  
#endif
	for (int it=0;it<(NumbAtoms*NumbTimes);it++)
    {
#ifndef WATERCLUSTER
		MCAngles[PHI][it] = 0.0;
		MCAngles[CTH][it] = 1.0;
		MCAngles[CHI][it] = 0.0;
#endif
		double phi  = MCAngles[PHI][it];
		double cost = MCAngles[CTH][it];
		double chi = MCAngles[CHI][it];
		double sint = sqrt(1.0 - cost*cost);
  
		MCCosine[AXIS_X][it] = sint*cos(phi);
		MCCosine[AXIS_Y][it] = sint*sin(phi);
		MCCosine[AXIS_Z][it] = cost;
		/*
		MCCosinex[AXIS_X][it] = cost*cos(phi)*cos(chi)-sin(phi)*sin(chi);
		MCCosinex[AXIS_Y][it] = cost*sin(phi)*cos(chi)+cos(phi)*sin(chi);
		MCCosinex[AXIS_Z][it] = -sint*cos(chi);
		MCCosiney[AXIS_X][it] = -cost*cos(phi)*sin(chi)-sin(phi)*cos(chi);
		MCCosiney[AXIS_Y][it] = -cost*sin(phi)*sin(chi)+cos(phi)*cos(chi);
		MCCosiney[AXIS_Z][it] = sint*sin(chi);
		*/
	}
}

void initLattice_config(double **pos)
// treat atoms and molecules separateley
// generate a cubic lattice if NumbAtoms = m^3, m - integer	
// nslices = Number of time slices	
{
	//cout<<"in initLattice"<<endl;
	const char *_proc_ = __func__;    // "initLattice_config";

	int natoms = 0;                   // number of atoms 
	int nmolcs = 0;                   // number of molecules

	for (int type=0;type<NumbTypes;type++)
    if ((MCAtom[type].molecule == 1) || (MCAtom[type].molecule == 2)|| (MCAtom[type].molecule == 3)||(MCAtom[type].molecule == 4)) nmolcs += MCAtom[type].numb;//Modified by Tapas Sahoo 
    else natoms += MCAtom[type].numb; 

	// ----- INITIAL CONFIGURATION FOR ATOMS ----------------------

#ifdef IOWRITE
	// box size per particle for atoms only: 
	double abox = BoxSize/pow((double)natoms,1.0/(double)NDIM); 
	double shift[NDIM];

	for (int id=0; id<NDIM; id++) shift[id] = 0.5*abox;  
    
	for (int type=0; type<NumbTypes; type++)   // count molecules only
    if (MCAtom[type].molecule == 0)
    {
		int offset = MCAtom[type].offset;         
		int maxnum = offset + MCAtom[type].numb*NumbTimes;         
      
		for (int atom=offset;atom<maxnum;atom+=NumbTimes)
		{
			for (int id=0;id<(NDIM-1);id++) 
				if (shift[id] > BoxSize) 
		      	{
					shift[id]    = 0.5*abox; 
					shift[id+1] += abox;
				}
 
			for (int id=0;id<NDIM;id++)   // set the center of the box at the origin
				pos[id][atom] = shift[id] - 0.5*BoxSize;

			shift[0] += abox;
		} 
	}

	// ----- INITIAL CONFIGURATION FOR MOLECULES -------------------

	// box size per particle for molecules only: 

	abox = BoxSize/pow((double)nmolcs,1.0/(double)NDIM); 

	cout<<"abox="<<abox<<" "<<BoxSize<<endl;
  
	for (int id=0;id<NDIM;id++)
	shift[id] = 0.5*abox;  
    
	for (int type=0;type<NumbTypes;type++)   // count molecules only
    if ((MCAtom[type].molecule == 1) || (MCAtom[type].molecule == 2)|| (MCAtom[type].molecule == 3)||(MCAtom[type].molecule == 4 ))//Modified by Tapas Sahoo
	{
		int offset = MCAtom[type].offset;         
		int maxnum = offset + MCAtom[type].numb*NumbTimes;         

		for (int atom=offset;atom<maxnum;atom+=NumbTimes)
		{  
			for (int id=0;id<(NDIM-1);id++) 
		    if (shift[id] > BoxSize) {shift[id] = 0.5*abox; shift[id+1] += abox;}
 
			for (int id=0;id<NDIM;id++) // set the center of the box at the origin
		    {
				pos[id][atom] = shift[id] - 0.5*BoxSize;

				if ((natoms == nmolcs))   // to avoid the overlap between particles
			    pos[id][atom] += BoxSize;

				cout<<atom<<" "<<id<<" "<<pos[id][atom]<<endl;

			}	 
 
			shift[0] += abox;
		}  // END loop over atoms
	}    // END loop over types
#endif
	if (InitMCCoords)
	{
		ifstream inputCoords("LatticeConfig.xyz");
		if (inputCoords.is_open ())
		{
			double LatticeCoords[NDIM][NumbAtoms];
			for (int type=0; type<NumbTypes; type++)   // count molecules only
			{
    			if (MCAtom[type].molecule == 0)
    			{
					int offset = MCAtom[type].offset;         
					int maxnum = offset + MCAtom[type].numb*NumbTimes;         
      
					for (int atom=offset; atom<maxnum; atom+=NumbTimes)
					{
						int gatom  = atom/NumbTimes;
						for (int id=0; id<NDIM; id++) 
						{
							inputCoords>>LatticeCoords[id][gatom];
							pos[id][atom] = LatticeCoords[id][gatom];
						}
						for (int id=0; id<NDIM; id++) inputCoords>>DipoleCoords[id][gatom];
					}	
				}	
			}

			for (int type=0; type<NumbTypes; type++)   // count molecules only
			{
    			if ((MCAtom[type].molecule == 1) || (MCAtom[type].molecule == 2)|| (MCAtom[type].molecule == 3)||(MCAtom[type].molecule == 4 ))
    			{
					int offset = MCAtom[type].offset;         
					int maxnum = offset + MCAtom[type].numb*NumbTimes;         
      
					for (int atom=offset; atom<maxnum; atom+=NumbTimes)
					{
						int gatom  = atom/NumbTimes;
						for (int id=0; id<NDIM; id++) 
						{
							inputCoords>>LatticeCoords[id][gatom];
							pos[id][atom] = LatticeCoords[id][gatom];
						}
						for (int id=0; id<NDIM; id++) inputCoords>>DipoleCoords[id][gatom];
					}	
				}	
			}
		} 
		inputCoords.close();
	}
}

void replInitial_config(double **pos)
// replicate configurations for all time slices
{
    for (int id = 0; id < NDIM; id++)	
    {
        for (int atom = 0; atom < NumbAtoms; atom++)
        {
            for (int it = 0; it < NumbTimes; it++)
            {
	            pos[id][atom*NumbTimes+it] = pos[id][atom*NumbTimes];
            }
        }
    }
}	

void initChain_config(double **pos)
{
	cout<<"in initChain"<<endl;
	const char *_proc_ = __func__;    // "initChain_config";

    double shift[NDIM];
    for (int id = 0; id < NDIM; id++)
    {
	    shift[id] = 0.0;
    }

	int NumbAtoms1;
	if (ENT_SIM) NumbAtoms1 = NumbAtoms/2;
	else NumbAtoms1 = NumbAtoms;

	double LatticeTheta = 0.0; //0.25*M_PI;
	double LatticePhi   = 0.25*M_PI;
    for (int atom = 0; atom < NumbAtoms1; atom++)
    {
        for (int it = 0; it < NumbTimes; it++)
        {
            int ii = it + atom*NumbTimes;
   		    for (int id = 0; id < NDIM; id++)   // set the center of the box at the origin
    	    {
	            pos[id][ii] = shift[id];
	        }
        }

		shift[0] += Distance*sin(LatticeTheta)*cos(LatticePhi);
		shift[1] += Distance*sin(LatticeTheta)*sin(LatticePhi);
		shift[2] += Distance*cos(LatticeTheta);
    }
	if (ENT_SIM)
	{
		for (int id = 0; id < NDIM; id++) shift[id] = 0.0;

		for (int atom = NumbAtoms1; atom < NumbAtoms; atom++)
		{
			for (int it = 0; it < NumbTimes; it++)
			{
				int ii = it + atom*NumbTimes;
				for (int id = 0; id < NDIM; id++)   // set the center of the box at the origin
				{
					pos[id][ii] = shift[id];
				}
			}
			shift[0] += Distance*sin(LatticeTheta)*cos(LatticePhi);
			shift[1] += Distance*sin(LatticeTheta)*sin(LatticePhi);
			shift[2] += Distance*cos(LatticeTheta);
		}

		shift[0] -= Distance*sin(LatticeTheta)*cos(LatticePhi);
		shift[1] -= Distance*sin(LatticeTheta)*sin(LatticePhi);
		shift[2] -= Distance*cos(LatticeTheta);
		for (int atom = 0; atom < NumbAtoms1; atom++)
		{
			for (int it = 0; it < NumbTimes; it++)
			{
				int ii = it + atom*NumbTimes;
				for (int id = 0; id < NDIM; id++)   // set the center of the box at the origin
				{
					pos[id][ii] = shift[id];
				}
			}
			shift[0] -= Distance*sin(LatticeTheta)*cos(LatticePhi);
			shift[1] -= Distance*sin(LatticeTheta)*sin(LatticePhi);
			shift[2] -= Distance*cos(LatticeTheta);
		}
	}	
}

void initCluster_config(double **pos, double **MCAngles)
{
	cout<<"in initCluster_config"<<endl;
	const char *_proc_ = __func__;    // "initCluster_config";

	string fname="initial_euler_angles_and_com.txt";
	cout<<"Initial configurations are read from - "<<endl;
	cout<<fname<<endl;

	ifstream fid(fname.c_str(),ios::in);

	double ** _comcoord = doubleMatrix(NumbAtoms,NDIM);  
	double ** _eulerangle = doubleMatrix(NumbAtoms,NDIM);  

	for (int iatom=0; iatom<NumbAtoms; iatom++) {
		for (int id=0; id<NDIM; id++) {
			fid >> _eulerangle[iatom][id];
			cout<<_eulerangle[iatom][id]<<"   ";
		}
		cout<<endl;
	}

	for (int iatom=0; iatom<NumbAtoms; iatom++) {
		for (int id=0; id<NDIM; id++) {
			fid >> _comcoord[iatom][id];
			cout<<_comcoord[iatom][id]<<"    ";
		}
		cout<<endl;
	}

    fid.close();

	for (int iatom = 0; iatom < NumbAtoms; iatom++)
	{
		for (int it = 0; it < NumbTimes; it++)
		{
			int ii = it + iatom*NumbTimes;
			for (int id = 0; id < NDIM; id++)   // set the center of the box at the origin
			{
				pos[id][ii] = _comcoord[iatom][id];
				MCAngles[id][ii] = _eulerangle[iatom][id];
			}
		}
	}

	free_doubleMatrix(_comcoord);
	free_doubleMatrix(_eulerangle);
}

#ifdef PROPOSED
void proposedGrid(void)
{
    dcost = 2.0/((double)NCOST - 1.0); 
    dphi  = 2.0*M_PI/((double)NPHI - 1.0);
    for (int ict = 0; ict < NCOST; ict++)
	{
 		for (int ip = 0; ip < NPHI; ip++)
		{
			int ii = ip + ict*NPHI;
        	costProposed[ii] = -1.0+ict*dcost;
			phiProposed[ii] = ip*dphi;
		}
	}
}
#endif
void ParamsPotential(void)
{
	if (DipoleMoment)
	{
		double DipoleMomentAU   = DipoleMoment/AuToDebye;
		DipoleMomentAU2= DipoleMomentAU*DipoleMomentAU;
	}
	if (Distance) RR   = Distance/BOHRRADIUS;
#ifdef EWALDSUM
	boxLength          = 20.0/BOHRRADIUS;
	double boxLength2  = boxLength*boxLength;
	double boxLength3  = boxLength2*boxLength;
	alpha              = 1.0*BOHRRADIUS;
	alpha2             = alpha*alpha;
	double alpha3      = alpha2*alpha;
	prefSelf           = 2.0*alpha3/(3.0*sqrt(M_PI));
	prefBfun           = 2.0*alpha/sqrt(M_PI);
	prefUk1            = 2.0*M_PI/boxLength;
	prefUk2            = pow(M_PI/(alpha*boxLength),2);
	prefUk3            = 0.5/boxLength3;
	KMAX               = 20;
#endif
}

//Tapas implemented the below section//
#ifdef GAUSSIANMOVE
void ProposedMCCoords()
{
	int atype = 0;
/*
	double convert         = 1.66054*2.998*pow(10,-3)/1.0546;
	double FrequencyHO     = 78.6; //in cm-1
	double FrequencyK      = FrequencyHO*CMRECIP2KL; //in Kelvin

    double prefacGauss1    = 0.5*MCAtom[atype].mass*FrequencyHO*convert;
	double prefacGauss2    = tanh(FrequencyK*MCTau);
    double prefacGauss3    = sinh(FrequencyK*MCTau);
*/
	double FrequencyK      = AuToKelvin;

    double prefacGauss1    = 0.5;
	double prefacGauss2    = tanh(FrequencyK*MCTau);
    double prefacGauss3    = sinh(FrequencyK*MCTau);

    double afacGauss       = prefacGauss1/prefacGauss2;
    double bfacGauss       = prefacGauss1/prefacGauss3;
	cout<<"afacGauss = "<<afacGauss<<" bfacGauss = "<<bfacGauss<<endl;
	
//Tridiagonal Matrix formation; starts

	int it0, it1, ii0, iip1, iim1;

	for (it0 = 0; it0 < NumbTimes; it0++)
	{
		for (it1 = 0; it1 < NumbTimes; it1++)
		{
			ii0 = it1 + it0*NumbTimes;
			AMat[ii0] = 0.0;
		}
	}

	for (it0 = 0; it0 < NumbTimes; it0++)
	{
		for (it1 = 0; it1 < NumbTimes; it1++)
		{
			if (it1 == it0)
			{
				ii0  = it1 + it0*NumbTimes;
				iip1 = (it1+1) + it0*NumbTimes;
				iim1 = (it1-1) + it0*NumbTimes;
				if ((it1+1) < NumbTimes)
   				{
					AMat[iip1]      = -bfacGauss;
   				}
				if ((it1 == 0) || (it1 == (NumbTimes-1)))
				{
					AMat[ii0]       = afacGauss;
				}
				else
				{
					AMat[ii0]       = 2.0*afacGauss;
				}
   				if ((it1-1) >= 0)
   				{
					AMat[iim1] 		= -bfacGauss;
   				}
			}
		}
	}
	for (it0 = 0; it0 < NumbTimes; it0++)
	{
		for (it1 = 0; it1 < NumbTimes; it1++)
		{
			int ii0 = it1 + it0*NumbTimes;
			cout<< AMat[ii0]<<"  ";
		}
		cout<<endl;
	}
	for (it0 = 0; it0 < NumbTimes; it0++)
	{
		Eigen[it0] = 0.0;
	}
	diag(AMat, Eigen, NumbTimes);
	for (it0 = 0; it0 < NumbTimes; it0++)
    {
      	cout<<Eigen[it0]<<"   ";
    }
}

void GetRandomCoords()
{
	int atom0, it0, it1, id0, ii0, t0;

	double sum;

	//RngStream Rng[1];
   	RngStream Rng[omp_get_num_procs()];     // initialize a parallel RNG named "Rng"
	double mu =0.0;
	double sigma;

	for (atom0 = 0; atom0 < NumbAtoms; atom0++)
	{
		for (id0 = 0; id0 < NDIM; id0++)
		{
			for (it0 = 0; it0 < NumbTimes; it0++)
			{
				sum = 0.0;
   				for (it1 = 0; it1 < NumbTimes; it1++)
   				{
					if (Eigen[it1] < 0.0) 
					{
						cout << "Eigen values must be positive"<<endl;
						exit(0);
					}
					sigma = sqrt(1.0/(2.0*Eigen[it1]));

					ii0 = it0 + it1*NumbTimes;
					sum += AMat[ii0]*rnorm(Rng, mu, sigma);
   				}
				t0 = it0 + atom0*NumbTimes;
				gausscoords[id0][t0]  = sum;//*BOHRRADIUS;
			}
		}
    }
}

void diag(double *a, double *w, int N)
{
	int LDA = N;
    int n = N, lda = LDA, info, lwork;
    double wkopt;
    double *work;
    lwork = -1;
    char u1[] = "Vectors";
    char u2[] = "Upper";
    dsyev_( u1, u2, &n, a, &lda, w, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = new double[lwork];
    dsyev_( u1, u2, &n, a, &lda, w, work, &lwork, &info );
    if( info > 0 )
    {
        cout<<"The algorithm failed to compute eigenvalues."<<endl;
        exit( 1 );
    }
    delete [] work;
} 

void print_matrix( char *desc, int m, int n, double* a, int lda )
{
    int i, j;
    cout<<desc<<endl;
    for( i = 0; i < m; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            cout<<a[i+j*lda]<<"  ";
        }
        cout<<endl;
    }
}
#endif

/*
#ifdef WATERCLUSTER
void GetInitialCoords()
{
	const char *_proc_=__func__;    // "init_rotdens"; 

	string fname="initial_euler_angles_and_com.txt";
	string pathr=PathToDensity;
	fname = pathr+fname;
	cout<<fname<<endl;

	ifstream fid(fname.c_str(),ios::in);

	if (!fid.good())
	_io_error(_proc_,IO_ERR_FOPEN,fname.c_str());

//---- read  grid information -------------------------------------
    fid >> THgrd>>PHgrd;

	cout<<"Rotational density matrix elements are read!"<<endl;
    cout<<"THgrd="<<THgrd<<" PHgrd="<<PHgrd<<endl;

	_cosg_indices[type] = doubleMatrix(THgrd*PHgrd,THgrd*PHgrd);  
	_rotden_indices[type] = doubleMatrix(THgrd*PHgrd,THgrd*PHgrd);  
	_rotderv_indices[type] = doubleMatrix(THgrd*PHgrd,THgrd*PHgrd);  
	_rotesqr_indices[type] = doubleMatrix(THgrd*PHgrd,THgrd*PHgrd);  

	int ii, jj;
	double cc, dd;
	for (int i=0; i<THgrd; i++) {
		for (int j=0; j<PHgrd; j++) {
			ii = j+i*PHgrd;
			for (int k=0; k<THgrd; k++) {
				for (int l=0; l<PHgrd; l++) {
					jj = l+k*PHgrd;
					fid >> _cosg_indices[type][ii][jj]>>_rotden_indices[type][ii][jj]>>_rotderv_indices[type][ii][jj]>>_rotesqr_indices[type][ii][jj]; 
					//cout<<_cosg_indices[type][ii][jj]<<" "<<_rotden_indices[type][ii][jj]<<_rotderv_indices[type][ii][jj]<<_rotesqr_indices[type][ii][jj]<<endl;
				}
			}
		}
	}

    fid.close();

	double cost,phi;
	double cstep=2.0/(double)(THgrd-1.0);
    double pstep=2.0*M_PI/(double)(PHgrd-1.0);
	cost_indices = new double [THgrd*PHgrd];
	phi_indices = new double [THgrd*PHgrd];

	for (int i=0; i<THgrd; i++) {
		cost = i*cstep-1.0;
		for (int j=0; j<PHgrd; j++) {
			phi = j*pstep;
			ii=j+i*PHgrd;
			cost_indices[ii] = cost;
			phi_indices[ii] = phi;
		}
	}

	MCIndices = new int [NumbTimes*NumbAtoms]; 
	for (int it=0; it<NumbTimes*NumbAtoms; it++) {
		MCIndices[it]=0;
	} 
}
*/
