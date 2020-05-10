//
//         main()
//

#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
//#include <stdlib>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "mc_confg.h"
#include "mc_setup.h"
#include "mc_input.h"
#include "mc_randg.h"
#include "mc_utils.h"
#include "mc_poten.h"
#include "mc_piqmc.h"
#include "mc_qworm.h"
#include "mc_estim.h"
#include "mc_const.h"

#include <stdio.h>
//#include <mpi.h>
#include <omp.h>
#include "rngstream.h"
#include "omprng.h"

using namespace std;

void MCWormAverage(void);
void MCWormAverageReset(void);

double avergCount;   // # of calls of get_estim inside a block
double avergCountENT;   // # of calls of get_estim inside a block
double totalCountENT;

void PIMCPass(int,int);

void MCGetAverage(void);
void MCResetBlockAverage(void);
void MCGetAveragePIMC(void);
void MCResetBlockAveragePIMC(void);
void MCGetAveragePIGS(void);
void MCResetBlockAveragePIGS(void);
void MCGetAveragePIGSENT(int);
void MCResetBlockAveragePIGSENT(void);

void MCSaveBlockAverages(long int);

void MCSaveAcceptRatio(long int,long int,long int);

#ifdef DEBUG_POT
extern "C" void caleng_(double *com_1, double *com_2, double *E_2H2O, double *Eulang_1, double *Eulang_2);
#endif

//--------------- BLOCK AVERAGE ------------

double _dbpot;       // potential energy differencies, block average  added by Hui Li
double _bpot;       // kinetic   energy, block average
double _btotal;
double _bkin;       // potential energy, block average
double _bCv;        // heat capacity, block average
double _bCv_trans;  // translational heat capacity, block average
double _bCv_rot;   //  rotational heat capacity, block average

double _dpot_total;  // potential energy differences, global average  added by Hui Li
double _bcostheta;
#
double _ucompx;
double _ucompy;
double _ucompz;
#
double _abs_ucompx;
double _abs_ucompy;
double _abs_ucompz;
#
vector<double> _cdipoleXYZ;
vector<double> _cdipoleX;
vector<double> _cdipoleY;
vector<double> _cdipoleZ;
vector<double> _cdipoleXY;
double _beiej[4];
double _bei[3];
//
double _bnm;
double _bdm;
double _trOfDensitySq;

double _brot;       // rotational kin energy, block average
double _brot1;       // rotational kin energy, block average
double _brotsq;     // rotational energy square, block average
double _rotsq_total; // rotational energy square, global average
double _Cv_total;    // heat capacity, global average
double _Cv_trans_total;    // translational heat capacity, global average
double _Cv_trans_1_total;    // translational heat capacity, global average
double _Cv_trans_2_total;    // translational heat capacity, global average
double _Cv_rot_total;    // rotational heat capacity, global average

fstream _feng;      // save accumulated energy
fstream _fdc;      // save accumulated energy
fstream _fangins;      // save accumulated energy
fstream _fengins;      // save accumulated energy
fstream _fentropy;
fstream _fswap;

//---------------- ESTIMATORS ---------------

void SaveEnergy    (const char [],double,long int); // block average
void SaveSumEnergy (double,double);                 // accumulated average

void SaveInstantEnergy ();                 // accumulated average

#ifdef DDCORR
void SaveDipoleCorr(const char [],double,long int);
#endif
#ifdef ORDERPARA
void SaveOrderCorr(const char fname [], double acount, long int blocknumb);
#endif
void SaveAngularDOF(const char [],double,long int);
void SaveTrReducedDens(const char [], double , long int );

void InitTotalAverage(void);
void DoneTotalAverage(void);
void SaveInstantConfig(const char [],long int); //Saving instataneous configurations

//---------------- INITIAL MCCOORDS -----------

extern "C" void initconf_(double *MCCooInit, double *MCAngInit, int *Ntotal, int *NBoson, int *PIndex,int *Rindex);

//---------------- PRINT PERMUTATION TABLE -----

extern "C" void prtper_(int *PIndex,int *NBoson,long int *blockCount);

//----------- subroutine that performs some initialization work for vh2h2.f-----

extern "C" void vinit_();
double MCAccepSwap;
double MCAccepUnSwap;
#ifdef PROPOSED
int iChooseOld = 0;
int iChooseNew;
int iChoose;
#endif

//-------------------------------------------

int main(int argc, char *argv[])
{
Distribution = "unSwap";
#ifdef PAIRDENSITY
	readPairDensity();
	cout<<"TAPAS"<<endl;
	exit(11);
#endif
   	randomseed(); //set seed according to clock
	if (ENT_ENSMBL == EXTENDED_ENSMBL) srand(time(NULL));
#ifdef POTH2
   vinit_();
#endif
//int numprocs, rank, namelen;
//char processor_name[MPI_MAX_PROCESSOR_NAME];


//MPI_Init(&argc, &argv);
//MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//MPI_Get_processor_name(processor_name, &namelen);
//MPI_Status mpistatus;

 
//-------------------------
	MPIsize = 1;            // this is for InitRandom()          
	MPIrank = MPI_MASTER;   // this is for InitRandom()          
//-----------------------------

	int restart = 0;
	IOReadParams(FINPUT,restart); // read input parameters 
// get the number of rotational steps treated by worker CPUs
// chunksize = NumbRotTimes / numprocs;

// pass numprocs to NProcs, which is a global variable
// NProcs = numprocs;

	MCInitParams();               // set system dependent parameters
	MCSetUnits(); 
	MCMemAlloc();

	MemAllocMCCounts();
	MemAllocQWCounts();

//head CPU comes in after memory declaration.
//if ( rank == MPI_MASTER ){// in the master cpu section

// only internal system of units after this point 
 
	MCInit();
#ifdef PROPOSED
proposedGrid();
#endif
ParamsPotential();
#ifdef GAUSSIANMOVE
	ProposedMCCoords();
#endif

	if (WORM)
	MCWormInit();

// clean the permutation sampling from the last run
   	if (FileExist(FPERMU))
   	_io_error("QMC ->",IO_ERR_FEXST,FPERMU);

   	if (!restart) // new run, generate new status, rnd() streams and config files     
   	{
// check the existence of the old files first
  
      	if (FileExist(FSTATUS))    
     	_io_error("QMC ->",IO_ERR_FEXST,FSTATUS);
       
      	if (FileExist(FCONFIG)) 
     	_io_error("QMC ->",IO_ERR_FEXST,FCONFIG);

      	if (FileExist(FRANDOM)) 
     	_io_error("QMC ->",IO_ERR_FEXST,FRANDOM);

      	if (FileExist(FSEED)) 
     	_io_error("QMC ->",IO_ERR_FEXST,FSEED);

//--------------------------------------------------------
//    	generate tables - potentals, configurations etc, status

	   	MCStartBlock = 0;  // PIMCRESTART //
//    	SEED for head CPU
      	SEED = 985456376;
      	RandomInit(MPIrank,MPIsize);

    	MCConfigInit();                // generate initial configurations
		
    	for (int it=0;it<NumbAtoms*NumbTimes;it++)
    	{
        	for (int id=0;id<NDIM;id++)
        	{
            	cout<<"it " << it<<" id "<< id<< " "<< MCCoords[id][it]<<endl;
        	}
    	} 
    	cout<<"  "<<endl;
    	cout<<"  "<<endl;
		for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
		{
        	int offset0 = NumbTimes*atom0;
			for (int atom1 = 0; atom1 < NumbAtoms; atom1++)
			{
				if (atom0 != atom1)
				{
            		int offset1 = NumbTimes*atom1;

            		int it = ((NumbRotTimes - 1)/2);
            		int t0 = offset0 + it;
            		int t1 = offset1 + it;

            		double dr;
            		double dr2 = 0.0;
            		for (int id=0;id<NDIM;id++)
            		{
                		dr   = (MCCoords[id][t0] - MCCoords[id][t1]);
                		dr2 += (dr*dr);
            		}
            		double r = sqrt(dr2);
					cout<<"atom0 "<<atom0<<" atom1 "<<atom1<<" RCOM "<<r<<endl;
				}
			}
		}
/*
    	cout<<"  "<<endl;
    	cout<<"  "<<endl;
    	for (int it=0;it<NumbAtoms*NumbTimes;it++)
    	{
        	for (int id=0;id<NDIM;id++)
        	{
            	cout<<"it " << it<<" id "<< id<< " "<< MCCosine[id][it]<<endl;
        	}
    	}	 
*/
//    	read in initial MCCoords and MCAngles
#ifdef IOWRITE
      	if(InitMCCoords)
      	{
         	cout<<"read in MCCoords and MCAngles from xyz.init file"<<endl;
         	int ntotal=NumbAtoms*NumbTimes;
         	int iperm,nboson;
         	nboson=0;
         	if(BOSONS)
         	{
            	nboson=MCAtom[BSTYPE].numb;
         	}
//       	initconf_(MCCooInit,MCAngInit,&ntotal,&MCAtom[BSTYPE].numb,PIndex,RIndex);
         	initconf_(MCCooInit,MCAngInit,&ntotal,&nboson,PIndex,RIndex);
         	if(BOSONS)
         	{
            	cout<<"read PIndex: "<<" ";
            	for(int id=0;id<MCAtom[BSTYPE].numb;id++)
            	cout<<PIndex[id]<<" ";

            	cout<<endl;

            	cout<<"read RIndex: "<<" ";
            	for(int id=0;id<MCAtom[BSTYPE].numb;id++)
            	cout<<RIndex[id]<<" ";

            	cout<<endl;
         	}

         	for (int it=0;it<NumbAtoms*NumbTimes;it++)
         	{
            	for(int id=0;id<NDIM;id++)
            	{
               		MCCoords[id][it]=MCCooInit[it*NDIM+id];
               		MCAngles[id][it]=MCAngInit[it*NDIM+id];
            	}
         	}
      	}
#endif
		delete [] MCCooInit;
		delete [] MCAngInit;

    	//string fname = MCFileName + IO_EXT_XYZ;
    	//IOxyz(IOWrite,fname.c_str());  // test output of initial config
      	//string fname = MCFileName;
      	//IOxyzAng(IOWrite,fname.c_str()); // test output of initial config
	
#ifdef DEBUG_POT
		//Testing the potential
		cout<<"Printing of values related to the caleng potential!"<<endl;
		int it=0;
		int atom0=0;
		double Eulang1[NDIM];

		int offset0 = atom0*NumbRotTimes;
		int tm0=offset0 + it/RotRatio;

		Eulang1[PHI]=1.0;
		Eulang1[CTH]=1.0;
		Eulang1[CHI]=1.0;

		cout<<" Eulang1[PHI] "<<Eulang1[PHI];
		cout<<" Eulang1[CTH] "<<Eulang1[CTH];
		cout<<" Eulang1[CHI] "<<Eulang1[CHI];
		cout<<endl;

		int atom1=1;
		double spot=0.0;
		if (atom1 != atom0)                    // skip "self-interaction"
		{	
			int offset1 = atom1*NumbRotTimes;
			if (MCType[atom1] == IMTYPE)
			{
				int t0 = offset0 + it;
				int t1 = offset1 + it;
				double com1[NDIM];
				double com2[NDIM];
				double Eulang2[NDIM];
				double E_2H2O;
				for (int id=0; id<NDIM; id++)
				{
					com1[id] = MCCoords[id][t0];
					com2[id] = MCCoords[id][t1];
					cout<<" com1 "<<com1[id]<<" com2 "<<com2[id]<<endl;
				}
				Eulang2[PHI]=1.0;
				Eulang2[CTH]=-1.0;
				Eulang2[CHI]=1.0;

				cout<<" Eulang2[PHI] "<<Eulang2[PHI];
				cout<<" Eulang2[CTH] "<<Eulang2[CTH];
				cout<<" Eulang2[CHI] "<<Eulang2[CHI];
				cout<<endl;
				caleng_(com1, com2, &E_2H2O, Eulang1, Eulang2);
				spot += E_2H2O;
			}
		}   // END sum over atoms
		cout<<" Caleng value "<<spot<<endl;
#endif

//--------------------------------------------------------

//PIMCRESTART begins here//
      	StatusIO(IOWrite,FSTATUS);
      	ConfigIO(IOWrite,FCONFIG);
      	TablesIO(IOWrite,FTABLES);
      	RandomIO(IOWrite,FRANDOM);
		SeedIO(IOWrite,FSEED); 
//PIMCRESTART ends here//

      	if (WORM)
      	QWormsIO(IOWrite,FQWORMS);
   	}

   	MCWormAverageReset();      // debug worm

   	InitTotalAverage();      // DUMP 

/*
#ifdef NORATIOTRICK
    ifstream fid("readConf.xyz");

	string sbuff;
    //string name;
	int offset;
	int type = 0;
	int atom = 0;  // first atom # will be 1, NOT 0 

	//fid>>MaxnTimes;
	//getline (fid,name);

	for (int type=0;type<NumbTypes;type++)
	for (int atom=0;atom<MCAtom[type].numb;atom++)
	for (int it=0;it<NumbTimes;it++)
	{
		fid>>sbuff;         // skip an atom type

		offset=MCAtom[type].offset+NumbTimes*atom;
		for (int id=0;id<NDIM;id++)
		{
			fid>>MCCoords[id][offset+it];
			fid>>MCAngles[id][offset+it];
		}
	}
    fid.close();
	//string fnamer = "readConf";
	//IOxyzAng(IORead,fnamer.c_str());
#endif
*/
// --- RESTART/START NEW RUN ----------------------------
//PIMCRESTART begins here//
   	if (restart) // new run, generate new status, rnd() streams and config files     
   	{
		StatusIO(IORead,FSTATUS);  // load MCStartBlock
   		ConfigIO(IORead,FCONFIG);  // load atoms/molecules positions
   		TablesIO(IORead,FTABLES);  // load permutation tables
   		RandomIO(IORead,FRANDOM);  // load rnd streams 
		SeedIO(IORead,FSEED);
      	string fname = MCFileName;
      	IOxyzAng(IOWrite,fname.c_str()); // test output of initial config
	}
//PIMCRESTART ends here//
    for (int it=0;it<NumbAtoms*NumbTimes;it++)
    {
        double phi  = MCAngles[PHI][it];
        double cost = MCAngles[CTH][it];
        double sint = sqrt(1.0 - cost*cost);
        double chi  = MCAngles[CHI][it];

        MCCosine[AXIS_X][it] = sint*cos(phi);
        MCCosine[AXIS_Y][it] = sint*sin(phi);
        MCCosine[AXIS_Z][it] = cost;
    }

	string fnamew = "writeConf";
	IOxyzAng(IOWrite,fnamew.c_str());

   	if (WORM)
   	QWormsIO(IORead,FQWORMS);

// BEGIN loop over blocks ------------------------

   	InitPotentials(); 

   	if (ROTATION)
    	InitRotDensity();

   	if((IREFLX == 1 || IREFLY == 1 || IREFLZ == 1) && MCAtom[IMTYPE].molecule != 2)
   	{
      	cout<<IREFLX<<" "<<IREFLY<<" "<<IREFLZ<<" "<<MCAtom[IMTYPE].molecule<<endl;
      	nrerror("main","SYMMETRIZATION IS ONLY IMPLEMENTED FOR TOPS");
   	}
// MPI BLOCK 1 in MASTER start
// send run information to all the worker CPUs
/*
   for(int iproc=1;iproc<numprocs;iproc++)
   {
      MPI_Send(&NumbTypes,1,MPI_INT,iproc,tagNumbTypes,MPI_COMM_WORLD);
      MPI_Send(&NumbAtoms,1,MPI_INT,iproc,tagNumbAtoms,MPI_COMM_WORLD);
      MPI_Send(&rhoprp[0],SizeRotDen,MPI_DOUBLE,iproc,tagrho,MPI_COMM_WORLD);
      MPI_Send(&erotpr[0],SizeRotDen,MPI_DOUBLE,iproc,tagero,MPI_COMM_WORLD);
      MPI_Send(&vtable[0],SizePotTab,MPI_DOUBLE,iproc,tagpot,MPI_COMM_WORLD);
      MPI_Send(&WORM,1,MPI_INT,iproc,tagWORM,MPI_COMM_WORLD);
      MPI_Send(&MCType[0],NumbAtoms,MPI_INT,iproc,tagMCType,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].molecule,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_molecule,MPI_COMM_WORLD);
      MPI_Send(&RotRatio,1,MPI_INT,iproc,tagRotRatio,MPI_COMM_WORLD);
      MPI_Send(&Worm.type,1,MPI_INT,iproc,tagWormtype,MPI_COMM_WORLD);
      MPI_Send(&NumbRotTimes,1,MPI_INT,iproc,tagNumbRotTimes,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].offset,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_offset,MPI_COMM_WORLD);
      MPI_Send(&IMTYPE,1,MPI_INT,iproc,tagIMTYPE,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].numb,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_numb,MPI_COMM_WORLD);
      MPI_Send(&MCTau,1,MPI_DOUBLE,iproc,tagMCTau,MPI_COMM_WORLD);
      MPI_Send(&ROTATION,1,MPI_INT,iproc,tagROTATION,MPI_COMM_WORLD);
      MPI_Send(&MCStartBlock,1,MPI_LONG,iproc,tagMCStartBlock,MPI_COMM_WORLD);
      MPI_Send(&NumberOfMCBlocks,1,MPI_LONG,iproc,tagNumberOfMCBlocks,MPI_COMM_WORLD);
      MPI_Send(&NumberOfMCPasses,1,MPI_LONG,iproc,tagNumberOfMCPasses,MPI_COMM_WORLD);
//    MPI_Send(&TZMAT[1][1],NDIM,MPI_DOUBLE,iproc,tagTZMAT,MPI_COMM_WORLD);
      MPI_Send(&RotDenType,1,MPI_INT,iproc,tagRotDenType,MPI_COMM_WORLD);
      MPI_Send(&RotOdEvn,1,MPI_INT,iproc,tagRotOdEvn,MPI_COMM_WORLD);
      MPI_Send(&RotEoff,1,MPI_DOUBLE,iproc,tagRotEoff,MPI_COMM_WORLD);
      MPI_Send(&X_Rot,1,MPI_DOUBLE,iproc,tagARot,MPI_COMM_WORLD);
      MPI_Send(&Y_Rot,1,MPI_DOUBLE,iproc,tagBRot,MPI_COMM_WORLD);
      MPI_Send(&Z_Rot,1,MPI_DOUBLE,iproc,tagCRot,MPI_COMM_WORLD);
      MPI_Send(&MCRotTau,1,MPI_DOUBLE,iproc,tagMCRotTau,MPI_COMM_WORLD);
   }
// tagrunning = (tagNumberOfMCPasses + 2) % tagUpper;
// tagrunning = (tagTZMAT + 2) % tagUpper;
   tagrunning = (tagMCRotTau + 2) % tagUpper;

// the following loop passes the same tagrunning to all workers, in order to guarantee all CPUs have the same running tag
   for(int iproc=1;iproc<numprocs;iproc++)
   {
      MPI_Send(&tagrunning,1,MPI_INT,iproc,tagTAGRUN,MPI_COMM_WORLD);
   }
// MPI BLOCK 1 in MASTER done
*/
    InitMCEstims();
   	//InitTotalAverage();      // DUMP 
   
   	ResetMCCounts();
   	ResetQWCounts();

/*
// 	try openmp for loop
   	int chunksize,nthrds;
   	#pragma omp parallel
   	{
   	int tid=omp_get_thread_num();
   	nthrds=omp_get_num_threads();
   	chunksize=NumbRotTimes/nthrds;
   	cout<<"chunksize="<<chunksize<<endl;
   int itini=chunksize*tid;
   int itfnl=itini+chunksize;
// #pragma omp for
   for (int itrot=itini;itrot<itfnl-1;itrot++)
   {
      cout<<"itrot="<<itrot<<" tid="<<tid<<endl;
   }
   }

   for (int itrot=chunksize-1;itrot<NumbRotTimes;itrot=itrot+chunksize)
   {
      cout<<"itrot="<<itrot<<endl;
   }
   for (int itrot=nthrds*chunksize;itrot<NumbRotTimes;itrot++)
   {
      cout<<"itrot="<<itrot<<endl;
   }
*/

// get the chunksize and # of threads for the parallel rotation
	#pragma omp parallel
	{
    	NThreads=omp_get_num_threads();
      	chunksize=NumbRotTimes/NThreads;
   	} // end omp

   	cout<<"NThreads="<<NThreads<<" chunksize="<<chunksize<<endl;

   	printf("OpenMP version: %d\n", _OPENMP);
   	//randomseed(); //set seed according to clock
	//RngStream Rng[omp_get_num_procs()];     // initialize a parallel RNG named "Rng"

	//Monte Carlo simulation begins here//
   	long int blockCount = MCStartBlock;  
   	while (blockCount<NumberOfMCBlocks) // START NEW BLOCK      
   	{      
    	blockCount++; 
       	if (PIMC_SIM) MCResetBlockAveragePIMC();
       	if (PIGS_SIM) MCResetBlockAveragePIGS();
       	if (ENT_SIM) MCResetBlockAveragePIGSENT();
     
       	long int passCount = 0;        // BEGIN NEW MC PASS
       	long int passTotal = 0;        // total number of Moves = passCount*time 

       	while (passCount++ < NumberOfMCPasses) 
       	for (int time=0; time<NumbTimes; time++)
       	{
        	passTotal++;  

           	for (int type=0;type<NumbTypes;type++)
           	if (WORM && (type == Worm.type)) 
           	{

           		MCWormMove();

           		if (!Worm.exists)
           		{
					/* reactive */
               		int rt = nrnd2(NumbTimes);   // select a time slice

               		if ((type == BSTYPE) || (type == FERMTYPE))
               			MCBisectionMoveExchange (type,rt);
                 	else
                 		MCBisectionMove (type,rt);
						/* reactive */
#ifdef MOVECENTROIDTEST
               		if  (time == 0)
           			if ((type == BSTYPE) || (type == FERMTYPE))
           				MCMolecularMoveExchange(type);        
               		else
               			MCMolecularMove(type);        
#endif
           		}
        	} 
       		else
			{
		    	PIMCPass(type,time);
			}

   			if (blockCount>NumberOfEQBlocks)        // skip equilibration steps
   			{
				// evaluate averages 
           		if (passTotal % MCSKIP_AVERG == 0)   // skip correlated configurations
           		{
               		if (WORM)
               		{
                		if (!Worm.exists)
                		{
                			//MCGetAverage();
	                        //prtper_(PIndex,&MCAtom[BSTYPE].numb,&blockCount);

                   			// print the instantaneous xyz and prl info for the closed path
                   			if(PrintXYZprl)
                   			{
                       			stringstream bc;                // convert block # to string
                       			bc.width(IO_BLOCKNUMB_WIDTH);
                       			bc.fill('0');
                       			bc<<blockCount;
                       			string fname = MCFileName + bc.str();  // file name prefix including block #
                       			IOxyzAng(IOWrite,fname.c_str());
                       			PrintXYZprl = 0;
                   			}

                   		}
                   		else  
                		MCWormAverage();  /// TEST ONLY
                	}	   
                	else
                 	{
                    	//MCGetAverage();
                    	if (PIMC_SIM) MCGetAveragePIMC();
                    	if (PIGS_SIM) MCGetAveragePIGS();
						if (ENT_ENSMBL == BROKENPATH) MCGetAveragePIGSENT(IMTYPE);
					    //omp_set_num_threads(1);
					    //MCGetAverage();
					    //omp_set_num_threads(NThreads);
                        //print the instantaneous xyz and prl info for the closed path
                    	if(PrintXYZprl)
                    	{
							stringstream bc;                // convert block # to string
							bc.width(IO_BLOCKNUMB_WIDTH);
							bc.fill('0');
							bc<<blockCount;
							string fname = MCFileName + bc.str();  // file name prefix including block #
							IOxyzAng(IOWrite,fname.c_str());
							PrintXYZprl = 0;
                    	}
                 	}
              	}
				if (ENT_ENSMBL == EXTENDED_ENSMBL) MCGetAveragePIGSENT(IMTYPE);

// DUMP, save global average. if avergCount = 0, then MCGetAverage is never called in this block.  All Saving steps are skipped.

              	if (passTotal % MCSKIP_TOTAL == 0 && avergCount)  
              	{
               		sumsCount += 1.0;                 
					//if (ENT_SIM) //SaveSumEnergy (totalCountENT,sumsCount);
					if (!ENT_SIM) SaveSumEnergy (totalCount,sumsCount);
            	}
			}  
          
        	if (passTotal % MCSKIP_RATIO == 0)
        	MCSaveAcceptRatio(passTotal,passCount,blockCount);
		}                               
       
		// END  loop over MC passes (time slices) -----------------

		if (blockCount>NumberOfEQBlocks && avergCount)   // skip equilibration steps
		{
   			MCSaveBlockAverages(blockCount);
			// save accumulated interatomic distribution

#ifdef HISTOGRAM
			if (blockCount > 1500)
			SaveGxyzSum(MCFileName.c_str(),totalCount);
#endif
#ifdef NEWDENSITY
			if (blockCount > 600)
			SaveGraSum(MCFileName.c_str(),totalCount);
#endif

#ifdef IOWRITE
			if(IMPURITY && MCAtom[IMTYPE].molecule == 1)
	    	SaveDensities2D(MCFileName.c_str(),totalCount,MC_TOTAL);

	   		if(IMPURITY && MCAtom[IMTYPE].molecule == 3)
	     	SaveDensities2D(MCFileName.c_str(),totalCount,MC_TOTAL);

           	if(IMPURITY && MCAtom[IMTYPE].molecule == 2)
           	{
            	SaveDensities3D(MCFileName.c_str(),totalCount,MC_TOTAL); // this step takes lots of space.  temporarily turnned off
              	SaveRho1D(MCFileName.c_str(),totalCount,MC_TOTAL);
				// SaveRhoThetaChi(MCFileName.c_str(),totalCount,MC_TOTAL); // we don't need 2d angular distribution comparison now
           	}

#endif
			if (ROTATION && PIMC_SIM)                  // DUMP  accumulated average
				SaveRCF(MCFileName.c_str(),totalCount,MC_TOTAL);
		}

		SaveInstantConfig(MCFileName.c_str(),blockCount);

//PIMCRESTART begins // 
		//  CHECKPOINT: save status, rnd streams and configs ------
	// The below segment will save the data at each 200 blocks interval. One may change the interval by changing blockCount%200 with blockCount%any number//

		MCStartBlock = blockCount; 
		if (blockCount % 200 == 0)
		{
			IOFileBackUp(FSTATUS); StatusIO(IOWrite,FSTATUS);
			IOFileBackUp(FCONFIG); ConfigIO(IOWrite,FCONFIG);
			IOFileBackUp(FTABLES); TablesIO(IOWrite,FTABLES);
			IOFileBackUp(FRANDOM); RandomIO(IOWrite,FRANDOM);      
			IOFileBackUp(FSEED); SeedIO(IOWrite,FSEED);
		}

       	if (WORM)
       	{
        	IOFileBackUp(FQWORMS);
        	QWormsIO(IOWrite,FQWORMS);
    	}
	}  
	// END  loop over blocks ------------------------
	// } // master cpu if block ends

	//MPI_Finalize();

	DoneTotalAverage();  // DUMP

	if (WORM)
		MCWormDone();

	DoneMCEstims();
	DonePotentials();

	if (ROTATION)
		DoneRotDensity();

	MCMemFree();
	RandomFree();

	MFreeMCCounts();
	MFreeQWCounts();

	return 1; 
}

void PIMCPass(int type,int time)
{
  // skip solvent and translation moves for rotations only
	if(TRANSLATION)
	{
#ifdef GAUSSIANMOVE
		MCMolecularMoveGauss(type);
#else
		if (PIGS_SIM || ENT_SIM)
		{	
			if (time == 0) { MCMolecularMovePIGS(type); }        
			MCBisectionMovePIGS(type,time);
			MCMolecularMoveNaiveNorm(type);
			//MCMolecularMoveNaive(type);
		}

		if (PIMC_SIM)
		{	
			if (time == 0) MCMolecularMove(type);        
// move the solvent particles
			MCBisectionMove(type,time);
		}	
#endif
	}

	if ((type == IMTYPE) && ROTATION && (MCAtom[type].molecule == 1))  // rotational degrees of freedom
    	MCRotationsMove(type);
   	if ((type == IMTYPE) && ROTATION && (MCAtom[type].molecule == 3))  // rotational degrees of freedom
    	MCRotationsMove(type);
	if ((type == IMTYPE) && ROTATION && (MCAtom[type].molecule == 2))  // non-linear rotor rotation added by Toby
    	MCRotations3D(type);
	if ((type == IMTYPE) && ROTATION && (MCAtom[type].molecule == 4))
    	MCRotationsMove(type);
}

void MCResetBlockAveragePIMC(void) 
{
	avergCount = 0.0;

	ResetMCEstims();

	ResetMCCounts();
	ResetQWCounts();

	_dbpot     = 0.0;  //added by Hui Li
	_bpot      = 0.0;
	_btotal    = 0.0;
	_bkin        = 0.0;

	_brot        = 0.0;
	_brotsq      = 0.0;
	_bCv         = 0.0;
	_bCv_trans   = 0.0;
	_bCv_rot     = 0.0;

	PrintXYZprl = 0;

	PrintYrfl   = 1;
	PrintXrfl   = 1;
	PrintZrfl   = 1;
#ifdef DDCORR
	if(MCAtom[IMTYPE].numb > 1)
	{
		int NDIMDP = NumbAtoms*(NumbAtoms+1)/2;
		_cdipoleXYZ.resize(NDIMDP);
		_cdipoleX.resize(NDIMDP);
		_cdipoleY.resize(NDIMDP);
		_cdipoleZ.resize(NDIMDP);
		_cdipoleXY.resize(NDIMDP);
		
		for (int idp = 0; idp < NDIMDP; idp++)
		{
			_cdipoleXYZ[idp] = 0.0;
			_cdipoleX[idp] = 0.0;
			_cdipoleY[idp] = 0.0;
			_cdipoleZ[idp] = 0.0;
			_cdipoleXY[idp] = 0.0;
		}
	}
#endif

#ifdef ORDERPARA	
	if(MCAtom[IMTYPE].numb > 1)
	{
		for (int id=0; id<=NDIM; id++){
			_beiej[id]=0.0;
		}	
		for (int id=0; id<NDIM; id++){
			_bei[id]=0.0;
		}	
	}
#endif	

	_bcostheta = 0.0;
	_ucompx     = 0.0;
	_ucompy     = 0.0;
	_ucompz     = 0.0;
}

void MCResetBlockAveragePIGS(void) 
{
	avergCount = 0.0;

	ResetMCEstims();

	ResetMCCounts();
	ResetQWCounts();

	_dbpot     = 0.0;  //added by Hui Li
	_bpot      = 0.0;
	_btotal    = 0.0;
	_bcostheta = 0.0;
	_ucompx     = 0.0;
	_ucompy     = 0.0;
	_ucompz     = 0.0;
	_abs_ucompx     = 0.0;
	_abs_ucompy     = 0.0;
	_abs_ucompz     = 0.0;
#ifdef ORDERPARA	
	if(MCAtom[IMTYPE].numb > 1)
	{
		for (int id=0; id<=NDIM; id++){
			_beiej[id]=0.0;
		}	
		for (int id=0; id<NDIM; id++){
			_bei[id]=0.0;
		}	
	}
#endif	
#ifdef DDCORR
	if(MCAtom[IMTYPE].numb > 1)
	{
		int NDIMDP = NumbAtoms*(NumbAtoms+1)/2;
		_cdipoleXYZ.resize(NDIMDP);
		_cdipoleX.resize(NDIMDP);
		_cdipoleY.resize(NDIMDP);
		_cdipoleZ.resize(NDIMDP);
		_cdipoleXY.resize(NDIMDP);
		
		for (int idp = 0; idp < NDIMDP; idp++)
		{
			_cdipoleXYZ[idp] = 0.0;
			_cdipoleX[idp] = 0.0;
			_cdipoleY[idp] = 0.0;
			_cdipoleZ[idp] = 0.0;
			_cdipoleXY[idp] = 0.0;
		}
	}
#endif
	_bkin        = 0.0;

	_brot        = 0.0;
	_brot1        = 0.0;
	_brotsq      = 0.0;
	_bCv         = 0.0;
	_bCv_trans   = 0.0;
	_bCv_rot     = 0.0;

	PrintXYZprl = 0;

	PrintYrfl   = 1;
	PrintXrfl   = 1;
	PrintZrfl   = 1;

}

void MCResetBlockAveragePIGSENT(void) 
{
	avergCount = 0.0;

	ResetMCEstims();
	ResetMCCounts();
	ResetQWCounts();

	if (ENT_ENSMBL == EXTENDED_ENSMBL) {
		MCAccepSwap = 0.0;
		MCAccepUnSwap = 0.0;
	}
    _bnm       = 0.0;
    _bdm       = 0.0;
	_trOfDensitySq = 0.0;

	if (ENT_SIM) avergCountENT = 0.0;
	PrintXYZprl = 0;
	PrintYrfl   = 1;
	PrintXrfl   = 1;
	PrintZrfl   = 1;
}

void MCGetAveragePIGSENT(int type) 
{
	avergCount       += 1.0;
	totalCount       += 1.0;  

	double snm, sdm;
	if (ENT_ENSMBL == EXTENDED_ENSMBL)
	{
		snm = MCAccepSwap/(MCAccepSwap+MCAccepUnSwap);
		sdm = MCAccepUnSwap/(MCAccepSwap+MCAccepUnSwap);
	}
	if (ENT_ENSMBL == BROKENPATH)
	{
		snm = GetEstimNM(type);
		sdm = GetEstimDM(type);
	}
    _bnm             += snm;
    _bdm             += sdm;
#ifdef SWAP
	_trOfDensitySq   += sdm/snm;
	_trOfDensitySq_total += sdm/snm;
#endif
#ifdef UNSWAP
	_trOfDensitySq   += snm/sdm;
	_trOfDensitySq_total += snm/sdm;
#endif

//  reflect for MF molecule
	if(IREFLY == 1)
	{
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_XZ();
        }
//      else
//      {
//      	cout<<"Not in Reflect_MF_XZ"<<endl;
//      }
	}

    if(IREFLX == 1)
    {
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_YZ();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_YZ"<<endl;
//      }
    }

    if(IREFLZ == 1)
    {
        double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
            Reflect_MF_XY();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_XY"<<endl;
//      }
    }

    if(IROTSYM == 1)
    {
        double rndrot=rnd7();

        if(rndrot < 0.5)
        RotSymConfig();
    }
}

void MCGetAveragePIGS(void) 
{
	avergCount       += 1.0;
	totalCount       += 1.0;  

	double spot       = GetPotEnergyPIGS(); // pot energy and density distributions
	_bpot            += spot;                     // block average for pot energy
	_pot_total       += spot;
	double stotal     = GetTotalEnergy();         // Total energy
	_btotal          += stotal;                   // kin+pot
	_total           += stotal;
	double srot1      = (stotal - spot);
	_brot1           += srot1;
	_rot_total1      += srot1;

	if(MCAtom[IMTYPE].numb > 1)
	{
#ifdef IOWRITE
		double cosTheta   = 0.0;
		double compxyz[NDIM];
		double abs_compxyz[NDIM];
		GetCosThetaPIGS(cosTheta, abs_compxyz, compxyz);
		double scostheta  = cosTheta;
		double scompx     = compxyz[0];
		double scompy     = compxyz[1];
		double scompz     = compxyz[2];
		double abs_scompx = abs_compxyz[0];
		double abs_scompy = abs_compxyz[1];
		double abs_scompz = abs_compxyz[2];

		_bcostheta       += scostheta; 

		_ucompx          += scompx;
		_ucompy          += scompy;
		_ucompz          += scompz;

		_abs_ucompx          += abs_scompx;
		_abs_ucompy          += abs_scompy;
		_abs_ucompz          += abs_scompz;
#endif		

#ifdef ORDERPARA		
		double eiej[4];
		double ei[3];
		for (int id=0; id<=NDIM; id++){
			eiej[id]=0.0;
		}
		for (int id=0; id<NDIM; id++){
			ei[id]=0.0;
		}
		GetOrderCorrPIGS(eiej, ei);
		for (int id=0; id<=NDIM; id++){
			_beiej[id]+=eiej[id];
		}
		for (int id=0; id<NDIM; id++){
			_bei[id]+=ei[id];
		}
#endif		

#ifdef DDCORR
		int NDIMDP = NumbAtoms*(NumbAtoms+1)/2;
		double DipoleCorrXYZ[NDIMDP];
		double DipoleCorrX[NDIMDP];
		double DipoleCorrY[NDIMDP];
		double DipoleCorrZ[NDIMDP];
		double DipoleCorrXY[NDIMDP];

		for (int idp = 0; idp < NDIMDP; idp++)
		{
			DipoleCorrXYZ[idp] = 0.0;
			DipoleCorrX[idp]   = 0.0;
			DipoleCorrY[idp]   = 0.0;
			DipoleCorrZ[idp]   = 0.0;
			DipoleCorrXY[idp]  = 0.0;
		}

		GetDipoleCorrelationPIGS(DipoleCorrXYZ, DipoleCorrX, DipoleCorrY, DipoleCorrZ, DipoleCorrXY);

		for (int idp = 0; idp < NDIMDP; idp++)
		{
			_cdipoleXYZ[idp] += DipoleCorrXYZ[idp];
			_cdipoleX[idp]   += DipoleCorrX[idp];
			_cdipoleY[idp]   += DipoleCorrY[idp];
			_cdipoleZ[idp]   += DipoleCorrZ[idp];
			_cdipoleXY[idp]  += DipoleCorrXY[idp];
			_cdipoleXYZ_total[idp] += DipoleCorrXYZ[idp];
			_cdipoleX_total[idp]   += DipoleCorrX[idp];
			_cdipoleY_total[idp]   += DipoleCorrY[idp];
			_cdipoleZ_total[idp]   += DipoleCorrZ[idp];
			_cdipoleXY_total[idp]  += DipoleCorrXY[idp];
		}
#endif
	}

	double srot;
	if (ROTATION)
	{
		if(MCAtom[IMTYPE].molecule == 1)
		{
			srot      = GetRotEnergy();           // kin energy
		}

		if(MCAtom[IMTYPE].molecule == 3)
		{
			srot      = GetRotPlanarEnergy();     // kin energy
		}

		if(MCAtom[IMTYPE].molecule == 2)
		{
			//srot          = GetRotE3D();
			srot          = 0.0;
		}
        
		if(MCAtom[IMTYPE].molecule == 4)
		{
			srot      = GetRotEnergyPIGS();           // kin energy
		}
		_brot        += srot;
		_rot_total   += srot;
	}

//  reflect for MF molecule
	if(IREFLY == 1)
	{
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_XZ();
        }
//      else
//      {
//      	cout<<"Not in Reflect_MF_XZ"<<endl;
//      }
	}

    if(IREFLX == 1)
    {
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_YZ();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_YZ"<<endl;
//      }
    }

    if(IREFLZ == 1)
    {
        double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
            Reflect_MF_XY();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_XY"<<endl;
//      }
    }

    if(IROTSYM == 1)
    {
        double rndrot=rnd7();

        if(rndrot < 0.5)
        RotSymConfig();
    }
}

void MCGetAveragePIMC(void) 
{
	avergCount       += 1.0;
	totalCount       += 1.0;  
#ifdef HISTOGRAM
	GetDensities();
#endif

    double skin       = 0.;
	if(TRANSLATION)
	{
		skin          = GetKinEnergy();           // kin energy
	}
	_bkin            += skin;                     // block average for kin energy
	_kin_total       += skin;                     // accumulated average 

	double spot       = GetPotEnergy_Densities(); // pot energy and density distributions
	_bpot            += spot;                     // block average for pot energy
	_pot_total       += spot;


	if(MCAtom[IMTYPE].numb > 1)
	{
#ifdef IOWRITE
		double cosTheta   = 0.0;
		double compxyz[3];
		compxyz[0] = 0.0;
		compxyz[1] = 0.0;
		compxyz[2] = 0.0;

		GetCosThetaPIMC(cosTheta, compxyz);
		double scostheta  = cosTheta;
		double scompx     = compxyz[0];
		double scompy     = compxyz[1];
		double scompz     = compxyz[2];

		_bcostheta       += scostheta;
		_costheta_total  += scostheta;

		_ucompx          += scompx;
		_ucompy          += scompy;
		_ucompz          += scompz;
		_ucompx_total    += scompx;
		_ucompy_total    += scompy;
		_ucompz_total    += scompz;
#endif

#ifdef ORDERPARA		
		double eiej[4];
		double ei[3];
		for (int id=0; id<=NDIM; id++){
			eiej[id]=0.0;
		}
		for (int id=0; id<NDIM; id++){
			ei[id]=0.0;
		}
		GetOrderCorrPIMC(eiej, ei);
		for (int id=0; id<=NDIM; id++){
			_beiej[id]+=eiej[id];
		}
		for (int id=0; id<NDIM; id++){
			_bei[id]+=ei[id];
		}

#endif		

#ifdef DDCORR
		int NDIMDP = NumbAtoms*(NumbAtoms+1)/2;
		double DipoleCorrXYZ[NDIMDP];
		double DipoleCorrX[NDIMDP];
		double DipoleCorrY[NDIMDP];
		double DipoleCorrZ[NDIMDP];
		double DipoleCorrXY[NDIMDP];

		for (int idp = 0; idp < NDIMDP; idp++)
		{
			DipoleCorrXYZ[idp] = 0.0;
			DipoleCorrX[idp]   = 0.0;
			DipoleCorrY[idp]   = 0.0;
			DipoleCorrZ[idp]   = 0.0;
			DipoleCorrXY[idp]  = 0.0;
		}

		GetDipoleCorrelationPIMC(DipoleCorrXYZ, DipoleCorrX, DipoleCorrY, DipoleCorrZ, DipoleCorrXY);

		for (int idp = 0; idp < NDIMDP; idp++)
		{
			_cdipoleXYZ[idp] += DipoleCorrXYZ[idp];
			_cdipoleX[idp]   += DipoleCorrX[idp];
			_cdipoleY[idp]   += DipoleCorrY[idp];
			_cdipoleZ[idp]   += DipoleCorrZ[idp];
			_cdipoleXY[idp]  += DipoleCorrXY[idp];
			_cdipoleXYZ_total[idp] += DipoleCorrXYZ[idp];
			_cdipoleX_total[idp]   += DipoleCorrX[idp];
			_cdipoleY_total[idp]   += DipoleCorrY[idp];
			_cdipoleZ_total[idp]   += DipoleCorrZ[idp];
			_cdipoleXY_total[idp]  += DipoleCorrXY[idp];
		}
#endif
	}

	double srot;
	if (ROTATION)
	{
		if(MCAtom[IMTYPE].molecule == 1)
		{
			srot      = GetRotEnergy();           // kin energy
		}

		if(MCAtom[IMTYPE].molecule == 3)
		{
			srot      = GetRotPlanarEnergy();     // kin energy
		}

		if(MCAtom[IMTYPE].molecule == 2)
		{
			srot          = GetRotE3D();
		}
        
		if(MCAtom[IMTYPE].molecule == 4)
		{
			srot      = GetRotPlanarEnergy();     // kin energy
		}
		_brot        += srot;
		_rot_total   += srot;

		_brotsq      += ErotSQ;
		_rotsq_total += ErotSQ;

	     GetRCF(); 

	}

	double stotal     = skin + srot + spot;
	_btotal          += stotal;                   // kin+pot
	_total           += stotal;

#ifdef IOWRITE
// accumulate terms for Cv
   double sCv;
   sCv = -0.5*(double)(NDIM*NumbAtoms)/(MCBeta*MCTau)
         -((double)(NDIM*NumbAtoms)*0.5/MCTau - skin - spot - srot)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin - spot - srot)
         + (2.0/MCBeta)*(0.5*(double)(NDIM*NumbAtoms)/MCTau - skin) + Erot_termSQ - ErotSQ;

   _bCv += sCv;
   _Cv_total += sCv;

// accumulate terms for translational Cv
   double sCv_trans;
   double sCv_trans_1=-((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin);
   double sCv_trans_2=(2.0/MCBeta)*0.5*(double)(NDIM*NumbAtoms)/MCTau - skin;
   sCv_trans = -0.5*(double)(NDIM*NumbAtoms)/(MCBeta*MCTau)
         -((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)
         + (2.0/MCBeta)*(0.5*(double)(NDIM*NumbAtoms)/MCTau - skin);

   _bCv_trans += sCv_trans;
   _Cv_trans_total += sCv_trans;
   _Cv_trans_1_total += sCv_trans_1;
   _Cv_trans_2_total += sCv_trans_2;

// accumulate terms for rotational Cv
   double sCv_rot;
   sCv_rot = -srot*srot + Erot_termSQ - ErotSQ;

   _bCv_rot += sCv_rot;
   _Cv_rot_total += sCv_rot;

// check whether there're bosons in the system

	if (BOSONS) 
	{
		GetExchangeLength();
		GetPermutation();
	}

#ifdef DEBUG_PIMC
	if (NDIM !=3) 
		nrerror("MCGetAverage()","Area estimators for 3D case only"); 
#endif
	if (BOSONS and MCAtom[IMTYPE].molecule == 1) 
	GetAreaEstimators();

	if (BOSONS)
	{
		int iframe = 0;
		GetAreaEstim3D(iframe);
	}

	if ((BOSONS and MCAtom[IMTYPE].molecule == 2) && ISPHER == 0) 
	{
		int iframe = 1;
		GetAreaEstim3D(iframe);
	}
#endif

//  reflect for MF molecule
	if(IREFLY == 1)
	{
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_XZ();
        }
//      else
//      {
//      	cout<<"Not in Reflect_MF_XZ"<<endl;
//      }
	}

    if(IREFLX == 1)
    {
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_YZ();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_YZ"<<endl;
//      }
    }

    if(IREFLZ == 1)
    {
        double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
            Reflect_MF_XY();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_XY"<<endl;
//      }
    }

    if(IROTSYM == 1)
    {
        double rndrot=rnd7();

        if(rndrot < 0.5)
        RotSymConfig();
    }
}

void MCWormAverageReset(void)
{
}

void MCWormAverage(void)
{
}

void MCSaveBlockAverages(long int blocknumb) 
{
  	const char *_proc_=__func__;    //   MCSaveBlockAverages()

  	stringstream bc;                // convert block # to string
  	bc.width(IO_BLOCKNUMB_WIDTH); 
  	bc.fill('0');
  	bc<<blocknumb;

  	string fname = MCFileName + bc.str();  // file name prefix including block #

  	//-----------------------------------------------------------------
  	// densities

  	// toby temporarily by-pass the density save for the non-linear rotor
#ifdef IOWRITE 
  	if( IMPURITY && MCAtom[IMTYPE].molecule == 1)
    {
      	SaveDensities1D    (fname.c_str(),avergCount);       
      	SaveDensities2D    (fname.c_str(),avergCount,MC_BLOCK);
    }
    if( IMPURITY && MCAtom[IMTYPE].molecule == 3)
    {
      	SaveDensities1D    (fname.c_str(),avergCount);       
      	SaveDensities2D    (fname.c_str(),avergCount,MC_BLOCK);
    }
  	if( IMPURITY && MCAtom[IMTYPE].molecule == 2)
    {
      	SaveDensities1D    (fname.c_str(),avergCount);

      	SaveRho1D          (fname.c_str(),avergCount,MC_BLOCK);

      	// SaveDensities3D    (fname.c_str(),avergCount,MC_BLOCK); // this step takes lots of space. temporarily turnned off

      	SaveRhoThetaChi    (fname.c_str(),avergCount,MC_BLOCK);

      	// IOxyzAng(IOWrite,fname.c_str());
    }

	if (ROTATION) 
	{
    	SaveRCF            (fname.c_str(),avergCount,MC_BLOCK); 
	}
#endif
	if ((PIGS_SIM) || (PIMC_SIM))
	{
		SaveEnergy(MCFileName.c_str(),avergCount,blocknumb);

		if(MCAtom[IMTYPE].numb > 1)
		{
			//SaveAngularDOF(MCFileName.c_str(),avergCount,blocknumb);
#ifdef DDCORR
			SaveDipoleCorr(MCFileName.c_str(),avergCount,blocknumb);
#endif
#ifdef ORDERPARA
			SaveOrderCorr(MCFileName.c_str(),avergCount,blocknumb);
#endif
		}
	}

	if (ENT_SIM)
	{
		//SaveEnergy(MCFileName.c_str(),avergCountENT,blocknumb);
		//SaveAngularDOF(MCFileName.c_str(),avergCountENT,blocknumb);
		//SaveDipoleCorr(MCFileName.c_str(),avergCountENT,blocknumb);
		SaveTrReducedDens(MCFileName.c_str(),avergCount,blocknumb);
	}	

#ifdef IOWRITE
	if (BOSONS) 
    SaveExchangeLength (MCFileName.c_str(),avergCount,blocknumb);

	if (BOSONS && MCAtom[IMTYPE].molecule ==1)  
    SaveAreaEstimators (MCFileName.c_str(),avergCount,blocknumb);

	if (BOSONS)
    {
    	int iframe = 0;
      	SaveAreaEstim3D (MCFileName.c_str(),avergCount,blocknumb,iframe);
    }

 	if ((BOSONS && MCAtom[IMTYPE].molecule ==2) && ISPHER == 0)
    {
    	int iframe = 1;
      	SaveAreaEstim3D (MCFileName.c_str(),avergCount,blocknumb,iframe);
    }
#endif

}

void SaveEnergy (const char fname [], double acount, long int blocknumb)
{
	const char *_proc_=__func__;    //  SaveEnergy()
 
	fstream fid;
	string fenergy;

	fenergy  = fname; 
	fenergy += IO_EXT_ENG; 

	fid.open(fenergy.c_str(),ios::app | ios::out);
	io_setout(fid);

	if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());

	if (PIMC_SIM)
	{	
		fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 
		fid << setw(IO_WIDTH) << _bkin*Units.energy/acount << BLANK;    // potential anergy 2
		fid << setw(IO_WIDTH) << _brot*Units.energy/acount << BLANK;    // rot energy 5  
		fid << setw(IO_WIDTH) << _bpot*Units.energy/acount << BLANK;    // potential anergy 2
		fid << setw(IO_WIDTH) << _btotal*Units.energy/acount << BLANK;  //total energy including rot energy 
		fid << endl;
	}	

	if (PIGS_SIM)
	{	
		fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 
		fid << setw(IO_WIDTH) << _brot*Units.energy/acount << BLANK;    // rot energy 5  
		fid << setw(IO_WIDTH) << _brot1*Units.energy/acount << BLANK;    // rot energy 5  
		fid << setw(IO_WIDTH) << _bpot*Units.energy/acount << BLANK;    // potential anergy 2
		fid << setw(IO_WIDTH) << _btotal*Units.energy/acount << BLANK;  //total energy including rot energy 
		fid << endl;
	}	
//
	/*
	if (acount != 0.0)
	{
		fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 
		fid << setw(IO_WIDTH) << _brot*Units.energy/acount << BLANK;    // rot energy 5  
		fid << setw(IO_WIDTH) << _brot1*Units.energy/acount << BLANK;    // rot energy 5  
		fid << setw(IO_WIDTH) << _bpot*Units.energy/acount << BLANK;    // potential anergy 2
		fid << setw(IO_WIDTH) << _btotal*Units.energy/acount << BLANK;  //total energy including rot energy 
	    fid << endl;
	}
	*/
//
#ifdef IOWRITE
    fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 
    fid << setw(IO_WIDTH) << _bkin*Units.energy/avergCount << BLANK;    // kinetic energy 2
    fid << setw(IO_WIDTH) << _bpot*Units.energy/avergCount << BLANK;    // potential anergy 3
    fid << setw(IO_WIDTH) <<(_bkin+_bpot)*Units.energy/avergCount << BLANK;  // 4
   
    // if (ROTATION)
    fid << setw(IO_WIDTH) << _brot*Units.energy/avergCount << BLANK;    // rot energy 5
    fid << setw(IO_WIDTH) << _brotsq*(Units.energy*Units.energy)/avergCount << BLANK;    // rot energy square
    // fid << setw(IO_WIDTH) << _dbpot*Units.energy/avergCount << BLANK;    // potential differencies added by Hui Li 
    fid << setw(IO_WIDTH) <<(_bkin+_bpot+_brot)*Units.energy/avergCount << BLANK;  //total energy including rot energy 
    fid << setw(IO_WIDTH) << _bCv/avergCount << BLANK; // heat capacity
    fid << setw(IO_WIDTH) << _bCv_trans/avergCount << BLANK; // translational heat capacity
    fid << setw(IO_WIDTH) << _bCv_rot/avergCount << BLANK; // rotational heat capacity
    fid << endl;
#endif
    fid.close();
}

void SaveInstantConfig(const char fname [], long int blocknumb)
{
	const char *_proc_=__func__;   
 
	fstream fid;
	string fconfig;

	fconfig  = fname; 
	fconfig += IO_EXT_XYZ; 

	fid.open(fconfig.c_str(),ios::app | ios::out);
	io_setout(fid);

	if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fconfig.c_str());

	fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 

	if (TRANSLATION)
	{
		for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
		{
			for (int it = 0; it < NumbTimes; it++) 
			{
				int offset0 = NumbTimes*atom0;
				int t0      = offset0 + it;

				fid << setw(IO_WIDTH) << MCCoords[AXIS_X][t0] << BLANK;
				fid << setw(IO_WIDTH) << MCCoords[AXIS_Y][t0] << BLANK;
				fid << setw(IO_WIDTH) << MCCoords[AXIS_Z][t0] << BLANK;
			}
		}
	}

	if (ROTATION)
	{
		for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
		{
			for (int it = 0; it < NumbRotTimes; it++) 
			{
				int offset0 = NumbTimes*atom0;
				int t0      = offset0 + it;
				fid << setw(IO_WIDTH) << MCAngles[CTH][t0] << BLANK;
				fid << setw(IO_WIDTH) << MCAngles[PHI][t0] << BLANK;
				fid << setw(IO_WIDTH) << MCAngles[CHI][t0] << BLANK;
			}
		}
	}
	fid << endl;
    fid.close();
}

#ifdef DDCORR
void SaveDipoleCorr(const char fname [], double acount, long int blocknumb)
{
    const char *_proc_=__func__;    

    fstream fid;
    string fenergy;

    fenergy  = fname;
    fenergy += "Dipole.corr";

    fid.open(fenergy.c_str(),ios::app | ios::out);
    io_setout(fid);

    if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());
	
	if (acount != 0)
	{
   		fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;   
		if (ENT_SIM) int NumbAtoms1 = NumbAtoms/2;
		else int NumbAtoms1 = NumbAtoms;
		int NDIMDP = NumbAtoms1*(NumbAtoms1+1)/2;

		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleXYZ[idp]/acount<< BLANK;
		}
		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleX[idp]/acount<< BLANK;
		}
		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleY[idp]/acount<<BLANK;
		}
		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleZ[idp]/acount<< BLANK;
		}
		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleXY[idp]/acount<< BLANK;
		}
    	fid << endl;
	}
    fid.close();
}
#endif

#ifdef ORDERPARA
void SaveOrderCorr(const char fname [], double acount, long int blocknumb)
{
    const char *_proc_=__func__;    

    fstream fid;
    string fenergy;

    fenergy  = fname;
    fenergy += "OrderPara.corr";

    fid.open(fenergy.c_str(),ios::app | ios::out);
    io_setout(fid);

    if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());
	
	if (acount != 0)
	{
   		fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;   

		for (int id=0; id<=NDIM; id++){ 
    		fid << setw(IO_WIDTH) << _beiej[id]/acount<< BLANK;
		}
		for (int id=0; id<NDIM; id++){ 
    		fid << setw(IO_WIDTH) << _bei[id]/acount<< BLANK;
		}
    	fid << endl;
	}
    fid.close();
}
#endif

void SaveAngularDOF(const char fname [], double acount, long int blocknumb)
{
    const char *_proc_=__func__;    

    fstream fid;
    string fenergy;

    fenergy  = fname;
    fenergy += ".dof";

    fid.open(fenergy.c_str(),ios::app | ios::out);
    io_setout(fid);

    if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());

	if (acount != 0)
	{
        fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;   
        fid << setw(IO_WIDTH) << _bcostheta/acount << BLANK;
        fid << setw(IO_WIDTH) << _ucompx/acount << BLANK;
        fid << setw(IO_WIDTH) << _ucompy/acount << BLANK;
        fid << setw(IO_WIDTH) << _ucompz/acount << BLANK;
        fid << setw(IO_WIDTH) << _abs_ucompx/acount << BLANK;
        fid << setw(IO_WIDTH) << _abs_ucompy/acount << BLANK;
        fid << setw(IO_WIDTH) << _abs_ucompz/acount << BLANK;
        fid << endl;
	}
    fid.close();
}

void SaveTrReducedDens(const char fname [], double acount, long int blocknumb)
{
    const char *_proc_=__func__;

    fstream fid;
    string fenergy;

    fenergy  = fname;
    fenergy += ".rden";

    fid.open(fenergy.c_str(),ios::app | ios::out);
    io_setout(fid);

    if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());

    fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;
	fid << setw(IO_WIDTH) << _bnm/avergCount << BLANK;
	fid << setw(IO_WIDTH) << _bdm/avergCount << BLANK;
	if ((ENT_ENSMBL != BROKENPATH) && (ENT_ENSMBL != EXTENDED_ENSMBL))
	{
		fid << setw(IO_WIDTH) << _trOfDensitySq/avergCount << BLANK;
	}
    fid << endl;
    fid.close();
}

void SaveSumEnergy (double acount, double numb)  // global average
{
	const char *_proc_=__func__;    //  SaveSumEnergy()
 
    if (PIMC_SIM)
	{	
		_feng << setw(IO_WIDTH_BLOCK) << numb <<BLANK;    
		_feng << setw(IO_WIDTH) << _kin_total*Units.energy/acount << BLANK;    
		_feng << setw(IO_WIDTH) << _rot_total*Units.energy/acount << BLANK;   
		_feng << setw(IO_WIDTH) << _pot_total*Units.energy/acount << BLANK;    
		_feng << setw(IO_WIDTH) << _total*Units.energy/acount << BLANK;   
	}
    if (PIGS_SIM)
	{	
		_feng << setw(IO_WIDTH_BLOCK) << numb <<BLANK;    
		_feng << setw(IO_WIDTH) << _rot_total*Units.energy/acount << BLANK;   
		_feng << setw(IO_WIDTH) << _rot_total1*Units.energy/acount << BLANK;   
		_feng << setw(IO_WIDTH) << _pot_total*Units.energy/acount << BLANK;    
		_feng << setw(IO_WIDTH) << _total*Units.energy/acount << BLANK;   
	}
//
#ifdef IOWRITE
	_feng << setw(IO_WIDTH_BLOCK) << numb <<BLANK;    
	_feng << setw(IO_WIDTH) << _kin_total*Units.energy/acount << BLANK;    
	_feng << setw(IO_WIDTH) << _pot_total*Units.energy/acount << BLANK;    
	_feng << setw(IO_WIDTH) <<(_kin_total+_pot_total)*Units.energy/acount << BLANK;   

// if (ROTATION)
	{
		_feng << setw(IO_WIDTH) <<_rot_total*Units.energy/acount << BLANK;   
		_feng << setw(IO_WIDTH) <<_rotsq_total*(Units.energy*Units.energy)/acount << BLANK;   
	}
//_feng << setw(IO_WIDTH) << _dpot_total*Units.energy/acount << BLANK;   //added by Hui Li 
	_feng << setw(IO_WIDTH) <<(_kin_total+_pot_total+_rot_total)*Units.energy/acount << BLANK;  //total energy including rot  
// Cv
	double Cv = 0.5*(double)(NDIM*NumbAtoms*NumbTimes*Temperature)-(_kin_total+_pot_total+_rot_total)*Units.energy/acount;
	Cv = Cv*Cv + _Cv_total*Units.energy*Units.energy/acount;
	Cv = -Cv*MCBeta/Temperature;
	_feng << setw(IO_WIDTH) <<Cv << BLANK;

	double Cv_trans = 0.5*(double)(NDIM*NumbAtoms*NumbTimes*Temperature)-_kin_total*Units.energy/acount;
	Cv_trans = Cv_trans*Cv_trans+ _Cv_trans_total*Units.energy*Units.energy/acount;
	Cv_trans = -Cv_trans*MCBeta/Temperature;
	_feng << setw(IO_WIDTH) <<Cv_trans << BLANK;

// if (ROTATION)
	{
		double Cv_rot = -_rot_total*Units.energy/acount;
		Cv_rot = Cv_rot*Cv_rot+ _Cv_rot_total*Units.energy*Units.energy/acount;
		Cv_rot = -Cv_rot*MCBeta/Temperature;
		_feng << setw(IO_WIDTH) <<Cv_rot << BLANK;
	}
#endif

//_feng << setw(IO_WIDTH) << _Cv_trans_1_total*Units.energy*Units.energy/acount << BLANK;
//_feng << setw(IO_WIDTH) << _Cv_trans_2_total*Units.energy*Units.energy/acount << BLANK;

	_feng << endl;
}

#ifdef DDCORR
void SaveSumDipoleCorr(double acount, double numb)
{
    const char *_proc_=__func__;

    _fdc << setw(IO_WIDTH_BLOCK) << numb <<BLANK;
	if (ENT_SIM) int NumbAtoms1 = NumbAtoms/2;
	else int NumbAtoms1 = NumbAtoms;
	int NDIMDP = NumbAtoms1*(NumbAtoms1+1)/2;

	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleXYZ_total[idp]/acount << BLANK;
	}
	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleX_total[idp]/acount << BLANK;
	}
	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleY_total[idp]/acount << BLANK;
	}
	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleZ_total[idp]/acount << BLANK;
	}
	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleXY_total[idp]/acount << BLANK;
	}
    _fdc << endl;
}
#endif

void SaveInstantEnergy()
{
    const char *_proc_=__func__;

    double srotinst, spotinst, stotalinst;
    if (MCAtom[IMTYPE].molecule == 4)
    {
		if (PIGS_SIM)
		{	
			srotinst   = GetRotEnergyPIGS();
			spotinst   = GetPotEnergyPIGS(); 
			stotalinst = GetTotalEnergy();
		}	
		if (PIMC_SIM)
		{	
			srotinst   = GetRotPlanarEnergy(); 
			spotinst   = GetPotEnergy_Densities(); 
			stotalinst = 0.0;
		}	
    }
    _fengins << setw(IO_WIDTH) << spotinst << BLANK;
    _fengins << setw(IO_WIDTH) << stotalinst << BLANK;
    _fengins << setw(IO_WIDTH) << srotinst << BLANK;
    _fengins << endl;
}

void InitTotalAverage(void)  // DUMP
{
	const char *_proc_=__func__;    

	totalCount = 0.0;  // need to save in the status file 
	if (ENT_SIM) totalCountENT = 0.0;  // need to save in the status file 
   	sumsCount = 0.0;  // counter for accum sums

	_kin_total = 0.0;   // need a function to reset all global average
	_pot_total = 0.0;
	_total = 0.0;

	_dpot_total = 0.0;  //added by Hui Li

	_rot_total = 0.0;
	_rot_total1 = 0.0;
	_rotsq_total = 0.0;
	_Cv_total = 0.0;
	_Cv_trans_total = 0.0;
	_Cv_trans_1_total = 0.0;
	_Cv_trans_2_total = 0.0;
	_Cv_rot_total = 0.0;

	// open files for output
	// ENERGY
//
	if ((PIGS_SIM) || (PIMC_SIM))
	{
		string fenergy;

		fenergy  = MCFileName + IO_SUM; 
		fenergy += IO_EXT_ENG; 
	 
		if (FileExist(fenergy.c_str()))   // backup the output of previous simulations 
		IOFileBackUp(fenergy.c_str());

		_feng.open(fenergy.c_str(), ios::out);
		io_setout(_feng);

		if (!_feng.is_open())
		_io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());
	}

//
#ifdef INSTANT
    string fangularins;

    fangularins  = MCFileName + "_instant";
    fangularins += ".dof";

    if (FileExist(fangularins.c_str()))   // backup the output of previous simulations 
    IOFileBackUp(fangularins.c_str());

#ifdef BINARY
    _fangins.open(fangularins.c_str(), ios::out | ios::binary);
#else
    _fangins.open(fangularins.c_str(), ios::out);
#endif
    io_setout(_fangins);

    if (!_fangins.is_open())
    _io_error(_proc_,IO_ERR_FOPEN,fangularins.c_str());
#endif
}

void DoneTotalAverage(void)
{
  _feng.close();
}

void MCSaveAcceptRatio(long int step,long int pass,long int block)
{
	int w = 8; 
 
	cout << "BLOCK:" << setw(w) << block << BLANK;
	cout << "PASS:"  << setw(w) << pass  << BLANK;
	cout << "STEP:"  << setw(w) << step  << BLANK;

	if(TRANSLATION)
	{
		for (int type=0;type<NumbTypes;type++)
		if (WORM && (type == Worm.type))
		{
			double ratio_open    = QWAccep[0][QW_OPEN] /QWTotal[0][QW_OPEN];
			double ratio_close   = QWAccep[0][QW_CLOSE]/QWTotal[0][QW_CLOSE];

			double ratio_advance = QWAccep[0][QW_ADVANCE]/QWTotal[0][QW_ADVANCE];
			double ratio_recede  = QWAccep[0][QW_RECEDE] /QWTotal[0][QW_RECEDE];

			double ratio_swap    = QWAccep[0][QW_SWAP]/QWTotal[0][QW_SWAP];

			cout<<setw(w)<<"open/close";
			cout<<" [ "; 
			cout<<QWTotal[0][QW_OPEN]/countQW<<"-"<<QWTotal[0][QW_CLOSE]/countQW; 
			cout<<" ] "; 
			cout<<BLANK;                     
 
			cout<<setw(w)<<ratio_open   << BLANK; // accept ratio for "open"   move
			cout<<setw(w)<<ratio_close  << BLANK; // accept ratio for "close"  move

			cout<<setw(w)<<"advance/recede";
			cout<<" [ "; 
			cout<<QWTotal[0][QW_ADVANCE]/countQW<<"-"<<QWTotal[0][QW_RECEDE]/countQW; 
			cout<<" ] "; 
			cout<<BLANK;                     
 
			cout<<setw(w)<<ratio_advance<< BLANK;  // accept ratio for "advance" move 
			cout<<setw(w)<<ratio_recede << BLANK;  // accept ratio for "recede"  move

			cout<<setw(w)<<"swap";
			cout<<" [ "; 
			cout<<QWTotal[0][QW_SWAP]/countQW; 
			cout<<" ] "; 
			cout<<BLANK;                     
			cout<<setw(w)<<ratio_swap << BLANK;  // accept ratio for "recede"  move
		}
		else
		{
#ifndef GAUSSIANMOVE
      		double ratio_molec = MCAccep[type][MCMOLEC]/MCTotal[type][MCMOLEC];
      		double ratio_multi = MCAccep[type][MCMULTI]/MCTotal[type][MCMULTI];
 
      		cout<<setw(w)<<MCAtom[type].type<<BLANK; // atom type

      		cout<<setw(w)<<ratio_molec<<BLANK;       // accept ratio for "molecular" move
      		cout<<setw(w)<<ratio_multi<<BLANK;       // accept ratio for multilevel move 
			if (PIGS_SIM) {
				double ratio_end_beads = MCAccepEndBeads[type][MCMOLEC]/MCTotalEndBeads[type][MCMOLEC];
				cout<<setw(w)<<ratio_end_beads<<BLANK;       // accept ratio for multilevel move 
			}
#endif
   		}
	}

   if (ROTATION)
   {
      double ratio_rotat = MCAccep[IMTYPE][MCROTAT]/MCTotal[IMTYPE][MCROTAT];
      cout<<BLANK<< "Rot: ";
      cout<<setw(w)<<ratio_rotat<<BLANK;       // accept ratio for "molecular" move
   }

   cout << endl;

}
