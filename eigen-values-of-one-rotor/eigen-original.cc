#include <stdio.h>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include "eigen.h"
#include "setup.h"
#include "input.h"
/* Complex datatype */
//struct _dcomplex { double re, im; };
//typedef struct _dcomplex dcomplex;

using namespace std;

int main()
{
	vinit_();

	IOReadParams(FINPUT);

    double cosTheta[nSize];
	double weightsTheta[nSize];
	thetagrid(nSize, cosTheta, weightsTheta);

	normalizationCheck(nSize, nSizeRot);

    double phi[2*nSize+1];
	double weightsPhi;
	phigrid(nSize, phi, weightsPhi);
	
	int nSizeTotal = sizeloop(nSizeRot, nSkip);

	double hMatrixRe[nSizeTotal][nSizeTotal];
	double hMatrixIm[nSizeTotal][nSizeTotal];

//Bra loop started here//

	int indexBra = 0;

	for (int i1Bra = 0; i1Bra<nSizeRot; i1Bra++)
	{
		if (i1Bra%nSkip) continue;

		for (int m1Bra = -i1Bra; m1Bra<=i1Bra; m1Bra++)
		{
			int phase1Bra;
			if (m1Bra>0) phase1Bra = pow((-1.),(m1Bra));
       		else phase1Bra=1.;

//Ket loop started here//

			int indexKet = 0;

			for (int i1Ket = 0; i1Ket<nSizeRot; i1Ket++)
			{		
				if (i1Ket%nSkip) continue;

				for (int m1Ket = -i1Ket; m1Ket<=i1Ket; m1Ket++)
				{
					int phase1Ket;
					if (m1Ket>0) phase1Ket = pow((-1.),(m1Ket));
   					else phase1Ket=1.;

//loop over coordinates start here//

                   	double sumRotor1Re 		 = 0.0;
                   	double sumRotor1Im 		 = 0.0;

					for (int iXRotor1 = 0; iXRotor1<nSize; iXRotor1++)
                    {
   				        double xRotor1        = cosTheta[iXRotor1];
						double theta1         = acos(xRotor1);

						double associatedLegBra1 = phase1Bra*normalizedPlm(i1Bra, m1Bra, xRotor1)*sqrt(weightsTheta[iXRotor1]);
						double associatedLegKet1 = phase1Ket*normalizedPlm(i1Ket, m1Ket, xRotor1)*sqrt(weightsTheta[iXRotor1]);

                   		double sumPhiRotor1Re = 0.0;
                   		double sumPhiRotor1Im = 0.0;

						for (int iPhiRotor1 = 0; iPhiRotor1<(2*nSize+1); iPhiRotor1++)
						{
							double spherHarmonicsReBra1  = associatedLegBra1*cos((double)m1Bra*phi[iPhiRotor1])*sqrt(weightsPhi);
							double spherHarmonicsImBra1  = associatedLegBra1*sin((double)m1Bra*phi[iPhiRotor1])*sqrt(weightsPhi);
							double spherHarmonicsReKet1  = associatedLegKet1*cos((double)m1Ket*phi[iPhiRotor1])*sqrt(weightsPhi);
							double spherHarmonicsImKet1  = associatedLegKet1*sin((double)m1Ket*phi[iPhiRotor1])*sqrt(weightsPhi);

//potential is called here
							double potl = 10.;
							double vpot 				= potl*cosTheta[iXRotor1];
//potential call is end here
							sumPhiRotor1Re += (spherHarmonicsReBra1*spherHarmonicsReKet1 + spherHarmonicsImBra1*spherHarmonicsImKet1)*vpot;
							sumPhiRotor1Im += (spherHarmonicsReBra1*spherHarmonicsImKet1 - spherHarmonicsImBra1*spherHarmonicsReKet1)*vpot;

						}
						sumRotor1Re          += sumPhiRotor1Re;
                        sumRotor1Im          += sumPhiRotor1Im;
					}

//loop over coordinates end here//


					double energyRotor;
					if (i1Bra == i1Ket && m1Bra == m1Ket) 
					{
						energyRotor = rotEnergy(i1Bra);
					}
					else
					{
						energyRotor = 0.0;
					}

#ifndef TESTNORM 
#ifdef EXCLUDEROTATION
					hMatrixRe[indexBra][indexKet] = sumRotor1Re;
					hMatrixIm[indexBra][indexKet] = sumRotor1Im;
#else
					hMatrixRe[indexBra][indexKet] = sumRotor1Re + energyRotor;
					hMatrixIm[indexBra][indexKet] = sumRotor1Im;
#endif
#else
					hMatrixRe[indexBra][indexKet] = sumRotor1Re;
					hMatrixIm[indexBra][indexKet] = sumRotor1Im;
#endif

#ifdef IOREAD
					if (sumRotor1Re > 1e-15)
					{
						cout<<" i1Bra, i1Ket, m1Bra, m1Ket "<<i1Bra<<"  "<<i1Ket<<"  "<<m1Bra<<"  "<<m1Ket<<"  "<< setw(10)<<sumRotor1Re<<"     "<<setw(10)<<sumRotor1Im<<endl;
					}
#endif

					indexKet += 1;
				}
			}

//Ket loop ended here//

			indexBra += 1;
		}
	}

//Bra loop ended here//
	
	double c[nSizeTotal*nSizeTotal];
	double d[nSizeTotal*nSizeTotal];

	int ii = 0;
	for (int i=0; i<nSizeTotal; i++)
	{
		for (int j=0; j<nSizeTotal; j++)
		{
			c[ii] = hMatrixRe[i][j];
			d[ii] = hMatrixIm[i][j];
			if(hMatrixRe[i][j] != hMatrixRe[j][i])
			{
				cout<<i<<"   "<<j<<"   "<<setw(20)<<hMatrixRe[i][j]<<"      "<<setw(20)<<hMatrixRe[j][i]<<endl;
			}
			ii++;
		}
	}
    matdiaghermitian(nSizeTotal, c, d);
#ifdef IOREAD
	normCheckGroundState(nSizeTotal, cosTheta, weightsTheta, phi, weightsPhi, c, d);
#endif
	normCheckGroundState(nSizeTotal, cosTheta, weightsTheta, phi, weightsPhi, c, d);
	densityCosTheta(nSizeTotal, cosTheta, weightsTheta, phi, weightsPhi, c, d);

	return 0;
}

void matdiagsymmetric(int nSizeTotal, double *a)
{
//Lapack starts here
    /* Locals */

	int N = nSizeTotal;
	int LDA = N;
	a[LDA*N];
    int n = N, lda = LDA, info, lwork;
    double wkopt;
    double* work;

    /* Local arrays */
    double w[N];

    /* Query and allocate the optimal workspace */
    lwork = -1;
    char msg1[] = "Vectors";
    char msg2[] = "Upper";
    dsyev_( msg1, msg2, &n, a, &lda, w, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );

    /* Solve eigenproblem */
    dsyev_( msg1, msg2, &n, a, &lda, w, work, &lwork, &info );

    /* Check for convergence */
    if(info>0)
    {
        cout<<"The algorithm failed to compute eigenvalues.\n";
        exit(1);
    }

    /* Print eigenvalues */
    char msg3[] = "Eigenvalues";
//    print_matrix(msg3, 1, n, w, 1);

    /* Print eigenvectors */
/*    char msg4[] = "Eigenvectors (stored columnwise)";
    print_matrix(msg4, n, n, a, lda );*/

    /* Free workspace */
    free((void*)work);
    exit(0);

// End of DSYEV Example 
}

void matdiaghermitian(int nSizeTotal, double *c, double *d)
{
//Lapack starts here
    /* Locals */

	c[nSizeTotal*nSizeTotal];
	d[nSizeTotal*nSizeTotal];

	int N = nSizeTotal;
	int LDA = N;
    int n = N, lda = LDA, info, lwork;
	dcomplex a[LDA*N];

	for (int i = 0; i < nSizeTotal*nSizeTotal; i++)
	{
		a[i] = {c[i],d[i]};
	}
	dcomplex wkopt;
	dcomplex* work;

    /* Local arrays */
    double w[N], rwork[3*N-2];

    /* Query and allocate the optimal workspace */
    lwork = -1;
    char msg1[] = "Vectors";
    char msg2[] = "Lower";
	zheev_( msg1, msg2, &n, a, &lda, w, &wkopt, &lwork, rwork, &info );
	lwork = (int)wkopt.re;
	work = (dcomplex*)malloc( lwork*sizeof(dcomplex) );

    /* Solve eigenproblem */
    zheev_( msg1, msg2, &n, a, &lda, w, work, &lwork, rwork, &info );
	for (int i = 0; i < nSizeTotal*nSizeTotal; i++)
    {
        c[i] = a[i].re;
        d[i] = a[i].im;
    }

    /* Check for convergence */
    if( info > 0 ) 
	{
		cout<< "The algorithm failed to compute eigenvalues."<<endl;
        exit( 1 );
    }

    /* Print eigenvalues */
    char msg3[] = "Eigenvalues";
    char msg4[] = "Eigenvectors (stored columnwise)";

	print_rmatrix( msg3, 1, n, w, 1 );
	/* Print eigenvectors */
//    print_matrix( msg4, n, n, a, lda );

    /* Free workspace */
    free( (void*)work );
//    exit( 0 );
} /* End of ZHEEV Example */


/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, double* a, int lda ) 
{
	int i, j;

	stringstream bc;
    bc.width(1);
    bc.fill('0');
    bc<<nSizeRot;
    string fname = "Eigen-values-nSizeRot" + bc.str()+".txt";

    ofstream myfile;
    myfile.open(fname.c_str());

    for( i = 0; i < m; i++ ) 
	{
    	for( j = 0; j < n; j++ )
		{
			myfile<<setw(14)<< a[i+j*lda] <<endl;
        }
	}
	myfile.close();
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, dcomplex* a, int lda ) 
{
	int i, j;
    for( i = 0; i < m; i++ ) 
	{
		double sum = 0.;
    	for( j = 0; j < n; j++ ) 
		{
			sum += a[i+j*lda].re*a[i+j*lda].re+a[i+j*lda].im*a[i+j*lda].im;
		}
		//cout <<sum<<endl;
	}
}

/* Auxiliary routine: printing a matrix */
/*void print_matrix_real( char* desc, int m, int n, double* a, int lda )
{
    int i, j;

    stringstream bc;
	bc.width(1);
	bc.fill('0');
    bc<<nSizeRot;
	string fname = "Eigen-vectors-nSizeRot" + bc.str()+".txt";

    ofstream myfile;
    myfile.open(fname.c_str());

    for(i=0; i<m; i++)
    {
        for(j = 0; j<n; j+( = 0; i < m; i++ )
    {
        double sum = 0.;
        for( j = 0; j < n; j++ )
        {
            sum += a[i+j*lda].re*a[i+j*lda].re+a[i+j*lda].im*a[i+j*lda].im;
        }
        cout <<sum<<endl;
    })
        {
            myfile<<setw(14)<<a[i+j*lda]<<endl;
        }
        cout<<endl;
    }
	myfile.close();
}
*/
void normalizationCheck(int nSize, int nSizeRot)
{
    double cosTheta[nSize];
	double weightsTheta[nSize];
	thetagrid(nSize, cosTheta, weightsTheta);

    double phi[2*nSize+1];
	double weightsPhi;
	phigrid(nSize, phi, weightsPhi);

    ofstream myfile;
    myfile.open("Norm-Legendre-Polynomial.txt");
	for (int i = 0; i<nSize; i++)
	{
		for (int j = 0; j<nSize; j++)
		{
			double sum=0.0;
			for (int k = 0; k<nSize; k++)
			{
				double x=cosTheta[k];
				sum += normalizedPl0(i, x)*normalizedPl0(j, x)*weightsTheta[k];
			}
			myfile<<" i "<<i<<" j "<<j<<"   "<<sum<<endl;
		}
	}
	myfile.close();

    myfile.open("Norm-Spherical-Harmonics.txt");

	for (int i = 0; i<nSizeRot; i++)
	{
		for (int j = 0; j<nSizeRot; j++)
		{
			for (int m = 0; m<i+1; m++)
			{
				for (int n = 0; n<j+1; n++)
				{
	            	double sumThetaPhi=0.0;
	        	    for (int k = 0; k<nSize; k++)
    	        	{
       	        		double x=cosTheta[k];
						double associatedLegBra = normalizedPlm(i, m, x)*sqrt(weightsTheta[k]);
						double associatedLegKet = normalizedPlm(j, n, x)*sqrt(weightsTheta[k]);
						
						double sumPhi = 0.0;
	        	    	for (int l = 0; l<(2*nSize+1); l++)
    	        		{
							double spherHarmonicsReBra = associatedLegBra*cos((double)m*phi[l])*sqrt(weightsPhi);
							double spherHarmonicsImBra = associatedLegBra*sin((double)m*phi[l])*sqrt(weightsPhi);
							double spherHarmonicsReKet = associatedLegKet*cos((double)n*phi[l])*sqrt(weightsPhi);
							double spherHarmonicsImKet = associatedLegKet*sin((double)n*phi[l])*sqrt(weightsPhi);
            	    		sumPhi += spherHarmonicsReBra*spherHarmonicsReKet+spherHarmonicsImBra*spherHarmonicsImKet;
            			}
						sumThetaPhi += sumPhi;
            		}
            		myfile<<" i "<<i<<" j "<<j<<" m "<<m<<" n "<<n<<"   "<<sumThetaPhi<<endl;
				}
			}
		}
	}
    myfile.close();

}	

double normalizedPlm(int j,int m, double x)
{
	if (j < m) cerr<<"j < m"<<endl;
	double jv=(double)j;
	double mv=fabs((double)m);
	return sqrt((jv+.5)*exp(lgamma(jv-mv+1.)-lgamma(jv+mv+1.)))*plgndr(j,abs(m),x); 

	// for the normalization constant sqrt(...) refer to 6.55 in Quantum Chemistry. 
	//lgamma is the natural logarithim of the gamma function: gamma(n) = (n-1)! 
}

double normalizedPl0(int j,double x)
{
	return sqrt((double)j+.5)*plgndr(j,0,x);
}

void thetagrid(int nSize, double *cosTheta, double *weightsTheta)
{
	//Gauss-Legendre quadrature points and wiights are generated 
	//by gauleg.f

    double x1=-1.;
    double x2=1.;
    double x[nSize], w[nSize];
    cosTheta[nSize];
	weightsTheta[nSize];

    gauleg(x1,x2,x,w,nSize);

    ofstream myfile;
    myfile.open("DVR-theta.txt");

    for (int i=0; i<nSize; i++)
    {
		cosTheta[i]     = x[i];
		weightsTheta[i] = w[i];

        myfile<<"DVR theta = "<<cosTheta[i]<<"  weights =  "<<weightsTheta[i]<<endl;
    }
    myfile.close();
}

void phigrid(int nSize, double *phi, double &weightsPhi)
{

	const double PI  =3.141592653589793238463;
    phi[2*nSize+1];

    ofstream myfile;
    myfile.open("DVR-phi.txt");

    for (int i=0; i<(2*nSize+1); i++)
    {
        phi[i]          = 2.0*PI*((double)i+1.0)/(2.0*(double)nSize+1.0);
        weightsPhi      = 1.0/(2.0*(double)nSize+1.0);

        myfile<<"DVR phi = "<<phi[i]<<"  weights =  "<<weightsPhi<<endl;
    }
    myfile.close();
}

int sizeloop(int nSizeRot, int nSkip)
{
    int indexBra = 0;

    for (int i1Bra = 0; i1Bra<nSizeRot; i1Bra++)
    {
        if (i1Bra%nSkip) continue;

        for (int m1Bra = -i1Bra; m1Bra<=i1Bra; m1Bra++)
        {
			indexBra++;
		}
	}
	return indexBra;
}


double rotEnergy(int l)
{
	double energyj0_v0    = -36117.5942855;  // unit cm^-1
	double energyj1_v0    = -35999.1009407;  // unit cm^-1
	//double rotConstant    = 0.5*(energyj1_v0 - energyj0_v0)/wavenumber;
	double rotConstant    = 1.;
	double energy         = rotConstant*(double)l*((double)l + 1.0);

	return energy;
}

void normCheckGroundState(int nSizeTotal, double *cosTheta, double *weightsTheta, double *phi, double weightsPhi, double *c, double *d)
{
	c[nSizeTotal*nSizeTotal];
    d[nSizeTotal*nSizeTotal];

    cosTheta[nSize];
    weightsTheta[nSize];

    phi[2*nSize+1];


	cout <<"Normalization test of wavefunctions"<<endl;
	for( int i = 0; i < nSizeTotal; i++ )
    {
        double sum = 0.;
        for( int j = 0; j < nSizeTotal; j++ )
        {   
            sum += c[i+j*nSizeTotal]*c[i+j*nSizeTotal]+d[i+j*nSizeTotal]*d[i+j*nSizeTotal];
        }
        cout <<" Psi "<<"["<<i<<"] = "<<sum<<endl;
    }


    vector<double> psiRe;
    vector<double> psiIm;
    for (int iXRotor1   = 0; iXRotor1<nSize; iXRotor1++)
    {
        double xRotor1  = cosTheta[iXRotor1];
        double theta1   = acos(xRotor1);

        for (int iPhiRotor1    = 0; iPhiRotor1<(2*nSize+1); iPhiRotor1++)
        {

  			int i = 0;
            //Bra loop started here//

            int indexBra = 0;
            double sumRe = 0.;
            double sumIm = 0.;

            for (int i1Bra = 0; i1Bra<nSizeRot; i1Bra++)
            {
           		if (i1Bra%nSkip) continue;

                for (int m1Bra = -i1Bra; m1Bra<=i1Bra; m1Bra++)
                {
                    int phase1Bra;
                    if (m1Bra>0) phase1Bra = pow((-1.),(m1Bra));
                    else phase1Bra=1.;

                    double associatedLegBra1     = phase1Bra*normalizedPlm(i1Bra, m1Bra, xRotor1)*sqrt(weightsTheta[iXRotor1]);
                    double spherHarmonicsReBra1  = associatedLegBra1*cos((double)m1Bra*phi[iPhiRotor1])*sqrt(weightsPhi);
                    double spherHarmonicsImBra1  = associatedLegBra1*sin((double)m1Bra*phi[iPhiRotor1])*sqrt(weightsPhi);

                    double psiGroundRe           = spherHarmonicsReBra1*c[i+indexBra*nSizeTotal] - spherHarmonicsImBra1*d[i+indexBra*nSizeTotal];
                    double psiGroundIm           = spherHarmonicsImBra1*c[i+indexBra*nSizeTotal] + spherHarmonicsReBra1*d[i+indexBra*nSizeTotal];

                    sumRe += psiGroundRe;
                    sumIm += psiGroundIm;
                    indexBra++;
                }
            }
            psiRe.push_back(sumRe);
            psiIm.push_back(sumIm);
        }
    }
    double sum = 0.;
    for (int i = 0; i < psiRe.size(); i++)
    {
        sum +=  psiRe[i]*psiRe[i]+psiIm[i]*psiIm[i];

    }
	cout <<"Normalization test of ground state wavefuntion :"<<endl;
    cout<<"Norm : "<<sum<<endl;
	psiRe.clear();
	psiIm.clear();
}


void densityCosTheta(int nSizeTotal, double *cosTheta, double *weightsTheta, double *phi, double weightsPhi, double *c, double *d)
{
	c[nSizeTotal*nSizeTotal];
    d[nSizeTotal*nSizeTotal];

    cosTheta[nSize];
    weightsTheta[nSize];

    phi[2*nSize+1];


	vector<double> psiRe;
    vector<double> psiIm; 
    for (int iXRotor1   = 0; iXRotor1<nSize; iXRotor1++)
    {
        double xRotor1  = cosTheta[iXRotor1];
        double theta1   = acos(xRotor1);

        for (int iPhiRotor1    = 0; iPhiRotor1<(2*nSize+1); iPhiRotor1++)
        {

            int i = 0;
            //Bra loop started here//

            int indexBra = 0;
            double sumRe = 0.;
            double sumIm = 0.;

            for (int i1Bra = 0; i1Bra<nSizeRot; i1Bra++)
            {
                if (i1Bra%nSkip) continue;

                for (int m1Bra = -i1Bra; m1Bra<=i1Bra; m1Bra++)
                {
                    int phase1Bra;
                    if (m1Bra>0) phase1Bra = pow((-1.),(m1Bra));
                    else phase1Bra=1.;

                    double associatedLegBra1     = phase1Bra*normalizedPlm(i1Bra, m1Bra, xRotor1);
                    double spherHarmonicsReBra1  = associatedLegBra1*cos((double)m1Bra*phi[iPhiRotor1]);
                    double spherHarmonicsImBra1  = associatedLegBra1*sin((double)m1Bra*phi[iPhiRotor1]);

                    double psiGroundRe           = spherHarmonicsReBra1*c[i+indexBra*nSizeTotal] - spherHarmonicsImBra1*d[i+indexBra*nSizeTotal];
                    double psiGroundIm           = spherHarmonicsImBra1*c[i+indexBra*nSizeTotal] + spherHarmonicsReBra1*d[i+indexBra*nSizeTotal];

                    sumRe += psiGroundRe;
                    sumIm += psiGroundIm;
                    indexBra++;
                }
            }
            psiRe.push_back(sumRe);
            psiIm.push_back(sumIm);
        }
    }

    string fname1 = "Density-CosTheta.txt";
    string fname2 = "Density-Theta.txt";

    ofstream myfile1, myfile2;
    myfile1.open(fname1.c_str());
    myfile2.open(fname2.c_str());

	int iii = 0;
    for (int iXRotor1   = 0; iXRotor1<nSize; iXRotor1++) 
	{
        double xRotor1  = cosTheta[iXRotor1];
        double theta1   = acos(xRotor1);

		double sum = 0.;
        for (int iPhiRotor1    = 0; iPhiRotor1<(2*nSize+1); iPhiRotor1++) 
		{
			sum += (psiRe[iii]*psiRe[iii] + psiIm[iii]*psiIm[iii])*weightsPhi;
			
			iii++;
        }
		myfile1<<xRotor1<<"     "<<sum<<endl;
		myfile2<<theta1<<"     "<<sum*sin(theta1)<<"     "<<endl;
    }
	myfile1.close();
	myfile2.close();
}

/*
//==========================bins============================//
    vector<double> cosThetaTld;
	for (int iXRotor1   = 0; iXRotor1<nSize; iXRotor1++) {
    	double xRotor1  = cosTheta[iXRotor1];
        double theta1   = acos(xRotor1);

		for (int iXRotor2    = 0; iXRotor2<nSize; iXRotor2++) {
        	double xRotor2   = cosTheta[iXRotor2];
            double theta2    = acos(xRotor2);

          	for (int iPhiRotor1    = 0; iPhiRotor1<(2*nSize+1); iPhiRotor1++) {
				for (int iPhiRotor2     = 0; iPhiRotor2<(2*nSize+1); iPhiRotor2++) {

					double xx = sin(theta1)*cos(phi[iPhiRotor1])*sin(theta2)*cos(phi[iPhiRotor2]) + sin(theta1)*sin(phi[iPhiRotor1])*sin(theta2)*sin(phi[iPhiRotor2]) + cos(theta1)*cos(theta2);
					cosThetaTld.push_back(xx);
				}
			}
		}
	}

	//cout << "extended vector size = " << cosThetaTld.size() << endl;
	for(int i = 0; i < cosThetaTld.size(); i++)
	{
    	//cout << "value of vec [" << i << "] = " << cosThetaTld[i] << endl;
    	cout << i<<"      "<<cosThetaTld[i] << endl;
   	}


	double ctmin = -1.1;
	double ctmax = 1.1;
	int numbins  = 101;
	double dbins = (ctmax - ctmin)/((double)numbins - 1.0);
	double hist[numbins - 1];

	double sum1 = 0.;
	for (int i = 0; i < numbins - 1; i++)
	{
		double binsi = ctmin + (double)i*dbins;
		double binsf = binsi + dbins;
		hist[i] = 0.;

		for (int j = 0; j < cosThetaTld.size(); j++)
		{
			if (cosThetaTld[j] >= binsi && cosThetaTld[j] <= binsf) hist[i]+=1.;
		}
		sum1 += hist[i]*dbins;
	}
	double sum = 0.;
	for (int i = 0; i < numbins - 1; i++) {
		//cout<<i<<"  "<<hist[i]/cosThetaTld.size()<<endl;
		//sum += (hist[i]/(double)cosThetaTld.size());
		sum += hist[i]/sum1*dbins;
		double binsi = ctmin + (double)i*dbins;
        double binsf = binsi + dbins;
		cout<<binsi<<"  "<<hist[i]/sum1<<endl;
	}
	cout<<sum<<endl;
*/
