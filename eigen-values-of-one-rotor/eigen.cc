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

	double hMatrix[nSizeTotal][nSizeTotal];

	int indexBra  = 0;
	for (int iBra = 0; iBra<nSizeRot; iBra++)
	{
		if (iBra%nSkip) continue;

		int indexKet  = 0;
		for (int iKet = 0; iKet<nSizeRot; iKet++)
		{		
			if (iKet%nSkip) continue;

           	double sum = 0.0;
			for (int iXRotor  = 0; iXRotor<nSize; iXRotor++)
            {
   		        double xRotor = cosTheta[iXRotor];
				double theta  = acos(xRotor);
				double plBra  = normalizedPl0(iBra, xRotor)*sqrt(weightsTheta[iXRotor]);
				double plKet  = normalizedPl0(iKet, xRotor)*sqrt(weightsTheta[iXRotor]);

				double vpot   = -10.0*xRotor;
				sum          += plBra*vpot*plKet;
			}
			double energyRotor;
			if (iBra == iKet) energyRotor = rotEnergy(iBra);
			else energyRotor = 0.0;
			hMatrix[indexBra][indexKet] = sum + energyRotor;

			indexKet += 1;
		}
		indexBra += 1;
	}
	
	int ii = 0;
	double c[nSizeTotal*nSizeTotal];
	for (int i=0; i<nSizeTotal; i++)
	{
		for (int j=0; j<nSizeTotal; j++)
		{
			if(hMatrix[i][j] != hMatrix[j][i])
			{
				cout<<i<<"   "<<j<<"   "<<setw(20)<<hMatrix[i][j]<<"      "<<setw(20)<<hMatrix[j][i]<<endl;
			}
			c[ii] = hMatrix[i][j];
			ii++;
		}
	}
	matdiagsymmetric(nSizeTotal, c);
#ifdef IOREAD
	normCheckGroundState(nSizeTotal, cosTheta, weightsTheta, c);
#endif
	densityCosTheta(nSizeTotal, cosTheta, weightsTheta, c);

	return (0);
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
    print_matrix(msg3, 1, n, w, 1);

    /* Print eigenvectors */
/*    char msg4[] = "Eigenvectors (stored columnwise)";
    print_matrix(msg4, n, n, a, lda );*/

    /* Free workspace */
    free((void*)work);
}

void print_matrix( char* desc, int m, int n, double* a, int lda )
{
    int i, j;

    stringstream bc;
    bc.width(1);
    bc.fill('0');
    bc<<nSizeRot;
    string fname = "Eigen-values-nSizeRot" + bc.str()+".txt";

    ofstream myfile;
    myfile.open(fname.c_str());

    for(i=0; i<m; i++)
    {
        for(j = 0; j<n; j++)
        {
            myfile<<setw(14)<<a[i+j*lda]<<endl;
        }
        cout<<endl;
    }
    myfile.close();
}

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
	        	    	for (int l = 0; l<2*nSize+1; l++)
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

    for (int i=0; i<2*nSize+1; i++)
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

    for (int iBra = 0; iBra<nSizeRot; iBra++)
    {
        if (iBra%nSkip) continue;
		indexBra++;
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

void normCheckGroundState(int nSizeTotal, double *cosTheta, double *weightsTheta, double *c)
{
	c[nSizeTotal*nSizeTotal];

    cosTheta[nSize];
    weightsTheta[nSize];

	cout <<"Normalization test of wavefunctions"<<endl;
	for( int i = 0; i < nSizeTotal; i++ )
    {
        double sum = 0.;
        for( int j = 0; j < nSizeTotal; j++ )
        {   
            sum += c[i+j*nSizeTotal]*c[i+j*nSizeTotal];
        }
        cout <<" Psi "<<"["<<i<<"] = "<<sum<<endl;
    }

    vector<double> psiRe;
    for (int iXRotor   = 0; iXRotor<nSize; iXRotor++)
    {
        double xRotor  = cosTheta[iXRotor];
        double theta   = acos(xRotor);

        int indexBra   = 0;
        double sum     = 0.;

        for (int iBra  = 0; iBra<nSizeRot; iBra++)
        {
      		if (iBra%nSkip) continue;
            sum += c[indexBra*nSizeTotal]*normalizedPl0(iBra, xRotor)*sqrt(weightsTheta[iXRotor]);
            indexBra++;
        }
        psiRe.push_back(sum);
    }
    double sum = 0.;
    for (int i = 0; i < psiRe.size(); i++) sum +=  psiRe[i]*psiRe[i];
	cout <<"Normalization test of ground state wavefuntion :"<<endl;
    cout<<"Norm : "<<sum<<endl;
	psiRe.clear();
}

void densityCosTheta(int nSizeTotal, double *cosTheta, double *weightsTheta, double *c)
{
	c[nSizeTotal*nSizeTotal];

    cosTheta[nSize];
    weightsTheta[nSize];

   	vector<double> psiRe;
    for (int iXRotor   = 0; iXRotor<nSize; iXRotor++)
    {
        double xRotor  = cosTheta[iXRotor];
        double theta   = acos(xRotor);

        int indexBra   = 0;
        double sum     = 0.;

        for (int iBra  = 0; iBra<nSizeRot; iBra++)
        {
            if (iBra%nSkip) continue;
            sum += c[indexBra*nSizeTotal]*normalizedPl0(iBra, xRotor);
            indexBra++;
        }
		sum *= sqrt(weightsTheta[iXRotor]);
        psiRe.push_back(sum);
    }

    string fname = "Density-CosTheta.txt";

    ofstream myfile;
    myfile.open(fname.c_str());

	double sum = 0.0;
    for (int iXRotor   = 0; iXRotor<psiRe.size(); iXRotor++) 
	{
        double xRotor  = cosTheta[iXRotor];
        double theta   = acos(xRotor);
		double density = psiRe[iXRotor]*psiRe[iXRotor];
		sum += density;
			
		myfile<<iXRotor<<"      "<<xRotor<<"     "<<theta<<"    "<<density<<endl;
    }
	cout<<"Norm  "<<sum<<endl;
	myfile.close();
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

          	for (int iPhiRotor1    = 0; iPhiRotor1<2*nSize+1; iPhiRotor1++) {
				for (int iPhiRotor2     = 0; iPhiRotor2<2*nSize+1; iPhiRotor2++) {

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
