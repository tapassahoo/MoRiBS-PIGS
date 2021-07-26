#include <cmath>

//#include <cblas.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "io-xyz.h"
#include "xyz-water-utils.h"

#include "mbpol.h"
#include "test-mbpol.h"
#include "mc_setup.h"

void InitMBpol(void)
{
	try {
		ifstream ifs("/home/tapas/MoRiBS-PIGS/dimer.xyz");

		if (!ifs)
			throw runtime_error("could not open the XYZ file");

		//string comment;
		kit::io::load_xyz(ifs, comment, elements, crd);
	} catch (const exception& e) {
		cerr << " ** Error ** : " << e.what() << endl;
		exit(11);
	}

	if (!kit::is_water(elements, crd, false)) {
		cerr << " ** Error ** : not water?" << endl;
		exit(12);
	}

    //body-fixed coordinates//
    O1_init[0]=0.0;
    O1_init[1]=0.0;
    O1_init[2]=0.065621658;

    H1_init[0]=-0.75739503834488564;
    H1_init[1]=0.0;
    H1_init[2]=-0.52073166675198157;

    H2_init[0]=0.75739503834488564;
    H2_init[1]=0.0;
    H2_init[2]=-0.52073166675198157;
}
void calmbpoleng(double *Eulang1, double *Eulang2, double *com1, double *com2, double &energy)
{
	const double pi = 4.0*atan(1.0);
    const size_t nw = 2;//elements.size()/3;

    x2o::mbpol pot;

	//Some parameters defined for blas routines//
    char trans='N';
    int Ncart=3;
    int INCX=1.0;
    int INCY=1.0;

    double one=1.0;
    double zero=0.0;

	double phi1 = Eulang1[0];
	double theta1 = Eulang1[1];
	double chi1 = Eulang1[2];

	double phi2 = Eulang2[0];
	double theta2 = Eulang2[1];
	double chi2 = Eulang2[2];

	rotation_matrix(phi1,theta1,chi1,rotmat1);
	rotation_matrix(phi2,theta2,chi2,rotmat2);

	dgemv_(&trans,&Ncart,&Ncart,&one,rotmat1,&Ncart,O1_init,&INCX,&zero,O1_rot,&INCY);
	dgemv_(&trans,&Ncart,&Ncart,&one,rotmat2,&Ncart,O1_init,&INCX,&zero,O2_rot,&INCY); 
	O2_rot[2]=O2_rot[2]+com2[2];

	dgemv_(&trans,&Ncart,&Ncart,&one,rotmat1,&Ncart,H1_init,&INCX,&zero,H1_rot,&INCY);
	dgemv_(&trans,&Ncart,&Ncart,&one,rotmat2,&Ncart,H1_init,&INCX,&zero,H3_rot,&INCY); 
	H3_rot[2]=H3_rot[2]+com2[2];

	dgemv_(&trans,&Ncart,&Ncart,&one,rotmat1,&Ncart,H2_init,&INCX,&zero,H2_rot,&INCY);
	dgemv_(&trans,&Ncart,&Ncart,&one,rotmat2,&Ncart,H2_init,&INCX,&zero,H4_rot,&INCY); 
	H4_rot[2]=H4_rot[2]+com2[2];

	crd[0]=O1_rot[0];
	crd[1]=O1_rot[1];
	crd[2]=O1_rot[2];
	crd[3]=H1_rot[0];
	crd[4]=H1_rot[1];
	crd[5]=H1_rot[2];
	crd[6]=H2_rot[0];
	crd[7]=H2_rot[1];
	crd[8]=H2_rot[2];

	crd[9]=O2_rot[0];
	crd[10]=O2_rot[1];
	crd[11]=O2_rot[2];
	crd[12]=H3_rot[0];
	crd[13]=H3_rot[1];
	crd[14]=H3_rot[2];
	crd[15]=H4_rot[0];
	crd[16]=H4_rot[1];
	crd[17]=H4_rot[2];

	energy = pot(nw, &(crd[0]));
}

void rotation_matrix(double phi, double theta, double chi,double* matrix){

	double cp=cos(phi);
	double sp=sin(phi);
	double ct=cos(theta);
	double st=sin(theta);
	double ck=cos(chi);
	double sk=sin(chi);

	matrix[0]=cp*ct*ck-sp*sk;
	matrix[3]=-cp*ct*sk-sp*ck;
	matrix[6]=cp*st;
	matrix[1]=sp*ct*ck+cp*sk;
	matrix[4]=-sp*ct*sk+cp*ck;
	matrix[7]=sp*st;
	matrix[2]=-st*ck;
	matrix[5]=st*sk;
	matrix[8]=ct;
};
