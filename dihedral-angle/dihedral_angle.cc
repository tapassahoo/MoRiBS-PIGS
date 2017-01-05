#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

double sqr(double x);
void VectorNormalisation(double *v);
double DotProduct(double *v, double *w);
void CrossProduct(double *v, double *w, double *cross);

#define PI 3.14159265

int main()
{
	//First point
	double x1   = 0.000;
	double y1   = 1.000; 
	double z1   = 0.000;

	//Second point
    double x2   = 0.000; 
    double y2   = 0.000;
    double z2   = 0.000;

	//Third point
    double x3   = 1.000; 
    double y3   = 0.000;
    double z3   = 0.000; 

	//Fourth point
    double x4   = 1.000;
    double y4   = 0.000;
    double z4   = -20000.000;

	//b1
	double b_a[3];
	b_a[0] 		= -(x1 - x2);
	b_a[1] 		= -(y1 - y2);
	b_a[2] 		= -(z1 - z2);

	//b2
	double b_c[3];
	b_c[0] 		= x2 - x3;
	b_c[1] 		= y2 - y3;
	b_c[2] 		= z2 - z3;

	//b3
	double c_d[3];
	c_d[0] 		= x4 - x3;
	c_d[1] 		= y4 - y3;
	c_d[2] 		= z4 - z3;

	double n1[3];
	double n2[3];
	double m[3];
//	VectorNormalisation(b_c);
//	VectorNormalisation(b_a);
//	VectorNormalisation(c_d);

	CrossProduct(b_a, b_c, n1);
	CrossProduct(b_c, c_d, n2);
	CrossProduct(n1, b_c, m);

	double x = DotProduct(n1, n2);
	double y = DotProduct(m, n2);

	double angle = 180.0 / PI * atan2(y, x);
    cout<<angle<<endl;
}

double sqr(double x){ return x*x; }

void VectorNormalisation(double *v) { double lenght = sqrt(sqr(v[0]) + sqr(v[1] + sqr(v[2]))); v[0] /= lenght; v[1] /= lenght; v[2] /= lenght; }

double DotProduct(double *v, double *w) { return (v[0] * w[0] + v[1] * w[1] + v[2] * w[2]); }

void CrossProduct(double *v, double *w, double *cross) 
{
	cross[0] = w[1] * v[2] - w[2] * v[1];
	cross[1] = w[2] * v[0] - w[0] * v[2];
	cross[2] = w[0] * v[1] - w[1] * v[0];
}
