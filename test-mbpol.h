#ifndef _TEST_MBPOL_H
#define _TEST_MBPOL_H 1

#include "io-xyz.h"
#include "xyz-water-utils.h"

#include "mbpol.h"

void rotation_matrix(double phi, double theta,double chi,double* matrix);
extern "C" {extern void dgemv_(char*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*); }
void InitMBpol(void);
void calmbpoleng(double *Eulang1, double *Eulang2, double *com1, double *com2, double &energy);

#endif  //mc_pimc.h
