/*******************************************************************************
*  Copyright (C) 2009-2015 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
*
********************************************************************************
*/
/*
   ZHEEV Example.
   ==============

   Program computes all eigenvalues and eigenvectors of a complex Hermitian
   matrix A:

   (  9.14,  0.00) ( -4.37, -9.22) ( -1.98, -1.72) ( -8.96, -9.50)
   ( -4.37,  9.22) ( -3.35,  0.00) (  2.25, -9.51) (  2.57,  2.40)
   ( -1.98,  1.72) (  2.25,  9.51) ( -4.82,  0.00) ( -3.24,  2.04)
   ( -8.96,  9.50) (  2.57, -2.40) ( -3.24, -2.04) (  8.44,  0.00)

   Description.
   ============

   The routine computes all eigenvalues and, optionally, eigenvectors of an
   n-by-n complex Hermitian matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.

   Example Program Results.
   ========================

 ZHEEV Example Program Results

 Eigenvalues
 -16.00  -6.76   6.67  25.51

 Eigenvectors (stored columnwise)
 (  0.34,  0.00) ( -0.55,  0.00) (  0.31,  0.00) ( -0.70,  0.00)
 (  0.44, -0.54) (  0.26,  0.18) (  0.45,  0.29) (  0.22, -0.28)
 ( -0.48, -0.37) ( -0.52, -0.02) ( -0.05,  0.57) (  0.15,  0.08)
 (  0.10, -0.12) ( -0.50,  0.28) ( -0.23, -0.48) (  0.34, -0.49)
*/
#include <stdlib.h>
#include <stdio.h>

/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

/* ZHEEV prototype */
extern void zheev( char* jobz, char* uplo, int* n, dcomplex* a, int* lda,
                double* w, dcomplex* work, int* lwork, double* rwork, int* info );
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, dcomplex* a, int lda );
extern void print_rmatrix( char* desc, int m, int n, double* a, int lda );

/* Parameters */
#define N 4
#define LDA N

/* Main program */
int main() {
        /* Locals */
        int n = N, lda = LDA, info, lwork;
        dcomplex wkopt;
        dcomplex* work;
        /* Local arrays */
        /* rwork dimension should be at least max(1,3*n-2) */
        double w[N], rwork[3*N-2];
        dcomplex a[LDA*N] = {
           { 9.14,  0.00}, {-4.37,  9.22}, {-1.98,  1.72}, {-8.96,  9.50},
           { 0.00,  0.00}, {-3.35,  0.00}, { 2.25,  9.51}, { 2.57, -2.40},
           { 0.00,  0.00}, { 0.00,  0.00}, {-4.82,  0.00}, {-3.24, -2.04},
           { 0.00,  0.00}, { 0.00,  0.00}, { 0.00,  0.00}, { 8.44,  0.00}
        };
        /* Executable statements */
        printf( " ZHEEV Example Program Results\n" );
        /* Query and allocate the optimal workspace */
        lwork = -1;
        zheev( "Vectors", "Lower", &n, a, &lda, w, &wkopt, &lwork, rwork, &info );
        lwork = (int)wkopt.re;
        work = (dcomplex*)malloc( lwork*sizeof(dcomplex) );
        /* Solve eigenproblem */
        zheev( "Vectors", "Lower", &n, a, &lda, w, work, &lwork, rwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Print eigenvalues */
        print_rmatrix( "Eigenvalues", 1, n, w, 1 );
        /* Print eigenvectors */
        print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
        /* Free workspace */
        free( (void*)work );
        exit( 0 );
} /* End of ZHEEV Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, dcomplex* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", a[i+j*lda].re, a[i+j*lda].im );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }
}
