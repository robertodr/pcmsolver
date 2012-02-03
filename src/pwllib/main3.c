/************
 *  main.c  *
 ************/


/*=================================================*
 *  Hauptprogramm zum Wavelet-Galerkin-Verfahren.  *
 *=================================================*/
 
 
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "sparse2.h"
#include "sparse.h"
#include "basis_pwl.h"
#include "dwt_pwl.h"
#include "WEM_pwl.h"
#include "WEMRHS_pwl.h"
#include "WEMPCG_pwl.h"
#include "WEMPGMRES_pwl.h"
#include "compression_pwl.h"
#include "interpolate_pwl.h"
#include "read_points.h"
#include "postproc_pwl.h"
#include "topology_pwl.h"
#include "precond_pwl.h"
#include "energy_pwl.h"
#include "volume.h"
#include "kern.h"


#if !defined pi
#define pi 3.1415926535897932385
#endif  


int main()
{
    sparse2         S_i, S_e;	/* komprimierte Steifigkeitsmatrix    */
    sparse		G;		/* Massenmatrix			      */
    vector3		*P;		/* Punkteliste der Einskalenbasis     */
    unsigned int	**F;		/* Elementliste der Einskalenbasis    */
    vector3		****T;		/* Oberflaecheninterpolation          */
    vector3		***U;		/* Knotenpunkte                       */
    unsigned int	nf;           	/* Laenge von F  	      	      */
    unsigned int	np;           	/* Laenge von P  	      	      */
    unsigned int	p;		/* Anzahl der Patches                 */
    unsigned int	M;	      	/* 2^M*2^M Elemente pro Patch         */
    element_pwl		*E;		/* hierarchische Elementliste         */
    wavelet_pwl		*W;		/* Liste der Wavelets                 */
    double		*u, *v;		/* approximierter Dichtevektor        */
    double		*rhs;		/* rechte Seite des Gleichungssystems */
    double		res;		/* berechnete Austauschenergien       */ 
    time_t		t1, t2, t3;	/* Sartzeit/Zwischenzeit/Endzeit      */
    unsigned int	i, j;		/* Laufindizes                        */
    double		det;		/* Determinante im anistropen Fall    */

    /*========================================================*/
    unsigned int	CASE = 3;	/* FLAG FOR CHOICE OF BIE */
    /*========================================================*/

    /* Ausgabe */
    switch (CASE)
        {  case 1:  printf("Pure Poisson Equation with 2nd kind integral equation:\n");
                printf("======================================================\n\n");
                break;
        case 2:  printf("Pure Poisson Equation with 1st kind integral equation:\n");
            printf("======================================================\n\n");
            break;
        case 3:  printf("Poisson-Boltzmann Equation:\n");
            printf("===========================\n\n");
            break;
        case 4:  printf("Anisotropic Dielectric:\n");
            printf("=======================\n\n");
            break;
        default: printf("ERROR: non-existent scheme\n"); return(0);
        }

    /* Initialisierung */
    read_points(&U,&p,&M);
    nf = p*(1<<M)*(1<<M);	        /* Anzahl der Patches */
    time(&t1);

    /* Topologie bestimmen */
    init_interpolate_pwl(&T,U,p,M);
    np = gennet_pwl(&P,&F,U,p,M);
    free_points(&U,p,M);
    volume(F,T,p,M);

    /* erstelle Element-/Waveletliste */
    printf("Number of levels:                %d \n",M);
    printf("Number of parameter domains:     %d \n",p);
    printf("Number of ansatz functions:      %d \n",np);
    printf("Computing the overhead:          ");
    generate_elementlist_pwl(&E,P,F,p,M);
    generate_waveletlist_pwl(&W,E,p,M,np);
    set_quadrature_level_pwl(W,E,p,M,np);
    simplify_waveletlist_pwl(W,E,p,M,np);
    complete_elementlist_pwl(W,E,p,M,np);
    time(&t2);
    printf("Computation time:                %g secs.\n\n",difftime(t2,t1));

    /* erstes Paar Steifigkeitsmatrizen aufstellen */
    printf("Computing the 1st pair of system matrices: \n");
    compression_pwl(&S_i,W,E,p,M,np);
    WEM_pwl(&S_i,W,P,E,T,p,M,SingleLayerInt,DoubleLayerInt,2*pi);
    postproc_pwl(&S_i,W,E,p,M);
    time(&t3);      /* Zwischenzeit */
    printf("Computation time:                %g secs.\n\n",difftime(t3,t2));

    /* zweites Paar Steifigkeitsmatrizen aufstellen */
    time(&t2);
    printf("Computing the 2nd pair of system matrices: \n");
    compression_pwl(&S_e,W,E,p,M,np);
    WEM_pwl(&S_e,W,P,E,T,p,M,SingleLayerExt,DoubleLayerExt,-2*pi);
    for (i=0; i<np; i++)  {  
        for (j=0; j<S_e.row_number[i]; j++) {
            S_e.value1[i][j] /= epsilon;
        }        
    }
    postproc_pwl(&S_e,W,E,p,M);
    time(&t3);      /* Zwischenzeit */
    printf("Computation time:                %g secs.\n\n",difftime(t3,t2));

    /* loese Gleichungssystem */
    u = (double*) calloc(np,sizeof(double));
    v = (double*) calloc(np,sizeof(double));
    WEMRHS_pwl2(&rhs,W,E,T,p,M,np);
    i = WEMPCG_pwl(&S_i,rhs,u,eps,W,F,p,M);		/* u = V_i^(-1)*N_f */
    printf("Solving the 1st linear system:   %d iterations\n",i);
    memset(rhs,0,np*sizeof(double));
    for (i=0; i<np; i++){			   	/* rhs = V_e*u */
        for (j=0; j<S_e.row_number[i]; j++) {
            rhs[i] += S_e.value1[i][j] * u[S_e.index[i][j]];
        }
    }
    i = WEMPGMRES_pwl3(&S_i,&S_e,rhs,v,eps,W,F,p,M);	/* solve complicated_system u = A^(-1)*rhs */ 
    printf("Solving the 2nd linear system:   %d iterations\n",i);   
    for (i=0; i<np; i++) u[i] -= 4*pi*v[i]; 	/* u = u - 4*pi*v */ 

    time(&t2);
    printf("Computation time:                %g secs.\n\n",difftime(t2,t3));

    /* Energie-Berechnung */
    tdwtLin(u,F,M,p,np);
    res = energy_pwl(u,F,T,p,M);
    time(&t2);
    printf("Over-all computation time:       %g secs.\n",difftime(t2,t1));

    /* Speicherplatz freigeben */
    printf("\n\n\n");
    if ((CASE == 3) || (CASE == 4)) free_sparse2(&S_e);
    free_waveletlist_pwl(&W,np);
    free_elementlist_pwl(&E,p,M);
    free_interpolate_pwl(&T,p,M);
    free_patchlist_pwl(&F,nf);
    free_sparse2(&S_i);
    free(rhs);
    free(P);
    free(u);
    free(v);
    return(0);
}