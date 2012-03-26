/****************
 *  Topology.c  *
 ****************/


/*====================================*
 *  Erstellt die topologischen Daten  *
 *====================================*/


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "vector2.h"
#include "vector3.h"
#include "constants.h"
#include "topology.h"


/*===================*
 *  Einfache Knoten  *
 *===================*/

unsigned int search_point(vector3 x, vector3 *P, unsigned int *pz);


void init_grid_SK(vector3 **P, unsigned int ***F, vector3 ***U, unsigned int p, unsigned int m, unsigned int *np, unsigned int *nf);


void refine_grid_SK(vector3 **P, unsigned int ***F, vector3 ***U, unsigned int p, unsigned int m, unsigned int M, unsigned int *np, unsigned int *nf);


unsigned int search_point(x,P,pz)
/* Suchfunktion zu gennet: Falls der Punkt x in
   der Punkteliste P vorhanden ist, liefert search_point den
   Index dieses Punktes, ansonsten wird x zur Liste hinzugefuegt
   und der entsprechende Index zurueckgegeben. */
vector3		x, *P;
unsigned int 	*pz;
{
unsigned int	k;

/* Suche in Punkteliste nach x */
for (k=0; (k < *pz) && (vector3_norm(vector3_sub(P[k],x)) > tol); k++);

if (k == *pz) P[(*pz)++] = x;  /* Punkt nicht in Punkteliste -> anhaengen */

return(k);	/* Index des Punktes x als Funktionsergebnis */
}


void init_grid_SK(P, F, U, p, m, np, nf)
/* Erstellt die Punkte- und Indexliste fuer Level 0. */
vector3 **P;                    /* Zeiger auf die Punkteliste          */
unsigned int ***F;              /* Zeiger auf die lokale Basisliste    */
vector3 ***U;                   /* Gitterpunkte                        */
unsigned int p;                 /* Anzahl der Parametergebiete         */
unsigned int m;                 /* 2^m*2^m Patches pro Parametergebiet */
unsigned int *np, *nf;          /* sizeof(P) bzw. sizeof(F)            */
{
    unsigned int n = 1 << m;    /* n*n Patches pro Parametergebiet     */
    unsigned int i;             /* Laufindizes fuer Parametergebiet    */

/* Speicherplatz allokieren: worst case */
    (*P) = (vector3 *) malloc(4 * p * sizeof(vector3));
    (*F) = (unsigned int **) malloc(p * sizeof(unsigned int *));

/* Eckpunkte der Parametergebiete bestimmen */
    *np = 0;
    *nf = p;
    for (i = 0; i < p; i++) {
        /* Bestimme die 4 Eckpunkte des Patches */
        (*F)[i] = (unsigned int *) malloc(4*sizeof(unsigned int));
	(*F)[i][0] = search_point(U[i][0][0],*P,np);
        (*F)[i][1] = search_point(U[i][0][n],*P,np);
        (*F)[i][2] = search_point(U[i][n][n],*P,np);
        (*F)[i][3] = search_point(U[i][n][0],*P,np);
    }
    return;
}


void refine_grid_SK(P, F, U, p, m, M, np, nf)
/* Erstellt die Punkte- und Indexliste fuer alle zusaetzlichen
   Gitterpunkte des Level m. */
vector3 **P;                    /* Zeiger auf die Punkteliste    */
unsigned int ***F;              /* Zeiger auf die Indexliste     */
vector3 ***U;                   /* Gitterpunkte                  */
unsigned int p;                 /* Anzahl der Parametergebiete   */
unsigned int m;                 /* aktuelles Level               */
unsigned int M;                 /* hoechstes Level               */
unsigned int *np, *nf;          /* sizeof(P) bzw. sizeof(F)      */
{
    unsigned int n = 1 << (m - 1);      /* n = 2^(m-1)                   */
    unsigned int S = 1 << (M - m);      /* Schrittweite zum nexten Punkt */
    signed int i1, i2, i3;      /* Laufindizes                   */
    unsigned int fz;            /* Patchzaehler                  */
    unsigned int kz;		/* Kantenzaehler        	       */
    unsigned int **K;           /* Kantenliste                         */

/* Speicherplatz allokieren: worst case */
    K = (unsigned int**) malloc(4*p*sizeof(unsigned int*));
    (*P) = (vector3 *) realloc(*P,p*(2*n+1)*(2*n+1)*sizeof(vector3));
    (*F) = (unsigned int **) realloc(*F, 4 * (*nf) * sizeof(unsigned int *));
    for (i1 = *nf; i1 < 4 * (*nf); i1++)
        (*F)[i1] = (unsigned int *) calloc(4, sizeof(unsigned int));

/* Kopiere altes Gitter in neues Gitter */
    fz = (p - 1) * n * n;
    for (i1 = p - 1; i1 >= 0; i1--) {
        for (i2 = n - 1; i2 >= 0; i2--) {
            for (i3 = n - 1; i3 >= 0; i3--) {
                (*F)[4 * fz + 2 * n * (2 * i2 + 1) + 2 * i3][3] = (*F)[fz + n * i2 + i3][3];
                (*F)[4 * fz + 2 * n * (2 * i2 + 1) + 2 * i3 + 1][2] = (*F)[fz + n * i2 + i3][2];
                (*F)[4 * fz + 4 * n * i2 + 2 * i3 + 1][1] = (*F)[fz + n * i2 + i3][1];
                (*F)[4 * fz + 4 * n * i2 + 2 * i3][0] = (*F)[fz + n * i2 + i3][0];
            }
        }
        fz -= n * n;
    }

/* Bestimme Verfeinerung */
*nf *= 4;
fz = kz = 0;
for (i1=0; i1<p; i1++)
{  
   /* Bestimme 1. Kante */
   for (i2=0; (i2 < kz) && (((*F)[fz][0] != K[i2][n+1]) || ((*F)[fz+2*n-1][1] != K[i2][0])); i2++);
   if (i2 == kz)	/* Kante noch nicht vorhanden -> unterteilen und neue Punkte merken */
   {  K[kz] = (unsigned int*) malloc((n+2)*sizeof(unsigned int));
      K[kz][0  ] = (*F)[fz      ][0];
      K[kz][n+1] = (*F)[fz+2*n-1][1];
      for (i3=0; i3<n; i3++)
      {  (*P)[*np] = U[i1][0][S*(2*i3+1)];
	 (*F)[fz+2*i3][1] = (*F)[fz+2*i3+1][0] = K[kz][i3+1] = (*np)++;
	 }
      kz++;
      }
   else 	/* Kante schon vorhanden */
   {  for (i3=0; i3<n; i3++) (*F)[fz+2*i3][1] = (*F)[fz+2*i3+1][0] = K[i2][n-i3];
      }

   /* Bestimme 2. Kante */
   for (i2=0; (i2 < kz) && (((*F)[fz+2*n-1][1] != K[i2][n+1]) || ((*F)[fz+4*n*n-1][2] != K[i2][0])); i2++);
   if (i2 == kz)	/* Kante noch nicht vorhanden -> unterteilen und neue Punkte merken */
   {  K[kz] = (unsigned int*) malloc((n+2)*sizeof(unsigned int));
      K[kz][0  ] = (*F)[fz+2*n  -1][1];
      K[kz][n+1] = (*F)[fz+4*n*n-1][2];
      for (i3=0; i3<n; i3++)
      {  (*P)[*np] = U[i1][S*(2*i3+1)][2*S*n];
	 (*F)[fz+2*n*(2*i3+1)-1][2] = (*F)[fz+4*n*(i3+1)-1][1] = K[kz][i3+1] = (*np)++;
	 }
      kz++;
      }
   else 	/* Kante schon vorhanden */
   {  for (i3=0; i3<n; i3++) (*F)[fz+2*n*(2*i3+1)-1][2] = (*F)[fz+4*n*(i3+1)-1][1] = K[i2][n-i3];
      }

   /* Bestimme 3. Kante */
   for (i2=0; (i2 < kz) && (((*F)[fz+4*n*n-1][2] != K[i2][n+1]) || ((*F)[fz+2*n*(2*n-1)][3] != K[i2][0])); i2++);
   if (i2 == kz)	/* Kante noch nicht vorhanden -> unterteilen und neue Punkte merken */
   {  K[kz] = (unsigned int*) malloc((n+2)*sizeof(unsigned int));
      K[kz][0  ] = (*F)[fz+4*n*n-1    ][2];
      K[kz][n+1] = (*F)[fz+2*n*(2*n-1)][3];
      for (i3=0; i3<n; i3++)
      {  (*P)[*np] = U[i1][2*S*n][S*(2*(n-i3)-1)];
         (*F)[fz+4*n*n-2*(i3+1)][2] = (*F)[fz+4*n*n-2*i3-1][3] = K[kz][i3+1] = (*np)++;
         }
      kz++;
      }
   else 	/* Kante schon vorhanden */
   {  for (i3=0; i3<n; i3++) (*F)[fz+4*n*n-2*(i3+1)][2] = (*F)[fz+4*n*n-2*i3-1][3] = K[i2][n-i3];
      }

   /* Bestimme 4. Kante */
   for (i2=0; (i2 < kz) && (((*F)[fz+2*n*(2*n-1)][3] != K[i2][n+1]) || ((*F)[fz][0] != K[i2][0])); i2++);
   if (i2 == kz)	/* Kante noch nicht vorhanden -> unterteilen und neue Punkte merken */
   {  K[kz] = (unsigned int*) malloc((n+2)*sizeof(unsigned int));
      K[kz][0  ] = (*F)[fz+2*n*(2*n-1)][3];
      K[kz][n+1] = (*F)[fz            ][0];
      for (i3=0; i3<n; i3++)
      {  (*P)[*np] = U[i1][S*(2*(n-i3)-1)][0];
         (*F)[fz+4*n*(n-i3)-2*n][0] = (*F)[fz+4*n*(n-i3-1)][3] = K[kz][i3+1] = (*np)++;
         }
      kz++;
      }
   else 	/* Kante schon vorhanden */
   {  for (i3=0; i3<n; i3++) (*F)[fz+4*n*(n-i3)-2*n][0] = (*F)[fz+4*n*(n-i3-1)][3] = K[i2][n-i3];
      }

   /* Verfeinere das Patch */
   for (i2=1; i2<2*n; i2++)
   {  for (i3=1; i3<2*n; i3++)
      {  if ((i2%2 == 1) || (i3%2 == 1))
         {  (*P)[*np] = U[i1][S*i2][S*i3];
	    (*F)[fz+2*n*i2+i3][0] = (*F)[fz+2*n*i2+i3-1][1] = (*F)[fz+2*n*(i2-1)+i3-1][2] = (*F)[fz+2*n*(i2-1)+i3][3] = (*np)++;
	    }
	 }
      }
   fz += 4*n*n;
   }

/* Speicherplatz wieder freigeben */
for (i1=0; i1<kz; i1++) free(K[i1]);
free(K);
return;
}


unsigned int gennet_SK(P, F, U, p, M)
/* Erstellt die Punkte- und Patchliste in hierarchischer Weise
   und liefert als Funktionsergebnis die Laenge der Punkteliste,
   die vom Geschlecht der Oberflaeche abhaengig ist. */
vector3 **P;                    /* Zeiger auf die Punkteliste          */
unsigned int ***F;              /* Zeiger auf die Patchliste           */
vector3 ***U;                   /* Gitterpunkte                        */
unsigned int p;                 /* Anzahl der Patches                  */
unsigned int M;                 /* 2^M*2^M Patches pro Parametergebiet */
{
    unsigned int m;             /* Laufindex fuer das Level            */
    unsigned int np;            /* Laenge von P                        */
    unsigned int nf;            /* Laenge von F                        */

    init_grid_SK(P, F, U, p, M, &np, &nf);
    for (m = 1; m <= M; m++)
        refine_grid_SK(P, F, U, p, m, M, &np, &nf);
    return (np);
}


/*===================*
 *  Doppelte Knoten  *
 *===================*/

void init_grid_DK(vector3 **P, unsigned int ***F, vector3 ***U, unsigned int p, unsigned int m, unsigned int *np, unsigned int *nf);


void refine_grid_DK(vector3 **P, unsigned int ***F, vector3 ***U, unsigned int p, unsigned int m, unsigned int M, unsigned int *np, unsigned int *nf);


void init_grid_DK(P, F, U, p, m, np, nf)
/* Erstellt die Punkte- und Indexliste fuer Level 0. */
vector3 **P;                    /* Zeiger auf die Punkteliste          */
unsigned int ***F;              /* Zeiger auf die lokale Basisliste    */
vector3 ***U;                   /* Gitterpunkte                        */
unsigned int p;                 /* Anzahl der Parametergebiete         */
unsigned int m;                 /* 2^m*2^m Patches pro Parametergebiet */
unsigned int *np, *nf;          /* sizeof(P) bzw. sizeof(F)            */
{
    unsigned int n = 1 << m;    /* n*n Patches pro Parametergebiet     */
    unsigned int i;             /* Laufindizes fuer Parametergebiet    */

/* Speicherplatz allokieren */
    (*P) = (vector3 *) malloc(4 * p * sizeof(vector3));
    (*F) = (unsigned int **) malloc(p * sizeof(unsigned int *));

/* Eckpunkte der Parametergebiete bestimmen */
    *np = 4 * p;
    *nf = p;
    for (i = 0; i < p; i++) {
        /* Bestimme die 4 Eckpunkte des Patches */
        (*P)[4 * i    ] = U[i][0][0];
        (*P)[4 * i + 1] = U[i][0][n];
        (*P)[4 * i + 2] = U[i][n][0];
        (*P)[4 * i + 3] = U[i][n][n];

        (*F)[i] = (unsigned int *) malloc(4 * sizeof(unsigned int));
        (*F)[i][0] = 4 * i;
        (*F)[i][1] = 4 * i + 1;
        (*F)[i][2] = 4 * i + 3;
        (*F)[i][3] = 4 * i + 2;
    }
    return;
}


void refine_grid_DK(P, F, U, p, m, M, np, nf)
/* Erstellt die Punkte- und Indexliste fuer alle zusaetzlichen
   Gitterpunkte des Level m. */
vector3 **P;                    /* Zeiger auf die Punkteliste    */
unsigned int ***F;              /* Zeiger auf die Indexliste     */
vector3 ***U;                   /* Gitterpunkte                  */
unsigned int p;                 /* Anzahl der Parametergebiete   */
unsigned int m;                 /* aktuelles Level               */
unsigned int M;                 /* hoechstes Level               */
unsigned int *np, *nf;          /* sizeof(P) bzw. sizeof(F)      */
{
    unsigned int n = 1 << (m - 1);      /* n = 2^(m-1)                   */
    unsigned int S = 1 << (M - m);      /* Schrittweite zum nexten Punkt */
    signed int i1, i2, i3;      /* Laufindizes                   */
    unsigned int fz;            /* Patchzaehler                  */

/* Speicherplatz allokieren */
    (*P) = (vector3 *) realloc(*P, 4 * (*np) * sizeof(vector3));
    (*F) = (unsigned int **) realloc(*F, 4 * (*nf) * sizeof(unsigned int *));
    for (i1 = *nf; i1 < 4 * (*nf); i1++)
        (*F)[i1] = (unsigned int *) calloc(4, sizeof(unsigned int));

/* Kopiere altes Gitter in neues Gitter */
    fz = (p - 1) * n * n;
    for (i1 = p - 1; i1 >= 0; i1--) {
        for (i2 = n - 1; i2 >= 0; i2--) {
            for (i3 = n - 1; i3 >= 0; i3--) {
                (*F)[4 * fz + 2 * n * (2 * i2 + 1) + 2 * i3][3] = (*F)[fz + n * i2 + i3][3];
                (*F)[4 * fz + 2 * n * (2 * i2 + 1) + 2 * i3 + 1][2] = (*F)[fz + n * i2 + i3][2];
                (*F)[4 * fz + 4 * n * i2 + 2 * i3 + 1][1] = (*F)[fz + n * i2 + i3][1];
                (*F)[4 * fz + 4 * n * i2 + 2 * i3][0] = (*F)[fz + n * i2 + i3][0];
            }
        }
        fz -= n * n;
    }

/* Bestimme Verfeinerung */
    fz = 0;
    *nf *= 4;
    for (i1 = 0; i1 < p; i1++) {
        for (i2 = 0; i2 <= 2 * n; i2++) {
            for (i3 = 0; i3 <= 2 * n; i3++) {
                if ((i2 % 2 == 1) || (i3 % 2 == 1)) {
                    (*P)[*np] = U[i1][i2 * S][i3 * S];
                    if ((i2 < 2 * n) && (i3 < 2 * n))
                        (*F)[fz + 2 * n * i2 + i3][0] = *np;
                    if ((i2 < 2 * n) && (0 < i3))
                        (*F)[fz + 2 * n * i2 + i3 - 1][1] = *np;
                    if ((0 < i2) && (0 < i3))
                        (*F)[fz + 2 * n * (i2 - 1) + i3 - 1][2] = *np;
                    if ((0 < i2) && (i3 < 2 * n))
                        (*F)[fz + 2 * n * (i2 - 1) + i3][3] = *np;
                    (*np)++;
                }
            }
        }
        fz += 4 * n * n;
    }

/* Speicherplatz wieder freigeben */
    return;
}


unsigned int gennet_DK(P, F, U, p, M)
/* Erstellt die Punkte- und Patchliste in hierarchischer Weise
   und liefert als Funktionsergebnis die Laenge der Punkteliste,
   die vom Geschlecht der Oberflaeche abhaengig ist. */
vector3 **P;                    /* Zeiger auf die Punkteliste          */
unsigned int ***F;              /* Zeiger auf die Patchliste           */
vector3 ***U;                   /* Gitterpunkte                        */
unsigned int p;                 /* Anzahl der Patches                  */
unsigned int M;                 /* 2^M*2^M Patches pro Parametergebiet */
{
    unsigned int m;             /* Laufindex fuer das Level            */
    unsigned int np;            /* Laenge von P                        */
    unsigned int nf;            /* Laenge von F                        */

    init_grid_DK(P, F, U, p, M, &np, &nf);
    for (m = 1; m <= M; m++)
        refine_grid_DK(P, F, U, p, m, M, &np, &nf);
    return (np);
}


/*=======================*
 *  Universelle Routine  *
 *=======================*/

void free_patchlist(F,nf)
/* gibt den Speicherplatz der (nf,4)-(unsigned int)-Patchliste F frei */
unsigned int	***F;
unsigned int	nf;
{
unsigned int	k;
for (k=0; k<nf; k++)  free((*F)[k]);
free(*F);
return;
}
