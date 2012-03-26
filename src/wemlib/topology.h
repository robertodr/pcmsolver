#ifndef TOPOLOGY
#define TOPOLOGY
/****************
 *  Topology.h  *
 ****************/
#ifdef __cplusplus
extern "C" {
#endif


/*====================================*
 *  Erstellt die topologischen Daten  *
 *====================================*/


    unsigned int gennet_DK(vector3 **P, unsigned int ***F, vector3 ***T, unsigned int p, unsigned int m);
    unsigned int gennet_SK(vector3 **P, unsigned int ***F, vector3 ***T, unsigned int p, unsigned int m);
/* berechnet Punkt- und Patchliste */


    void free_patchlist(unsigned int ***F, unsigned int nf);
/* gibt den Speicherplatz der Patchliste frei */

#ifdef __cplusplus
}
#endif
#endif
