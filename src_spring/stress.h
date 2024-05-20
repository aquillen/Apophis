/**
 * @file        stress.h
 * @brief       springs 
 * @author      Alice Quilen 
 */


#ifndef _STRESS_H
#define _STRESS_H

#include"particle.h"
#include"spring.h"

struct stresstensor {  // for each node
        double sigxx; //stress tensor
        double sigyy;
        double sigzz;
        double sigxy;
        double sigyz;
        double sigxz;
        double eig1;  // eigenvalues of stress tensor eig1 is biggest
        double eig2;
        double eig3;
        double maxF;  // max Force of a spring
        int s_index;  // index of spring giving max Force
        int fail;  // is there material failure 
};


// in stress.c
extern struct stresstensor *stressvec; // global so can be reached by routines here
void update_stresstensor(struct reb_simulation* const r,double nodevol);
int markfailure(struct reb_simulation* const r, int npert,
       double pmass_div, double sigt_I, double sigt_S);   
void sfilename(struct reb_simulation* const r,char *root, double tp, char *fname);
void killsprings(struct reb_simulation* const r, int npert);
void spring_force_one(struct reb_simulation* const r,int i, 
   double *Fx, double *Fy, double *Fz, double *Lx, double *Ly, double *Lz);

#endif // _STRESS_H


