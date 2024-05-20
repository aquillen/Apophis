/**
 * @file        quat.h
 * @brief       quat
 * @author      Alice Quilen 
 */


#ifndef _QUAT_H
#define _QUAT_H

#include"particle.h"
#include"spring.h"


// a structure to hold quaturnion 
struct quaturnion {
   double w;
   double x;
   double y;
   double z; 
};

// in quat.c

void rotate_vec_quaturnion(double x,double y,double z,
   struct quaturnion q, double *xr, double *yr, double *zr);
void q_normalize(struct quaturnion* q);
struct quaturnion uv_rotate(double ux, double uy, double uz, double vx, double vy, double vz);
void q_rotate_body(struct reb_simulation* const r, int il, int ih,struct  quaturnion q);


#endif // _QUAT_H


