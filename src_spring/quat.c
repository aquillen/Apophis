
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "rebound.h"
#include "tools.h"
#include "output.h"

#include"particle.h"
#include"spring.h"
#include"quat.h"


// rotate a vector via a quaturnion q
void rotate_vec_quaturnion(double x,double y,double z, 
   struct quaturnion q, double *xr, double *yr, double *zr)
{
   
   *xr =(1-2*q.y*q.y - 2*q.z*q.z)*x +   (2*q.x*q.y - 2*q.z*q.w)*y + (2*q.x*q.z + 2*q.y*q.w)*z;
   *yr =  (2*q.x*q.y + 2*q.z*q.w)*x + (1-2*q.x*q.x - 2*q.z*q.z)*y + (2*q.y*q.z - 2*q.x*q.w)*z;
   *zr =  (2*q.x*q.z - 2*q.y*q.w)*x +   (2*q.y*q.z + 2*q.x*q.w)*y + (1-2*q.x*q.x-2*q.y*q.y)*z;

}

// normalize a quaturnion 
void q_normalize(struct quaturnion* q)
{
   double u = sqrt(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
   q->w /= u;
   q->x /= u;
   q->y /= u;
   q->z /= u;
}

// get a quaturnion to describe the shortest arc rotation
// between a vector u and a vector v 
// https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
struct quaturnion uv_rotate(double ux, double uy, double uz, double vx, double vy, double vz)
{
   double udotv = ux*vx + uy*vy + uz*vz;
   double u2    = ux*ux + uy*uy + uz*uz;
   double v2    = vx*vx + vy*vy + vz*vz;
   struct quaturnion q;
   q.w = udotv + sqrt(u2*v2);
   q.x = uy*vz - uz*vy; // cross product
   q.y = uz*vx - ux*vz;
   q.z = ux*vy - uy*vx;
   q_normalize(&q);
   return q;
 
  // recipe:
  // q.w   = dot(u, v) + sqrt(length_2(u) * length_2(v))
  // q.xyz = cross(u, v)
  // then normalize 
  // does not yet take into account possibility of u = -v !!!!!!!!
  //   in this case q.w is 0 and you need to find any vector orthogonal to u to put in q.xyz 
}


// rotate extended body using a quaturnion q
//     particles with indices in range [il,ih) 
// rotate positions and velocity using a quaturnion q
// don't change center of mass position or center of mass velocity 
void q_rotate_body(struct reb_simulation* const r, int il, int ih, struct quaturnion q)
{
   struct reb_particle* particles = r->particles;
   double xc,yc,zc;
   compute_com(r,il, ih, &xc, &yc, &zc);
   double vxc,vyc,vzc;
   compute_cov(r,il, ih, &vxc, &vyc, &vzc);
   double xr,yr,zr;
   double vxr,vyr,vzr;
   for(int i=il;i<ih;i++){
      double x0 = particles[i].x - xc;
      double y0 = particles[i].y - yc;
      double z0 = particles[i].z - zc;
      rotate_vec_quaturnion(x0,y0,z0, q, &xr,&yr,&zr);
      particles[i].x = xr + xc;
      particles[i].y = yr + yc;
      particles[i].z = zr + zc;

      double vx0 = particles[i].vx - vxc;
      double vy0 = particles[i].vy - vyc;
      double vz0 = particles[i].vz - vzc;
      rotate_vec_quaturnion(vx0,vy0,vz0, q, &vxr,&vyr,&vzr);
      particles[i].vx = vxr + vxc;
      particles[i].vy = vyr + vyc;
      particles[i].vz = vzr + vzc;
   }
}

