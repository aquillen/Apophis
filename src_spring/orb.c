
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"

#include "spring.h"
#include "kepcart.h"
#include "orb.h"


// extern int NS; // number of springs
// extern int NPERT;  // number of external perturbing bodies


// make a PlutoCharon binary with two masses point m1,m2 
// add those two masses 
// m1 = mp and m2 = mratio*mp
// masses are separated by distance sep, circular orbit for binary
// center of mass of PC and velocity set to be in orbit with given orbital elements
// return mean motion around PC
// r1 is radius for display of Pluto
// orbital elements are aa,ee,ii, longnode, argperi, meananom
// these refer to orbit of resolved body [0,N) about binary
// the resolved body is assumed to be alone and at origin
double add_pluto_charon_kep(struct reb_simulation* const r,double mp, 
    double mratio, double sep,
    double aa,double ee, double ii,double longnode,double argperi,double meananom, 
    double r1)
{
   // const int i1= r->N;
   // const int i2= i1+1;
   struct reb_particle pt;
   double m1 = mp; // 
   pt.ax = 0.0; pt.ay = 0.0; pt.az = 0.0;
   pt.vx = 0.0; pt.vy = 0.0; pt.vz = 0.0;
   pt.x  = 0.0; pt.y  = 0.0; pt.z  = 0.0;
   double m2 = mratio*m1;
   double mtot  = m1 + m2;
   double vc_bin = sqrt(r->G*mtot/sep);
   double x1 = sep*m2/mtot;  // initial orientation of binary is here
   double x2 = -sep*m1/mtot;
   double vy1 =  vc_bin*m2/mtot;
   double vy2 = -vc_bin*m1/mtot;
   double r2 = r1*pow(m2/m1,0.333333); // display radius
   if (r2 < 1.0) r2=1.0;
   //
   double GM = r->G*(mtot+1.0); // +1 here for extended body
   OrbitalElements orbel;
   orbel.a = aa;
   orbel.e = ee;
   orbel.i = ii;
   orbel.argperi = argperi;
   orbel.longnode= longnode;
   orbel.meananom= meananom;
   PhaseState ps;
   cartesian(GM, orbel, &ps);
   double om_orb = sqrt(GM/aa)/aa;

   // ratio of reduced mass to binary mass
   double muBratio = mratio/pow(1.0 + mratio,2.0); 
   double fac = 1.0 + (3.0/8.0)*muBratio*pow(sep/aa,2.0);
   om_orb *= fac; // correct mean motion by binary quad moment
   
   pt.x = ps.x + x2; pt.y = ps.y + 0.0; pt.z = ps.z + 0.0;
   pt.vx = ps.xd + 0.0; pt.vy = ps.yd + vy2; pt.vz = ps.zd + 0.0; 
   pt.m = m2; pt.r = r2;
   reb_simulation_add(r,pt);   // add Charon

   pt.x = ps.x + x1; pt.y = ps.y + 0.0; pt.z = ps.z + 0.0;
   pt.vx = ps.xd + 0.0; pt.vy = ps.yd + vy1; pt.vz = ps.zd + 0.0; 
   pt.m = m1; pt.r = r1;
   reb_simulation_add(r,pt);   // Pluto is last
   return om_orb;
}

/////////////////////////////////////////////////////
// add a point mass 
// the point mass added has mass m1, radius r1 (for display)
// if ip>=0 put point mass in orbit about another point mass with index ip
// if ip==-1 
//    the extended mass has indices [il,ih)
//    add a single point mass with extended body orbiting it
//    the pt mass added is at origin
//    if extended body is rotating then it still rotates
//       however if orbit is tilted then obliquity will not be the same
// return mean motion
// positions and velocity computed with orbelements
// orbels computed with GM taking into account both masses
/////////////////////////////////////////////////////
double add_pt_mass_kep(struct reb_simulation* const r, 
    int il, int ih, int ip, double m1, double r1,
    double aa,double ee, double ii,
    double longnode,double argperi,double meananom)  
{
   struct reb_particle* particles = r->particles;
   struct reb_particle pt;
   pt.m = m1;
   pt.r = r1;
   pt.ax = 0; pt.ay = 0; pt.az = 0;
   double GM,m0;
   double  x0=0.0;  double y0=0.0; double  z0=0.0;
   double vx0=0.0; double vy0=0.0; double vz0=0.0;
   if (ip<0){ // moving resolved body
     m0 = sum_mass(r,il,ih);
     subtractcom(r,il,ih);  // zero center of mass position of extended body
     subtractcov(r,il,ih);  // subtract center of mass velocity of extended body
   }
   else { // new particle has motion w.r.t to particle at ip
     m0 = particles[ip].m;
     x0 = particles[ip].x;    y0 = particles[ip].y;   z0 = particles[ip].z;
     vx0 = particles[ip].vx; vy0 = particles[ip].vy; vz0 = particles[ip].vz;
   }
   GM = r->G*(m1+m0); 
   OrbitalElements orbel;
   orbel.a = aa;
   orbel.e = ee;
   orbel.i = ii;
   orbel.argperi = argperi;
   orbel.longnode= longnode;
   orbel.meananom= meananom;
   PhaseState ps;
   cartesian(GM, orbel, &ps);
   double om_orb = sqrt(GM/aa)/aa;

   if (ip<0){  // move position of extended body
     pt.x   = 0.0; pt.y   = 0.0; pt.z   = 0.0;
     pt.vx  = 0.0; pt.vy  = 0.0; pt.vz  = 0.0; // new particle is at origin
     move_resolved(r,ps.x,ps.y,ps.z,ps.xd,ps.yd,ps.zd, il,ih);
   }
   else { // position and velocity of new particle with respected to mass at index ip
     pt.x       = ps.x + x0;
     pt.y       = ps.y + y0;
     pt.z       = ps.z + z0;
     pt.vx      = ps.xd + vx0;
     pt.vy      = ps.yd + vy0;
     pt.vz      = ps.zd + vz0;
   }
   reb_simulation_add(r,pt); // add the new particle
   return om_orb; // return mean motion
}

////////////////////////////////////////////////////////////////
// return mean motion, set up a single point mass perturber in orbit
// m1 is mass of perturber
// mball is mass of extended body  [0,N)
// use orbital elements to add m1 to particle list
// origin is extended body
double add_one_mass_kep(struct reb_simulation* const r, double m1,
    double aa,double ee, double ii,
    double longnode,double argperi,double meananom, 
    double mball, double r1)
{
   struct reb_particle pt;
   double GM = r->G*(m1+mball); // mball here is extended body
   OrbitalElements orbel;
   orbel.a = aa;
   orbel.e = ee;
   orbel.i = ii;
   orbel.argperi = argperi;
   orbel.longnode= longnode;
   orbel.meananom= meananom;
   PhaseState ps;
   cartesian(GM, orbel, &ps);
   double om_orb = sqrt(GM/aa)/aa;

   pt.m       = m1;
   pt.x       = ps.x;
   pt.y       = ps.y;
   pt.z       = ps.z;
   pt.vx      = ps.xd;
   pt.vy      = ps.yd;
   pt.vz      = ps.zd;
   pt.r       = r1;
   pt.ax = 0; pt.ay = 0; pt.az = 0;
   reb_simulation_add(r,pt);
   return om_orb;   
}


//////////////////////////////////////////////
// helper routine for drifting two masses (a binary)
//     (a binary) using Beauge et al. 06's
// formulae, Beauge, Michtchenko, Ferraz-Mello 2006, MNRAS 365, 1160
// we can write evolution as dr^2/dt^2 = -v/taua - (v-v_c)/taue
// inv_taua and inv_taue are in units of 1/time 
// inv_taua = (1/a) da/dt, inv_taue = (1/e) de/dt
//  if inv_tau_a >0 drifting apart !!!!!!!!!!!!!!
//
// cc>0 and alpha=0 corresponds to no eccentricity change
// with taua^-1 = 2C(1-alpha), taue^-1 = C*alpha from Beauge et al.
// kk = taua/taue
// eqns of motion can be written
// dv/dt = - C*(1-alpha) v - E(v-vc)
//       = -invtaua*v/2 - invtaue(v-vc)
// 
// This routine drifts a binary!!!!
// The drift is done trying to maintain momentum of center of mass of binary
// this routine will probably work for three or more body system
// and will probably work for binaries with various mass ratios
// this routine returns velocity changes dvx,dvy,dvz
// this routine needs x,y,z,vx,vy,vz which are the differences for the two masses
//    double x  = particles[im2].x - particles[im1].x; 
// 
//    particles[im1].vx -=  m2*dvx/(m1+m2);   etc
//    particles[im2].vx +=  m1*dvx/(m1+m2); 
// ACQ May 2021
/////////////////////////////////////////////////////////////
void dodrift_dv(double tstep,double inv_taua, double inv_taue,  double GMM, double x, double y, double z, double vx, double vy, double vz, double *dvx, double *dvy, double *dvz)
{

      double rdotv = x*vx + y*vy + z*vz;
      double rad = sqrt(x*x+y*y+z*z);
      double r2 = rad*rad;
      double vc = sqrt(GMM/rad); // circular velocity
      //  rxL is  vector in plane of orbit, perp to r, 
         //  direction of rotation
      // r x v x r = r x -L
      // vector identity axbxc = (adotc)b-(adotb)c
      double rcrossl_x = r2*vx - rdotv*x; 
      double rcrossl_y = r2*vy - rdotv*y; 
      double rcrossl_z = r2*vz - rdotv*z;
      double vl = sqrt(rcrossl_x*rcrossl_x+ rcrossl_y*rcrossl_y+ 
                       rcrossl_z*rcrossl_z); // length of rcrossl
      rcrossl_x /= vl; rcrossl_y /= vl; rcrossl_z /= vl; // unit vector now
      double vcvec_x = vc*rcrossl_x; 
      double vcvec_y = vc*rcrossl_y; 
      double vcvec_z = vc*rcrossl_z;
      // difference between velocity and vc
      double dd_vc_x = vx - vcvec_x;
      double dd_vc_y = vy - vcvec_y;
      double dd_vc_z = vz - vcvec_z;

// compute changes in velocity 
      *dvx =  tstep*(vx*inv_taua/2.0 + dd_vc_x*inv_taue);
      *dvy =  tstep*(vy*inv_taua/2.0 + dd_vc_y*inv_taue);
      *dvz =  tstep*(vz*inv_taua/2.0 + dd_vc_z*inv_taue);

}

//////////////////////////////////////////////
// drift the semi-major axis of two masses
//     (a binary) using Beauge et al. 06's recipe
// This routine drifts a binary!!!!
// if il1 != -1:
//    resolved body 1 contains particles [il1,ih1)
// else ih1 is index of point mass 1
// if il2 != -1:
//    resolved body 2 contains particles [il2,ih2)
// else ih2 is index of point mass 2
//////////////////////////////////////////////
void dodrift_twobodies(struct reb_simulation* const r, double tstep, double inv_taua, double inv_taue, int il1, int ih1, int il2, int ih2)  
{
   printf("hello\n");
}


//////////////////////////////////////////////
// drift the semi-major axis of two point masses
//     (a binary) using Beauge et al. 06's
// This routine drifts a binary!!!!
// The two point masses have indexes: im1, im2 
/////////////////////////////////////////////////////////////
void dodrift_bin(struct reb_simulation* const r, double tstep, double inv_taua, double inv_taue, int im1, int im2)  
{
   struct reb_particle* particles = r->particles;

   // int im1 = r->N -1; // index for primary perturber
   // int im2 = r->N -2; // index for member of binary 
   double m1 = particles[im1].m;
   double m2 = particles[im2].m;
   double GMM = r->G *(m1+m2);
   
      double x  = particles[im2].x - particles[im1].x; 
      double y  = particles[im2].y - particles[im1].y; 
      double z  = particles[im2].z - particles[im1].z;
      double vx = particles[im2].vx -particles[im1].vx; 
      double vy = particles[im2].vy -particles[im1].vy; 
      double vz = particles[im2].vz -particles[im1].vz;

      double dvx,dvy,dvz; // compute changes in velocity
      dodrift_dv(tstep,inv_taua, inv_taue, GMM, x,y,z,vx,vy,vz, 
              &dvx, &dvy, &dvz);

// update velocities , in such a way as to conserve
// momentum of binary, lower mass is affected much more than higher mass
      particles[im1].vx -=  m2*dvx/(m1+m2); 
      particles[im1].vy -=  m2*dvy/(m1+m2);
      particles[im1].vz -=  m2*dvz/(m1+m2);
      particles[im2].vx +=  m1*dvx/(m1+m2); 
      particles[im2].vy +=  m1*dvy/(m1+m2);
      particles[im2].vz +=  m1*dvz/(m1+m2);
    
}


//////////////////////////////////////////////
// just like the previous routine: dodrift_bin
//    but one mass is a resolved body!  the other is a point mass, indexed with im1
//    resolved body contains particles [il,ih)
/////////////////////////////////////////////////////////////
void dodrift_res(struct reb_simulation* const r, double tstep, double inv_taua, double inv_taue, int im1, int il, int ih)
{
   struct reb_particle* particles = r->particles;

   double m1 = particles[im1].m;
   double m2 = sum_mass(r,il,ih);
   double xc =0.0; double yc =0.0; double zc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc);
   double vxc =0.0; double vyc =0.0; double vzc =0.0;
   compute_cov(r,il, ih, &vxc, &vyc, &vzc);


   double GMM = r->G *(m1+m2);
   
      double x  = xc               - particles[im1].x; 
      double y  = yc               - particles[im1].y; 
      double z  = zc               - particles[im1].z;
      double vx = vxc              - particles[im1].vx; 
      double vy = vyc              - particles[im1].vy; 
      double vz = vzc              - particles[im1].vz;

      double dvx,dvy,dvz; // compute changes in velocity
      dodrift_dv(tstep,inv_taua, inv_taue, GMM, x,y,z,vx,vy,vz, 
              &dvx, &dvy, &dvz);

// update velocities , in such a way as to conserve
// momentum of binary, lower mass is affected much more than higher mass
      particles[im1].vx -= m2*dvx/(m1+m2); 
      particles[im1].vy -= m2*dvy/(m1+m2);
      particles[im1].vz -= m2*dvz/(m1+m2);
      double cvx         = m1*dvx/(m1+m2);  // changes to velocity of resolved body
      double cvy         = m1*dvy/(m1+m2); 
      double cvz         = m1*dvz/(m1+m2); 
      move_resolved(r,0.0,0.0,0.0,cvx,cvy,cvz, il,ih); // change only velocities
    
}


//////////////////////////////////////////////
// just like the previous routine: dodrift_bin and dodrift_res
//    but both masses are resolved bodes! 
//    resolved body 1 contains particles [il1,ih1)
//    resolved body 2 contains particles [il2,ih2)
/////////////////////////////////////////////////////////////
void dodrift_twores(struct reb_simulation* const r, double tstep, double inv_taua, double inv_taue,  int il1, int ih1, int il2, int ih2)
{
   // struct reb_particle* particles = r->particles;

   double m1 = sum_mass(r,il1,ih1);
   double m2 = sum_mass(r,il2,ih2);
   double xc1 =0.0; double yc1 =0.0; double zc1 =0.0;
   double xc2 =0.0; double yc2 =0.0; double zc2 =0.0;
   compute_com(r,il1, ih1, &xc1, &yc1, &zc1);
   compute_com(r,il2, ih2, &xc2, &yc2, &zc2);
   double vxc1 =0.0; double vyc1 =0.0; double vzc1 =0.0;
   double vxc2 =0.0; double vyc2 =0.0; double vzc2 =0.0;
   compute_cov(r,il1, ih1, &vxc1, &vyc1, &vzc1);
   compute_cov(r,il2, ih2, &vxc2, &vyc2, &vzc2);

   double GMM = r->G *(m1+m2);
   
      double x  = xc2 - xc1; 
      double y  = yc2 - yc1; 
      double z  = zc2 - zc1; 
      double vx = vxc2 - vxc1; 
      double vy = vyc2 - vyc1; 
      double vz = vzc2 - vzc1; 

      double dvx,dvy,dvz; // compute changes in velocity
      dodrift_dv(tstep,inv_taua, inv_taue, GMM, x,y,z,vx,vy,vz, 
              &dvx, &dvy, &dvz);

// update velocities , in such a way as to conserve
// momentum of binary, lower mass is affected much more than higher mass
      double cvx1 = -1.0*m2*dvx/(m1+m2);  // changes to velocity of resolved body1
      double cvy1 = -1.0*m2*dvy/(m1+m2);
      double cvz1 = -1.0*m2*dvz/(m1+m2);
      double cvx2 =      m1*dvx/(m1+m2);  // changes to velocity of resolved body2
      double cvy2 =      m1*dvy/(m1+m2); 
      double cvz2 =      m1*dvz/(m1+m2); 
      move_resolved(r,0.0,0.0,0.0,cvx1,cvy1,cvz1, il1,ih1); // change only velocities
      move_resolved(r,0.0,0.0,0.0,cvx2,cvy2,cvz2, il2,ih2); // change only velocities
    
}

///////////////////////////////////////////////////////////////

// add a single extra mass  of mass m1
// use cartesian coordinates to add m1 to particle list
void add_one_mass_cartesian(struct reb_simulation* const r, double m1, double r1, double x,double y, double z, double vx, double vy, double vz)
{
   struct reb_particle pt;
   pt.m       = m1;
   pt.x       = x;
   pt.y       = y;
   pt.z       = z;
   pt.vx      = vx;
   pt.vy      = vy;
   pt.vz      = vz;
   pt.r       = r1;
   pt.ax = 0; pt.ay = 0; pt.az = 0;
   reb_simulation_add(r,pt);
}


#define RSOFT 0.01

// apply quadrupole force of particle with index ic onto all particles 
// by adding to particle accelerations
// this takes into account quadrupole of the planet, using its mass and radius 
// J2p is unitless 
// Rplusp is radius of planet
// ic is the array index of the planet with the J2 (and where you find its mass)
// see https://en.wikipedia.org/wiki/Geopotential_model for formulas, I checked, seems ok
// taking into account orientation of planet's north pole
// pole at (sin thetap cos phip, sin thetap sin phip, cos thetap)
//            with thetap=0 pole is along + z axis
void quadJ2pole(struct reb_simulation* const r,double J2p, double Rplusp,double phip, double thetap, int ic)
{
   if (ic >= r->N) return;
   if (J2p == 0.0) return;
   double Gmc = r->G * r->particles[ic].m;
   double cst = Gmc*J2p*pow(Rplusp,2.0);
   double xhat = sin(thetap)*cos(phip);
   double yhat = sin(thetap)*sin(phip);
   double zhat = cos(thetap);
   for(int i =0;i < r->N;i++){
     if (i != ic){
        double dx = r->particles[i].x - r->particles[ic].x;
        double dy = r->particles[i].y - r->particles[ic].y;
        double dz = r->particles[i].z - r->particles[ic].z;
        double rr = sqrt(dx*dx + dy*dy + dz*dz);
        double r5 = pow(rr + RSOFT,5.0) ; // arbitrary softening here
        double cos2 = pow((dx *xhat + dy*yhat + dz*zhat)/rr,2.0);
        double Fx = cst*dx*1.5/r5 *(5.0*cos2 -1.0); // these are accelerations not forces
        double Fy = cst*dy*1.5/r5 *(5.0*cos2 -1.0); 
        double Fz = cst*dz*1.5/r5 *(5.0*cos2 -3.0);
        // double m = r->particles[i].m;
        r->particles[i].ax += Fx; // 
        r->particles[i].ay += Fy;
        r->particles[i].az += Fz;
        double mfac = r->particles[i].m/r->particles[ic].m;
        r->particles[ic].ax -= Fx*mfac; // opposite force
        r->particles[ic].ay -= Fy*mfac;
        r->particles[ic].az -= Fz*mfac;
     }
   } 
}  

// forces can also be written as
// Fx = GM J2 x 3/2  (5 z^2/r^2 - 1)
// Fy = GM J2 y 3/2  (5 z^2/r^2 - 1)
// Fz = GM J2 z 3/2  (5 z^2/r^2 - 3)
// here I replaced z^2 with dot product of pole direction with vector

//

// put in a C22,S22 quadrupole term, no rotation angle specified , so w.r.t underlying coord system
void quadC22S22(struct reb_simulation* const r,double C22, double S22, double Rplusp, int ic)
{
   if (ic >= r->N) return;
   if (C22*S22 == 0.0) return;
   double Gmc = r->G * r->particles[ic].m;
   double Cst = Gmc*C22*pow(Rplusp,2.0)*3.0;
   double Sst = Gmc*S22*pow(Rplusp,2.0)*6.0;
   for(int i =0;i < r->N;i++){
     if (i != ic){
        double dx = r->particles[i].x - r->particles[ic].x;
        double dy = r->particles[i].y - r->particles[ic].y;
        double dz = r->particles[i].z - r->particles[ic].z;
        double x2 = dx*dx; double y2 = dy*dy; double z2 = dz*dz;
        double rr = sqrt(x2 + y2 + z2 + RSOFT*RSOFT);// softening here
        double r7 = pow(rr,7.0); 
        double Fx =(Cst*dx/r7)*(-3.0*x2 +7.0*y2 +2*z2); // these are accelerations not forces
        Fx +=  (Sst*dy/r7)*(-4*x2 + y2 + z2);
        double Fy =(Cst*dy/r7)*(-7.0*x2 +3.0*y2 -2*z2); 
        Fy +=  (Sst*dx/r7)*(x2 - 4*y2 + z2);
        double Fz =(Cst*dz/r7)*-5*(x2 -y2); 
        Fz +=  (Sst*dz/r7)*(-5*dx*dy);
        // double m = r->particles[i].m;
        r->particles[i].ax += Fx; // 
        r->particles[i].ay += Fy;
        r->particles[i].az += Fz;
        double mfac = r->particles[i].m/r->particles[ic].m;
        r->particles[ic].ax -= Fx*mfac; // opposite force
        r->particles[ic].ay -= Fy*mfac;
        r->particles[ic].az -= Fz*mfac;
      }
   }

}

// for a planetary potential 
// U = -GM/r (1 + C20 + C22 + S22)
void quadforce(struct reb_simulation* const r, double C20, double C22, double S22, double Rplusp, double Mars_omega, int ic)
{
   if (ic >= r->N) return;
   if (fabs(C22) + fabs(S22) + fabs(C20)== 0.0) return;
   double Gmc = r->G * r->particles[ic].m *Rplusp*Rplusp; // GmR^2
   double phi_Mars = fmod(Mars_omega*r->t,2.0*M_PI);
   for(int i =0;i < r->N;i++){
     if (i != ic){
        double dx = r->particles[i].x - r->particles[ic].x;
        double dy = r->particles[i].y - r->particles[ic].y;
        double dz = r->particles[i].z - r->particles[ic].z;
        double x2 = dx*dx; double y2 = dy*dy; double z2 = dz*dz;
        double w = sqrt(x2 + y2 + RSOFT*RSOFT)  ;
        double rr = sqrt(x2 + y2 + z2 + RSOFT*RSOFT); // softening here
        double r4 = pow(rr,4.0); 
        double phi = atan2(dy,dx) - phi_Mars;
        // double theta = atan2(w,dz);
        double ctheta = dz/rr; 
        double stheta = w/rr; // = w/rr
        double cphi2 = cos(2*phi);
        double sphi2 = sin(2*phi);
        double dUdr   = 3*Gmc/r4*(C20*(3*ctheta*ctheta - 1.0)/2 
                                + (C22*cphi2 + S22*sphi2)*3*stheta*stheta);
        double dUdtor = -Gmc/r4*(-C20 + 2*C22*cphi2 + 2*S22*sphi2)*3*stheta*ctheta;
        double dUdpors= -Gmc/r4*(-C22*sphi2 + S22*cphi2)*6*stheta;
        double Fx = -dUdr*dx/rr;
        double Fy = -dUdr*dy/rr;
        double Fz = -dUdr*dz/rr;
        Fx -= dUdtor*dx*dz/(rr*w);
        Fy -= dUdtor*dy*dz/(rr*w);
        Fz += dUdtor*w/rr;
        Fx += dUdpors*dy/w;
        Fy -= dUdpors*dx/w;
        r->particles[i].ax += Fx; //  accelerations on nodes, mult by m to get total force
        r->particles[i].ay += Fy;
        r->particles[i].az += Fz;
        double mfac = r->particles[i].m/r->particles[ic].m;
        r->particles[ic].ax -= Fx*mfac; // opposite force
        r->particles[ic].ay -= Fy*mfac;
        r->particles[ic].az -= Fz*mfac;
     }
   }
}


/////////////////////////////////////////////
// compute orbital elements for two bodies
// if il1 != -1:
//    resolved body 1 contains particles [il1,ih1)
// else ih1 is index of point mass 1
// if il2 != -1:
//    resolved body 2 contains particles [il2,ih2)
// else ih2 is index of point mass 2
/////////////////////////////////////////////
void orbels_twobodies(struct reb_simulation* const r, int il1, int ih1, int il2, int ih2, double *aa, double *ee, double *ii, double *meananom, double *longnode, double *argperi)
{
   struct reb_particle* particles = r->particles;

   double x1 =0.0; double y1 =0.0; double z1 =0.0;
   double x2 =0.0; double y2 =0.0; double z2 =0.0;
   double vx1 =0.0; double vy1 =0.0; double vz1 =0.0;
   double vx2 =0.0; double vy2 =0.0; double vz2 =0.0;
   double m1; double m2;

   if (il1>=0){ // body 1 is resolved
     m1 = sum_mass(r,il1,ih1);
     compute_com(r,il1, ih1, &x1, &y1, &z1); // compute center of mass
     compute_cov(r,il1, ih1, &vx1,&vy1,&vz1);
   }
   else { // point mass 1
     m1 = particles[ih1].m; 
     x1 =particles[ih1].x;   y1 =particles[ih1].y;  z1 =particles[ih1].z; 
     vx1=particles[ih1].vx; vy1 =particles[ih1].vy; vz1=particles[ih1].vz; 
   }

   if (il2>=0){ // body 2 is resolved
     m2 = sum_mass(r,il2,ih2);
     compute_com(r,il2, ih2, &x2, &y2, &z2);
     compute_cov(r,il2, ih2, &vx2,&vy2,&vz2);
   }
   else { // point mass 2
     m2 = particles[ih2].m;
     x2 =particles[ih2].x;  y2 =particles[ih2].y;  z2 =particles[ih2].z; 
     vx2=particles[ih2].vx; vy2=particles[ih2].vy; vz2=particles[ih2].vz; 
   }

   double GMM = r->G *(m1+m2);
   
   PhaseState ps;
   OrbitalElements orbel;
   ps.x  = x2-x1; 
   ps.y  = y2-y1; 
   ps.z  = z2-z1; 
   ps.xd = vx2-vx1; 
   ps.yd = vy2-vy1; 
   ps.zd = vz2-vz1; 
   keplerian(GMM, ps, &orbel);
   *aa = orbel.a;
   *ee = orbel.e;
   *ii = orbel.i;
   *argperi = orbel.argperi;
   *longnode= orbel.longnode;
   *meananom= orbel.meananom;

}


/////////////////////////////////////////////////////
// put two resolved bodies in a binary orbit 
// return mean motion
// positions and velocity computed with orbelements
// orbels computed with GM taking into account both masses
// resolved bodies are [il1,ih1) and [il2,ih2)
/////////////////////////////////////////////////////
double add_twores_bin(struct reb_simulation* const r, 
    int il1, int ih1, int il2, int ih2,
    double aa,double ee, double ii,
    double longnode,double argperi,double meananom)  
{
   // struct reb_particle* particles = r->particles;
   double m1 = sum_mass(r,il1,ih1);
   double m2 = sum_mass(r,il2,ih2);
   double GMM = r->G*(m1+m2);

   OrbitalElements orbel;
   PhaseState ps;
   orbel.a =aa;
   orbel.e =ee;
   orbel.i =ii;
   orbel.argperi =argperi;
   orbel.longnode=longnode;
   orbel.meananom=meananom;
   cartesian(GMM, orbel, &ps);
 
// move both bodies to origin
   subtractcom(r,il1,ih1);  // subtract center of mass only from resolved body 1
   subtractcom(r,il2,ih2);  // ditto for body 2
   subtractcov(r,il1,ih1);  // center of velocity subtracted only from resolved body 1
   subtractcov(r,il2,ih2);  // ditto for body 2
   // move the second body
   move_resolved(r,ps.x,ps.y,ps.z,ps.xd,ps.yd,ps.zd,il2,ih2);
 
// move to center of mass frame, whole system
   subtractcom(r,0,r->N);  
   subtractcov(r,0,r->N); 
   
   double meanmotion = sqrt(GMM)/pow(aa,1.5);
   return meanmotion;
}


// the resolved body with indexes in [il,ih)
// exert an angular acceleration on it about center of mass
// alphax,alphay,alphaz, is the angular acceleration 
// does not change center of mass coordinates or velocity
// for each particle dv = (r cross alpha) *dt
void drift_spin(struct reb_simulation* const r, double tstep, int il, int ih, double alphax, double alphay, double alphaz)
{
   double xc,yc,zc;
   compute_com(r,il, ih, &xc, &yc, &zc);// compute center of mass
   // double vxc,vyc,vzc;
   // compute_cov(r,il, ih, &vxc, &vyc, &vzc);

   double alpha = sqrt(alphax*alphax + alphay*alphay + alphaz*alphaz);
   // size limit
   if (alpha > 1e-5){
     for(int i=il;i<ih;i++){
       double dx =  r->particles[i].x - xc;
       double dy =  r->particles[i].y - yc;
       double dz =  r->particles[i].z - zc;
       double rcrosso_x = -dy*alphaz + dz*alphay;   // r cross alpha
       double rcrosso_y = -dz*alphax + dx*alphaz;
       double rcrosso_z = -dx*alphay + dy*alphax;
                // change the spin with respect to center of mass
       r->particles[i].vx +=  rcrosso_x*tstep;
       r->particles[i].vy +=  rcrosso_y*tstep;
       r->particles[i].vz +=  rcrosso_z*tstep;
     }
   }
}


// open an ephemeris file and save it in an ephemeris structure
// mu is GM of the particle in the ephemeris
// jd0 is the julian date where t=0 in our simulation
// assuming xyz of ephemeris in km
// store in m and seconds
void read_eph(char* filename, struct ephem* new_ephem){
    new_ephem->Npt = 0;
    double scalefac = new_ephem->scalefac; // conversion factor
    FILE *fp;
    fp = fopen(filename,"r");
    // printf("read_eph file openned %s\n",filename);
    int iis = 0;
    char ss[300];
    while(iis ==0){  // find start of data
        fgets(ss,300,fp);
        if (ss[0] == '$'){
            iis = 1;
        }
    }
    struct ephem_pt pt;
    double jd,x,y,z;
    const double day_to_sec = 24*60*60.0;
    char* pch;
    while (iis==1){
        fgets(ss,300,fp);
        if (ss[0] == '$'){
            iis = 2;  // check that we are not at the end of data
        }
        else{
            pch = strtok (ss,","); // separate string into comment separated tokens
            sscanf(pch,"%lf",&jd);
            pch = strtok (NULL, ","); // should be the date as a string
            pch = strtok (NULL, ",");
            sscanf(pch,"%lf",&x);
            pch = strtok (NULL, ",");
            sscanf(pch,"%lf",&y);
            pch = strtok (NULL, ",");
            sscanf(pch,"%lf",&z);
            pt.jd = jd;
            pt.x = x*scalefac; // notice sign flip!
            pt.y = y*scalefac; // km to meters
            pt.z = z*scalefac;
            double dday = jd - new_ephem->jd0;
            pt.t = dday*day_to_sec; // time in seconds
            ephem_pt_add(new_ephem,pt);
            // printf("%.2e %.2e %.2e %.2e\n",pt.jd,pt.x,pt.y,pt.z);
        }
    }
    fclose(fp);
    // new_ephem->t0 = new_ephem->ephem_arr[0].t;
    // printf("read_eph file read in %s\n",filename);
}

int Nptmax = 0;
// add a time point to ephemeris! bare routine
void ephem_pt_add(struct ephem* ephem1, struct ephem_pt pt)
{
    int Npt = ephem1->Npt;
    while (Nptmax<= Npt){
            Nptmax += 128;
            ephem1->ephem_arr = realloc(ephem1->ephem_arr,sizeof(struct ephem_pt)*Nptmax);
    }
    ephem1->ephem_arr[Npt] = pt;
    ephem1->Npt++;
    
}

// figure out xyz position of ephemeris mass at time t
// via linear interpolation
// return 0 if outside ephemeris time window
int ephem_pos(struct ephem* ephem1, double t, double *xE, double *yE, double *zE){
    int Npt = ephem1->Npt;
    double t0 = ephem1->ephem_arr[0].t;
    double tmax = ephem1->ephem_arr[Npt-1].t;
    double dt_eph = ephem1->dt_eph;
    *xE=0.0; *yE=0.0; *zE=0.0;
    if (t < t0){
        return 0;
    }
    if (t >= tmax){
        return 0;
    }
    int index = (int)((t-t0)/dt_eph);
    double lam = fmod(t-t0,dt_eph)/dt_eph;
    // printf("ephempos index %d lam %lf\n",index,lam);
    // linear interpolation
    *xE = (1.-lam)*ephem1->ephem_arr[index].x + lam*ephem1->ephem_arr[index+1].x;
    *yE = (1.-lam)*ephem1->ephem_arr[index].y + lam*ephem1->ephem_arr[index+1].y;
    *zE = (1.-lam)*ephem1->ephem_arr[index].z + lam*ephem1->ephem_arr[index+1].z;
    return 1;
}

// add an ephemeris particle to particle list so that it can be displayed
// and so that gravitational force is automatically added
// position is at t=0
void add_particle_ephem(struct reb_simulation* r, struct ephem* ephem1)
{
    double x,y,z;
    double t=0.0;
    ephem_pos(ephem1, t, &x, &y, &z);  // get position
    printf("add particle ephem %lf %lf %lf\n",x,y,z);
    struct reb_particle pt;
    pt.x = x;
    pt.y = y;
    pt.z = z;
    pt.r = ephem1->R;
    pt.m = ephem1->M;
    ephem1->pi = r->N; // store particle index in ephemeris structure
    reb_simulation_add(r,pt);
}

// update ephemeris particle position
int update_ephem_pos(struct reb_simulation* r,struct ephem* ephem1)
{
    double x,y,z;
    int status=0;
    status = ephem_pos(ephem1, r->t, &x, &y, &z);  // get position
    if (status >0){
        r->particles[ephem1->pi].x =x; // update particle position
        r->particles[ephem1->pi].y =y;
        r->particles[ephem1->pi].z =z;
        // printf("emphempos: %lf %lf %lf\n",x,y,z);
    }
    return status;

}
