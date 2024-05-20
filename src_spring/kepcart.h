// used by kepcart

#ifndef _KEPCART_
#define _KEPCART_

// phasespace and orbital element structures

typedef struct {
  double x, y, z, xd, yd, zd;
} PhaseState;

// orbital elements: 
//   a semi-major axis
//   e eccentricity 
//   i inclination
//   longnode  longitude of the ascending node 
//   argperi argument of pericenter
//   meananom  mean anomaly
typedef struct {
  double a, e, i, longnode, argperi, meananom;
} OrbitalElements;


#define GG 6.6732e-8      // gravitational constant cgs 
#define CC 2.997924562e10 // speed of light, cm/s 
#define Msol 1.989e33     // mass of sun,  g 
#define AU 1.49597892e13  // astronomical unit, cm 

/* in kepcart.cpp */

// handing probability distributions randoms
double rayleigh(double sigma);
double powerlaw(double xmin, double xmax, double gamma);

// kepcart conversions
double ecc_ano(double e,double l);  // solves kepler's eqn
double ecc_anohyp(double e,double l);  // solves kepler's eqn
void keplerian(double GM, PhaseState state, OrbitalElements *orbel);
void cartesian(double GM, OrbitalElements orbel, PhaseState *state);
double kepler(double ecc, double mean_anom);

// for integrating with f,g functions
void kepstep(double t, double M1,
   PhaseState state, PhaseState *newstate);
double solvex(double r0dotv0, double alpha,
                double M1, double r0, double dt);
double C_prussing(double y);
double S_prussing(double y);

/*
void kepstep(double t, double M1,
   double x,  double y,  double z,
   double vx, double vy, double vz,
   double *xnew, double *ynew, double *znew,
   double *vxnew, double *vynew, double *vznew);
*/


#endif

