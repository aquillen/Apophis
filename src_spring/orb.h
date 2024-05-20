/**
 * @file        orb.h
 * @brief       orb 
 * @author      Alice Quilen 
 */


#ifndef _ORB_H
#define _ORB_H

#include"particle.h"
#include"spring.h"


// a structure to hold ephemerus times
struct ephem_pt {
   double jd;  // Julian date of ephemerus data point
   double x;
   double y;
   double z;
   double t; 
};

// a structure to hold a list of ephemerus points for a point mass
struct ephem {
    int Npt;  // number of points in ephemerus
    struct ephem_pt* ephem_arr; // array of positions and times
    double M;  // mass kg
    double R;  // radius m
    double dt_eph; // time interval between points
    double jd0;  // time in ephem_arr is w.r.t to this julian date in days
    int pi;   // index in particle list
    double scalefac; // to convert to km, and includes a sign that can be adjusted
};

// in orb.c
double add_pluto_charon_kep(struct reb_simulation* const r,double mp,
    double mratio, double sep,
    double aa,double ee, double ii,double longnode,double argperi,double meananom,
    double r1);
double add_pt_mass_kep(struct reb_simulation* const r,
    int il, int ih, int ip, double m1, double r1,
    double aa,double ee, double ii,
    double longnode,double argperi,double meananom);
double add_one_mass_kep(struct reb_simulation* const r, double m1,
    double aa,double ee, double ii,
    double longnode,double argperi,double meananom,
    double mball, double r1);
double add_twores_bin(struct reb_simulation* const r,
    int il1, int ih1, int il2, int ih2,
    double aa,double ee, double ii,
    double longnode,double argperi,double meananom);
void dodrift_dv(double tstep,double inv_taua, double inv_taue,  double GMM, double x, double y, double z, double vx, double vy, double vz, double *dvx, double *dvy, double *dvz);
void dodrift_twobodies(struct reb_simulation* const r, double tstep, double inv_taua, double inv_taue, int il1, int ih1, int il2, int ih2);
void dodrift_bin(struct reb_simulation* const r, double tstep, double inv_taua, double inv_taue, int im1, int im2);
void dodrift_res(struct reb_simulation* const r, double tstep, double inv_taua, double inv_taue, int im1, int il, int ih);
void dodrift_twores(struct reb_simulation* const r, double tstep, double inv_taua, double inv_taue,  int il1, int ih1, int il2, int ih2);
void add_one_mass_cartesian(struct reb_simulation* const r, double m1, double r1, double x,double y, double z, double vx, double vy, double vz);
void quadJ2pole(struct reb_simulation* const r,double J2p, double Rplusp,double phip, double thetap, int ic);
void quadC22S22(struct reb_simulation* const r,double C22, double S22, double Rplusp, int ic);
void quadforce(struct reb_simulation* const r, double C20, double C22, double S22, double Rplusp, double Mars_omega, int ic);
void orbels_twobodies(struct reb_simulation* const r, int il1, int ih1, int il2, int ih2, double *aa, double *ee, double *ii, double *meananom, double *longnode, double *argperi);
void drift_spin(struct reb_simulation* const r, double tstep, int il, int ih, double alphax, double alphay, double alphaz);
void read_eph(char* filename,  struct ephem* ephem1);
void ephem_pt_add(struct ephem* ephem1, struct ephem_pt pt);
void add_particle_ephem(struct reb_simulation* r, struct ephem* ephem1);
int update_ephem_pos(struct reb_simulation* r,struct ephem* ephem1);


#endif // _ORB_H


