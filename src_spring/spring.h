/**
 * @file        spring.h
 * @brief       springs 
 * @author      Alice Quilen 
 */


#ifndef _SPRING_H
#define _SPRING_H


#include"particle.h"

// spring structure 
struct spring  {
	double ks;   // spring constant
	double rs0;  // distance of no force
	double gamma; // damping coefficient
	double k_heat;  //  heat diffusion coefficient!
	int    i;    // vertex 1 referring to a particle
	int    j;    // vertex 2 referring to a particle
};

// coordination structure for each node
struct coord {
       int ncoord;  // number of springs connected to this node
       int coordlist[50];  // a list of spring indices ** note number here is hard wired, giving a max!
       double lx;
       double ly;
       double lz;
};

struct node {
        int surf;  // is 1 if is a surface node
        double temp;  // temperature
        double cv;    // specific heat, probably integrated over mass of node?
};


// globals
extern struct spring* springs;
extern int NS; // numbers of springs
extern int NPERT; // number of point masses
// extern double b_distance; // mush formation
// extern double mush_distance; // mush spring connection 
// extern double t_reform; // formation springs timescale
// extern double gamma_all; // for gamma  of all springs
// extern struct reb_particle* init_particles;  // initial positions of particles 

// list of springs related subroutines

// in spring.c

double mean_radius(struct reb_simulation* const r,int il, int ih);
void set_gamma_fac(double gamma_fac, int il, int ih);
void set_gamma(double new_gamma, int il, int ih);
int *marksurface_cone(struct reb_simulation* const r, double mind,
      double rb,double slope);
int *marksurface(struct reb_simulation* const r, int N_bennu, double mind);
int *marksurface_football(struct reb_simulation* const r, double mind,
      double ax,double by,double cz);
double detI(double Ixx,double Iyy,double Izz,double Ixy,double Iyz,double Ixz);
double dEdt_total(struct reb_simulation* const r);
double dEdt(struct reb_simulation* const r,struct spring spr);
double max_radius(struct reb_simulation* r, int il, int ih);
double min_radius(struct reb_simulation* r,int il, int ih);
double sum_mass(struct reb_simulation* const r,int il, int ih);
double mindist(struct reb_simulation* const r,int imin, int imax);
double compute_rot_kin(struct reb_simulation* const r, int il, int ih);
int add_spring_i(struct reb_simulation* const r,int i1, int i2,  struct spring spring_vals);
double spring_potential_energy(struct reb_simulation* const r);
double grav_potential_energy(struct reb_simulation* const r,int il, int ih);
double strain(struct reb_simulation* const r,struct spring spr);
double spring_length(struct reb_simulation* const r, struct spring spr);
int within_shape(struct reb_simulation* r,int il, int ih,
         double minr, double maxr, double x,double y, double z);
int nearest_to_shape(struct reb_simulation* r,int il, int ih, double x,double y, double z);
void del_spring(struct reb_simulation* const r, int i);
void spr_ang_mid(struct reb_simulation* const r, struct spring spr, double xc, double yc, double zc, double *rmid, double *thetamid, double *phimid);
void spr_xyz_mid(struct reb_simulation* const r, struct spring spr, double xc, double yc, double zc, double *xmid, double *ymid, double *zmid);
void normalize(struct reb_simulation* const r, struct reb_particle *pt);
void springs_add(struct reb_simulation* const r,struct spring spr);
void compute_com(struct reb_simulation* const r,int il, int ih, double *xc, double *yc, double *zc);
void compute_cov(struct reb_simulation* const r,int il, int ih, double *vxc, double *vyc, double *vzc);
void centerbody(struct reb_simulation* const r,int il, int ih);
void subtractcov(struct reb_simulation* const r,int il, int ih);
void subtractcom(struct reb_simulation* const r,int il, int ih);
void spring_forces(struct reb_simulation* const r);
void zero_accel(struct reb_simulation* const r);
void connect_springs_dist(struct reb_simulation* const r, double h_dist, int i0, int imax, struct spring spring_vals);
void connect_springs_dist_p(struct reb_simulation* const r, double h_dist, int i0, int imax, int j0, struct spring spring_vals);
void rand_football_from_sphere(struct reb_simulation* r, double dist, double ax, double by, double cz, double total_mass);
void rand_football(struct reb_simulation* const r, double dist, double ax, double by, double cz, double total_mass);
void rand_rectangle(struct reb_simulation* const r, double dist, double ax, double by, double cz, double total_mass);
void rand_rectangle_2d(struct reb_simulation* const r, double dist, double ax, double by, double total_mass);
void rand_cone(struct reb_simulation* const r, double dist, double rb, double slope,  double total_mass);
// void set_gamma(struct reb_simulation* const r,double new_gamma, int il, int ih)
// void set_gamma_fac(struct reb_simulation* const r,double gamma_fac, int il, int ih)
void spin(struct reb_simulation* const r,int il, int ih, double omegax, double omegay, double omegaz);
void make_binary_spring(struct reb_simulation* const r,double m1, double m2, double sep, double omegax, double omegay, double omegaz, struct spring spring_vals);
void measure_L(struct reb_simulation* const r,int il, int ih, double *llx, double *lly, double *llz);
void measure_L_origin(struct reb_simulation* const r,int il, int ih, double *llx, double *lly, double *llz);
void mom_inertia(struct reb_simulation* const r,int il, int ih, double *Ixx, double *Iyy, double *Izz, double *Ixy, double *Iyz, double *Ixz);
void compute_semi(struct reb_simulation* const r, int il,int ih, int im1, double *aa, double *meanmo, double *ee, double *ii, double *LL);
void compute_Lorb(struct reb_simulation* const r, int il,int ih, int npert, double *llx, double *lly, double *llz);
void compute_semi_bin(struct reb_simulation* const r, int il, int ih, int npert, double *aa, double *meanmo, double *ee, double *ii, double *LL);
void total_mom(struct reb_simulation* const r,int il, int ih, double *ppx, double *ppy, double *ppz);
void rotate_origin(struct reb_simulation* const r, int il, int ih, double alpha, double beta, double ggamma);
void rotate_body(struct reb_simulation* const r, int il, int ih, double alpha, double beta, double ggamma);
void rotate_vector(double x, double y, double z, double *xr, double *yr, double *zr, double alpha, double beta, double ggamma);
void invI(double Ixx,double Iyy,double Izz,double Ixy,double Iyz,double Ixz, double *Axx,double *Ayy,double *Azz,double *Axy,double *Ayz,double *Axz);
void eigenvalues(double Ixx, double Iyy,  double Izz, double Ixy, double Iyz, double Ixz, double *eig1, double *eig2, double *eig3);
void body_spin(struct reb_simulation* const r, int il, int ih, double *omx, double *omy, double *omz, double *big, double *middle, double *small);
void adjust_ks(struct reb_simulation* const r, int npert, double ksnew, double gammanew, double kheatnew, double rmin, double rmax);
void adjust_ks_abc(struct reb_simulation* const r, int npert, double ksnew, double gammanew, double kheatnew,  double a, double b, double c, double x0,double y0, double z0, int inside);
void adjust_ks_abc_fac(struct reb_simulation* const r, int npert, double ks_fac, double gamma_fac, double kheat_fac, int mm, double phi0, double a, double b, double c, double x0,double y0, double z0, int inside);
void adjust_mass_abc(struct reb_simulation* const r, int npert, double mfac, double a, double b, double c, double x0,double y0, double z0, int inside);
void adjust_mass_side(struct reb_simulation* const r, int npert, double mfac, double xmin);
void move_resolved(struct reb_simulation* r, double dx,double dy,double dz, double dvx,double dvy,double dvz, int il,int ih);
void read_vertex_file(struct reb_simulation* r, char *filename);
void rmshape_vertices(struct reb_simulation* const r, int N_bennu);
void rand_bennu(struct reb_simulation* const r, double dist, double total_mass);
void rescale_xyz(struct reb_simulation* const r, int il, int ih, double scalefac);
void connect_springs_dist_nodemax(struct reb_simulation* const r, double h_dist_min,  double h_dist_max, int i0, int imax, struct spring spring_vals, int nodemax);
void spring_force_k_one(struct reb_simulation* const r,int i_s, double *Fx, double *Fy, double *Fz, double *Lx, double *Ly, double *Lz);
void stretch(struct reb_simulation* const r,int il,int ih, double scalex,double scaley,double scalez);
void eigenvecs(double Cxx,double Cyy,double Czz,double Cxy, double Cyz, double Cxz, double eigv,double *eigvecx,double *eigvecy, double *eigvecz);
void rotate_to_principal(struct reb_simulation* const r, int il, int ih, double *ee1,double *ee2, double *ee3);
void vec_mul(double Axx, double Axy, double Axz, double Ayx, double Ayy, double Ayz, double Azx, double Azy, double Azz, double bx,  double by,  double bz, double *cx, double *cy, double *cz);
double mean_L(struct reb_simulation* const r);
double fill_hcp(struct reb_simulation* r, double dd, double ax, double by, double cz, double total_mass);
double fill_cubic(struct reb_simulation* r, double dd, double ax, double by, double cz, double total_mass);
double Young_mush(struct reb_simulation* const r, int il, int ih, double rmin, double rmax);
double Young_mush_big(struct reb_simulation* const r, int il, int ih);
void store_xyz0(struct reb_simulation* const r, struct reb_particle* init_particles); 

// in heat.c 

struct node *mknodevec(struct reb_simulation* const r, double cv,double T0);
void print_heat(struct reb_simulation* const r, int npert, char* filename, int ndt, double powerfac);
void addto_heatvec(struct reb_simulation* const r);
void norm_heatvec(int ndt);
void hfilename(struct reb_simulation* const r,char *root, double tp, char *fname);
void nfilename(struct reb_simulation* const r,char *root, double tp, char *fname);
void surface_nodes(int *surfp, struct node *nodevec, int ntot,  double Tsurf);
void print_node(struct reb_simulation* const r, int npert, struct node *nodevec, char* filename);
void transport_heat(struct reb_simulation* const r, int npert, struct node *nodevec, double dt);
void heat_nodes_tidal(struct reb_simulation* const r, int npert, struct node *nodevec, double dt);
void adjust_spring_temp_lin(struct reb_simulation* r, struct node *nodevec, struct spring mush, double dksdT, double dgammadT, double dkheatdT);
void adjust_spring_temp_ts(struct reb_simulation* r, struct node *nodevec, struct spring mush_hot, struct spring mush_cold, double Ttrans);
void heat_nodes_radiogenic(struct reb_simulation* r, struct node *nodevec, double dote_radg, int npert);
void adjust_nodes_cp(struct reb_simulation* const r, int npert, struct node *nodevec, double cpnew, double Ttrans, int above);
double Kthermal_mush(struct reb_simulation* const r, int il, int ih, double rmin, double rmax);


#endif // _SPRING_H


