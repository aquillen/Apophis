/**
 * @file        heat.h
 * @brief       springs 
 * @author      Alice Quilen 
 */


#ifndef _HEAT_H
#define _HEAT_H


#include"particle.h"
#include"spring.h"

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


#endif // _HEAT_H


