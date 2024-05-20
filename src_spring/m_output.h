/**
 * @file        m_output.h
 * @brief       springs 
 * @author      Alice Quilen 
 */


#ifndef _M_OUTPUT_H
#define _M_OUTPUT_H

#include"particle.h"
#include"spring.h"

// in m_output.c

void print_surf(struct reb_simulation* const r, int il, int ih, int *surfarr, char* filename);
void print_extended(struct reb_simulation* const r, int il, int ih, char* filename,int isfirst);
void print_extended_simp(struct reb_simulation* const r, int il, int ih, char* filename);
void print_extended_2nodes(struct reb_simulation* const r, int il, int ih, char* filename);
void print_pm(struct reb_simulation* const r, int ip, int ipert, char* filename);
void print_tab(struct reb_simulation* const r, int npert, char* filename);
void print_bin(struct reb_simulation* const r, int npert, char* filename);
void write_springs(struct reb_simulation* const r,char *fileroot, int index);
void read_springs(struct reb_simulation* const r,char *fileroot, int index);
void read_particles(struct reb_simulation* const r,char *fileroot, int index);
void write_particles(struct reb_simulation* const r,char *fileroot, int index);
void toistring(char *istring, int i);
void print_covar(struct reb_simulation* const r, int il, int ih, struct reb_particle* init_particles, 
   char* filename,int isfirst);
void print_stress(struct reb_simulation* const r, int il, int ih, char* fileroot, int index);


#endif // _M_OUTPUT_H

