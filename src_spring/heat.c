
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"  
#include "tools.h"  // from rebound 
#include "output.h" // from rebound?
#include "spring.h"
#include "kepcart.h"


extern int NS; // number of springs
extern int NPERT;  // number of external perturbing bodies
double *heatvec; // global so can be reached by all routines here
// heatvec is an array for each spring


// print out heat file
// npert is number of point mass perturbers
// ndt is number of times heatvec is added to 
//   (need to normalize so units of power)
// powerfac is how much of heat to distribute at center of spring 
//    rather than at nodes
// also outputs rotated xy, after rotates assuming xy plane for orbit and toward pertuber
void print_heat(struct reb_simulation* const r, int npert, char* filename, int ndt, double powerfac)
{
   struct reb_particle* particles = r->particles;

   norm_heatvec(ndt);  // normalize power so is average power

   FILE *fpo;
   fpo = fopen(filename, "w");
   fprintf(fpo,"#%.2e xyz dedt xrot yrot \n",r->t);
   int il =0;
   int ih =r->N-npert; // NPERT?
   double xc,yc,zc;
   compute_com(r,il, ih, &xc, &yc, &zc);
   int im1 = r->N -1; // index for primary perturber
   // double m1 = particles[im1].m;
   double x1 = particles[im1].x; double y1 = particles[im1].y; double z1 = particles[im1].z;
   double theta = atan2(y1-yc,x1-xc);
   double ct = cos(theta); double st = sin(-theta);  // note sign!
   fprintf(fpo,"#%.3f %.3f %.3f %.3f %.3f %.3f %.3f \n",x1,y1,z1,xc,yc,zc,theta);
// theta =0 when m1 is at +x direction compared to resolved body
   double xmid,ymid,zmid;

   for(int i=0;i<NS;i++){
      spr_xyz_mid(r, springs[i], xc, yc, zc, &xmid, &ymid, &zmid);
       // rotate around center of body in xy plane
   // after rotation:
   // +x is toward perturber
   // +y is direction of rotation of perturber w.r.t to body
   // so -y is headwind direction for body in body
   // and +y is tailwind surface on body?????
      double xrot = xmid*ct - ymid*st;  
      double yrot = xmid*st + ymid*ct;
      // double power = dEdt(r,springs[i]);
      
      double power = heatvec[i];
      // double powerfac = 0.5; // fraction of heat goes to center of spring
      // heat on center of spring
      fprintf(fpo,"%d %.3f %.3f %.3f %.5e %.3f %.3f\n"
          ,i,xmid,ymid,zmid,power*powerfac,xrot,yrot);
      // put some heat on nodes as well as center of spring
      if (powerfac < 1.0){ 
         int ii = springs[i].i; 
         int jj = springs[i].j;
         double xpi = particles[ii].x - xc;
         double ypi = particles[ii].y - yc;
         double zpi = particles[ii].z - zc;
         double xpj = particles[jj].x - xc;
         double ypj = particles[jj].y - yc;
         double zpj = particles[jj].z - zc;
         fprintf(fpo,"%d %.3f %.3f %.3f %.5e %.3f %.3f\n"
           ,i,xpi,ypi,zpi,power*(1.0-powerfac)*0.5,xrot,yrot);
         fprintf(fpo,"%d %.3f %.3f %.3f %.5e %.3f %.3f\n"
           ,i,xpj,ypj,zpj,power*(1.0-powerfac)*0.5,xrot,yrot);
      }
   }
   fclose(fpo);
   for (int i=0;i<NS;i++) heatvec[i] = 0; // set back to zero
}

// add to heat store vector 
// stores power for each spring
// uses a heatvec array to store this
// adds dE/dt to array -- this is power  (not dependent on dt)
// so power is added together
// actually later on we should normalize this
// by the number of timesteps we did it
void addto_heatvec(struct reb_simulation* const r)
{
   static int first = 0;
   if (first==0){
      first=1;
      heatvec = malloc(NS*sizeof(double));
      for(int i=0;i<NS;i++) heatvec[i]=0.0;
   }
   if (r->dt >0){
      for(int i=0;i<NS;i++){
         double power = dEdt(r,springs[i]);
         heatvec[i] += power;
      }
   }
}

// normalize heatvec
// ndt is the number of timesteps used to store power!
void norm_heatvec(int ndt)
{
   for(int i=0;i<NS;i++){
      heatvec[i] /= ndt;
   }
}


// make a heat file name depending on numbers of tp
void hfilename(struct reb_simulation* const r,char *root, double tp, char *fname){
   int xd = (int)(r->t/tp);
   char junks[20];
   sprintf(junks,"%d",xd);
   sprintf(fname,"%s_",root);
   if (xd < 100000) strcat(fname,"0");
   if (xd < 10000)  strcat(fname,"0");
   if (xd < 1000)   strcat(fname,"0");
   if (xd < 100)    strcat(fname,"0");
   if (xd < 10)     strcat(fname,"0");
   strcat(fname,junks);
   strcat(fname,"_heat.txt");
}

// make a node file name depending on numbers of tp
void nfilename(struct reb_simulation* const r,char *root, double tp, char *fname){
   int xd = (int)(r->t/tp);
   char junks[20];
   sprintf(junks,"%d",xd);
   sprintf(fname,"%s_",root);
   if (xd < 100000) strcat(fname,"0");
   if (xd < 10000)  strcat(fname,"0");
   if (xd < 1000)   strcat(fname,"0");
   if (xd < 100)    strcat(fname,"0");
   if (xd < 10)     strcat(fname,"0");
   strcat(fname,junks);
   strcat(fname,"_node.txt");
}

// allocate a node vector, index per each mass node (including point masses)
// initialize each node with a specific heat and temperature
struct node *mknodevec(struct reb_simulation* const r, double cv,double T0){
   struct node *nodevec;
   nodevec = malloc(sizeof(struct node)*(r->N+10));
   for (int i=0;i< r->N+10;i++){
        nodevec[i].surf = 0;  // initialize as if in interior
        nodevec[i].temp = T0; // temperature
        nodevec[i].cv   = cv; // initialize specific heat      
   }
   return nodevec;
}

// make sure nodevec.surf is 1 or zero if node is on surface
// intialize surface temperature on each surface node
// array surfp tells whether node is on surface
// array of nodes nodevec stores node information
void surface_nodes(int *surfp, struct node *nodevec, int ntot,  double Tsurf){
   for (int i=0;i<ntot;i++){
      nodevec[i].surf = surfp[i];  // transfer marked surface to nodevec
      if (surfp[i] == 1){
        nodevec[i].temp = Tsurf; // initialize surface temperature
      }
   }
}

// print out node temperature file
// npert is number of point mass perturbers
void print_node(struct reb_simulation* const r, int npert, struct node *nodevec, char* filename)
{
   struct reb_particle* particles = r->particles;

   FILE *fpo;
   fpo = fopen(filename, "w");
   fprintf(fpo,"#%.2e i xyz vxvyvz T cv surf m xrot yrot \n",r->t);
   
   double xc,yc,zc;
   int il = 0;
   int ih =  r->N-1;
   compute_com(r,il, ih, &xc, &yc, &zc);
   int im1 = r->N -1; // index for primary perturber
   double x1 = particles[im1].x; double y1 = particles[im1].y; double z1 = particles[im1].z;
   double theta = atan2(y1-yc,x1-xc);
   double ct = cos(theta); double st = sin(-theta);  // note sign!
   fprintf(fpo,"#%.3f %.3f %.3f %.3f %.3f %.3f %.3f \n",x1,y1,z1,xc,yc,zc,theta);
   
   for (int i=il;i<ih;i++) {
      double x = particles[i].x;
      double y = particles[i].y;
      double z = particles[i].z;
      double vx = particles[i].vx;
      double vy = particles[i].vy;
      double vz = particles[i].vz;
      double m = particles[i].m;
      double xrot = x*ct - y*st;  
      double yrot = x*st + y*ct;
      fprintf(fpo,"%4d  ",i);
      fprintf(fpo,"%10.6f %10.6f %10.6f ",x,y,z);
      fprintf(fpo,"%10.6f %10.6f %10.6f   ",vx,vy,vz);
      fprintf(fpo,"%10.6e %10.6f %3d %10.6f "
           ,nodevec[i].temp,nodevec[i].cv,nodevec[i].surf,m);
      fprintf(fpo,"%10.6e %10.6f ",xrot,yrot);
      fprintf(fpo,"\n");
   }
   fclose(fpo); 
      
}

// update temperatures on each node for the timestep with size dt
// conductive heat transport over spring network
double *newTvec;
void transport_heat(struct reb_simulation* const r, int npert, struct node *nodevec, double dt)
{
   static int first=0;
   if (first==0) { 
      first =1;
      newTvec = malloc(sizeof(double)*r->N);
   }
   for(int k=0;k<r->N;k++) newTvec[k] = 0.0;
   
   for(int i=0;i<NS;i++){ //loop over springs
      int ii = springs[i].i;
      int jj = springs[i].j;
      double qij = (nodevec[jj].temp - nodevec[ii].temp)*springs[i].k_heat;
      newTvec[ii] +=  qij*dt/(r->particles[ii].m*nodevec[ii].cv);
      newTvec[jj] -=  qij*dt/(r->particles[jj].m*nodevec[jj].cv);
   }
  
   // only update temperatures after all fluxes have been added up
   for(int k=0;k<r->N-npert;k++){
      if (nodevec[k].surf ==0) { // only do this to interior nodes 
          nodevec[k].temp += newTvec[k];
      }
   }
}

// apply tidal heat to internal nodes, raising their temperature
void heat_nodes_tidal(struct reb_simulation* const r, int npert, struct node *nodevec, double dt)
{
   for(int k=0;k<NS;k++){
      double power = dEdt(r,springs[k]); // heating from springs
      int ii = springs[k].i; 
      int jj = springs[k].j;
      double dT_ii = 0.5*power*dt/(r->particles[ii].m*nodevec[ii].cv);
      double dT_jj = 0.5*power*dt/(r->particles[jj].m*nodevec[jj].cv);
      if (nodevec[ii].surf ==0){  // only change temp of interior nodes 
          nodevec[ii].temp += dT_ii;
      }
      if (nodevec[jj].surf ==0){  // only change temp of interior nodes 
          nodevec[jj].temp += dT_jj;
      }
   }

}


// compute thermal conductivity 
// using midpoints of springs in radial range
//    from center of mass [rmin,rmax]
// uses rest lengths
// only computes center of mass using particles index range [il,ih)
double Kthermal_mush(struct reb_simulation* const r, int il, int ih, double rmin, double rmax){
  double sum=0.0;
  double xc,yc,zc;
  compute_com(r,il, ih, &xc, &yc, &zc); // center of mass coords for particles in range
  double rmid,thetamid,phimid;
  for (int i=0;i<NS;i++){
       spr_ang_mid(r, springs[i],xc,yc,zc, &rmid, &thetamid, &phimid);

       double rc = rmid; // center of spring
       if ((rc<rmax) && (rc > rmin)){
         double kappa = springs[i].k_heat;
         double Li = springs[i].rs0;
         sum += kappa*Li*Li;
       }
  }
  double volume = (4.0*M_PI/3.0)*(pow(rmax,3.0) - pow(rmin,3.0)); // in shell
  double K_T = sum/(4.0*volume);  // my guess for the thermal conductivity 
  return K_T; // return conductivity
}


// addjust spring properties using temperatures of nodes, linear model
// heat capacity not yet adjusted!
void adjust_spring_temp_lin(struct reb_simulation* r, struct node *nodevec, struct spring mush, double dksdT, double dgammadT, double dkheatdT)
{
   for (int k=0;k<NS;k++){
      int ii = springs[k].i; 
      int jj = springs[k].j;
      double Ti = nodevec[ii].temp;
      double Tj = nodevec[jj].temp;
      double Tave = 0.5*(Ti + Tj);
      springs[k].k_heat = mush.k_heat*(1. + Tave*dkheatdT);
      springs[k].ks     = mush.ks    *(1. + Tave*dksdT);
      springs[k].gamma  = mush.gamma *(1. + Tave*dgammadT);
   }
}

// addjust spring properties using temperatures of nodes, two state model
// heat capacity not yet adjusted!
void adjust_spring_temp_ts(struct reb_simulation* r, struct node *nodevec, struct spring mush_hot, struct spring mush_cold, double Ttrans)
{
   for (int k=0;k<NS;k++){
      int ii = springs[k].i; 
      int jj = springs[k].j;
      double Ti = nodevec[ii].temp;
      double Tj = nodevec[jj].temp;
      double Tave = 0.5*(Ti + Tj);
      if (Tave >= Ttrans){
         springs[k].k_heat = mush_hot.k_heat;
         springs[k].ks     = mush_hot.ks    ;
         springs[k].gamma  = mush_hot.gamma ;
      }
      else {
         springs[k].k_heat = mush_cold.k_heat;
         springs[k].ks     = mush_cold.ks    ;
         springs[k].gamma  = mush_cold.gamma ;
      }
   }
}


// add in a constant heat rate on each node, radiogenic
// dote_radg is the energy per unit mass per unit time
void heat_nodes_radiogenic(struct reb_simulation* r, struct node *nodevec, double dote_radg, int npert)
{
   for (int i=0;i<r->N -npert;i++){
      double dT_i = dote_radg*r->dt/nodevec[i].cv;
      if (nodevec[i].surf ==0){  // only change temp of interior nodes 
          nodevec[i].temp += dT_i;
      }
   }
}


void adjust_nodes_cp(struct reb_simulation* const r, int npert, struct node *nodevec, double cpnew, double Ttrans, int above) 
{
   int il =0;
   int ih =r->N - npert;
   for(int i=il;i<ih;i++){
      if ((nodevec[i].temp > Ttrans) && (above ==1)){
          nodevec[i].cv   = cpnew; // change specific heat if hot
      }
      if ((nodevec[i].temp < Ttrans) && (above ==0)){  // if cold
          nodevec[i].cv   = cpnew; 
      }
   }
      
}

