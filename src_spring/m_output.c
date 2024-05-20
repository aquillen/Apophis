
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
#include "m_output.h"
#include "stress.h"


// Contains routines:
// print_surf prints out surface particles from a list 
// print_extended  prints out positions of particles in an extended body
//      also prints out some rotation information 
// print_extended_simp similar to print_extended but less stuff is printed out
// print_pm print out information about a point mass
// write_springs, read_springs  write and read spring information (all of them) 
// write_particles, read_particles write and read particles (all of them)
// print_covar print out covariance matrix for an extended body, needed for Kabsh algorithm
// a couple more are used here!

extern int NS;  // numbers of springs

// printout out positions and velocities of only surface particles
void print_surf(struct reb_simulation* const r, int il, int ih, int *surfarr, char* filename)
{
   FILE *fpo;
   static int tstart =0;
   static double toff = 0.0;
   if (tstart==0){
     tstart=1;
     toff = r->t;
     printf("print_surf toff=%.6f\n",toff);
   }
   fpo = fopen(filename, "w");
   fprintf(fpo,"# t=%.8f #txyzvxvyvzaxayazrm\n",r->t); // actual time
   for(int i=0;i< r->N;i++){
      if (surfarr[i] ==1)
        fprintf(fpo,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6e\n",
          r->t - toff, 
          r->particles[i].x, r->particles[i].y, r->particles[i].z,
          r->particles[i].vx, r->particles[i].vy, r->particles[i].vz,
          r->particles[i].ax, r->particles[i].ay, r->particles[i].az,
          r->particles[i].r, r->particles[i].m); 
   }
   fclose(fpo);
   printf("\n print_surf: N=%d %s\n",r->N,filename);
}


// print out information for an extended body indices [il,ih)
// into filename
// first time this is called please set isfirst to 1 (prints first line)
// afterwards set isfirst to 0 (so does not print first line)
// this is not done automatically so that this routine can be called for more than 1 extended body
void print_extended(struct reb_simulation* const r, int il, int ih, char* filename,int isfirst)
{
   FILE *fpo;
   if (isfirst==1){
     fpo = fopen(filename, "w"); // first time create file
     fprintf(fpo,"# t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz KErot PEspr PEgrav Etot dEdtnow\n");
   }
   else {
     fpo = fopen(filename, "a");
   }
   fprintf(fpo,"%.3f ",r->t);
   double xc =0.0; double yc =0.0; double zc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc); // compute center of mass 
   double vxc =0.0; double vyc =0.0; double vzc =0.0;
   compute_cov(r,il, ih, &vxc, &vyc, &vzc); // compute center of velocity 
   fprintf(fpo,"%.6f %.6f %.6f ",xc,yc,zc);
   fprintf(fpo,"%.6f %.6f %.6f ",vxc,vyc,vzc);

   double omx,omy,omz, Ibig, Imid, Ismall, llx,lly,llz;
   // computes spin vector omega using inverse of moment of inertia matrix and angular momentum
   body_spin(r, il, ih, &omx, &omy, &omz, &Ibig, &Imid, &Ismall);
   measure_L(r,il, ih, &llx, &lly, &llz); // computes spin angular momentum of spinning body
   fprintf(fpo,"%.6e %.6e %.6e ",omx, omy, omz);
   fprintf(fpo,"%.6e %.6e %.6e ",llx,lly,llz);
   // fprintf(fpo,"%.6e %.6e %.6e ",Ibig, Imid, Ismall);
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   // compute moment of inertia matrix
   mom_inertia(r,il,ih, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   fprintf(fpo,"%.6e %.6e %.6e %.6e %.6e %.6e ",Ixx,Iyy,Izz,Ixy,Iyz,Ixz);

   double KErot = compute_rot_kin(r,il,ih); // total kinetic energy (including rotational)
   double pe_springs = spring_potential_energy(r); // potential energy in springs
   double pe_grav = grav_potential_energy(r,il,ih); // potential energy in gravity 
   double E_tot = KErot+pe_springs+pe_grav; 
   double dEdtnow = dEdt_total(r); // current dissipation from velocities of particles
                                  // that are connected by springs
   fprintf(fpo,"%.5e %.5e %.10e %.10e %.10e ",KErot,pe_springs,pe_grav,E_tot,dEdtnow);
   fprintf(fpo,"\n");
   fclose(fpo);
}

// similar to last routine except without dEdt stuff
// print out information for an extended body indices [il,ih)
void print_extended_simp(struct reb_simulation* const r, int il, int ih, char* filename)
{
   static int first=0;
   FILE *fpo;
   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"# t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz \n");
   }
   else {
     fpo = fopen(filename, "a");
   }
   fprintf(fpo,"%.3f ",r->t);
   double xc =0.0; double yc =0.0; double zc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc);
   double vxc =0.0; double vyc =0.0; double vzc =0.0;
   compute_cov(r,il, ih, &vxc, &vyc, &vzc);
   fprintf(fpo,"%.5f %.5f %.5f ",xc,yc,zc);
   fprintf(fpo,"%.5f %.5f %.5f ",vxc,vyc,vzc);

   double omx,omy,omz, Ibig, Imid, Ismall, llx,lly,llz;
   // computes omega using inverse of moment of inertia matrix and angular momentum
   body_spin(r, il, ih, &omx, &omy, &omz, &Ibig, &Imid, &Ismall);
   measure_L(r,il, ih, &llx, &lly, &llz); // computes spin angular momentum of spining body
   fprintf(fpo,"%.4e %.4e %.4e ",omx, omy, omz);
   fprintf(fpo,"%.4e %.4e %.4e ",llx,lly,llz);
   // fprintf(fpo,"%.4e %.4e %.4e ",Ibig, Imid, Ismall);
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   mom_inertia(r,il,ih, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",Ixx,Iyy,Izz,Ixy,Iyz,Ixz);

   fprintf(fpo,"\n");
   fclose(fpo);
}

// print out information for two specific nodes in an extended body indices [il,ih)
// the first node has largest x value at initial time
// the second node has largest x value at initial time
// same nodes are printed each time after they are chosen
void print_extended_2nodes(struct reb_simulation* const r, int il, int ih, char* filename)
    // struct node* nodevec, 
{
   static int first=0;
   static int iz=0; // index of particle with largest z value
   static int ix=0; // index of particle with largest x value
   FILE *fpo;
   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"#t xyzvxvyvz(xnode) xyzvxvyvz(znode) \n");
     iz = nearest_to_shape(r,il,ih, 0.0,0.0,10.0);
     ix = nearest_to_shape(r,il,ih,10.0,0.0, 0.0);
   }
   else {
     fpo = fopen(filename, "a");
   }
   fprintf(fpo,"%.3f ",r->t);
   fprintf(fpo,"%.5f %.5f %.5f ",
        r->particles[ix].x,r->particles[ix].y,r->particles[ix].z);
   fprintf(fpo,"%.5f %.5f %.5f ",
        r->particles[ix].vx,r->particles[ix].vy,r->particles[ix].vz);
   fprintf(fpo,"%.5f %.5f %.5f ",
        r->particles[iz].x,r->particles[iz].y,r->particles[iz].z);
   fprintf(fpo,"%.5f %.5f %.5f ",
        r->particles[iz].vx,r->particles[iz].vy,r->particles[iz].vz);
   fprintf(fpo,"\n");
   fclose(fpo);
}

#define NPMAXE 100
// print out information for a point mass with index ip
// which point mass it is ipert (0 being the central one)
//  ip is its index in the particle list
// ipert is only used to make each point mass output be separate from the others
// you can use this routine to output non-pert pointmasses as long as they
// have separate ipert numbers called here
void print_pm(struct reb_simulation* const r, int ip, int ipert, char* filename)
{
   struct reb_particle* particles = r->particles;
   static int firstarr[NPMAXE]; // this is just so we can keep a record
   // of which bodies we have been outputing
   static int first=0;
   if (first==0){
      first=1;
      for (int i=0;i<NPMAXE;i++) firstarr[i]=0;
   }
   FILE *fpo;
   if (firstarr[ipert]==0){  
     firstarr[ipert]=1;
     fpo = fopen(filename, "w");                     
     fprintf(fpo,"# t x y z vx vy vz m\n");
     printf("%s w openned\n",filename);
   }
   else {
     fpo = fopen(filename, "a");
   }
   double x = particles[ip].x;
   double y = particles[ip].y;
   double z = particles[ip].z;
   double vx = particles[ip].vx;
   double vy = particles[ip].vy;
   double vz = particles[ip].vz;
   fprintf(fpo,"%.3f ",r->t);
   fprintf(fpo,"%.5f %.5f %.5f ",x,y,z);
   fprintf(fpo,"%.5f %.5f %.5f ",vx,vy,vz);
   fprintf(fpo,"%.3e ",particles[ip].m);
   fprintf(fpo,"\n");
   fclose(fpo);
}

// prints various info for an extended body in orbit about a pt mass or about
// a binary 
void print_tab(struct reb_simulation* const r, int npert, char* filename)
{
   // struct reb_particle* particles = r->particles;
   static int first=0;
   FILE *fpo;

   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"# t a n e i omx omy omz A B C E lx ly lz ang lox loy loz Ixx Iyy Izz Ixy Iyz Ixz\n");
   }
   else {
     fpo = fopen(filename, "a");
   }
   int il = 0;
   int ih = r->N - npert;
   double a,meanmo,ecc,incl,Lorb;
   if (npert ==1) { // orbit about a pt mass
      compute_semi(r, il, ih, r->N-1, &a, &meanmo, &ecc, &incl, &Lorb);
   }
   else   // orbit about a binary
      compute_semi_bin(r, il, ih, npert, &a, &meanmo, &ecc, &incl, &Lorb);
   fprintf(fpo,"%.3f %.6e %.4e %.4f %.4e ",r->t, a, meanmo, ecc, incl);
   double omx,omy,omz, Ibig, Imid, Ismall;
   // computes spin using inverse of moment of inertia matrix and angular momentum vec
   body_spin(r, il, ih, &omx, &omy, &omz, &Ibig, &Imid, &Ismall);
   fprintf(fpo,"%.4e %.4e %.4e %.3f %.3f %.3f ",omx, omy, omz, Ibig, Imid, Ismall);
   double E = Young_mush(r, il, ih, 0.0, 0.5);
   fprintf(fpo,"%.3f ",E);
   double llx,lly,llz;
   measure_L(r,il, ih, &llx, &lly, &llz); // spin angular momentum of spining body
   fprintf(fpo,"%.5f %.5f %.5f ",llx,lly,llz);
   double llx_o,lly_o,llz_o;
   compute_Lorb(r, il,ih, npert, &llx_o, &lly_o, &llz_o);  // total orbital angular momentum
   double L_spin = sqrt(llx*llx + lly*lly + llz*llz);
   double L_orb = sqrt(llx_o*llx_o + lly_o*lly_o + llz_o*llz_o);
   double LsdotLo = llx*llx_o + lly*lly_o + llz*llz_o; 
   double angle = acos(LsdotLo/(L_spin*L_orb));  // angular difference  between 
            // spin angular momentum and orbital angular momentum
   fprintf(fpo,"%.5f ",angle);
   fprintf(fpo,"%.5f %.5f %.5f ",llx_o,lly_o,llz_o); // orbital angular momentum
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   mom_inertia(r,il,ih, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",Ixx,Iyy,Izz,Ixy,Iyz,Ixz);


   fprintf(fpo,"\n");
   fclose(fpo);
}

void print_bin(struct reb_simulation* const r, int npert, char* filename)
{
   struct reb_particle* particles = r->particles;
   static int mark =-1;
   static int first=0;
   FILE *fpo;
   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"# t xyz,vxyz(res) xyz,vxyz(PT) xyz,vxyz(PC)  ");
     fprintf(fpo," lx ly lz Ixx Iyy Izz Ixy Iyz Ixz\n  ");
     double r2max = 0.0;
     for(int i=0;i<r->N-npert;i++){
       double rc2 = particles[i].x*particles[i].x + particles[i].y*particles[i].y
                  + particles[i].z*particles[i].z;
       if (rc2 > r2max){   // mark particle with largest initial r2 value 
         mark = i;
         r2max = rc2;
       }
     }
   }
   else {
     fpo = fopen(filename, "a");
   }

   fprintf(fpo,"%.3f ",r->t);

   int il = 0; // index range for resolved body
   int ih = r->N - npert;
   double xc =0.0; double yc =0.0; double zc =0.0;
   double vxc =0.0; double vyc =0.0; double vzc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc); // center of mass of resolved body
   compute_cov(r,il, ih, &vxc, &vyc, &vzc); // center of velocity of resolved body
   int iml = r->N -npert; // index range for perturbing masses, binary
   int imh = r->N; 
   double xb=0.0; double yb=0.0; double zb=0.0; 
   double vxb=0.0; double vyb=0.0; double vzb = 0.0;
   compute_com(r,iml, imh, &xb, &yb, &zb); // center of mass of binary 
   compute_cov(r,iml, imh, &vxb, &vyb, &vzb); // center of velocity of binary  

   double xr = xc-xb; double vxr = vxc-vxb;  // of resolved w.r.t binary center
   double yr = yc-yb; double vyr = vyc-vyb;
   double zr = zc-zb; double vzr = vzc-vzb;
   // store resolved center vs bin center
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",xr,yr,zr,vxr,vyr,vzr);

   double xs = particles[mark].x - xc;   // marked w.r.t to resolved center
   double ys = particles[mark].y - yc; 
   double zs = particles[mark].z - zc; 
   double vxs = particles[mark].vx - vxc; 
   double vys = particles[mark].vy - vyc; 
   double vzs = particles[mark].vz - vzc; 
   // store marked w.r.t. to resolved center 
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",xs,ys,zs,vxs,vys,vzs);

   if (npert==2){
     double x_CP  = particles[iml].x  - particles[iml+1].x;  // charon w.r.t pluto
     double y_CP  = particles[iml].y  - particles[iml+1].y; 
     double z_CP  = particles[iml].z  - particles[iml+1].z; 
     double vx_CP = particles[iml].vx - particles[iml+1].vx; 
     double vy_CP = particles[iml].vy - particles[iml+1].vy; 
     double vz_CP = particles[iml].vz - particles[iml+1].vz; 
     // store Charon w.r.t Pluto
     fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",x_CP,y_CP,z_CP,vx_CP,vy_CP,vz_CP);
   }
   else {
     fprintf(fpo,"0 0 0 0 0 0 ");
   }

   double llx,lly,llz;
   measure_L(r,il, ih, &llx, &lly, &llz); // spin angular momentum of spining body
   fprintf(fpo,"%.5f %.5f %.5f ",llx,lly,llz);

   //compute moment of inertia matrix
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   mom_inertia(r,il,ih, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",Ixx,Iyy,Izz,Ixy,Iyz,Ixz);

   fprintf(fpo,"\n");
   fclose(fpo);
   //  xxxxxx
}


// write out springs to a file
void write_springs(struct reb_simulation* const r,char *fileroot, int index)
{
   FILE *fpo;
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot);
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_springs.txt");
   fpo = fopen(filename,"w");
   for(int i=0;i<NS;i++){
      double dr = spring_length(r,springs[i]);
      fprintf(fpo,"%d %d %.10e %.10e %.10e %.10e ",
        springs[i].i, springs[i].j, 
        springs[i].ks, springs[i].rs0, 
        springs[i].gamma, springs[i].k_heat);
      fprintf(fpo,"%.10e %.10e\n",
         strain(r,springs[i]), dr);
   }
   fclose(fpo);
   printf("\n write_springs: NS=%d %s\n",NS,filename);
}


// read springs file
void read_springs(struct reb_simulation* const r,char *fileroot, int index)
{
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot);
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_springs.txt");
   printf("\n reading in springs %s\n",filename);
   FILE *fpi;
   fpi = fopen(filename,"r");
   char string[300];
   struct spring spr;
   int  i,j; 
   double ks,rs0,gamma,k_heat;
   while(fgets(string,300,fpi) != NULL){
      sscanf(string,"%d %d %lf %lf %lf %lf",
        &i,&j,&ks,&rs0,&gamma,&k_heat);
      spr.i=i; spr.j=j;
      spr.ks=ks; spr.rs0=rs0; spr.gamma=gamma; spr.k_heat=k_heat; 
      springs_add(r,spr);
   }
   fclose(fpi);
   printf("read_springs: NS=%d\n",NS);
}

// read particle information from a file, all particles!
void  read_particles(struct reb_simulation* const r,char *fileroot, int index)
{
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot);
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_particles.txt");
   printf("\n reading in particles %s\n",filename);
   FILE *fpi;
   fpi = fopen(filename,"r");
   char string[300];
   struct reb_particle pt;
   pt.ax = 0.0; pt.ay = 0.0; pt.az = 0.0;

   fgets(string,300,fpi); // read in time
   double t_dump;
   sscanf(string,"%lf",&t_dump);
   r->t = t_dump;  // set time

   double x,y,z,vx,vy,vz,m,rad;
   while(fgets(string,300,fpi) != NULL){
      sscanf(string,"%lf %lf %lf %lf %lf %lf %lf %lf\n",
        &x, &y, &z, &vx, &vy, &vz, &rad, &m);
      pt.x = x;   pt.y = y;   pt.z = z;
      pt.vx = vx; pt.vy = vy; pt.vz = vz;
      pt.r = rad; pt.m = m; 
      reb_simulation_add(r,pt);
   }
   fclose(fpi);
   printf("read_particles: N=%d\n",r->N);

}

// print out particle information to a file, all particles!
// index only determines filename number
void write_particles(struct reb_simulation* const r,char *fileroot, int index)
{
   FILE *fpo;
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot);
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_particles.txt");
   fpo = fopen(filename,"w");
   fprintf(fpo,"%.6f\n",r->t); // print time
   for(int i=0;i< r->N;i++){
      fprintf(fpo,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
        r->particles[i].x, r->particles[i].y, r->particles[i].z,
        r->particles[i].vx, r->particles[i].vy, r->particles[i].vz,
        r->particles[i].r, r->particles[i].m); 
   }
   fclose(fpo);
   printf("\n write_particles: N=%d %s\n",r->N,filename);
}

// give me an integer string
void toistring(char *istring, int i)
{
   char junks[20];
   sprintf(junks,"%d",i);
   if (i < 100000) strcpy(istring,"0");
   if (i < 10000)  strcat(istring,"0");
   if (i < 1000)   strcat(istring,"0");
   if (i < 100)    strcat(istring,"0");
   if (i < 10)     strcat(istring,"0");
   strcat(istring,junks);
}



// print out covariance matrix for an extended body indices [il,ih) 
// into filename
// first time this is called please set isfirst to 1 (that way there is a first line with labels)
// afterwards set isfirst to 0
// this is not done automatically so that this routine can be called for more than 1 extended body
// the nodes positions to compare to are given in arrays x0,y0,z0 
//  these are particle positions w.r.t to center of mass
// the covar matrix is relevant for the Kabsch algorithm to find the body rotation matrix
// https://en.wikipedia.org/wiki/Kabsch_algorithm
void print_covar(struct reb_simulation* const r, int il, int ih, 
    struct reb_particle* init_particles, char* filename,int isfirst)
{
   FILE *fpo;
   if (isfirst==1){
     fpo = fopen(filename, "w"); // first time create file
     fprintf(fpo,"# t C_xx C_xy C_xz C_yx C_yy C_yz C_zx C_zy C_zz \n");
     double xc0=0.0; double yc0=0.0; double zc0=0.0;
      // first time this is run substract off center of mass from initial particle list in [il,ih)
     for (int i=il;i<ih;i++){
        xc0 += init_particles[i].x; 
        yc0 += init_particles[i].y; 
        zc0 += init_particles[i].z; 
     }
     xc0 /= (ih-il); // center of mass position
     yc0 /= (ih-il);
     zc0 /= (ih-il);
     for (int i=il;i<ih;i++){
       init_particles[i].x -= xc0; 
       init_particles[i].y -= yc0; 
       init_particles[i].z -= zc0; 
     }
   }
   else {
     fpo = fopen(filename, "a");
   }
   fprintf(fpo,"%.3f ",r->t);
   double xc =0.0; double yc =0.0; double zc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc); // current center of mass
   double cxx =0.0; double cxy =0.0; double cxz =0.0;
   double cyx =0.0; double cyy =0.0; double cyz =0.0;
   double czx =0.0; double czy =0.0; double czz =0.0;
   struct reb_particle* particles = r->particles;
   for (int i=il;i<ih;i++){
      double x = particles[i].x - xc;
      double y = particles[i].y - yc;
      double z = particles[i].z - zc;
      double x0 = init_particles[i].x;
      double y0 = init_particles[i].y;
      double z0 = init_particles[i].z;
// notice convention here cxy is sum x_current_position *  y_original_position!
      cxx += x*x0;
      cxy += x*y0;
      cxz += x*z0;
      cyx += y*x0;
      cyy += y*y0;
      cyz += y*z0;
      czx += z*x0;
      czy += z*y0;
      czz += z*z0;
   }
   fprintf(fpo,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",cxx,cxy,cxz,cyx,cyy,cyz,czx,czy,czz);
   fclose(fpo);

}

// print out stress file for a resolved body 
// requires global stressvec
void print_stress(struct reb_simulation* const r, int il, int ih, char* fileroot, int index)
{
   // struct reb_particle* particles = r->particles;

   FILE *fpo;
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot); 
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_stress.txt");
   fpo = fopen(filename,"w");
   fprintf(fpo,"#%.6f\n",r->t); // print time
   fprintf(fpo,"#i s_xx s_yy s_zz s_xy s_yz s_xz eig1 eig2 eig3 fail\n"); 

   double xc,yc,zc;
   compute_com(r,il, ih, &xc, &yc, &zc);
   for(int i=il;i<ih;i++){
      // double x = particles[i].x - xc;
      // double y = particles[i].y - yc;
      // double z = particles[i].z - zc;
      // double m = particles[i].m;
      // fprintf(fpo,"%d %.4f %.3f %.3f %.3f " ,i,m,x,y,z);
      fprintf(fpo,"%d " ,i);
      double sigxx = stressvec[i].sigxx;
      double sigyy = stressvec[i].sigyy;
      double sigzz = stressvec[i].sigzz;
      double sigxy = stressvec[i].sigxy;
      double sigyz = stressvec[i].sigyz;
      double sigxz = stressvec[i].sigxz;
      fprintf(fpo," %.4e %.4e %.4e %.4e %.4e %.4e",sigxx,sigyy,sigzz,sigxy,sigyz,sigxz);
      double eig1  = stressvec[i].eig1;
      double eig2  = stressvec[i].eig2;
      double eig3  = stressvec[i].eig3;
      fprintf(fpo," %.4e %.4e %.4e " ,eig1,eig2,eig3);
      fprintf(fpo," %d " ,stressvec[i].fail);
      fprintf(fpo,"\n");
   }

}


