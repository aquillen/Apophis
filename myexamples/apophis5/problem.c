

/**
 * resolved mass spring model
 * using the leap frog integrator. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h" // from rebound
#include "tools.h"   // from rebound
#include "output.h"  // from rebound

#include "spring.h"  // for springs 
#include "orb.h"     // for orbit stuff
#include "m_output.h"     // for io
#include "quat.h"    // for quaturnions
// #include "stress.h"  // stress tensor  

int NS;  // global numbers of springs
struct spring* springs;  // global for springs list

double gamma_all;  // for gamma of all springs
double t_damp;     // end faster damping, relaxation
double t_print;    // for table printout 
double t_datadump; // for extended body outputs
char froot[30];    // needed for output files
int npert=0;       // number of perturbative point masses
int neph =0;       // number of ephemeris point masses
struct ephem ephem1; // ephemeris info
int stype=0;  // type of simulation: 
// =0, particle generation and relaxation 
// =1, run tidal encounter
// double nodevol;    // used by stress tensor computation
struct stresstensor *stressvec;  // stress tensor?, for some reason does not compile without

void heartbeat(struct reb_simulation* const r);

void additional_forces(struct reb_simulation* r){
   spring_forces(r); // spring forces
}

int main(int argc, char* argv[]){
        struct reb_simulation* r = reb_simulation_create();
     // This allows you to connect to the simulation display using
     // a web browser. To see it go to http://localhost:1234
        reb_simulation_start_server(r, 1234);
     // sim display starts with x to right, y up, z torward viewer?

        struct spring spring_mush; // spring parameters for mush
	// Setup constants
	r->integrator	= REB_INTEGRATOR_LEAPFROG;
	r->gravity	= REB_GRAVITY_BASIC;
	r->boundary	= REB_BOUNDARY_NONE;
	r->G 		= 6.67408E-11; // mks gravitational force
        r->additional_forces = additional_forces;  // setup callback function for additional forces
        double tmax = 0.0;  // if 0 then integrate forever

// things to set!  can be read in with parameter file
        double mball;          //  mass of resolved body 
        double dt; // timestep
        double b_distance,ks,mush_fac,gamma_fac; // spring parameters 
        double omegax,omegay,omegaz; // spin vector
        // double phi,theta,psi; // initial Euler angles   (radians)
        //    w.r.t to spin angular momentum direction 
        double lambda_L_deg; // ecliptic longitude, in degrees of angular momentum vector
        double beta_L_deg;   // ecliptic lattitude, in degrees
        unsigned int seed=1; // for random number generator
        char objfname[300]; //  vertex shape model filename
        double distcor; // multiplication factor to get shape model in m
        double shapevol; // volume of shape model in m^3 computed outside 
        char ephfname[300]; //  ephemeris filename
        double ME, RE, dt_eph, jd0,scalefac;  // relevant for ephemeris, 
           // mass, radius, timestep, julian date of t=0
        int sindex; // for reading in particles/springs for resolved body
        char saveroot[100]; // root for saved particles/springs filename
        stype = 0;
        // phi = 0.0; psi=0.0;  theta = 0.0;
        lambda_L_deg = 0.0; beta_L_deg = 0.0;

    if (argc ==1){
        strcpy(froot,"t1");   // to make output files
	dt	   = 100.;    // Timestep in s
        tmax       = 0.0;     // max integration time
        t_print    = 1000*dt;     // printouts for table
        t_datadump = 1e6*dt;     // printouts of  springs and particles
        stype = 0;   //  vanilla generate particles run
        mball      = 6e10;    // mass of extended body
        sprintf(objfname,"../shape/ApophisModel1.obj"); // is not in km 
        distcor = 300;        // 
        shapevol = 6e7;
        b_distance = 30.0;    // min separation between particles in m
        mush_fac    = 2.3;    // ratio of smallest spring distance to minimum interparticle dist
        ks          = 4e3;    // spring constant
        // spring damping
        gamma_all   = 1.0;    // final damping coeff
        gamma_fac   = 5.0;    // factor initial gamma is higher that gamma_all
        t_damp      = 1.0;    // gamma from initial gamma 
                              // to gamma_all for all springs at this time
        omegax      = 0.0;
        omegay      = 0.0;
        omegaz      = 0.2;    // initial spin
        seed = 1;
        neph = 0; // no external pt masses
        stype = 0;

     }
     else{
        FILE *fpi;
        fpi = fopen(argv[1],"r");
        char line[300];
        strcpy(froot,argv[1]); // fileroot for outputs
        // fgets(line,300,fpi);  sscanf(line,"%s",froot);  
        fgets(line,300,fpi);  sscanf(line,"%lf",&dt);     // timestep in s
        fgets(line,300,fpi);  sscanf(line,"%lf",&tmax);   // integrate to this time in s
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_print); // output timestep in s
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_datadump); // output particles and springs timestep in s
        fgets(line,300,fpi);  sscanf(line,"%d",&stype);   // what type of simulation! 
           // if s=0 then generate particles, otherwise read in from file
        fgets(line,300,fpi);  sscanf(line,"%s",saveroot); // root of filename for reading in particles 
        fgets(line,300,fpi);  sscanf(line,"%d",&sindex);  // index of saved particles/springs 
        fgets(line,300,fpi);  sscanf(line,"%lf",&mball);  // mass of extended body in kg 
        fgets(line,300,fpi);  sscanf(line,"%s",objfname); // filename for shape model .obj file
        // sprintf(objfname,"ApophisModel1.obj"); // is probably not in km 
        fgets(line,300,fpi);  sscanf(line,"%lf",&distcor);    // multiply shape model by this to get in m
        fgets(line,300,fpi);  sscanf(line,"%lf",&shapevol);    // shape model volume in m3
        fgets(line,300,fpi);  sscanf(line,"%lf",&b_distance); // min interparticle distance in m
        fgets(line,300,fpi);  sscanf(line,"%lf",&mush_fac);   // sets max spring length
        fgets(line,300,fpi);  sscanf(line,"%lf",&ks);         // spring constant
        fgets(line,300,fpi);  sscanf(line,"%lf",&gamma_fac);  // factor initial gamma is higher
        fgets(line,300,fpi);  sscanf(line,"%lf",&gamma_all);  // damping final
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_damp);     // time to switch
        fgets(line,300,fpi);  sscanf(line,"%lf",&omegax);     // initial body spin in s-1
        fgets(line,300,fpi);  sscanf(line,"%lf",&omegay);     // initial body spin
        fgets(line,300,fpi);  sscanf(line,"%lf",&omegaz);     // initial body spin
        fgets(line,300,fpi);  sscanf(line,"%lf",&lambda_L_deg); // ecliptic longitude (deg)
        fgets(line,300,fpi);  sscanf(line,"%lf",&beta_L_deg);   // ecliptic latitude (deg)
        fgets(line,300,fpi);  sscanf(line,"%u" ,&seed);       // random number seed, if 1 then use time of day
        fgets(line,300,fpi);  sscanf(line,"%d",&neph);        // number of ephemeris point masses
        if (neph>0){
          fgets(line,300,fpi);  sscanf(line,"%s",ephfname);   // ephemeris filename
          fgets(line,300,fpi);  sscanf(line,"%lf",&ME);       // mass of Earth in kg
          fgets(line,300,fpi);  sscanf(line,"%lf",&RE);       // equatorial radius of Earth in m
          fgets(line,300,fpi);  sscanf(line,"%lf",&jd0);      // time of t=0 in julian days
          fgets(line,300,fpi);  sscanf(line,"%lf",&dt_eph);   // time between eph points
          fgets(line,300,fpi);  sscanf(line,"%lf",&scalefac); // for units km to m
          ephem1.M  = ME;
          ephem1.dt_eph  = dt_eph;
          ephem1.jd0  = jd0;
          ephem1.R  = RE;
          ephem1.scalefac  = scalefac;
          read_eph(ephfname, &ephem1); // read in from ephem file
          printf("read in ephemeris from %s\n",ephfname);
        }
        

     }
     npert = 0; 
     

/// end parameters of things to set from parm file /////////////////////////

    r->dt=dt;                     // set integration timestep
    r->softening      = b_distance/100.0;	// Gravitational softening length


   // set properties of springs
   spring_mush.gamma     = gamma_fac*gamma_all; // initial damping coefficient
   spring_mush.ks        = ks; // spring constant
   double mush_distance  = b_distance*mush_fac; 
       // distance for connecting springs

   FILE *fpr;
   char fname[300];
   sprintf(fname,"%s_run.txt",froot);  // output filename, 
   // for storing some information about the simulation
   fpr = fopen(fname,"w");

   NS=0; // start with no springs 

   if (seed != 1)  // seed random number generator
      srand(seed); // specify if not 1
   int il,ih;

   // if we are reading in a previously made set of particles
   // then we only need the shape mode to find the mean radius
   read_vertex_file(r, objfname); // read in shape model
   // result is a bunch of points near surface, sizescale is not yet adjusted!
   int N_vertex = r->N; // numbers of vertex particles 
   rescale_xyz(r, 0, N_vertex, distcor);  // put shape model in m  (hopefully)
   // scales by distcor
   double mean_rad =  mean_radius(r,0,N_vertex); // mean radius of shape model
   printf("mean_radius of shape model %.2e\n",mean_rad);
   if (stype >0){
      rmshape_vertices(r, N_vertex); // remove shape model vertices
   }

// create resolved body particle distribution using shape model (which is probably in km)
   if (stype ==0){
      rand_bennu(r, b_distance, mball);  // fill shape model with internal points,
          //  separated by at least b_distance, make total mass be mball 
  
      rmshape_vertices(r, N_vertex); // remove shape model vertices

      il=0;  // index range for the single resolved body
      ih=r->N;
      subtractcom(r,il,ih); // move reference frame to resolved body
      subtractcov(r,il,ih); // center of velocity subtracted 
      printf("body made \n");
      // spin it
      spin(r,il, ih, omegax, omegay, omegaz);
      subtractcov(r,il,ih); // center of velocity subtracted 
      printf("spinning now \n");
      // make springs, all pairs connected within interparticle distance mush_distance
      connect_springs_dist(r,mush_distance, il, ih, spring_mush);
      printf("springs made\n");
      // rotate body via Euler angles zxz convention
      // rotate_body(r, il, ih,-psi,-theta,-phi); // rotate body by Euler angles 
      // reverse order and minus signs to be consistent with 
      // convention to go from body to inertial frame 
      // printf("Euler angle rotation done \n");
   }
   else { // read in a previous set of saved particles and springs
      read_particles(r,saveroot, sindex);
      read_springs(r,saveroot, sindex);
      il = 0;
      ih = r->N;
      r->t = 0;  // reset time 
   }
   // would be done in either case of simulation

   // set up display
   double rball = max_radius(r, il, ih); // find maximum radius
   const double boxsize = 2.0*rball;    // display
   // const double boxsize = 1e8;    // display
   reb_simulation_configure_box(r, boxsize, 1, 1, 1); // simulation display set up

   // add ephemeris particle to the end of the particle list for display
    if (neph>0){
        add_particle_ephem(r, &ephem1);
        printf("ephem particle added to particles list \n");
    }
    // you have to do this after resolved body is set up 
    // because it is last particle in particle array 


   double llx,lly,llz;
   measure_L(r,il, ih, &llx, &lly, &llz); // measure angular momentum of body
   // printf("L %e %e %e\n",llx,lly,llz);
   // we want to rotate body so that L is in direction set by longitude,lattitude
   double xr,yr,zr;
   // construct a unit vector for our desired L orientation
   xr = cos(lambda_L_deg*M_PI/180.0) *cos(beta_L_deg*M_PI/180.0);
   yr = sin(lambda_L_deg*M_PI/180.0) *cos(beta_L_deg*M_PI/180.0);
   zr = sin(beta_L_deg*M_PI/180.0);
   struct quaturnion q = uv_rotate(llx,lly,llz,xr,yr,zr); // get a quaturnion that will rotate u to v
   q_rotate_body(r, il, ih, q);// do the rotation!
   printf("body rotated for desired L orientation\n");
   // I checked that L was in the z direction  if r was in z direction
   // note that ephemeris and Earth particle are not rotated 
   
   // nodevol = shapevol/(ih-il);
   // create a stress tensor
   // update_stresstensor(r, nodevol);

   ///////////////////////////////////
   // print out stuff to a file and to terminal
   double ddr = 0.4*rball; // mini radius for computing Young modulus
   printf("ddr = %.3e mush_distance =%.3e\n",ddr,mush_distance);
   fprintf(fpr,"mush_distance (m) %.3e\n",mush_distance);

   double Emush = Young_mush(r,il,ih, 0.0, ddr); // compute from springs out to ddr
   printf("Youngs_modulus (Pa) %.3e\n",Emush);
   fprintf(fpr,"Youngs_modulus (Pa) %.3e\n",Emush);
    // ks in units of force/distance

   double Egrav = r->G*mball*mball/pow(mean_rad,4.0);
   printf("Egrav (Pa) %.3e\n",Egrav);
   fprintf(fpr,"Egrav (Pa) %.3e\n",Egrav);

   double LL = mean_L(r);  // mean spring length
   printf("mean_L  (m) %.3e\n",LL);
   fprintf(fpr,"mean_L (m) %.3e\n",LL);

   // relaxation timescale, gamma is a unit of 1/time 
   // note no 2.5 here!
   double tau_relax = 1.0*gamma_all*0.5*(mball/(r->N -npert -neph))/spring_mush.ks; 
                 // Kelvin Voigt relaxation time with final gamma
   // this is m gamma/k , units ok, gamma_all is t-1 ks is m/t^2
   printf("relaxation_time (s) %.3e\n",tau_relax);
   fprintf(fpr,"relaxation_time (s) %.3e\n",tau_relax);

   double t_grav = sqrt(r->G*mball/pow(mean_rad,3.0));
   printf("t_grav (s) %.3e\n",t_grav);
   fprintf(fpr,"t_grav (s) %.3e\n",t_grav);

   double rho_a = mball*3.0/(M_PI*4.0)/pow(mean_rad,3.0);
   printf("rho_a (kg/m3) %.3e\n",rho_a);
   fprintf(fpr,"rho_a (kg/m3) %.3e\n",rho_a);

   double c_s = sqrt(Emush/rho_a);
   printf("c_s (m/s) %.3e\n",c_s);
   fprintf(fpr,"c_s (m/s) %.3e\n",c_s);

   double vgrav = sqrt(r->G*mball/mean_rad);
   printf("vgrav (m/s) %.3e\n",vgrav);
   fprintf(fpr,"vgrav (m/s) %.3e\n",vgrav);

   double dtr = b_distance/c_s;
   printf("dtr (s) %.3e\n",dtr);
   fprintf(fpr,"dtr (s) %.3e\n",dtr);

   // ratio of numbers of particles to numbers of springs for resolved body
   double Nratio = (double)NS/(ih-il);
   printf("N=%d NS=%d NS/N=%.1f\n", r->N, NS, Nratio);
   fprintf(fpr,"N %d\n", r->N);
   fprintf(fpr,"NS %d\n",  NS);
   fprintf(fpr,"NS/N %.1f\n", Nratio);

   // compute spin vector, angular momentum, and moments of inertia of extended body
   double omx,omy,omz,I_big,I_middle,I_small;
   body_spin(r, il, ih, &omx, &omy, &omz, &I_big,&I_middle,&I_small);
   measure_L(r,il, ih, &llx, &lly, &llz);
   printf("Omega (s-1) (%.3e, %.3e, %.3e)\n",omx,omy,omz);
   fprintf(fpr,"Omega (s-1) (%.3e, %.3e, %.3e)\n",omx,omy,omz);
   printf("L (kg m^2/s) (%.3e, %.3e, %.3e)\n",llx,lly,llz);
   fprintf(fpr,"L (kg m^2/s) (%.3e, %.3e, %.3e)\n",llx,lly,llz);
   printf("I_mid/I_big %.3f, I_small/I_big %.3f\n",I_middle/I_big,I_small/I_big);
   fprintf(fpr,"I_mid/I_big %.3f, I_small/I_big %.3f\n",I_middle/I_big,I_small/I_big);

   fclose(fpr);


   r->heartbeat = heartbeat;
   // centerbody(r,il,ih);  // move reference frame to resolved body
   // this substracts resolved body com from all particles! don't do it because
   // of ephemeris particle!

   // max integration time
   if (tmax ==0.0)
      reb_simulation_integrate(r, INFINITY);
   else
      reb_simulation_integrate(r, tmax);

   reb_simulation_free(r);  // free up
}


void heartbeat(struct reb_simulation* const r){
    int il = 0; // resolved body index range
    int ih = r->N -npert - neph;

    // stuff to do every timestep
    if (reb_simulation_output_check(r,r->dt)) {
        subtractcom(r,il,ih); // move resolved body to origin
        subtractcov(r,il,ih); // move center of velocity of resolved body to orgin
        if (neph>0){
            int estatus = update_ephem_pos(r,&ephem1); // move ephemeris body position
            if (estatus == 0) // stop simulation if we are outside ephemeris time window 
                reb_simulation_stop(r); 
        }
    }

    if (reb_simulation_output_check(r,10.0*r->dt)){
            reb_simulation_output_timing(r, 0);  // print time
    }
    if (fabs(r->t - t_damp) < 0.9*r->dt) set_gamma(gamma_all,0,r->N-npert-neph); 
            // damp initial bounce only 
            // reset gamma only at t near t_damp
	
    // saves for restarts and for analyzing particles
    if (reb_simulation_output_check(r,t_datadump)) {
        int ii = (int)(r->t/t_datadump);
        write_springs(r,froot, ii);
        write_particles(r,froot, ii);
        // update_stresstensor(r, nodevol);
        // print_stress(r,il,ih,froot, ii);
    }

    static int first=0;
    static char extendedfile[50]; 
    static char covarfile[50];
    static char pmfile[50];
    static struct reb_particle* init_particles;  // for storing initial positions of particles 
    if (first==0){
        first=1;
        sprintf(extendedfile,"%s_ext.txt",froot);
        print_extended(r,il,ih,extendedfile,1);
        
        if (stype==1){
           sprintf(pmfile,"%s_E.txt",froot); // emphem position file
           print_pm(r, r->N-1, 0, pmfile);
           sprintf(covarfile,"%s_cov.txt",froot); // covariance matrix file
           init_particles = malloc(r->N*sizeof(struct reb_particle));
           store_xyz0(r,init_particles); // store positions in arrays
           // struct reb_particle  pt = init_particles[0];
           // printf("%lf",pt.x);
           print_covar(r, il,ih,init_particles, covarfile, 1);
        }
    }
    else{
        if (reb_simulation_output_check(r,t_print)) {
           print_extended(r,il,ih,extendedfile,0);
           if (stype==1){
              print_pm(r, r->N-1, 0, pmfile);
              print_covar(r, il,ih,init_particles, covarfile, 0);
           }
        }
    }


}






