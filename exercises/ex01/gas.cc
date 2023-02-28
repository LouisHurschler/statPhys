#include <cmath>
#include <cstdlib>
#include <cstdio>
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
// SIMULATION PARAMETERS
//
const int    num_atoms = 200;    // number of atoms
const double box_edge  = 50.0;   // edge of the square box [nm]
const double ini_temp  = 300.0;  // initial temperature [K]
//
const int    num_steps = 5000;   // number of time steps
const double time_step = 0.5;    // time-step size [ps]
//
const double mas_atom  = 40.0;   // mass of the atom [g/mol] -> about 40 for Argon
const double rad_atom  = 0.2;    // radius of the atom [nm]  -> about 0.2 for Argon 
//
const int    rand_ini  = 230570; // random-number generator seed
//
const int print_frequency = 100; // frequency of printout
//
// DERIVED QUANTITIES
//
const double cst_pi        = 4.0*atan(1.0); // circular constant
const double cst_avogadro  = 6.02214076E23; // Avogadro constant [1/mol]
const double cst_boltzmann = 1.380649E-23;  // Boltzmann constant [J/K]
const double cst_idealgas  = cst_boltzmann*cst_avogadro/1000.0; // ideal-gas constant [kJ/(mol K)]
//
const double box_area      = box_edge*box_edge; // area of the square box [nm^2]
const double ini_kine      = num_atoms*cst_idealgas*ini_temp; // initial kinetic energy [kJ/mol]
const double ini_vel       = sqrt(2.0*ini_kine/(num_atoms*mas_atom)); // initial velocity norm given to all atoms [nm/ps]
//
// Notes:
//
// (1) With the above units, we can calculate the kinetic energy directly
//     in [kJ/mol] as kin = 0.5*mas*vel^2, with mas in [g/mol] and vel in [nm/ps],
//     as we have: x 10^6 for nm^2/ps^2 => m^2/s^2, x 10^-3 for g/mol => kg/mol,
//     and x 10^-3 for J/mol => kJ/mol
//
// (2) The conversion between kinetic energy in [kJ/mol] and temperature in [K]
//     is given by the equipartition theorem (with two degrees of freedom per
//     atom in 2D), as kin = num_atoms*cst_idealgas*temp with cst_idealgas
//     in [kJ/(mol K)]
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
// INCLUDE VISUALIZATION VARIABLES & FUNCTIONS (ACTIVE WITH -DVISUAL COMPILATION)
//
#include "visual.h"
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
// ALREADY WRITTEN BY PHIL
void   init_config( double crd[], double vel[] );
double process_wall_collisions( double crd[], double vel[]);
int    process_atom_collisions( double crd[], double vel[] );
void   propagate_in_time( double crd[], double vel[]) ;
//
// TO BE WRITTEN BY THE STUDENTS
double calculate_temperature( double vel[] );
double calculate_pressure ( double kick );
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
int main() {
  //
  // allocate configuration arrays for the num_atoms atoms
  double crd[2*num_atoms]; // array of atomic coordinates [nm]
  double vel[2*num_atoms]; // array of atomic velocities [nm/ps]
  //
  // initialize the random-number generator and the visualization
  srand(rand_ini);
  init_visual();
  // assign random initial coordinates and velocities
  // (for the velocities, we give all atoms the same velocity norm
  // compatible with the initial kinetic energy, but random velocity
  // orientations
  init_config(crd, vel);
  //
  // prepare for the cumulative averaging of the pressure
  double pres_cumul = 0.0;
  //
  bool early_quit = false;
  //
  // loop over time-steps
  for ( int istep = 0; istep < num_steps; istep++ ) {
    //
    // update time value
    double time = istep * time_step;
    //
    // update visualization
    draw_visual( crd );
    //
    // process possible elastic collisions between the atoms
    int num_collisions = process_atom_collisions(crd, vel);
    //    
    // process possible elastic collisions of atoms with the box wall
    // (and keep track of the "kick", i.e. the sum of the
    // normal outwards-directed velocities)
    double kick = process_wall_collisions(crd, vel);   
    //
    // calculate the instantaneous temperature and pressure
    // => for now, these functions return zero
    // => the functions are to be implemented by the students!!!
    double temp = calculate_temperature(vel);
    double pres = calculate_pressure(kick);
    //
    // cummulate the pressure for averaging
    pres_cumul += pres;
    //
    // report results ever so many steps
    if ( istep != 0 && !(istep % print_frequency) )
      printf ("t [ps] %-10.0f / Col %-2d / T [K] %-10.1f / P [kJ/(mol nm^2)] cur %-10.5f ave %-10.5f\n",
	      time, num_collisions, temp, pres, pres_cumul/(istep+1));
    //
    // propagate the coordinates using the velocities
    propagate_in_time(crd, vel) ;
    //
    // check for a possible termination by clicking off the window
    if ( early_quit = wait_visual() ) break;
    //
  }
  //
  // if not already done, wait for the user to click off the window
  if ( !early_quit ) wait_for_close();
  //
  // terminate visualization
  kill_visual();
  //
  // THAT'S ALL FOLKS !!!
  return 0;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
void init_config( double crd[], double vel[] ) {
  //
  // Generate initial atomic coordinates at random within the rectangular box
  // Note that atoms may be too close to a wall (i.e. within less
  // than one atomic radius) and/or to each other (i.e. within less
  // than two atomic radii), but these issues will be solved by the
  // first calls to process_atom_collisions and process_wall_collisions
  //
  for ( int n = 0; n < num_atoms; n++ ) {
    crd[2*n  ] =  box_edge * (double) rand() / RAND_MAX;
    crd[2*n+1] =  box_edge * (double) rand() / RAND_MAX;
  }
  //
  // Generate initial atom velocities using the same velocity norm for all atoms
  // (appropriate for the chosen initial kinetic energy) along with random directions
  // Note that this does not math the realistic (Maxwell-Boltzmann) distribution
  // of velocities, but this incorrect distribution ill relax to the correct
  // one very rapidly when the atoms start to collide
  //
  for ( int n = 0; n < num_atoms; n++ ) {
    double theta = 2.0*cst_pi * (double) rand() / RAND_MAX;
    vel[2*n  ] = ini_vel * cos(theta);
    vel[2*n+1] = ini_vel * sin(theta);      
  }
  //
  return;
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

double process_wall_collisions( double crd[], double vel[] ) {
  //
  // Process collisions of atoms with any of the four boxes of the wall.
  // A collision happens if the atom is with in less than its radius from the wall.
  // In this case, we reset the position to be exactly the radius, and we
  // revert the normal component of the atomic velocity. This corresponds
  // to a hard-sphere with hard-wall elastic shock, where the
  // kinetic energy is conserved and a momentum amounting to 2 m v_n
  // is transferred to the wall, with m the mass and v_n the normal velocity
  // directed towards the outside of the box. She sum of all v_n
  // (i.e. the "kick") is returned by the function.
  //
  double kick = 0.0;
  //
  // loop over all atoms n
  for ( int n = 0; n < num_atoms; n++ ) {
    //    
    double x = crd[2*n];
    double y = crd[2*n+1];
    //
    if ( x < rad_atom ) {
      // collision with the left wall at x = 0 
      double dx = rad_atom - x;
      crd[2*n] = rad_atom;
      vel[2*n] = -vel[2*n];
      kick += vel[2*n];
    }
    if ( x > box_edge - rad_atom ) {
      // collision with the right wall at x = box_edge
      double dx = x - ( box_edge - rad_atom );
      crd[2*n] = box_edge - rad_atom;
      kick += vel[2*n];
      vel[2*n] = -vel[2*n];
    }
    if ( y < rad_atom ) {
      // collision with the top wall at y = 0
      double dy = rad_atom - y;
      crd[2*n+1] = rad_atom;
      vel[2*n+1] = -vel[2*n+1];
      kick += vel[2*n+1];
    }
    if ( y > box_edge - rad_atom ) {
      // collision with the bottom wall at y = box_edge
      double dy = y - ( box_edge - rad_atom );
      crd[2*n+1] = box_edge - rad_atom;
      kick += vel[2*n+1];
      vel[2*n+1] = -vel[2*n+1];
    }
  }
  //
  return  kick;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
int process_atom_collisions( double crd[], double vel[]) {
  //
  // Process collisions between atoms. A collision happens when two atoms are
  // within less than two radii away from each other.
  // In this case, we reset the position along the interatomic vector to
  // be exactly the radius, and we swap the velocity component of the two
  // atoms along the interatomic vector. This corresponds hard-sphere
  // elastic shock, where the kinetic energy and the momentum are both
  // conserved. The function returns the number of collisions that occurred.
  //  
  int num_atom_collisions = 0;
  //
  // loop over all atom pairs (m,n)
  for ( int m = 0; m < num_atoms - 1; m++ ) {
    for ( int n = m + 1; n < num_atoms; n++ ) {
      //
      double dx = crd[2*n  ] - crd[2*m  ];
      double dy = crd[2*n+1] - crd[2*m+1];
      double d2 = dx*dx + dy*dy;
      double d = sqrt(d2);
      //
      if ( d < 2.0*rad_atom ) {
	//
	// there is an atom-atom collision
	//
	// reset position to contact
	double fac = ( 2*rad_atom - d ) / ( 2.0*d );
	crd[2*m  ] -= fac * dx;
	crd[2*m+1] -= fac * dy;
	crd[2*n  ] += fac * dx;
	crd[2*n+1] += fac * dy;
	//
	// perform elastic collision
	double prjm = vel[2*m]*dx/d + vel[2*m+1]*dy/d;
	double prjn = vel[2*n]*dx/d + vel[2*n+1]*dy/d;
	vel[2*m  ] += (prjn - prjm) * dx/d;
	vel[2*m+1] += (prjn - prjm) * dy/d;
	vel[2*n  ] += (prjm - prjn) * dx/d;
	vel[2*n+1] += (prjm - prjn) * dy/d;
	//
	num_atom_collisions++;
	//	
      }
    }
  }
  //
  return  num_atom_collisions;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
void propagate_in_time( double crd[], double vel[] ) {
  //
  // propagate the atomic coordinates forward by one time-step
  // in free-flight based on the atomic velocities
  //
  for ( int n = 0; n < num_atoms; n++ ) {
    crd[2*n  ] += vel[2*n  ] * time_step;
    crd[2*n+1] += vel[2*n+1] * time_step;
  }
  //
  return;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
double calculate_temperature ( double vel[] ) {
  //
  // calculate the instantaneous temperature
  // => for now, this function returns zero
  // => the function is to be implemented by the students!!!
  //
  return 0.0;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
double calculate_pressure ( double kick ) { 
  //
  // calculate the instantaneous pressure
  // => for now, this function returns zero
  // => the function is to be implemented by the students!!!
  //
  return 0.0;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
