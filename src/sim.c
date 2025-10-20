#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"


double standard_error(double x, double x2, int nblock)
{
  double X = x/nblock;
	double X2 = x2/nblock;

	double stderr = sqrt((X2 - X*X)/(nblock-1));

	return stderr;
}

void print_standard_error(char *str, double x, double x2, int nblock)
{
  double s;

  s = standard_error(x, x2, nblock);
  printf("%s%g +/- %g\n", str, x / nblock, s);
  return;
}

void print_flucation(char *str, double x, double x2, int nblock){
	double X = x/nblock;
	double X2 = x2/nblock;

	double fluc = sqrt((X2 - X*X));

	printf("%s%g", str, fluc);
}

// Construct file names based on the parameters of the run, e.g.
// 0064_r0.500_T1.000_alpha0.10_dt010
// Note that sprintf returns the number of characters written to the string.
// When the format specifier "%4d" gives "  64" the format "%4.4d" instead
// gives "0064". (We don't want spaces in the file names.)
char *get_filename(Par *par)
{
  char *fname = malloc(FNAMESIZE);
  char *f = fname;

  f += sprintf(f, "%4.4d", par->n);
  f += sprintf(f, "_r%5.3f", par->rho);
  f += sprintf(f, "_T%5.3f", par->T);
#ifndef MC		// neither alpha nor deltat in an MC simulation
  f += sprintf(f, "_alpha%4.2f", par->alpha);
  f += sprintf(f, "_dt%3.3d", (int) rint(1000.0 * par->deltat));
#endif
  
  if (f - fname >= FNAMESIZE) {		// f - fname is the size of the string
    fprintf(stderr, "Error: too long file name: %s...\n", fname);
    exit(EXIT_FAILURE);
  }

  return fname;
}



// Determine size of the simulation cell in D dimensions
// from number of particles and density.
void size_from_rho(Par *par)
{
  int d;

  par->vol = par->n / par->rho;
  for (d = 0; d < D; d++)
    par->L[d] = exp(log(par->vol) / D);
}


// Concatenate two strings and return a pointer to this new string
char *add_strings(char *str1, char *str2)
{
  char *str = malloc(strlen(str1) + strlen(str2) + 1);
  strcpy(str, str1);
  strcat(str, str2);
  return str;
}



double *run_simulation(Par *par, double *atoms)
{
  int i, itherm, isamp, iblock;
  char *filename;
  FILE *estream = NULL;
  double *force = NULL;
  Measured *val = NULL, *vblock, *v1sum, *v2sum;
  double *dval, *dvblock, *dv1sum, *dv2sum;
  int iv, nval = sizeof(Measured) / sizeof(double);
#ifdef MC
  int naccept = 0;
#else
  int istep;
  int nstep = rint(1.0 / par->deltat);
#endif

  /*** 1. Write out some text to show what we are doing ***
#ifdef MC
    printf("Monte Carlo for %d particles at T = %g, rho = %g, L = %5.3f,  b = %g\n", par->n, par->T,
	   par->rho, par->L[0], par->b);
#else
    if (par->alpha > 0.0)
      printf("\n--- Langevin dynamics of a Lennard-Jones gas ---\n\n");
    else
      printf("\n--- Molecular dynamics of a Lennard-Jones gas ---\n\n");

    printf("Gas with %d particles at rho = %g, T = %g, alpha = %g, deltat = %g, L = %5.3f\n",
	   par->n, par->rho, par->T, par->alpha, par->deltat, par->L[0]);
#endif

    
  /*** 2. Initialize ***/
  init_ran(par->seed);
  if (!force)
    force = malloc(D * par->n * sizeof(double));

  // To store measured values
  // The same value may be accessed as val->epot or dval[0],
  // the second variable is convenient for looping over the measured quantities.
  if (!val) {
    val  = calloc(nval, sizeof(double));	// Values from a single measurement
    vblock  = calloc(nval, sizeof(double));	// Accumulate values from a block with nsamp measurements
    v1sum = calloc(nval, sizeof(double));	// Accumulate block averages
    v2sum = calloc(nval, sizeof(double));	// Accumulate block averages squared (for error estimates)

    dval = (double *) val;
    dvblock = (double *) vblock;
    dv1sum = (double *) v1sum;
    dv2sum = (double *) v2sum;
  }

  /*** 3. Open file for writing out the total energy ***/
  filename = get_filename(par);		  // Get file name from parameters
  estream = fopen(add_strings("efile/", filename), "w");
  if (!estream) {
    fprintf(stderr, "*** Cannot open efile/%s for write\n", filename);
    return atoms;
  }

  /*** 4. We need some configuration to start from ***/
  if (!atoms) {
    // Attempt to read from a file for the same parameter values
    atoms = read_conf(par, atoms, filename, 1);
    if (atoms)		// If successful there is no need to equilibrate
      par->ntherm = 0;
    else
      atoms = init_conf(par);	// Start from a simple configuration
  }

  double *pos = atoms;
#ifdef VEL
  double *vel = atoms + D * par->n;
#endif

  /*** This section is only for the initial steps of the lab ***/
  // Measure properties of the initial configuration (see measure.c)
  measure(par, atoms, val);
  printf("Potential energy = %g\n", val->epot);
#ifdef VEL
  printf("Kinetic energy = %g\n", val->ekin);
#endif

  step(par, atoms, force);	
	// // Fix this (3): remove printf and exit statements
  // printf("x, y, vx, vy = %g  %g  %g  %g\n", pos[0], pos[1], vel[0], vel[1]);
  // printf("...exiting\n");
  // exit(EXIT_SUCCESS);

  
  /*** 5. Equilibrate ***/
  if (par->ntherm) {
    printf("Equilibrate: %d...", par->ntherm);
    fflush(stdout);

    for (itherm = 0; itherm < par->ntherm; itherm++) {
#ifdef MC
      monte_carlo_sweep(par, atoms);
#else    
      // advance one unit of time (since nstep = 1/deltat)
      for (istep = 0; istep < nstep; istep++)
		step(par, atoms, force);
#endif
    }
    printf("done\n");
  // Write configuration to the named file in "conf/"
    write_conf(par, atoms, filename);
  }


  // Initialize variables
  for (iv = 0; iv < nval; iv++)
    dval[iv] = 0.0;

  /*** 6. Production part ***/

  printf("\nSimulate %d blocks x %d samples each: ", par->nblock, par->nsamp);
  fflush(stdout);
  // This is the main loop: calculate values for nblock blocks with nsamp samples each.
  for (iblock = 0; iblock < par->nblock; iblock++) {
		for (iv = 0; iv < nval; iv++) dvblock[iv] = 0.0;

    for (isamp = 0; isamp < par->nsamp; isamp++) {
#ifdef MC
      naccept += monte_carlo_sweep(par, pos);	// One MC sweep
#else
      for (istep = 0; istep < nstep; istep++)	// This advances one unit of time
	step(par, atoms, force);
#endif

      // Measure all kinds of things and accumulate results (file: measure.c).
      measure(par, atoms, val);
      
		if (isamp % 100 == 0)	// Write energy from a single measurement to file, but not too often
		fprintf(estream, "%g\n", val->etot);

      for (iv = 0; iv < nval; iv++)
		dvblock[iv] += dval[iv];
    }

    // These values should be average of par->nsamp measurements
    for (iv = 0; iv < nval; iv++) dvblock[iv] /= par->nsamp;
    
    for (iv = 0; iv < nval; iv++) dv1sum[iv] += dvblock[iv];
    for (iv = 0; iv < nval; iv++) dv2sum[iv] += dvblock[iv] * dvblock[iv];

    // Print out numbers to show that things are progressing
    printf("%d ", iblock + 1);	fflush(stdout);
  }	// End of loop "for (iblock..."
  printf("\n");

  if (estream) fclose(estream);

  print_standard_error("Potential E: ", v1sum->epot, v2sum->epot, par->nblock);
  print_standard_error("Pressure   : ", v1sum->pressure, v2sum->pressure, par->nblock);
#ifdef VEL
  print_standard_error("Kinetic E  : ", v1sum->ekin, v2sum->ekin, par->nblock);
  print_standard_error("Total E    : ", v1sum->etot, v2sum->etot, par->nblock);
#endif
  // Fix this (4): calculate and print out \delta E = sqrt(<E^2> - <E>^2)
	print_flucation("Sigma E    : ", v1sum->etot, v1sum->etot2, par->nblock);

#ifdef MC
  printf("Acceptance ratio = %g\n", naccept / (par->n * 1.0 * par->nsamp * par->nblock));
#endif

  // Write configuration to the named file in "conf/"
  write_conf(par, atoms, filename);

  if (force) free(force);

  return atoms;
}

// The function read_args mostly sets values to variables in the par struct.
// The string arg should contain strings like "N=50" and "T=1.5".
// The exception is "run" which makes the simulation start.
// Return value: 1 for success, 0 for failure.
int read_args(Par *par, char *arg)
{
  static double *atoms = NULL;
  char *s;

  
  // This is where the simulation is started  
  if (!strcmp(arg, "run")) {
    atoms = run_simulation(par, atoms);
    return 1;
  }

  // Most arguments have an '=' sign, if so then
  // (1) make s point to the part after the '=',
  // (2) change the '=' to the NULL character to mark end of string for arg.
  s = strchr(arg, '=');
  if (s)		// Let s point at string after '=', if that char is found
    *s++ = '\0';

  
  // Read a config from directory "conf", with file name in the string s.
  // "read=filename" tries to open the specified file,
  // "read" (with no file name) constructs a file name from the parameters.
  if (!strcmp(arg, "read")) {
    atoms = read_conf(par, atoms, s, 0);	// Defined in config.c

    if (atoms)	// Successful?
      return 1;
    else
      return 0;
  }

  // If there is no '=' in arg (and arg was neither "run" nor "read", handled above)
  // then something is wrong. Write message and return.
  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'\n", arg);
    return 0;
  }

  // Number of particles
  if (!strcmp(arg, "N")) {
    par->n = strtol(s, NULL, 10);
    return 1;
  }

  // Density
  if (!strcmp(arg, "rho")) {
    double oldrho = par->rho;
    par->rho = strtod(s, NULL);
    size_from_rho(par);

    if (atoms) {		// Rescale existing configuration
      double *pos = atoms;
      double fact = exp(log(par->rho/oldrho) / D);
      int id;
      for (id = 0; id < D * par->n; id++)
	pos[id] /= fact;
    }
    return 1;
  }

  // Temperature
  if (!strcmp(arg, "T")) {
    par->T = strtod(s, NULL);
    return 1;
  }

  
  // Number of blocks (length of the run)
  if (!strcmp(arg, "nblock")) {
    par->nblock = strtol(s, NULL, 10);
    return 1;
  }

  // Length of the termalization
  if (!strcmp(arg, "ntherm")) {
    par->ntherm = strtol(s, NULL, 10);
    return 1;
  }

  // Number of samples in each block
  if (!strcmp(arg, "nsamp")) {
    par->nsamp = strtol(s, NULL, 10);
    return 1;
  }

  // Seed for the random number generator (which is otherwise taken from the clock).
  if (!strcmp(arg, "seed")) {
    par->seed = strtol(s, NULL, 10);
    return 1;
  }

#ifdef MC
  // Maximum trial displacement of the particle (in each direction).
  if (!strcmp(arg, "b")) {
    par->b = strtod(s, NULL);
    return 1;
  }
#else
  // Friction parameter for Langevin dynamics.
  if (!strcmp(arg, "alpha")) {
    par->alpha = strtod(s, NULL);
    return 1;
  }
  // Size of the time step.
  if (!strcmp(arg, "deltat")) {
    par->deltat = strtod(s, NULL);
    return 1;
  }
#endif

  fprintf(stderr, "No such variable name: '%s'.\n", arg);
  return 0;
}

// This is where the program starts.
int main(int argc, char *argv[])
{
  int iarg;
  Par par;

  // Number of degrees of freedom per particle = df.
  // For both pos and vel: 2D, for pos only: D.
#ifdef VEL
  par.df = 2 * D;
#else
  par.df = D;
#endif
  par.n = 64;		// Initial values for some parameters.
  par.nblock = 10;
  par.seed = 0;
  par.ntherm = 1000;
  par.nsamp = 1000;
  par.mass = 1;
  par.rho = 0.6;
  size_from_rho(&par);	// To set par->L
#ifdef MC
  par.b = 0.4;
#else
  par.alpha = 0.0;
  par.deltat = 0.01;
#endif

  // Print information text if the program is started without arguments (i.e. argc==1).
  if (argc == 1) {
#ifdef MC
    printf("Usage: %s N=50 rho=0.5 T=0.1 alpha=1.0\n", argv[0]);
    printf("Optional arguments (with defaults) b=%g nblock=%d nsamp=%d ntherm=%d seed=%d\n",
	   par.b, par.nblock, par.nsamp, par.ntherm, par.seed);
#else
    printf("Usage: %s N=50 L=14.0 T=0.1 alpha=1.0\n", argv[0]);
    printf("Optional arguments (with defaults) deltat=%g nblock=%d nsamp=%d ntherm=%d seed=%d\n",
	   par.deltat, par.nblock, par.nsamp, par.ntherm, par.seed);
#endif
    exit(EXIT_SUCCESS);
  }

  // Interpret the commands given in the argument list.
  // The argument "run" will start the simulation (see "run_simulation" above).
  for (iarg = 1; iarg < argc; iarg++)
    if (!read_args(&par, argv[iarg]))
      exit(EXIT_FAILURE);

  printf("\n");
  exit(EXIT_SUCCESS);
}

