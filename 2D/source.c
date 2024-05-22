#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

#include <stdio.h>

#include "paramsDict.h"
#include "output_htg_MPI.h"

float u_in;
float u_out;
scalar f0[];

u.n[left]  = dirichlet(f0[]*(u_in) + (1-f0[])*(u_out));
u.t[left]  = dirichlet(0);
p[left]    = neumann(0);
f[left]    = f0[];

u.n[right] = neumann(0);
p[right]   = dirichlet(0);

int main (int argc, char * argv[])
{

  if (!readParams())
  	exit(1);
  
  init_grid (NCells);
  size (size_domm);

  u_in = (q_in/(pi*sq(radius_in)));
  u_out = q_out/(pi*sq(radius_out));
  
  rho1 = rho_in;	    // in
  rho2 = rho_out;		// out
  mu1 = mu_in;
  mu2 = mu_out;
  f.sigma = SIGMA;

  DT = maxDT;
  
  run();
}


event init (t = 0) {
  if (!restore (file = "restart")) {
    refine (x < 2*radius_in && sq(y) + sq(z) < 2.*sq(radius_in) && level < maxlevel);
    fraction (f0, sq(radius_in) - sq(y) - sq(z));
    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0});
	
    foreach() {
      f[] = f0[]*(x < radius_in);
      u.x[] = f[];
    }
    boundary ({f,u.x});
  }
}

event logfile (i++) {
  if (i == 0)
    fprintf (stderr,
	     "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf (stderr, "%g %g %d %d %d %ld %g %g\n", 
	   t, dt, mgp.i, mgpf.i, mgu.i,
	   grid->tn, perf.t, perf.speed);
}

//event snapshot (t = 1e-5; t += 1e-5; t <= 0.1) {
//  char name[80];
//  sprintf (name, "snapshot-%g", t);
//  dump (name);
//}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){fTol,uemax,uemax,uemax}, maxlevel);
}

int fileN = 0;
event output (t = 0.000001; t += 0.0001; t <= tmax) {
	  char path[]="./output";	char prefix[80];
	  sprintf(prefix, "fields-%i", fileN+7000000);

	  vector mpdl[];  vector mpdu[];
	  scalar mpdl_x[];  scalar mpdl_y[];  scalar mpdl_z[];
	  foreach(){
		mpdl_x[] = u.x[];	mpdl_y[] = u.y[];	mpdl_z[] = u.z[];
	  }

	  vector * vlist;
	  output_htg((scalar *) {f, p, mpdl_x, mpdl_y, mpdl_z}, (vector *) {vlist}, path, prefix, i, t);
	  printf("Fields written! (%i)\n", fileN);
	  fileN = fileN + 1;
}