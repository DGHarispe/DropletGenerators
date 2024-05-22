#include "distance.h"
#include "embed.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"

#include "tension.h"

#include "tag.h"

#include "output_htg_MPI.h"

#include "paramsDict.h"

scalar f0[];
scalar d[];

void fraction_from_stl (FILE * fp)
{

  bool restored = false;
  if(restore (file = "./output/dump-origIm", list = {d})){
	  fprintf(stderr, "Dumped distance read!\n");
	  restored = true;
  }
  
  coord * p = input_stl (fp);

  if (!restored){
  	distance (d, p);
  	while (adapt_wavelet ({d}, (double[]){5e-4*L0}, maxRlevel).nf);
  }

  vertex scalar phi[];
  foreach_vertex(){
    	phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
				d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  }
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
}

float u1;
float u2;

float lower_x = 0.0;
float lower_z = 0.0;
float higher_x = 0.0;
float higher_z = 0.0;

bool new_tt = false;

int main(int argc, char *argv[]) {
	
  if (!readParams())	//reads all the parameters and constants needed from parameters.txt (# for comments)
	  exit(1);
  
  coord * p = input_stl (fopen ("model.stl", "r"));
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  
  foreach_dimension()
  if (max.x - min.x > maxl)
  maxl = max.x - min.x;
  
  L0 = ownL0;

  size (size_domm);

  printf("The domain size is: %g!\n", size_domm);
    
  lower_x = oxc;
  higher_x = oxc + size_domm;
  
  origin (oxc, oyc, ozc);
	
  higher_z = ozc + size_domm;
	
  printf("The coordinates are: x:%g y:%g z:%g!\n", oxc, oyc, ozc);
  
  N = NCells;
	
  //variables with 1 is dispersed phase
  //variables with 2 is continuous phase
  u1 = (q1/((cwidth*cheight)));
  u2 = q2/((cwidth*cheight));

  mu1 = mu1_val;
  mu2 = mu2_val;
  rho1 = rho1_val;
  rho2 = rho2_val;

  f.sigma = SIGMA;

  DT = DTi;

  printf("Starting!\n");

  run();
}


u.n[front]  = dirichlet( -u1 * cs[]);
u.t[front]  = dirichlet(0.);
u.r[front]  = dirichlet(0.);
p[front]    = neumann(0.);
pf[front]   = neumann(0.);

f[front] = 1.;

u.n[left]  = dirichlet( u2 * cs[]);
u.t[left] = dirichlet(0);
u.r[left] = dirichlet(0);
p[left] = neumann(0);
pf[left] = neumann(0);

f[left] = 0.;

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
u.r[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[back] = neumann(0.);
u.t[back] = neumann(0.);
u.r[back] = neumann(0.);
p[back]   = dirichlet(0.);
pf[back]  = dirichlet(0.);


u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);
//p[embed] = neumann(0.);
f[embed] = 0.;


event init (t = 0)
{
  FILE * fp = fopen ("model.stl", "r");

  fraction_from_stl (fp);
  fclose (fp);

  refine ((z>higher_z-(size_domm/pow(2,maxRlevel))*4) && cs[] && level < maxRlevel);
	
  foreach() {
	  foreach_dimension()
	  	u.x[] = 0.0;
  }
}

int fileN = 0;
event output (t = 0.00000001; t += 0.0001; t <= tmax) {
	char path[]="./output";	char prefix[80];
	sprintf(prefix, "fields-%i", fileN+7000000);

	vector mpdl[];  vector mpdu[];
	scalar mpdl_x[];  scalar mpdl_y[];  scalar mpdl_z[];
	foreach(){
	mpdl_x[] = u.x[];	mpdl_y[] = u.y[];	mpdl_z[] = u.z[];
	}

	output_htg((scalar *) {f, cs, p, mpdl_x, mpdl_y, mpdl_z}, (vector *) {fs}, path, prefix, i, t);
	printf("Fields written! (%i)\n", fileN);
	fileN = fileN + 1;
}


event timePost (i++)
{
  
  fprintf(stderr, "%d %g %g\n", i, t, dt);
  
  if ((i > startSteps) && !new_tt){
	  DT = nMaxDT;
	  new_tt = true;
  }

}

#define IN_OUT_LET ((x<lower_x+(size_domm/pow(2,maxRlevel))*4) || (x>higher_x-(size_domm/pow(2,maxRlevel))*4) || (z>higher_z-(size_domm/pow(2,maxRlevel))*4))


event adapt (i++) {
	scalar region[];
	foreach(){
		if (cs[]){
			region[] = ((IN_OUT_LET) || fabs(d[])<((size_domm/pow(2,maxRlevel)))/2)*noise();
		}
	}

	adapt_wavelet ({cs,region,u,f}, (double[]){csTol,0.05,uemax,uemax,uemax,fTol}, maxRlevel);
}

event droplets (t += 0.0005)
{
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  int n = tag (m);
  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach (serial)
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f[]*p.x;
    }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  for (int j = 0; j < n; j++)
    printf ("drops data: %d %d %g %g %g\n", i,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
}