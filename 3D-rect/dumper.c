#include "distance.h"

#include "output_htg_MPI.h"

#include <stdio.h>

#include "paramsDict.h"

#define IN_OUT_LET ((x<lower_x+(size_domm/pow(2,maxRlevel))*4) || (x>higher_x-(size_domm/pow(2,maxRlevel))*4) || (z>higher_z-(size_domm/pow(2,maxRlevel))*4))

int main()
{
	
	if (!readParams())
		exit(1);

	//DT = 1e-14;
	
	L0 = ownL0;
	
	init_grid (16);
	
	size (size_domm);
	
	origin (oxc, oyc, ozc);
	
	float lower_x = oxc;
    float higher_x = oxc + size_domm;
	float higher_z = ozc + size_domm;

	scalar d[];
	
	coord * p = input_stl (fopen ("model.stl", "r"));
	
	distance (d, p);

	while (adapt_wavelet ({d}, (double[]){5e-4*L0}, maxlevel = maxRlevel).nf);
	
	refine(fabs(d[])<(size_domm/pow(2,maxRlevel))*2 && level < maxRlevel);
	refine (IN_OUT_LET && (d[]>-(size_domm/pow(2,maxRlevel))) && level < maxRlevel);
	
	dump (file = "./output/dump-origIm", list={d});
}
