#include"main.h"


// ionized imprity scattering calculation
double nu_ii_p_funct(double k, int counter, double beta_constant, double epsilon_s)
{
	// equation 3.29 of Ramu Thesis
	double A= 0.5*(2 + (beta_constant*beta_constant/(k*k)));
	double B= abs((A+1)/(A-1));
	
	double iiA= (pow(e,4)*m_h* N_ii)/(32*pi*pow(k,3)*epsilon_s*epsilon_s*epsilon_0*epsilon_0*pow(h_bar,3));
	double iiB= ((3*A -1)*(3*A -1)*log(B) - 18*A +12 -8/(A+1));
	
	double ii_p= iiA*iiB ;
	
	FILE *fid1;
	fid1= fopen("nu_iip.dat","w");
	for(int i=0;i<points;i++){
		fprintf(fid1,"%e     %e\n", energy_p[i], ii_p);
	}
	fclose(fid1);
		
	
	return ii_p;
	
}

