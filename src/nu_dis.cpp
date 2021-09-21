
#include"main.h"

void nu_dis(double T, double beta_constant, double epsilon_s)
{
    // According to equations in Ager's Mathematica code (sources are referenced there):

	double k, v;

	for (int counter = 0;counter<points;counter++)
	{
	        k = k_grid[counter];					
    		v = v_n[counter];

		nu_dislocation[counter] = (N_dis*e*e*e*e*k)/(h_bar*h_bar*epsilon_0*epsilon_0*epsilon_s*epsilon_s*c_lattice*c_lattice*v)
    *1/(pow(beta_constant,4)*pow((1+(4*k*k)/pow(beta_constant,2)),1.5))*2.43146974985767e42*1.60217657/1e8;
	
		//cout<<"nu_dislocation[counter] =  "<<nu_dislocation[counter]<<endl;

	}
    // The coefficient for conversion of unit to get 1/s:
}
