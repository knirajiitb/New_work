#include"main.h"


// Acoustic scattering calculation  
void nu_de_p_funct(int T_loop)
{
	
	double T=T_array[T_loop];
	double v;
	double k_dum;
	
	//cout<< "temp  "<<T<< endl;
		
	for (int counter = 0;counter<points;counter++)
	    {
		k_dum = k_grid[counter];
		v= v_p[counter];
			        
	        nu_deformation_p[counter][0][0] = (k_B*T*pow(E_deformation,2)*k_dum*k_dum)/(2*pi*h_bar*h_bar*C_long*v)*1e10*1.60217657/1e8; // equation taken from ramu thesis equation 3.32a 
		//cout<<"nu_deformation_p[counter] =  "<<nu_deform[counter]<<endl;
		// for unit conversion  *1e10*1.60217657/1e8
	    }
			
	// saving results  scattreing rate and energy
	FILE *fid1;
	fid1= fopen("ADPscatt.dat","w");
	fprintf(fid1,"#energy               nu_de_p\n");
	for(int i=0;i<points;i++)
	{
		fprintf(fid1,"%e      %e\n", energy_p[i], nu_deformation_p[i][0][0]);
	}
	fclose(fid1);
			
		
	
}
