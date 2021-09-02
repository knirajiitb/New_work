#include"main.h"


// Acoustic scattering calculation  
void nu_de_p_funct(int T_loop)
{
	
	double T=T_array[T_loop];
	double v;
	double k_dum;
	
	//cout<< "temp  "<<T<< endl;
	//cout<<"C_long = "<<C_long<<" dyme/cm^2"<<endl;	
	for (int counter = 0;counter<points;counter++)
	    {
		k_dum = k_grid[counter]*1e9;
		v = v_p[counter]/100;   // unit m-s
			        
	        nu_deformation_p[counter][0][0] = ((k_B*e)*T*pow(E_deformation,2)*(k_dum*k_dum))/(2*pi*h_bar*h_bar*(C_long/10.0)*v); 
	        // equation taken from ramu thesis equation 3.32a 
		
		/*
		cout<<"counter = "<<counter<<endl;
		cout<<"v = "<<v<<endl;
		cout<<"nu_deformation_p[counter][0][0] =  "<<nu_deformation_p[counter][0][0]<<endl;
		getchar();
		*/
	    }
	
	/*		
	// saving results  scattreing rate and energy
	FILE *fid1;
	fid1= fopen("de.dat","w");
	fprintf(fid1,"#energy               nu_de_p\n");
	for(int i=0;i<points;i++)
	{
		fprintf(fid1,"%e      %e\n", k_grid[i], nu_deformation_p[i][0][0]);
	}
	fclose(fid1);
	//*/			
	
}
