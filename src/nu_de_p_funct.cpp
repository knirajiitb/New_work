#include"main.h"


// Acoustic scattering calculation
void nu_de_p_funct()
{
	
	
	// from ramu thesis equation 3.32a 
	//double adp= (k_B*T*pow(E_deformation,2)*m_h*k)/(2*pi*pow(h_bar,3)*C_long);
	double T;
	double v;
	double k_dum;
	double nu_deform[limit2]={0};
	cout<<"trial,I'M getting compiled"<<endl;
	

	for (int counter = 0;counter<points;counter++)
	    		{
			        k_dum = k_grid[counter];
			        v= v_p[counter];					
	            		//nu_deform[counter] = nu_de_p_funct(k_dum,T);
	            		nu_deform[counter] = (k_B*T*pow(E_deformation,2)*k_dum*k_dum)/(2*pi*h_bar*h_bar*C_long*v)*1e10*1.60217657/1e8;
				//cout<<"nu_deformation_p[counter] =  "<<nu_deform[counter]<<endl;
			}
	
		FILE *fid1;
		fid1= fopen("ADPscatt.dat","w");
			for(int i=0;i<points;i++){
			fprintf(fid1,"%e     %e\n", energy_p[i], nu_deform[i]);
			}
		fclose(fid1);
			
		
	
}
