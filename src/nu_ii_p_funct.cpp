#include"main.h"


// ionized imprity scattering calculation
void nu_ii_p_funct()
{
	
	double v;
	double k_dum;
	double ii_p[limit2]={0};
	double A;
	double B;
	double iiB;
	double iiA;
	double epsilon_lf= epsilon_s[0];
	double beta_constant =  0.140607315390570; // used this from nu_iip.cpp for n-type
	//double beta_constant;
	// equation 3.29 of Ramu Thesis
	cout<<"screening length is ="<<beta_constant<<endl;
	//cout<<"ionized impurity = "<<N_ii<<endl;
		for (int counter = 0;counter<points;counter++)
	    		{
			        k_dum = k_grid[counter];
			        v= v_p[counter];
				A= 0.5*(2 + (beta_constant*beta_constant/(k_dum*k_dum)));
				
				B= abs((A+1)/(A-1));
				//cout<<"B = "<<B<<endl;
				//cout<<"A = "<<A<<endl;
				//cout<<"B = "<<B<<endl;
				
				iiA= (pow(e,4)* N_ii)/(32*pi*k_dum*k_dum*epsilon_0*epsilon_lf*epsilon_lf*epsilon_0*v*pow(h_bar,2));
				iiB = ((3*A -1)*(3*A -1)*log(B) - 18*A +12 -8/(A+1));
				
				ii_p[counter]= iiA*iiB ;
				
				
				
				}
				
				
				
				FILE *fid1;
				fid1= fopen("iipScatt.dat","w");
				for(int i=0;i<points;i++){
					fprintf(fid1,"%e     %e\n", energy_p[i], ii_p[i]);
				}
				
								
				fclose(fid1);
		
				
				
	
}
