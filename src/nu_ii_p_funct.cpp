#include"main.h"


// ionized imprity scattering calculation //---- equation 3.29 of Ramu Thesis
void nu_ii_p_funct()
{
	double T;
	int T_loop;
	double v;
	double k_dum;
	double ii_p[limit2]={0};
	double A;
	double B;
	double iiB;
	double iiA;
	double epsilon_lf= epsilon_s[0];
	//double beta_constant =  0.140607315390570;
	
	double beta_constant;
	double k , de, k_step;
    	double integral = 0;

	    for(int counter = 0;counter <= points-2; counter++)
	    {
		k = k_grid[counter];
		de = energy_p[counter+1]-energy_p[counter];
		k_step = k_grid[counter+1]-k_grid[counter];

		if (free_e ==0)  // DOSCAR is used
		{
		   integral = integral+ de*((1e6*Ds_p[counter])/volume1)*f0(energy_p[counter],efef_p,T)*(1-f0(energy_p[counter],efef_p,T));
		    
		    // unit is (1/m)^3; 1e6 is multiplied to convert volume from (1/cm)^3 to (1/m)^3
		}
		else
		{
		    integral = integral+k_step*(k/pi)*(k/pi)*f0(energy_p[counter],efef_p,T)*(1-f0(energy_p[counter],efef_p,T))*1e27;
		    // unit is (1/m)^3
		    // Part of equation (70) or Rode's book (book8)
		}
		
	    }
    
		    double bet = (e*e/(epsilon_s[T_loop]*epsilon_0*(k_B*e)*T)*integral);   // unit (1/m)^2
		    											// 	
		    bet = pow(bet,0.5);   // unit (1/m) 
		    
		    beta_constant = bet*1e-9;
		    //converted from 1/m to 1/nm
		    cout<<"screening length is ="<<beta_constant<<endl;
		    
		for (int counter = 0;counter<points;counter++)
	    		{
			        k_dum = k_grid[counter];
			        v= v_p[counter];
				A= 0.5*(2 + (beta_constant*beta_constant/(k_dum*k_dum)));
				
				B= abs((A+1)/(A-1));
				
				//cout<<"A = "<<A<<endl;
				//cout<<"B = "<<B<<endl;
				
				iiA= (pow(e,4)* abs(N_ii))/(32*pi*k_dum*k_dum*epsilon_0*epsilon_lf*epsilon_lf*epsilon_0*v*pow(h_bar,2));
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
