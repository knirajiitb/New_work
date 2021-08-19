#include"main.h"


// POP scattering calculation
void nu_So_p_funct(double T, int T_loop,double omega_LO)
{	
	double k_pp;
	double k_mm;
	
	double C_minus;
	double C_plus;
	double B_minus;
	double B_plus;
	
	
	double k_dum; 
	double N_op= N_poph(omega_LO,T);  // number of optical phonon at temperature T
	double v;
	
	double arr[points];
	double pop1,pop2,pop3;

	double kplus_grid_pop_p[limit2], kminus_grid_pop_p[limit2];
	int plus_index_pop_p[limit2], minus_index_pop_p[limit2];
	int minus_index, plus_index;
	
		
	for (int counter = 0;counter < points;counter++)
	{

		kplus_grid_pop_p[counter] = kplus(counter,omega_LO,points,energy_p);
		kminus_grid_pop_p[counter] = kminus(counter,omega_LO,points,energy_p);

		for (int i=0;i<points;i++)
		arr[i] = abs(k_grid[i] - kminus_grid_pop_p[counter]);
		minus_index =FindMinInd(arr,points);

		for (int i=0;i<points;i++)
		arr[i] = abs(k_grid[i] - kplus_grid_pop_p[counter]);
		plus_index =FindMinInd(arr,points);

		plus_index_pop_p[counter] = plus_index;
		minus_index_pop_p[counter] = minus_index;
			    
			
	}
		
	for(int counter=0;counter<points;counter++)
	{	
		k_dum=k_grid[counter];
		v=v_p[counter];
		k_mm= kminus_grid_pop_p[counter];
		k_pp= kplus_grid_pop_p[counter];
		C_minus= (k_dum*k_dum+k_mm*k_mm)/(2*k_mm*k_dum);	// c_minus from equation 3.22 of ramu thesis
 		C_plus= (k_dum*k_dum+k_pp*k_pp)/(2*k_pp*k_dum);	// c_minus from equation 3.22 of ramu thesis	
	 	 		
		B_minus= (((1+3*C_minus*C_minus)/2)*log(abs((1+C_minus)/(1-C_minus))) - 3*C_minus);		//from equation 3.22 of ramu thesis
		B_plus= (((1+3*C_plus*C_plus)/2)*log(abs((1+C_plus)/(1-C_plus))) - 3*C_plus);		//from equation 3.22 of ramu thesis
		 
		pop1= (e*e*omega_LO/(16*pi*h_bar*v*epsilon_0))*(1/epsilon_inf[T_loop] - 1/epsilon_s[T_loop]);
		pop1= pop1/e; // for unit conversion of h_bar
		pop2= B_plus*(N_op*(1-f0(energy_p[plus_index_pop_p[counter]],efef_p,T)) + (N_op+1)*f0(energy_p[plus_index_pop_p[counter]],efef_p,T));
		pop3= B_minus*(N_op*f0(energy_p[minus_index_pop_p[counter]],efef_p,T) + (N_op+1)*(1-f0(energy_p[minus_index_pop_p[counter]],efef_p,T)));
		
		nu_So_p[counter][0][0]= pop1*(pop2 + pop2);
		
	}
	
	FILE *fid1;
	fid1 = fopen("nu_So_pop.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e               %e  \n", energy_p[i], nu_So_p[i][0][0]);
	}
	fclose(fid1);
		
	
}
