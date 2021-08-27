#include"main.h"
// NPOP scattering calculation
void nu_npop_p_funct(double T)			
{
	double k_pp;
	double k_mm;
	double k_dum; 
	double arr[points];
	double v;
	double npop1,npop2,npop3,npop4;
	double N_op[limit5];
	double c_bar = (1.0/3.0)*C_long + (2.0/3)*C_trans; // unit  dyne/cm^2
	double E_npop[limit5];
	double omega_npop[limit5];
	double nu_npop_p_all[limit2][limit5]={0};	// scattering rate at each frequency
	double nu_npop_p_total[limit2]={0};        // sum of all scattering rate ... it is also saved as nu_npop_p[limit2][0][0]...
	
			
	cout<< "getting compiled niceley"<<endl;	
	 	
	
	double kplus_grid_npop[limit2], kminus_grid_npop[limit2];
	int plus_index_npop[limit2], minus_index_npop[limit2];
	int minus_index, plus_index;
	
	for (int counter = 0;counter < points;counter++)
	{	
		for(int loop=0;loop<npop_number;loop++)
		{
			E_npop[loop] = De_npop[loop]; 
			//cout<<"%e" <<E_npop <<endl;
			omega_npop[loop] = we_npop[loop]; 
			
			kplus_grid_npop[counter] = kplus(counter,omega_npop[loop],points,energy_p);
			kminus_grid_npop[counter] = kminus(counter,omega_npop[loop],points,energy_p);

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - kminus_grid_npop[counter]);
			minus_index =FindMinInd(arr,points);

			for (int i=0;i<points;i++)
			arr[i] = abs(k_grid[i] - kplus_grid_npop[counter]);
			plus_index =FindMinInd(arr,points);

			plus_index_npop[counter] = plus_index;
			minus_index_npop[counter] = minus_index;
		}
	}
	
	int len = sizeof(nu_npop_p_total)/sizeof(nu_npop_p_total[0]);	
	for (int counter = 0;counter<len;counter++)
		nu_npop_p_total[counter] = 0;
	
	for(int counter=0; counter<points;counter++)
	{	
		k_dum=k_grid[counter];   // unit is 1/nm
		v=v_p[counter]*1e-2;    // conveterd from unit  cm/s to m/s
		k_mm= kminus_grid_npop[counter];  // unit is 1/nm
		k_pp= kplus_grid_npop[counter];   // unit is 1/nm
		
		for(int loop=0;loop<npop_number;loop++)
		{	
			E_npop[loop] = De_npop[loop]; 
			//cout<<"%e" <<E_npop <<endl;
			omega_npop[loop] = we_npop[loop]; 
			
			N_op[loop] = N_poph(omega_npop[loop],T);
			
			npop1= E_npop[loop] *E_npop[loop] *omega_npop[loop]*k_dum*k_pp/(2*pi*h_bar*c_bar*v);
			cout<<"npop1 is = "<<npop1<<endl;
			npop2=  E_npop[loop]*E_npop[loop] *omega_npop[loop] *k_dum*k_mm/(2*pi*h_bar*c_bar*v);
			npop3= (N_op[loop]*(1-f0(energy_p[plus_index_npop[counter]],efef_p,T)) + (N_op[loop]+1)*f0(energy_p[plus_index_npop[counter]],efef_p,T));
			npop4=	(N_op[loop]*f0(energy_p[plus_index_npop[counter]],efef_p,T) + (N_op[loop]+1)*(1-f0(energy_p[plus_index_npop[counter]],efef_p,T)));
			
			if (energy_p[counter] < h_bar*omega_npop[loop])
			    {
				    nu_npop_p_all[counter][loop] = 0;
			    }
			    else
			    {    
				nu_npop_p_all[counter][loop] = npop1*npop3 +npop2*npop4;
			    }
			    
			    
		 	nu_npop_p_total[counter] = nu_npop_p_total[counter] + nu_npop_p_all[counter][loop];  
		 	nu_npop_p[counter][0][0]= nu_npop_p_total[counter];
		 }
			
				
	 }	
		  
	
	
	
	
	
	FILE *fid1;
	fid1 = fopen("npopscatt.dat","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%e          %e  \n", energy_p[i], 	nu_npop_p[i][0][0] );
	}
	fclose(fid1);
	
	
			
}
