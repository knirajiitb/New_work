
#include"main.h"

void nu_de(double T)
//deformation potential scattering rate according to equation (112) or Rode's book
{

	double k, v, nu[limit2][de_number]={0};

	for (int counter = 0;counter<points;counter++)
	{
	        k = k_grid[counter];					
    		v = v_n[counter];
    		nu[counter][0] = (k_B*T*pow(E_deformation[0],2)*k*k)/(3*pi*h_bar*h_bar*v*C_long)*(3-8*pow(c_n[counter],2)
+6*pow(c_n[counter],4))*1e10*1.60217657/1e8;
		// From equation (112) of Rode's book (book8):
		
		// 1e10 is coming from unit conversion (take a look at OneNote notes in Deformation potential //
		// section) and *1.60217657/1e8 is to get from cm/s (v) to hk/(md)


		if(de_number==2)		
		{
		nu[counter][1] = (k_B*T*pow(E_deformation[1],2)*k*k)/(3*pi*h_bar*h_bar*v*C_trans)*(3-8*pow(c_n[counter],2) 
		+ 6*pow(c_n[counter],4))*1e10*1.60217657/1e8;
		}
		
		if(de_number==3)		
		{

		nu[counter][1] = (k_B*T*pow(E_deformation[1],2)*k*k)/(3*pi*h_bar*h_bar*v*C_trans)*(3-8*pow(c_n[counter],2) 
		+ 6*pow(c_n[counter],4))*1e10*1.60217657/1e8;

		nu[counter][2] = (k_B*T*pow(E_deformation[2],2)*k*k)/(3*pi*h_bar*h_bar*v*C_za)*(3-8*pow(c_n[counter],2) 
		+ 6*pow(c_n[counter],4))*1e10*1.60217657/1e8;
		}
		
		nu_deformation[counter] = nu[counter][0] + nu[counter][1] + nu[counter][2];
		//cout<<"nu_deformation[counter] =  "<<nu_deformation[counter]<<endl;
	}
			
	
	FILE *fid1;
	fid1 = fopen("acoustic_scattering_rate.dat","w");
	
	if(de_number==1)
		fprintf(fid1,"# energy           nu_LA \n");
	else if(de_number==2)
		fprintf(fid1,"# energy           nu_LA 	nu_TA	\n");
	else
		fprintf(fid1,"# energy           nu_LA 	nu_TA		nu_ZA \n");
	
	if(de_number==1)
	{
		for (int i = 0; i < points; i++)		
			fprintf(fid1,"  %e    	%e   \n", energy_n[i], nu[i][0] );
	}

	if(de_number==2)
	{
		for (int i = 0; i < points; i++)		
			fprintf(fid1,"  %e    	%e    %e \n", energy_n[i], nu[i][0], nu[i][1] );
	}

	if(de_number==3)
	{
		for (int i = 0; i < points; i++)		
			fprintf(fid1,"  %e    	%e   	%e	%e \n", energy_n[i], nu[i][0], nu[i][1], nu[i][2] );
	}
	fclose(fid1);
	
}


