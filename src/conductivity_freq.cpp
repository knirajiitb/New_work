#include"main.h"

void conductivity_freq()
{	
	double t, dt;
	dt = 1/omega_s;
	//initial = 0;
	
	for(int j=0;j<len_freq;j++)
	{
		//cout<<"j = "<<j<<endl;			
		double tau, t, initial;
		tau = 75e-15;
		initial = -0.8e-12;

		for(int i=0;i<time_limit;i++)
		{
			//t = initial + i*dt + dt;
			t = initial + i*dt ;
					
			//cout<<"i = "<<i<<endl;			
			sigma_freqr[j] = sigma_freqr[j] + sigma_time[i]*cos(2*pi*freq[j]*t)*dt;
			sigma_freqi[j] = sigma_freqi[j] + sigma_time[i]*sin(2*pi*freq[j]*t)*dt;
			// unit 
			
			J_freqr[j] = J_freqr[j] + J_time[i]*cos(2*pi*freq[j]*t)*dt;
			J_freqi[j] = J_freqi[j] + J_time[i]*sin(2*pi*freq[j]*t)*dt;

			Efield_freqr[j] = Efield_freqr[j] + Efield_time[i]*cos(2*pi*freq[j]*t)*dt;
			Efield_freqi[j] = Efield_freqi[j] + Efield_time[i]*sin(2*pi*freq[j]*t)*dt;			
		}
		
		//sigma_freqr[j] = sigma_freqr[j]/(dt*time_limit);
		//sigma_freqi[j] = sigma_freqi[j]/(dt*time_limit);
	}
	
	
	
//-------------------------------------------- save sigma with freq --------------------------------------------------	
	FILE *fid1;
	fid1 = fopen("sigma_freq.dat","w");

	fprintf(fid1,"#index	 Frequency	Real conductivity  Imaginary conductivity\n" );
	for (int i = 0; i < len_freq; i++)
	{
		fprintf(fid1,"%d    %e		%e	%e \n", i+1, freq[i], sigma_freqr[i], sigma_freqi[i]);
	}
	fclose(fid1);
//----------------------------------------------------------------------------------------------------------------------
			    	      
//-------------------------------------------- save current density with freq --------------------------------------------------	

	fid1 = fopen("J_freq.dat","w");

	fprintf(fid1,"#index	 Frequency	Real Current Density  		Imaginary Current Density\n" );
	for (int i = 0; i < len_freq; i++)
	{
		fprintf(fid1,"%d    %e		%e	%e \n", i+1, freq[i], J_freqr[i], J_freqi[i]);
	}
	fclose(fid1);
//----------------------------------------------------------------------------------------------------------------------


//-------------------------------------------- save Efield with freq --------------------------------------------------	
	fid1 = fopen("Efield_freq.dat","w");

	fprintf(fid1,"#index	 Frequency	Real Efield     Imaginary Efield \n" );
	for (int i = 0; i < len_freq; i++)
	{
		fprintf(fid1,"%d    %e		%e	%e \n", i+1, freq[i], Efield_freqr[i], Efield_freqi[i]);
	}
	fclose(fid1);
//----------------------------------------------------------------------------------------------------------------------

}




