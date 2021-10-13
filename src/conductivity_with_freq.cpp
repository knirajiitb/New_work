#include"main.h"

void conductivity_with_freq(double T)
{	
	cout<<"conductivity with freq "<<endl;
	double Sa[points]={0}, Se[points]={0}, m1[points]={0}, m2[points]={0};

	double fir[points]={0}, fii[points]={0}, fir_old[points]={0}, fii_old[points]={0};	

	// Sa and Se calculated for in scattering rate	
	// If POP scattering is included
	if (scattering_mechanisms[1] == 1)
	{
		for (int i = 0;i < points;i++)
		{
                        Sa[i] = N_poph_atT *lambda_i_minus_grid[i]*(1-f_dist[minus_index_pop[0][i]])/(1-f_dist[i]) 						        *(kminus_grid_pop[0][i]/k_grid[i]);

                        Se[i] = (N_poph_atT+1) *lambda_i_plus_grid[i]*(1-f_dist[plus_index_pop[0][i]])/(1-f_dist[i]) 						        *(kplus_grid_pop[0][i]/k_grid[i]);
		}
	}

	int p,m;

	double de, integral_numerator1=0, integral_numerator2=0, integral_denominator=0;
	
	double fi[points]={0}, fi_old[points]={0};
	
	
//----------------------------------------------------------------------------------------------------	
	// Drude mobilility
	for (int i= 0;i < points;i++)
	{
		fi_old[i] = 1/(denom[i]) ;
	}

	for (int iteration = 1;iteration<iterations;iteration++)
	{            
		for (int i= 0;i < points;i++)
		{
			p = plus_index_pop[0][i];
			m = minus_index_pop[0][i];
			fi[i] = 1/(denom[i])*(1 + Sa[i]*fi_old[m] + Se[i]*fi_old[p]) ;
		}
		
		for (int i= 0;i < points;i++)
		{
			fi_old[i] = fi[i] ;
		}
	
	}
	integral_numerator1 = 0;
	integral_denominator = 0;
	
	for (int i = 0;i<=points-2;i++)
	{
		
	    de = (energy_n[i+1] - energy_n[i]);
	    integral_numerator1 = integral_numerator1 + de*(Ds_n[i]/volume1)*v_n[i]*v_n[i]*(-f0x1_f0[i]/(k_B*e*T))*fi[i];
		    // = int[fi(En)*DOS(En)*v(En)*v(En)*dEn*df0/de]

	    integral_denominator = integral_denominator + de*(Ds_n[i]/volume1)*f0(energy_n[i],efef_n,T);
		    // =int[f(En)*DOS(En)*dEn]
	}

	double mobility, tau1;
	
   	mobility = (-e/3.0)*integral_numerator1/integral_denominator;
	
	tau1 = mobility*(m*m_e)/e;
	
	//cout<<"mobility for drude model = "<<mobility<<endl;
	//cout<<"tau1 = "<<tau1<<endl;
	
	for(int j=0;j<len_freq;j++)
	{
		
	   	mobility_drude_freqr2[j] = mobility/(1+freq[j]*2*pi*freq[j]*2*pi*tau1*tau1);

	   	mobility_drude_freqi2[j] = mobility*tau1*2*pi*freq[j]/(1+freq[j]*2*pi*freq[j]*2*pi*tau1*tau1);

		sigma_drude_freqr2[j] = n_e*e*mobility_drude_freqr2[j];		
		sigma_drude_freqi2[j] = n_e*e*mobility_drude_freqi2[j];	
	}	
	
//----------------------------------------------------------------------------------------------------
	



//------------------------ loop for different frequency --------------------------------------	
	for(int j=0;j<len_freq;j++)
	{
		integral_numerator1 = 0;
		integral_numerator2 = 0;
		integral_denominator = 0;
				
		for (int i= 0;i < points;i++)
		{
			m1[i] = denom[i]/(denom[i]*denom[i] + freq[j]*freq[j]*4*pi*pi) ;
			m2[i] = (freq[j]*2*pi)/(denom[i]*denom[i] + freq[j]*freq[j]*4*pi*pi) ;
			
			fir_old[i] = m1[i];
			fii_old[i] = -m2[i];	
		}
		
		
		for (int iteration = 1;iteration<iterations;iteration++)
		{            
			for (int i=0;i<points;i++)
			{
				p = plus_index_pop[0][i];
				m = minus_index_pop[0][i];
				
				fir[i] = m1[i]*(1 + Sa[i]*fir_old[m] + Se[i]*fir_old[p]) + m2[i]*(Sa[i]*fii_old[m] + Se[i]*fii_old[p]);
				fii[i] = -m2[i]*(1 + Sa[i]*fir_old[m] + Se[i]*fir_old[p]) + m1[i]*(Sa[i]*fii_old[m] + Se[i]*fii_old[p]);
			}
			for (int i=0;i<points;i++)
			{
				fir_old[i] = fir[i];
				fii_old[i] = fii[i];				
			}				
		}
		
		
		for (int i = 0;i<=points-2;i++)
		{
		    de = (energy_n[i+1] - energy_n[i]);
		    integral_numerator1 = integral_numerator1 + de*(Ds_n[i]/volume1)*v_n[i]*v_n[i]*(-f0x1_f0[i]/(k_B*e*T))*fir[i];
			    // = int[fi(En)*DOS(En)*v(En)*v(En)*dEn*df0/de]

		    integral_numerator2 = integral_numerator2 + de*(Ds_n[i]/volume1)*v_n[i]*v_n[i]*(-f0x1_f0[i]/(k_B*e*T))*fii[i];
			    // = int[fi(En)*DOS(En)*v(En)*v(En)*dEn*df0/de]

		    integral_denominator = integral_denominator + de*(Ds_n[i]/volume1)*f0(energy_n[i],efef_n,T);
			    // =int[f(En)*DOS(En)*dEn]
		}

	   	mobility_freqr2[j] = (-e/3.0)*integral_numerator1/integral_denominator;

	   	mobility_freqi2[j] = (-e/3.0)*integral_numerator2/integral_denominator;

		sigma_freqr2[j] = n_e*e*mobility_freqr2[j];		
		sigma_freqi2[j] = n_e*e*mobility_freqi2[j];	
		/*
		cout<<"Frequency  =   "<<freq[j]<<endl;
		cout<<"Real Mobility = "<<mobility_freqr2[j]<<" cm^2/V-s "<<endl;
		cout<<"Imaginary Mobility = "<<mobility_freqi2[j]<<" cm^2/V-s "<<endl;
		cout<<"Real Mobility with drude model = "<<mobility_drude_freqr2[j]<<" cm^2/V-s "<<endl;
		cout<<"Imaginary Mobility with drude model = "<<mobility_drude_freqi2[j]<<" cm^2/V-s "<<endl;
		cout<<"Real Conductivity = "<<sigma_freqr2[j]<<endl;
		cout<<"Imaginary Conductivity = "<<sigma_freqi2[j]<<endl;			
		cout<<"Real Conductivity  with drude model = "<<sigma_drude_freqr2[j]<<endl;
		cout<<"Imaginary Conductivity  with drude model = "<<sigma_drude_freqi2[j]<<endl;	
		*/		
	}  // for loop for different frequency completed here	
//------------------------ loop for different frequency completed--------------------------------------


/*
// with Drude Model 	
//-------------------------------------------- save sigma with freq --------------------------------------------------	
	FILE *fid1;
	fid1 = fopen("sigma_freq_2.dat","w");

	fprintf(fid1,"#index	 Freq.	Real_Cond (S/cm)  Imaginary_Cond (S/cm)  Real_Cond_Drude (S/cm) Imaginary_Cond_Drude (S/cm)\n" );
	
	for (int i = 0; i < len_freq; i++)
	{
		fprintf(fid1,"%d    %e	%e	%e 	%e	%e \n", i+1, freq[i], sigma_freqr2[i], sigma_freqi2[i],
		 sigma_drude_freqr2[i], sigma_drude_freqi2[i]);
	}
	fclose(fid1);
//----------------------------------------------------------------------------------------------------------------------
			    	      
//-------------------------------------------- save mobility with freq --------------------------------------------------	

	fid1 = fopen("mobility_freq_2.dat","w");

	fprintf(fid1,"#index	 Freq.	 Real Mobility(cm^2/V-s)   Imaginary_Mobility    Real_Mobility_Drude    Imaginary_Mobility_Drude \n" );
	for (int i = 0; i < len_freq; i++)
	{
		fprintf(fid1,"%d    %e	%e	%e  	%e	%e \n", i+1, freq[i], mobility_freqr2[i], mobility_freqi2[i],
		mobility_drude_freqr2[i], mobility_drude_freqi2[i]);
	}
	fclose(fid1);
//----------------------------------------------------------------------------------------------------------------------
*/


// without Drude Model 
//-------------------------------------------- save sigma with freq --------------------------------------------------	
	FILE *fid1;
	fid1 = fopen("sigma_freq_2.dat","w");

	fprintf(fid1,"#index	 Freq.	Real_Cond (S/cm)  Imaginary_Cond (S/cm) \n" );
	
	for (int i = 0; i < len_freq; i++)
	{
		fprintf(fid1,"%d    %e	  %e	%e \n", i+1, freq[i], sigma_freqr2[i], sigma_freqi2[i]);
	}
	fclose(fid1);
//----------------------------------------------------------------------------------------------------------------------
			    	      
//-------------------------------------------- save mobility with freq --------------------------------------------------	

	fid1 = fopen("mobility_freq_2.dat","w");

	fprintf(fid1,"#index	 Freq.	 Real Mobility(cm^2/V-s)   Imaginary_Mobility\n" );
	for (int i = 0; i < len_freq; i++)
	{
		fprintf(fid1,"%d    %e	  %e	%e \n", i+1, freq[i], mobility_freqr2[i], mobility_freqi2[i]);
	}
	fclose(fid1);
//----------------------------------------------------------------------------------------------------------------------

}




