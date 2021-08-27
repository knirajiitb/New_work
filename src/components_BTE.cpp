#include"main.h"

void components_BTE(double T, int T_loop, double efefn, double efefp, int ii)
{

//---------------------------- components for BTE ------------------------------------------------------------------
	if(type == "n")
	{	
                   double k_dum;
		    beta_constant = beta(T, T_loop);
		    // unit 1/nm

		    cout<< "Inverse screening length, beta = "<<beta_constant<<" (1/nm)"<<endl;
		    double integral_numerator_n = 0;
		    double integral_denominator_n = 0;

		    //--------------------------- df0dz_integral_n -----------------------------------------
		    int factr = 10;
		    for(int counter = 0;counter<=points-2;counter++)
		    {
		        double dk = (k_grid[counter+1]-k_grid[counter])/factr;
		        for (int ss = 0;ss<=factr-1;ss++)
		        {
		            integral_numerator_n = integral_numerator_n + dk*(pow(((k_grid[counter]+ss*dk)/pi),2))
		                                *f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T))*energy_n[counter]/(k_B*T);
		            // Part of equation (54) of Rode's book
		            integral_denominator_n = integral_denominator_n + dk*pow(((k_grid[counter]+ss*dk)/pi),2)
		            *f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T));
		            // Part of equation (54) of Rode's book
		        }
		    }

		    df0dz_integral_n = integral_numerator_n/integral_denominator_n;
		    //cout<<"df0dz_integral_n = "<<df0dz_integral_n<<endl;
		    //--------------------------- df0dz_integral_n calculated -----------------------------------------
		    		    
//--------------------------------------common terms calculated -------------------------------------------------------

		for (int counter = 0;counter<points;counter++)
		{
			k_dum = k_grid[counter];
			// unit 1/nm
			//cout<<"counter+1 = "<<counter+1<<endl;
			//cout<<"k_dum = "<<k_dum<<endl;

		        df0dk_grid[counter] = df0dk(k_dum, T, efefn, coefficients_cond, kindex_cond, a11);
			 // unit (nm)

		        f_dist[counter] = f0(energy_n[counter],efefn,T);
		        //cout<<"In between "<<endl;
		        //cout<<"energy_n[counter]  = "<<energy_n[counter]<<endl;
		        //cout<<"efefn = "<<efefn<<endl;
		        //cout<<"T = "<<T<<endl;

		        thermal_driving_force[counter] = -1*v_n[counter]*df0dz(k_dum, efefn, T, df0dz_integral_n,coefficients_cond, kindex_cond, a11);

		        f0x1_f0[counter] = f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T));

		        electric_driving_force[counter] = -(1*E/h_bar)*df0dk_grid[counter]*1e-7;
			// unit is 1/s , hbar unit is eV-s so e i s not multiplied in numerator

		        //cout<<"df0dk_grid[counter] =  "<<df0dk_grid[counter]<<endl;
		        //cout<<"f_dist[counter]  =  "<<f_dist[counter]<<endl;
		        //cout<<"thermal_driving_force[counter] =  "<<thermal_driving_force[counter]<<endl;
		        //cout<<"f0x1_f0[counter] =  "<<f0x1_f0[counter]<<endl;
		        //cout<<"electric_driving_force[counter]  = "<<electric_driving_force[counter]<<endl;
		}

//--------------------------------------common terms calculated completed  --------------------------------------------

//----------------------POP scattering rate --------------------------------------------------------------------
			    
		 // polar optical phonon scattering 
		if (scattering_mechanisms[1]==1)
		{
		    N_poph_atT = N_poph(omega_LO,T);
		    //cout<<"N_poph_atT = "<<N_poph_atT<<endl;

			int minus_index, plus_index;
			double arr[points];

			for (int counter1 = 0;counter1 < points;counter1++)
			{

				kplus_grid_pop[counter1] = kplus(counter1,omega_LO,points,energy_n);
				kminus_grid_pop[counter1] = kminus(counter1,omega_LO,points,energy_n);

				for (int i=0;i<points;i++)
				arr[i] = abs(k_grid[i] - kminus_grid_pop[counter1]);
				minus_index =FindMinInd(arr,points);

				for (int i=0;i<points;i++)
				arr[i] = abs(k_grid[i] - kplus_grid_pop[counter1]);
				plus_index =FindMinInd(arr,points);

				plus_index_pop[counter1] = plus_index;
				minus_index_pop[counter1] = minus_index;
			    
			    //cout<<"counter1 = "<<counter1<<endl;
			    //cout<<"plus_index_pop[counter1] = "<<plus_index_pop[counter1]<<endl;
			    //cout<<"minus_index_pop[counter1] = "<<minus_index_pop[counter1]<<endl;	
			    //cout<<"kplus_grid_pop[counter1] = "<<kplus_grid_pop[counter1]<<endl;
			    //cout<<"kminus_grid_pop[counter1] = "<<kminus_grid_pop[counter1]<<endl;
			    //getchar();	
			}
		}
								
		// POP scattering
		if (scattering_mechanisms[1]==1)
		{
			for (int counter = 0;counter < points;counter++)
			{
			    Aminus_grid[counter] = Aminus(counter,omega_LO,points,energy_n);
			    Aplus_grid[counter] = Aplus(counter,omega_LO,points,energy_n);


			    betaplus_grid[counter] = betaplus(counter,omega_LO,epsilon_s[T_loop],epsilon_inf[T_loop],points);
			    betaminus_grid[counter] = betaminus(counter,omega_LO,epsilon_s[T_loop],epsilon_inf[T_loop],points);

			    lambda_i_plus_grid[counter] = abs(lambda_i_plus(counter,omega_LO,Aplus_grid[counter], epsilon_s[T_loop],epsilon_inf[T_loop], points));

			    lambda_i_minus_grid[counter] = abs(lambda_i_minus(counter,omega_LO,Aminus_grid[counter],                epsilon_s[T_loop],epsilon_inf[T_loop],points));

			    lambda_o_plus_grid[counter] = abs(lambda_o_plus(counter,omega_LO,Aplus_grid[counter],
			    epsilon_s[T_loop],epsilon_inf[T_loop],points));

			    lambda_o_minus_grid[counter] = abs(lambda_o_minus(counter, omega_LO,Aminus_grid[counter],
			    epsilon_s[T_loop],epsilon_inf[T_loop],points));

				//---------------------------- code to debug -------------------------------------------------------------
				//cout<<"counter = "<<counter<<endl;
				//cout<<"Aplus_grid[counter] =   "<<Aplus_grid[counter]<<endl;
				//cout<<"Aminus_grid[counter] =   "<<Aminus_grid[counter]<<endl;

				//cout<<"betaplus_grid[counter] =  "<<betaplus_grid[counter]<<endl;
				//cout<<"betaminus_grid[counter] =  "<<betaminus_grid[counter]<<endl;

				//cout<<"lambda_i_plus_grid[counter] =  "<<lambda_i_plus_grid[counter]<<endl;
				//cout<<"lambda_i_minus_grid[counter] =  "<<lambda_i_minus_grid[counter]<<endl;
				//cout<<"lambda_o_plus_grid[counter] =  "<<lambda_o_plus_grid[counter]<<endl;
				//cout<<"lambda_o_minus_grid[counter] =  "<<lambda_o_minus_grid[counter]<<endl;

				//getchar();
			}
		}

		// pop scattering
		if(scattering_mechanisms[1]==1)
			pop_So(T, efefn, ii);
//-------------------------------------------------------- pop_So calculated---------------------------------------

//----------------------POP scattering rate completed--------------------------------------------------------------------


//----------------------Ionized Impurity scattering rate --------------------------------------------------------------------
			        			
		// ionized impourity scattering
	        if (scattering_mechanisms[0]==1)
	        {
			for (int counter = 0;counter < points;counter++)
			{
				k_dum = k_grid[counter];
			    //cout<<"k_dum = "<<k_dum<<endl;
			    //cout<<"beta_constant =  "<<beta_constant<<endl;
			    //cout<<"c_n[counter] =  "<<c_n[counter]<<endl;
			    //cout<<"counter =  "<<counter+1<<endl;

			    B_ii = (4*(k_dum*k_dum)/(beta_constant*beta_constant))/(1+4*k_dum*k_dum/(beta_constant*beta_constant))
			    +8*(beta_constant*beta_constant+2*k_dum*k_dum)/(beta_constant*beta_constant+4*k_dum*k_dum)*(pow(c_n[counter],2))
			    +(3*pow(beta_constant,4)+
			      6*pow(beta_constant,2)*k_dum*k_dum-8*k_dum*k_dum*k_dum*k_dum)/((beta_constant*beta_constant+4*k_dum*k_dum)*k_dum*k_dum)*pow(c_n[counter],4);
			    // According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)


			    D_ii = 1+(2*pow(beta_constant,2)*(pow(c_n[counter],2)/pow(k_dum,2))+
			              (3*pow(beta_constant,4)*(pow(c_n[counter],4))/(4*pow(k_dum,4))));
			              // According to equation (92) of Semiconductors and Semimetals, volume 10 (Rode's chapter)

			    nu_ionizedimpurity[counter] = nu_ii(k_dum,counter,beta_constant,v_n[counter],epsilon_s[T_loop]);
			     // unit 1/second	
			    //cout<<"B_ii = "<<B_ii<<endl;
			    //cout<<"D_ii = "<<D_ii<<endl;
			    //cout<<"N_ii = "<<N_ii<<endl;
			    //cout<<"nu_ionizedimpurity[counter] = "<<nu_ionizedimpurity[counter]<<endl;
			}
	        }
//----------------------II scattering rate completed--------------------------------------------------------------------

//----------------------Acoustic scattering rate --------------------------------------------------------------------
		        
		// acoustic scattering 
	        if (scattering_mechanisms[3]==1)
	        {
	        	for (int counter = 0;counter<points;counter++)
	    		{
			        k_dum = k_grid[counter];					
	            		nu_deformation[counter] = nu_de(k_dum,counter,T,v_n[counter]);
				//cout<<"nu_deformation[counter] =  "<<nu_deformation[counter]<<endl;
			}
		}	
//----------------------Acoustic scattering rate completed--------------------------------------------------------------------

//----------------------PZ scattering rate --------------------------------------------------------------------
		        
		// Piezoelectric scattering
	        if (scattering_mechanisms[4]==1)
	        {
	        	for (int counter = 0;counter<points;counter++)
	    		{
			        k_dum = k_grid[counter];					
	         		nu_piezoelectric[counter] = nu_pe(k_dum,counter,T,P_piezo[T_loop],epsilon_s[T_loop],v_n[counter]);
				//cout<<"nu_piezoelectric[counter] =  "<<nu_piezoelectric[counter]<<endl;
			}
		}
//----------------------PZ scattering rate --------------------------------------------------------------------
			
//----------------------Disloaction scattering rate -------------------------------------------------------			
		// Dislocation scattering
	        if (scattering_mechanisms[6]==1)
	        {
	        	for (int counter = 0;counter<points;counter++)
	    		{
			        k_dum = k_grid[counter];					
			        nu_dislocation[counter] = nu_dis(k_dum,counter,T,beta_constant,epsilon_s[T_loop], v_n[counter]);
				//cout<<"nu_dislocation[counter] =  "<<nu_dislocation[counter]<<endl;
			}
		}
//----------------------Disloaction scattering rate completed ------------------------------------------			

//----------------------Alloy scattering rate ----------------------------------------------------			
			
		// Alloy scattering
	        if (scattering_mechanisms[7]==1)
	        {
	        	for (int counter = 0;counter<points;counter++)
	    		{
			        k_dum = k_grid[counter];					
	     	                nu_alloy[counter] = nu_alloy1(k_dum, v_n[counter]);
				//cout<<"nu_alloy[counter] =  "<<nu_alloy[counter]<<endl;
			}
		}
//----------------------Alloy scattering rate completed ------------------------------------------			

			
//---------------------- Neutral impurity scattering rate  ---- ------------------------------------------			

		// neutral impurity
	        if (scattering_mechanisms[9]==1)  
		{
	        	for (int counter = 0;counter<points;counter++)
	    		{
			        k_dum = k_grid[counter];					
						
				nu_neutralimpurity[counter] = nu_im(k_dum,counter,epsilon_s[T_loop],N_im[ii],v_n[counter]);
				//cout<<"nu_neutralimpurity[counter] =  "<<nu_neutralimpurity[counter]<<endl;
			}
		}
			
//---------------------- Neutral impurity scattering rate  completed ------------------------------------------			


//---------------------- Intervalley scattering rate  ---- ------------------------------------------			

		// intravalley scattering	
		if (scattering_mechanisms[8]==1)   
		{
			for (int aa=0;aa<iv_number;aa++)
			{
				N_e[aa] = N_poph(we[aa],T);
				//cout<<" N_e[aa] = "<<N_e[aa]<<endl;
			}

	            	//cout<<"For intervalley scattering   =   "<<endl;	
	        	for (int counter = 0;counter<points;counter++)
	    		{
			    k_dum = k_grid[counter];					
			    //cout<<"counter   =   "<<counter<<endl;
			    for (int aa = 0;aa<iv_number;aa++)
			    {
			        lambda_e_plus_grid[counter][aa] = abs(lambda_e_plus(counter,we[aa],rho,De[aa],nfv[aa],points));
			        lambda_e_minus_grid[counter][aa] = abs(lambda_e_minus(counter,we[aa],rho,De[aa],nfv[aa],points));
				//cout<<"counter = "<<counter<<endl;
			        //cout<<"aa   =   "<<aa<<endl;
			        //cout<<"lambda_e_plus_grid   =   "<<lambda_e_plus_grid[counter][aa]<<endl;
			        //cout<<"lambda_e_minus_grid   =   "<<lambda_e_minus_grid[counter][aa]<<endl;
			        
			        // Equation number 129 of rode book
			    }
			    //getchar();
			}
	        }

	 	// intravalley scattering
		if (scattering_mechanisms[8]==1)  
		{
			int len = sizeof(nu_iv_total)/sizeof(nu_iv_total[0]);	
			for (int counter = 0;counter<len;counter++)
				nu_iv_total[counter] = 0;
			
			//cout<<endl<<"Inside intravalley scattering"<<endl;
			//cout<<" iv_number = "<<iv_number<<endl;
		    for (int counter = 0;counter<points;counter++)
		    {
		    	//cout<<"counter = "<<counter<<endl;
			for (int aa = 0;aa<iv_number;aa++)
			{
			    //cout<<"aa = "<<aa<<endl;
			    double k_minus = kminus(counter,we[aa],points,energy_n);
			    //cout<<"k_minus = "<<k_minus<<endl;

			    double arr[points];
			    for (int i=0;i<points;i++)
				arr[i] = abs(k_grid[i] - k_minus);
			    int minus_index =FindMinInd(arr,points);


			    double k_plus = kplus(counter,we[aa],points,energy_n);
			    //cout<<"k_plus = "<<k_plus<<endl;

			    for (int i=0;i<points;i++)
				arr[i] = abs(k_grid[i] - k_plus);
			    int plus_index =FindMinInd(arr,points);

			    //cout<<"plus_index = "<<plus_index<<endl;

			    double f_negative = f0(energy_n[minus_index],efefn,T);
			    double f_positive =  f0(energy_n[plus_index],efefn,T);

			    //cout<<"f_negative = "<<f_negative<<endl;
			    //cout<<"f_positive = "<<f_positive<<endl;

			    if (energy_n[counter] < h_bar*we[aa])
			    {
				    nu_iv[counter][aa] = (N_e[aa] + f_positive) *lambda_e_plus_grid[plus_index][aa];
			    }
			    else
			    {    
				nu_iv[counter][aa] = (N_e[aa] + 1 - f_negative) * lambda_e_minus_grid[minus_index][aa] + 
							(N_e[aa] + f_positive)*lambda_e_plus_grid[plus_index][aa];
			    }
			    
			    
			    nu_iv_total[counter] = nu_iv_total[counter] + nu_iv[counter][aa];            
			    
			    //cout<<"nu_iv[counter][aa] = "<<nu_iv[counter][aa]<<endl;
			    //cout<<"nu_iv_total[counter] = "<<nu_iv_total[counter]<<endl;
			    //getchar(); 	
			     	
			}  
		    }
		}
//------------------interavalley scattering rate completed ----------------------------
		
//---------------------- NPOP scattering rate  ---- ------------------------------------------			

		// npop scattering
	        if (scattering_mechanisms[2]==1)  
	        {
			for (int aa=0;aa<npop_number;aa++)
			{
				N_npop[aa] = N_poph(we_npop[aa],T);
				//cout<<" N_npop[aa] = "<<N_npop[aa]<<endl;
			}   
			//getchar();	

			//cout<<"For npop scattering   =   "<<endl;	
			for (int counter = 0;counter<points;counter++)
	    		{
			    k_dum = k_grid[counter];					
	            	
			    for (int aa = 0;aa<npop_number;aa++)
			    {
			        lambda_e_plus_grid_npop[counter][aa] = abs(lambda_e_plus(counter,we_npop[aa],rho,De_npop[aa],1,points));
			        lambda_e_minus_grid_npop[counter][aa] = abs(lambda_e_minus(counter,we_npop[aa],rho,De_npop[aa],1,points));
			        // Equation number 129 of rode book  for npop nfv = 1 always
			        //cout<<"aa   =   "<<aa<<endl;
			        //cout<<"lambda_e_plus_grid_npop[counter][aa]   =   "<<lambda_e_plus_grid_npop[counter][aa]<<endl;
			        //cout<<"lambda_e_minus_grid_npop[counter][aa]   =   "<<lambda_e_minus_grid_npop[counter][aa]<<endl; 
			    }
			    //getchar();
			}
	        }

	 	// npop scattering
		if (scattering_mechanisms[2]==1)  
		{
			int len = sizeof(nu_npop_total)/sizeof(nu_npop_total[0]);	
			for (int counter = 0;counter<len;counter++)
				nu_npop_total[counter] = 0;

			//cout<<endl<<"Inside npop scattering"<<endl;
			//cout<<" npop_number = "<<npop_number<<endl;
		    for (int counter = 0;counter<points;counter++)
		    {
		    	//cout<<"counter = "<<counter<<endl;
			for (int aa = 0;aa<npop_number;aa++)
			{
			    //cout<<"aa = "<<aa<<endl;
			    double k_minus = kminus(counter,we_npop[aa],points,energy_n);
			    //cout<<"k_minus = "<<k_minus<<endl;

			    double arr[points];
			    for (int i=0;i<points;i++)
				arr[i] = abs(k_grid[i] - k_minus);
			    int minus_index =FindMinInd(arr,points);
			    //cout<<"minus_index = "<<minus_index<<endl;	

			    double k_plus = kplus(counter,we_npop[aa],points,energy_n);
			    //cout<<"k_plus = "<<k_plus<<endl;

			    for (int i=0;i<points;i++)
				arr[i] = abs(k_grid[i] - k_plus);
			    int plus_index =FindMinInd(arr,points);

			    //cout<<"plus_index = "<<plus_index<<endl;

			    double f_negative = f0(energy_n[minus_index],efefn,T);
			    double f_positive =  f0(energy_n[plus_index],efefn,T);

			    //cout<<"f_negative = "<<f_negative<<endl;
			    //cout<<"f_positive = "<<f_positive<<endl;

			    if (energy_n[counter] < h_bar*we_npop[aa])
			    {
				    nu_npop[counter][aa] = (N_npop[aa] + f_positive) *lambda_e_plus_grid_npop[plus_index][aa];
			    }
			    else
			    {    
				nu_npop[counter][aa] = (N_npop[aa] + 1 - f_negative) * lambda_e_minus_grid_npop[minus_index][aa] + 
							(N_npop[aa] + f_positive)*lambda_e_plus_grid_npop[plus_index][aa];
			    }
			    
			    
			    nu_npop_total[counter] = nu_npop_total[counter] + nu_npop[counter][aa];            
			    
			    //cout<<"nu_npop[counter][aa] = "<<nu_npop[counter][aa]<<endl;
			    //cout<<"nu_npop_total[counter] = "<<nu_npop_total[counter]<<endl;
			    //getchar(); 	
			}
		    }
		}
		
		
//--------------------------------------npop completed -------------------------------------------------------

		// total scattering rates calculated here
		for(int counter=0;counter<points;counter++)
		{
		        nu_el[counter] = nu_ionizedimpurity[counter]*scattering_mechanisms[0] + 
		        	nu_npop_total[counter] * scattering_mechanisms[2] + 
		        	nu_deformation[counter]*scattering_mechanisms[3] + 
		               nu_piezoelectric[counter]*scattering_mechanisms[4] + 
		         	nu_dislocation[counter] * scattering_mechanisms[6] + 
		         	nu_alloy[counter] * scattering_mechanisms[7] + 
		         	nu_iv_total[counter] * scattering_mechanisms[8] +
		         	nu_neutralimpurity[counter]*scattering_mechanisms[9];

			denom[counter] = (S_o_grid_total[counter]*scattering_mechanisms[1] + nu_el[counter]);	
			//cout<<"nu_el[counter] = "<<nu_el[counter]<<endl;
		}			

	 	
//-----------------------------------------------------------------------------------------------------------------

				
	//------------------------------------ components for BTE END --------------------------------------------------------------
	// ----------------------------saving data ---------------------------------------------------
		    /*
		    FILE *fid1;
		    fid1 = fopen("Aplus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, Aplus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("Aminus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, Aminus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("betaplus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, betaplus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("betaminus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, betaminus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("lambda_i_plus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, lambda_i_plus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("lambda_i_minus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, lambda_i_minus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("lambda_o_plus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, lambda_o_plus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("lambda_o_minus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, lambda_o_minus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("lambda_e_plus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, lambda_e_plus_grid[i]);
		fclose(fid1);

		    fid1 = fopen("lambda_e_minus.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, lambda_e_minus_grid[i]);
		fclose(fid1);
		    */
		    
		    /* 			
		    fid1 = fopen("nu_deformation.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_deformation[i]);
		fclose(fid1);

		    fid1 = fopen("nu_piezoelectric.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_piezoelectric[i]);
		fclose(fid1);

		    fid1 = fopen("nu_ionizedimpurity.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_ionizedimpurity[i]);
		fclose(fid1);
		    
		    fid1 = fopen("nu_dislocation.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_dislocation[i]);
		fclose(fid1);    

		    fid1 = fopen("nu_alloy.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_alloy[i]);
		fclose(fid1);

		    fid1 = fopen("nu_neutralimpurity.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_neutralimpurity[i]);
		fclose(fid1);
		    
		    //FILE *fid1;
		    fid1 = fopen("nu_el.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, nu_el[i]);
		fclose(fid1);

		    fid1 = fopen("df0dk_grid.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, df0dk_grid[i]);
		fclose(fid1);

		    fid1 = fopen("f_dist.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, f_dist[i]);
		fclose(fid1);

		    fid1 = fopen("thermal_driving_force.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, thermal_driving_force[i]);
		fclose(fid1);

		    fid1 = fopen("f0x1_f0.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, f0x1_f0[i]);
		fclose(fid1);

		    fid1 = fopen("elecgtric_driving_force.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %e\n", i+1, electric_driving_force[i]);
		fclose(fid1);
		    */
		    //getchar();

	// ----------------------------saving data ---------------------------------------------------
	//------------------------------ reading data -----------------------------------------------
		    /*
		    fid1 = fopen("nu_deformation.txt","r");
		    for (int i = 0; i < points; i++)
		    {
		        fgets(line, 1000, fid1);
		        sscanf(line, "%lf", &nu_deformation[i]);
		    }
		fclose(fid1);

		   fid1 = fopen("nu_ionizedimpurity.txt","r");
		    for (int i = 0; i < points; i++)
		    {
		        fgets(line, 1000, fid1);
		        sscanf(line, "%lf", &nu_ionizedimpurity[i]);
		    }
		fclose(fid1);

		   fid1 = fopen("nu_piezoelectric.txt","r");
		    for (int i = 0; i < points; i++)
		    {
		        fgets(line, 1000, fid1);
		        sscanf(line, "%lf", &nu_piezoelectric[i]);
		    }
		fclose(fid1);
		    */

	//------------------------------------------------------------------------------------------------------------------
	}  // if condition for type 'n'
	else 
	{
				
		// ionized impourity scattering
	        if (scattering_mechanisms[0]==1)
	        {
	        	nu_ii_p_funct(T_loop);
		}
		
		// POP scattering
	        if (scattering_mechanisms[1]==1)
		{
			nu_So_p_funct(T, T_loop, omega_LO);
		}
			
		// npop scattering
		if (scattering_mechanisms[2]==1)
		{
			nu_npop_p_funct(T);
		}

		// acoustic scattering
		if (scattering_mechanisms[3]==1)
		{
			nu_de_p_funct(T_loop);
		}
		//*/
			
	}
}


