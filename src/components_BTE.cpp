#include"main.h"

void components_BTE(double T, int T_loop, double efefn, double efefp, int ii)
{

	double integral_numerator = 0;
	double integral_denominator = 0;
	double k_dum;

//---------------------------- components for BTE ------------------------------------------------------------------
	if(type == "n")
	{	
		    beta_constant = beta(T, T_loop);
		    // unit 1/nm

		    cout<< "Inverse screening length, beta = "<<beta_constant<<" (1/nm)"<<endl;

		    //--------------------------- df0dz_integral -----------------------------------------
		    int factr = 10;
		    for(int counter = 0;counter<=points-2;counter++)
		    {
		        double dk = (k_grid[counter+1]-k_grid[counter])/factr;
		        for (int ss = 0;ss<=factr-1;ss++)
		        {
		            integral_numerator = integral_numerator + dk*(pow(((k_grid[counter]+ss*dk)/pi),2))
		                                *f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T))*energy_n[counter]/(k_B*T);
		            // Part of equation (54) of Rode's book
		            integral_denominator = integral_denominator + dk*pow(((k_grid[counter]+ss*dk)/pi),2)
		            *f0(energy_n[counter],efefn,T)*(1-f0(energy_n[counter],efefn,T));
		            // Part of equation (54) of Rode's book
		        }
		    }

		    df0dz_integral = integral_numerator/integral_denominator;
		    //cout<<"df0dz_integral for n type = "<<df0dz_integral<<endl;
		    //--------------------------- df0dz_integral calculated -----------------------------------------
		    		    
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

		        thermal_driving_force[counter] = -1*v_n[counter]*df0dz(k_dum, efefn, T, df0dz_integral,coefficients_cond, kindex_cond, a11);

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
			    
		// pop scattering
		if(scattering_mechanisms[1]==1)
			pop_So(T, efefn, ii, T_loop);
//-------------------------------------------------------- pop_So calculated---------------------------------------

//----------------------POP scattering rate completed--------------------------------------------------------------------


//----------------------Ionized Impurity scattering rate --------------------------------------------------------------------
			        			
		// ionized impourity scattering
	        if (scattering_mechanisms[0]==1)
	        {
	        	nu_ii(epsilon_s[T_loop]);
	        }
//----------------------II scattering rate completed--------------------------------------------------------------------

//----------------------Acoustic scattering rate --------------------------------------------------------------------
		        
		// acoustic scattering 
	        if (scattering_mechanisms[3]==1)
	        {
            		nu_de(T);
		}	
//----------------------Acoustic scattering rate completed--------------------------------------------------------------------

//----------------------PZ scattering rate --------------------------------------------------------------------
		        
		// Piezoelectric scattering
	        if (scattering_mechanisms[4]==1)
	        {
         		nu_pe(T,P_piezo[T_loop],epsilon_s[T_loop]);
		}
//----------------------PZ scattering rate --------------------------------------------------------------------
			
//----------------------Disloaction scattering rate -------------------------------------------------------			
		// Dislocation scattering
	        if (scattering_mechanisms[6]==1)
	        {
			nu_dis(T,beta_constant,epsilon_s[T_loop]);
		}
//----------------------Disloaction scattering rate completed ------------------------------------------			

//----------------------Alloy scattering rate ----------------------------------------------------			
			
		// Alloy scattering
	        if (scattering_mechanisms[7]==1)
	        {
		        nu_alloy1();
		}
//----------------------Alloy scattering rate completed ------------------------------------------			

			
//---------------------- Neutral impurity scattering rate  ---- ------------------------------------------			

		// neutral impurity
	        if (scattering_mechanisms[9]==1)  
		{						
			nu_im(epsilon_s[T_loop],N_im[ii]);
		}
			
//---------------------- Neutral impurity scattering rate  completed ------------------------------------------			


//---------------------- Intervalley scattering rate  ---- ------------------------------------------			

		// intravalley scattering	
		if (scattering_mechanisms[8]==1)   
		{
			nu_iv_n(T);
		}
//------------------interavalley scattering rate completed ----------------------------
		
//---------------------- NPOP scattering rate  ---- ------------------------------------------			

	 	// npop scattering
		if (scattering_mechanisms[2]==1)  
		{
			nu_npop_n(T);						
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
				
		//--------------------------- df0dz_integral -----------------------------------------
		int factr = 10;
		for(int counter = 0;counter<=points-2;counter++)
		{
			double dk = (k_grid[counter+1]-k_grid[counter])/factr;
			for (int ss = 0;ss<=factr-1;ss++)
			{
			    integral_numerator = integral_numerator + dk*(pow(((k_grid[counter]+ss*dk)/pi),2))
						*f0(energy_p[counter],efefp,T)*(1-f0(energy_p[counter],efefp,T))*energy_p[counter]/(k_B*T);
			    // Part of equation (54) of Rode's book
			    integral_denominator = integral_denominator + dk*pow(((k_grid[counter]+ss*dk)/pi),2)
			    *f0(energy_p[counter],efefp,T)*(1-f0(energy_p[counter],efefp,T));
			    // Part of equation (54) of Rode's book
			}
		}

		df0dz_integral = integral_numerator/integral_denominator;
		//cout<<"df0dz_integral for p type = "<<df0dz_integral<<endl;
		//--------------------------- df0dz_integral calculated -----------------------------------------		    		    
//--------------------------------------common terms calculated -------------------------------------------------------

		for (int counter = 0;counter<points;counter++)
		{
			k_dum = k_grid[counter];
			// unit 1/nm
			//cout<<"counter+1 = "<<counter+1<<endl;
			//cout<<"k_dum = "<<k_dum<<endl;

		        df0dk_grid[counter] = df0dk(k_dum, T, efefp, coefficients_val, kindex_val, b11);
			 // unit (nm)

		        f_dist[counter] = f0(energy_p[counter],efefp,T);
		        //cout<<"In between "<<endl;
		        //cout<<"energy_p[counter]  = "<<energy_p[counter]<<endl;
		        //cout<<"efefp = "<<efefp<<endl;
		        //cout<<"T = "<<T<<endl;

		        thermal_driving_force[counter] = -1*v_p[counter]*df0dz(k_dum, efefp, T, df0dz_integral,coefficients_val, kindex_val, b11);

		        f0x1_f0[counter] = f0(energy_p[counter],efefp,T)*(1-f0(energy_p[counter],efefp,T));

		        electric_driving_force[counter] = -(1*E/h_bar)*df0dk_grid[counter]*1e-7;
			// unit is 1/s , hbar unit is eV-s so e is not multiplied in numerator

		        //cout<<"df0dk_grid[counter] =  "<<df0dk_grid[counter]<<endl;
		        //cout<<"f_dist[counter]  =  "<<f_dist[counter]<<endl;
		        //cout<<"thermal_driving_force[counter] =  "<<thermal_driving_force[counter]<<endl;
		        //cout<<"f0x1_f0[counter] =  "<<f0x1_f0[counter]<<endl;
		        //cout<<"electric_driving_force[counter]  = "<<electric_driving_force[counter]<<endl;
		}

//--------------------------------------common terms calculated completed  --------------------------------------------

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

		// total scattering rates calculated here
		for(int counter=0;counter<points;counter++)
		{
		        nu_el[counter] = nu_ionizedimpurity_p[counter][0][0]*scattering_mechanisms[0] + 
		        	nu_npop_p[counter][0][0] * scattering_mechanisms[2] + 
		        	nu_deformation_p[counter][0][0]*scattering_mechanisms[3]; 
		       
			denom[counter] = (nu_So_p[counter][0][0]*scattering_mechanisms[1] + nu_el[counter]);	
			//cout<<"nu_el[counter] = "<<nu_el[counter]<<endl;
		}			


//----------------------in scattering terms are calulated here ----------------------------------------------------------------
			    
		 // polar optical phonon scattering 
		if (scattering_mechanisms[1]==1)
		{
		    N_poph_atT = N_poph(omega_LO,T);
		    //cout<<"N_poph_atT = "<<N_poph_atT<<endl;

			int minus_index, plus_index;
			double arr[points];
								
//------------------------- In scattering term for POP scattering rate calculation  --------------------------------
			double const1, c_plus[points]={0}, c_minus[points]={0}, C_plus[points]={0}, C_minus[points]={0};
			double A, B, C, const2[points]={0}, const3;
			
			// Eq no. 20 of paper coupled band Ramu paper
			const1 = e*e*omega_LO*(1/epsilon_inf[T_loop] - 1/epsilon_s[T_loop])*(1/epsilon_0)/(16*pi*(h_bar*e));
			
			for (int counter = 0;counter < points;counter++)
			{    
			    const2[counter]= const1/(v_p[counter]/100);		
			    
			    const3 = const2[counter];
			    
			    c_plus[counter] = (k_grid[counter]*k_grid[counter] + kplus_grid_pop[counter]*kplus_grid_pop[counter])/(2*k_grid[counter]*kplus_grid_pop[counter]);      // 
			    
			    c_minus[counter] = (k_grid[counter]*k_grid[counter] + kminus_grid_pop[counter]*kminus_grid_pop[counter])/(2*k_grid[counter]*kminus_grid_pop[counter]);     // 

			    A = abs((1.0+c_plus[counter])/(1.0-c_plus[counter]));
			    B = (c_plus[counter]+3.0*c_plus[counter]*c_plus[counter]*c_plus[counter])/2.0;
			    C = 2.0 + 3.0 * (c_plus[counter]*c_plus[counter]);
			    
			    C_plus[counter] = abs(B*log(A) - C);    //
			     
			    A = abs((1.0+c_minus[counter])/(1.0-c_minus[counter]));
			    B = (c_minus[counter]+3.0*c_minus[counter]*c_minus[counter]*c_minus[counter])/2.0;
			    C = 2.0 + 3.0 * (c_minus[counter]*c_minus[counter]);

			    C_minus[counter] = abs(B*log(A) - C);   // 

			    lambda_i_plus_grid[counter] = abs(const3*C_plus[counter]*((N_poph_atT+1)*(1 - f0(energy_p[counter],efefp,T))
			    + (N_poph_atT)*f0(energy_p[counter],efefp,T)));

			    lambda_i_minus_grid[counter] = abs(const3*C_minus[counter]*((N_poph_atT)*(1 - f0(energy_p[counter],efefp,T))
			    + (N_poph_atT+1)*f0(energy_p[counter],efefp,T)));

				//---------------------------- code to debug -------------------------------------------------------------
				//cout<<"counter = "<<counter<<endl;

				//cout<<"lambda_i_plus_grid[counter] =  "<<lambda_i_plus_grid[counter]<<endl;
				//cout<<"lambda_i_minus_grid[counter] =  "<<lambda_i_minus_grid[counter]<<endl;
				//getchar();
			}
			
			/*
			// save results
			
			FILE *fid1;
			fid1 = fopen("electric_driving_force_p.txt","w");
			
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", electric_driving_force[i]);	
			fclose(fid1);


			fid1 = fopen("df0dk_p.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", df0dk_grid[i]);	
			fclose(fid1);

			fid1 = fopen("thermal_driving_force_p.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", thermal_driving_force[i]);	
			fclose(fid1);
		
			fid1 = fopen("const2.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", const2[i]);	
			fclose(fid1);

			fid1 = fopen("lambda_i_plus_p.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", lambda_i_plus_grid[i]);	
			fclose(fid1);

			fid1 = fopen("lambda_i_minus_p.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", lambda_i_minus_grid[i]);	
			fclose(fid1);

			fid1 = fopen("C_plus_in_p.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", C_plus[i]);	
			fclose(fid1);

			fid1 = fopen("C_minus_in_p.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", C_minus[i]);	
			fclose(fid1);
			
			fid1 = fopen("k_minus.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", kminus_grid_pop[i]);	
			fclose(fid1);
			
			fid1 = fopen("k_plus.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", kplus_grid_pop[i]);	
			fclose(fid1);
			
			fid1 = fopen("plus_index.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%d \n", plus_index_pop[i]);	
			fclose(fid1);

			fid1 = fopen("minus_index.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%d \n", minus_index_pop[i]);	
			fclose(fid1);

			fid1 = fopen("denom.txt","w");		
			for (int i = 0; i < points; i++)
				fprintf(fid1,"%e \n", denom[i]);	
			fclose(fid1);
			*/
			
		}
		
		
//------------------------- In scattering term for POP scattering rate calculated  --------------------------------

//-------------------------------------------------------- pop_So calculated---------------------------------------

//----------------------POP scattering rate completed--------------------------------------------------------------------
			
	}
}


