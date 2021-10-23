#include"main.h"

void solve_g(double T)
{

	FILE *fid;	
	char line[1000];	
	
	double Si_grid_all[points][iterations+1]={0};
	double Si_grid_all_so_pop[points][iterations+1]={0};
	double Si_grid_all_pop[points][iterations+1]={0};

	double Si_th_grid_all[points][iterations+1]={0};
	double Si_th_grid_all_so_pop[points][iterations+1]={0};
	double Si_th_grid_all_pop[points][iterations+1]={0};

	if(geometry==1)
	{	

	//---------------------------------------- solve BTE for g ----------------------------------------------------------
		for (int i=0;i<points;i++)
		{
			g[i]=0;

			g_rta[i]=0;

			g_old[i]=0;
			g_pop[i]=0;

			g_th[i]=0;
			g_th_old[i]=0;
			g_th_pop[i]=0;

			Si_grid[i]=0;
			Si_pop_grid[i]=0;
			Si_th_grid[i]=0;
			Si_th_pop_grid[i]=0;

		}

		for(int i=0;i<points;i++)
		{
			Si_grid_all[i][0] = k_grid[i];  	
			result_g[i][0] = k_grid[i];
			result_g_pop[i][0] = k_grid[i];
			result_g_th[i][0] = k_grid[i];

		}


		for (int iteration = 0;iteration<iterations;iteration++)
		{

			for (int i = 0;i < points;i++)
			{
				Si_grid[i] = 0;
				Si_th_grid[i]=0;
			}
			
			// If POP scattering is included
			if (scattering_mechanisms[1] == 1)
			{
				for (int i = 0;i < points;i++)
				{
					Si_pop_grid[i] = 0;
					Si_th_pop_grid[i]=0;
				}
			}			

			// If POP scattering is included
			if (scattering_mechanisms[1] == 1)
			{
				for (int i = 0;i < points;i++)
				{

					for (int m3 = 0;m3 < pop_number;m3++)
					{	            

						Si_grid[i] = Si_grid[i] + Sa_pop[m3][i]*g[minus_index_pop[m3][i]] +
						Se_pop[m3][i]*g[plus_index_pop[m3][i]];

						Si_pop_grid[i] = Si_pop_grid[i] + Sa_pop[m3][i]*g_pop[minus_index_pop[m3][i]] +
						Se_pop[m3][i]*g_pop[plus_index_pop[m3][i]];


						Si_th_grid[i] = Si_th_grid[i] + Sa_pop[m3][i]*g_th[minus_index_pop[m3][i]] +
						Se_pop[m3][i]*g_th[plus_index_pop[m3][i]];


						Si_th_pop_grid[i] = Si_th_pop_grid[i] + Sa_pop[m3][i]*g_th_pop[minus_index_pop[m3][i]] +
						Se_pop[m3][i]*g_th_pop[plus_index_pop[m3][i]];
					}
				}
			}


			for (int i=0;i<points;i++)
			{
				Si_grid_all[i][iteration+1] = Si_grid[i];  	
				Si_th_grid_all[i][iteration+1] = Si_th_grid[i];  	
				
				g[i] = (Si_grid[i]+electric_driving_force[i])/(denom[i]);
				g_th[i] = (Si_th_grid[i] + thermal_driving_force[i])/(denom[i]);
			}


			// If POP scattering is included
			if (scattering_mechanisms[1] == 1)
			{
			    for (int i=0;i<points;i++)
			    {
				Si_grid_all_pop[i][iteration+1] = Si_pop_grid[i];  	
				Si_th_grid_all_pop[i][iteration+1] = Si_th_pop_grid[i];  	
			    }
			    
				if(type=="n")
				{
					for (int i=0;i<points;i++)
					{
						g_pop[i] = (Si_pop_grid[i] + electric_driving_force[i])/nu_pop_total[i];
						g_th_pop[i] = (Si_th_pop_grid[i] + thermal_driving_force[i])/nu_pop_total[i];
					}
				}	
				else
				{
					for (int i=0;i<points;i++)
					{
						g_pop[i] = (Si_pop_grid[i] + electric_driving_force[i])/nu_So_p[i][0][0];
						g_th_pop[i] = (Si_th_pop_grid[i] + thermal_driving_force[i])/nu_So_p[i][0][0];
					}
				}

				for (int i=0;i<points;i++)
				{
					result_g_pop[i][iteration+1] = g_pop[i];
        			    	result_g_th_pop[i][iteration+1] = g_th_pop[i];
				}
			}

			for (int i=0;i<points;i++)
			{
				result_g[i][iteration+1] = g[i];

				result_g_th[i][iteration+1] = g_th[i];

				//fprintf('Iteration %d in BTE: at T = %5.2f K. Average change in g = %e \n',
				//iteration,T,sum(g-g_old)/points);

				g_old[i] = g[i];
				g_th_old[i] = g_th[i];

				if (iteration==0)
					g_rta[i] = g[i] ;
			}
		}   // iteration loop end here

		//---------------------------------------- solve BTE for g finished---------------------------------------

		//---------------------------------------- Solve BTE for g and h for hall mobility --------------------------------

		// common for both n type and ptype
		if(Bfield!=0)
		{
			//cout<<"inside bfield condiction"<<endl;
			//getchar();
			
		    for(int i=0;i<points;i++)
		    {
			beta1[i]=0;
			gH[i]=0;
			hH[i]=0;

			gH_rta[i]=0;
			hH_rta[i]=0;

			gH_pop[i]=0;
			hH_pop[i]=0;

			Si_grid_g[i]=0;
			Si_grid_h[i]=0;

			Si_pop_grid_g[i]=0;
			Si_pop_grid_h[i]=0;
		   }

		    for(int i=0;i<points;i++)
		    {
			result_gH[i][0] = k_grid[i];
			result_hH[i][0] = k_grid[i];
		    }	

		    for (int iteration = 0;iteration<iterations;iteration++)
		    {

			double beta1_pop[points];
		     		
			for(int i=0;i<points;i++)
			{
				Si_grid_g[i]=0;
				Si_grid_h[i]=0;
			}
		    	
		    // If POP scattering is included
		    if (scattering_mechanisms[1] == 1)
		    {
			for (int i = 0;i < points;i++)
			{

				for (int m3 = 0;m3 < pop_number;m3++)
				{	            
					Si_grid_g[i] = Si_grid_g[i] + Sa_pop[m3][i]*gH[minus_index_pop[m3][i]] +
					Se_pop[m3][i]*gH[plus_index_pop[m3][i]];

					Si_pop_grid_g[i] = Sa_pop[m3][i]*gH_pop[minus_index_pop[m3][i]] +
					Se_pop[m3][i]*gH_pop[plus_index_pop[m3][i]];

					Si_grid_h[i] = Si_grid_h[i] + Sa_pop[m3][i]*hH[minus_index_pop[m3][i]] +
					Se_pop[m3][i]*hH[plus_index_pop[m3][i]];

					Si_pop_grid_h[i] = Sa_pop[m3][i]*hH_pop[minus_index_pop[m3][i]] +
					Se_pop[m3][i]*hH_pop[plus_index_pop[m3][i]];

					/*		
					cout<<"i = "<<i<<endl;
					getchar();
					//*/	
				}                            
			    }
			}
			


			for (int i=0;i<points;i++)
			{

				beta1[i] = e*(v_n[i]*0.01)*Bfield/((h_bar*e)*(k_grid[i]*pow(10,9))*(denom[i]));
				
				beta1_pop[i] = e*(v_n[i]*0.01)*Bfield/((h_bar*e)*(k_grid[i]*pow(10,9))* (nu_pop_total[i]));

			    	gH[i] = (Si_grid_g[i] + electric_driving_force[i] + 
			    	beta1[i] * Si_grid_h[i])/((denom[i])*(1 + beta1[i]*beta1[i]));
			    	
			    	hH[i] = (Si_grid_h[i] - beta1[i] * electric_driving_force[i] -
			    	beta1[i] * Si_grid_g[i])/((denom[i]) * (1 + beta1[i]*beta1[i]));

			    	gH_pop[i] = (Si_pop_grid_g[i] + electric_driving_force[i] + 
			    	beta1_pop[i] * Si_pop_grid_h[i])/((nu_pop_total[i] )*(1 + beta1_pop[i]*beta1_pop[i]));
			    	
			    	hH_pop[i] = (Si_pop_grid_h[i] - beta1_pop[i] * electric_driving_force[i] 
			    	- beta1_pop[i] * Si_pop_grid_g[i])/((nu_pop_total[i]) * (1 + beta1_pop[i]*beta1_pop[i]));
			}
			
			if (iteration==0)
			{
				for (int i=0;i<points;i++)
				{
					gH_rta[i] = gH[i] ;
					hH_rta[i] = hH[i] ;
				}	        	
			}

			for (int i=0;i<points;i++)
			{
				result_gH[i][iteration+1] = gH[i];
				result_hH[i][iteration+1] = hH[i];
			}
		    }  // iteration loop ends here
		    	
	 	}  // if condition !Bfield terminated

		//----------------------- Solve BTE for g and h for hall mobility finished-------------
		
	} // if condition for geometry==1 finished here
	
	else if(geometry==2)
	{
		
	//---------------------------------------- solve BTE for g ----------------------------------------------------------
		    for (int i=0;i<points;i++)
		    {
			g[i]=0;

			g_rta[i]=0;

			g_old[i]=0;
			g_pop[i]=0;
			g_so_pop[i]=0;

			g_th[i]=0;
			g_th_old[i]=0;
			g_th_pop[i]=0;
			g_th_so_pop[i]=0;

			Si_grid[i]=0;
			Si_pop_grid[i]=0;
			Si_so_pop_grid[i]=0;
			
			Si_th_grid[i]=0;
			Si_th_pop_grid[i]=0;
			Si_th_so_pop_grid[i]=0;
			
			if(Bfield!=0)
			{
				// magnetic field related arrays
				beta1[i]=0;
				gH[i]=0;
				hH[i]=0;

				gH_rta[i]=0;
				hH_rta[i]=0;

				gH_pop[i]=0;
				hH_pop[i]=0;

				gH_so_pop[i]=0;
				hH_so_pop[i]=0;

				Si_grid_g[i]=0;
				Si_grid_h[i]=0;

				Si_pop_grid_g[i]=0;
				Si_pop_grid_h[i]=0;

				Si_so_pop_grid_g[i]=0;
				Si_so_pop_grid_h[i]=0;
			}
		     }
		     	
		    for(int i=0;i<points;i++)
		    {
			Si_grid_all[i][0] = k_grid[i];  	
			Si_grid_all_pop[i][0] = k_grid[i];  	
			Si_grid_all_so_pop[i][0] = k_grid[i];  	

			Si_th_grid_all[i][0] = k_grid[i];  	
			Si_th_grid_all_pop[i][0] = k_grid[i];  	
			Si_th_grid_all_so_pop[i][0] = k_grid[i];  	

			result_g[i][0] = k_grid[i];
			result_g_pop[i][0] = k_grid[i];
			result_g_so_pop[i][0] = k_grid[i];
			
			result_g_th[i][0] = k_grid[i];
			result_g_th_pop[i][0] = k_grid[i];
			result_g_th_so_pop[i][0] = k_grid[i];
		    }


		    for (int iteration = 0;iteration<iterations;iteration++)
		    {
					
			/*
			for (int i=0;i<points;i++)
			{
			    cout<<"energy_n[i] = "<<energy_n[i]<<endl;
			    cout<<"g[i] = "<<g[i]<<endl;
			    getchar();
			}
			//*/
					
			for (int i = 0;i < points;i++)
			{
				Si_grid[i] = 0;
				Si_th_grid[i]=0;
			}
			
			// If POP scattering is included
			if (scattering_mechanisms[1] == 1)
			{
				for (int i = 0;i < points;i++)
				{
					Si_pop_grid[i] = 0;
					Si_th_pop_grid[i]=0;
				}
			}			
			
			// If SO POP scattering is included
			if (scattering_mechanisms[10] == 1)
			{
				for (int i = 0;i < points;i++)
				{
					Si_so_pop_grid[i] = 0;
					Si_th_so_pop_grid[i]=0;
				}
			}		


			    // If POP scattering is included
			    if (scattering_mechanisms[1] == 1)
			    {
				for (int i = 0;i < points;i++)
				{
					for (int m3 = 0;m3 < pop_number;m3++)
					{	            
											
						Si_grid[i] = Si_grid[i] + Sa_pop[m3][i]*g[minus_index_pop[m3][i]]+
						Se_pop[m3][i]*g[plus_index_pop[m3][i]];

						Si_pop_grid[i] = Si_pop_grid[i] + Sa_pop[m3][i]*g_pop[minus_index_pop[m3][i]]+
						Se_pop[m3][i]*g_pop[plus_index_pop[m3][i]];

						Si_th_grid[i] = Si_th_grid[i] + Sa_pop[m3][i]*g_th[minus_index_pop[m3][i]]+
						Se_pop[m3][i]*g_th[plus_index_pop[m3][i]];

						Si_th_pop_grid[i] = Si_th_pop_grid[i] + Sa_pop[m3][i]*g_th_pop[minus_index_pop[m3][i]]+
						Se_pop[m3][i]*g_th_pop[plus_index_pop[m3][i]];

					}
							            				     			         
			    }  // end of i loop for points
			}  // end of pop condition

			// If SO POP scattering is included
			if (scattering_mechanisms[10] == 1)
			{
				for (int i = 0;i < points;i++)
				{							
					for (int m3 = 0;m3 < so_pop_number;m3++)
					{	            

						Si_grid[i] = Si_grid[i] + Sa_so_pop[m3][i]*g[minus_index_so_pop[m3][i]] +
						Se_so_pop[m3][i]*g[plus_index_so_pop[m3][i]];

						Si_so_pop_grid[i] = Si_so_pop_grid[i] + 
						Sa_so_pop[m3][i]*g_so_pop[minus_index_so_pop[m3][i]] +
						Se_so_pop[m3][i]*g_so_pop[plus_index_so_pop[m3][i]];

						Si_th_grid[i] = Si_th_grid[i] + Sa_so_pop[m3][i]*g_th[minus_index_so_pop[m3][i]] +
						Se_so_pop[m3][i]*g_th[plus_index_so_pop[m3][i]];

						Si_th_so_pop_grid[i] = Si_th_so_pop_grid[i] + 
						Sa_so_pop[m3][i]*g_th_so_pop[minus_index_so_pop[m3][i]] +
						Se_so_pop[m3][i]*g_th_so_pop[plus_index_so_pop[m3][i]];
					}
				}  // end of i loop for points

			} // end of SO pop condiction 
	 			 			
			for (int i=0;i<points;i++)
			{
				Si_grid_all[i][iteration+1] = Si_grid[i];  	
				Si_th_grid_all[i][iteration+1] = Si_th_grid[i];  	
				
				g[i] = (Si_grid[i]+electric_driving_force[i])/(denom[i]);
				g_th[i] = (Si_th_grid[i] + thermal_driving_force[i])/(denom[i]);
			}

			// If POP scattering is included
			if (scattering_mechanisms[1] == 1)
			{
			    for (int i=0;i<points;i++)
			    {
				Si_grid_all_pop[i][iteration+1] = Si_pop_grid[i];  	
				Si_th_grid_all_pop[i][iteration+1] = Si_th_pop_grid[i];  	

				g_pop[i] = (Si_pop_grid[i] + electric_driving_force[i])/nu_pop_total[i];
				g_th_pop[i] = (Si_th_pop_grid[i] + thermal_driving_force[i])/nu_pop_total[i];

			    }    
			}

			// If So POP scattering is included
			if (scattering_mechanisms[10] == 1)
			{
			    for (int i=0;i<points;i++)
			    {
			    
				Si_grid_all_so_pop[i][iteration+1] = Si_so_pop_grid[i];  	
				Si_th_grid_all_so_pop[i][iteration+1] = Si_th_so_pop_grid[i];
				
				g_so_pop[i] = (Si_so_pop_grid[i] + electric_driving_force[i])/nu_so_pop_total[i];
				g_th_so_pop[i] = (Si_th_so_pop_grid[i] + thermal_driving_force[i])/nu_so_pop_total[i];
			    	 
			    } 	
			}


			for (int i=0;i<points;i++)
			{
				result_g[i][iteration+1] = g[i];
				result_g_th[i][iteration+1] = g_th[i];

			//fprintf('Iteration %d in BTE: at T = %5.2f K. Average change in g = %e \n',iteration,T,sum(g-g_old)/points);

				g_old[i] = g[i];
				g_th_old[i] = g_th[i];
			}
			
			// pop result
			if (scattering_mechanisms[1] == 1)
			{
			    for (int i=0;i<points;i++)
			    {
			    	result_g_pop[i][iteration+1] = g_pop[i];
			    	result_g_th_pop[i][iteration+1] = g_th_pop[i];
			    }
			}

			// so pop result
			if (scattering_mechanisms[10] == 1)
			{
			    for (int i=0;i<points;i++)
			    {
			    	result_g_so_pop[i][iteration+1] = g_so_pop[i];
			    	result_g_th_so_pop[i][iteration+1] = g_th_so_pop[i];
			    }
			}
			
			
			if (iteration==0)
			{
				for (int i=0;i<points;i++)
					g_rta[i] = g[i] ;
			}	            
		}   // iteration loop end here

	//---------------------------------------- solve BTE for g finished----------------------------------------------------------

	//---------------------------------------- Solve BTE for g and h for hall mobility---------------------------------------------

		if(Bfield!=0)
		{
		    for(int i=0;i<points;i++)
		    {
			result_gH[i][0] = k_grid[i];
			result_hH[i][0] = k_grid[i];
		    }	

		    for (int iteration = 0;iteration<iterations;iteration++)
		    {

			double beta1_pop[points];
			double beta1_so_pop[points];

			for (int i = 0;i < points;i++)
			{
				Si_grid_g[i] = 0;
				Si_grid_h[i] = 0;
			}		     		

		    // If POP scattering is included
		    if (scattering_mechanisms[1] == 1)
		    {
			for (int i = 0;i < points;i++)
			{

				for (int m3 = 0;m3 < pop_number;m3++)
				{	            
					Si_grid_g[i] = Si_grid_g[i] + Sa_pop[m3][i]*gH[minus_index_pop[m3][i]] +
					Se_pop[m3][i]*gH[plus_index_pop[m3][i]];

					Si_pop_grid_g[i] = Sa_pop[m3][i]*gH_pop[minus_index_pop[m3][i]] +
					Se_pop[m3][i]*gH_pop[plus_index_pop[m3][i]];

					Si_grid_h[i] = Si_grid_h[i] + Sa_pop[m3][i]*hH[minus_index_pop[m3][i]] +
					Se_pop[m3][i]*hH[plus_index_pop[m3][i]];

					Si_pop_grid_h[i] = Sa_pop[m3][i]*hH_pop[minus_index_pop[m3][i]] +
					Se_pop[m3][i]*hH_pop[plus_index_pop[m3][i]];

					/*		
					cout<<"i = "<<i<<endl;
					getchar();
					//*/
				}	                            
			    }
			}
			
		    // If SO POP scattering is included
		    if (scattering_mechanisms[10] == 1)
		    {
			for (int i = 0;i < points;i++)
			{

				for (int m3 = 0;m3 < so_pop_number;m3++)
				{	            
					Si_grid_g[i] = Si_grid_g[i] + Sa_so_pop[m3][i]*gH[minus_index_pop[m3][i]] +
					Se_so_pop[m3][i]*gH[plus_index_pop[m3][i]];

					Si_so_pop_grid_g[i] = Sa_so_pop[m3][i]*gH_so_pop[minus_index_pop[m3][i]] +
					Se_so_pop[m3][i]*gH_so_pop[plus_index_pop[m3][i]];

					Si_grid_h[i] = Si_grid_h[i] + Sa_so_pop[m3][i]*hH[minus_index_pop[m3][i]] +
					Se_so_pop[m3][i]*hH[plus_index_pop[m3][i]];

					Si_so_pop_grid_h[i] = Sa_so_pop[m3][i]*hH_so_pop[minus_index_pop[m3][i]] +
					Se_so_pop[m3][i]*hH_so_pop[plus_index_pop[m3][i]];

					/*		
					cout<<"i = "<<i<<endl;
					getchar();
					//*/	                           
				} 
			    }
			}


			for (int i=0;i<points;i++)
			{
				beta1[i] = e*(v_n[i]*0.01)*Bfield/((h_bar*e)*(k_grid[i]*pow(10,9)) *	(denom[i]));
				
			    	gH[i] = (Si_grid_g[i] + electric_driving_force[i] + 
			    	beta1[i] * Si_grid_h[i])/((denom[i])*(1 + beta1[i]*beta1[i]));
			    	
			    	hH[i] = (Si_grid_h[i] - beta1[i] * electric_driving_force[i] - 
			    	beta1[i] * Si_grid_g[i])/((denom[i]) * (1 + beta1[i]*beta1[i]));

			}
			
			// If POP scattering is included
			if (scattering_mechanisms[1] == 1)
			{
			
				for (int i=0;i<points;i++)
				{
					beta1_pop[i] = e*(v_n[i]*0.01)*Bfield/((h_bar*e)*(k_grid[i]*pow(10,9)) * (nu_pop_total[i]));
				    	gH_pop[i] = (Si_pop_grid_g[i] + electric_driving_force[i] + 
				    	beta1_pop[i] * Si_pop_grid_h[i])/((nu_pop_total[i] )*(1 + beta1_pop[i]*beta1_pop[i]));
				    	
				    	hH_pop[i] = (Si_pop_grid_h[i] - beta1_pop[i] * electric_driving_force[i] -
				    	beta1_pop[i] * Si_pop_grid_g[i])/((nu_pop_total[i]) * (1 + beta1_pop[i]*beta1_pop[i]));
				}
			}
						
			// If SO POP scattering is included
			if (scattering_mechanisms[10] == 1)
			{
			
				for (int i=0;i<points;i++)
				{
				beta1_so_pop[i] = e*(v_n[i]*0.01)*Bfield/((h_bar*e)*(k_grid[i]*pow(10,9)) * (nu_so_pop_total[i]));
			
			    	gH_so_pop[i] = (Si_so_pop_grid_g[i] + electric_driving_force[i] + 
			    	beta1_so_pop[i] * Si_so_pop_grid_h[i])/((nu_so_pop_total[i] )*(1 + beta1_so_pop[i]*beta1_so_pop[i]));
			    	
			    	hH_so_pop[i] = (Si_so_pop_grid_h[i] - beta1_so_pop[i] * electric_driving_force[i] -
			    	beta1_so_pop[i] * Si_so_pop_grid_g[i])/((nu_so_pop_total[i]) * (1 + beta1_so_pop[i]*beta1_so_pop[i]));
				}
			}

			if (iteration==0)
			{
				for (int i=0;i<points;i++)
				{
					gH_rta[i] = gH[i] ;
					hH_rta[i] = hH[i] ;
				}	        	
			}
			for (int i=0;i<points;i++)
			{
				result_gH[i][iteration+1] = gH[i];
				result_hH[i][iteration+1] = hH[i];
			}

		    } // end of iteration loop
		    	
	 	}  // if condition !Bfield terminated
		
	}  // end of else if condition geometry==2	    


	//------------- save data ---------------------------------------------------	    
	
	/*		
	FILE *fid1;

	fid1 = fopen("electric_driving_force.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%e \n", electric_driving_force[i]);	
	fclose(fid1);

	fid1 = fopen("thermal_driving_force_n.txt","w");
	for (int i = 0; i < points; i++)
		fprintf(fid1,"%e \n", thermal_driving_force[i]);	
	fclose(fid1);


	fid1 = fopen("g_all.txt","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%d \t", i+1);
		for(int j=0;j<iterations+1;j++)
		fprintf(fid1,"%e ", result_g[i][j]);
		fprintf(fid1," \n");
	}
	fclose(fid1);

	// pop scattering 
	if (scattering_mechanisms[1] == 1)
	{
		fid1 = fopen("g_all_pop.txt","w");
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
			fprintf(fid1,"%e ", result_g_pop[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);
	}

	fid1 = fopen("Si_all.txt","w");
	for (int i = 0; i < points; i++)
	{
		fprintf(fid1,"%d \t", i+1);
		for(int j=0;j<iterations+1;j++)
		fprintf(fid1,"%e ", Si_grid_all[i][j]);
		fprintf(fid1," \n");
	}
	fclose(fid1);


	// save pop result
	if (scattering_mechanisms[1] == 1)
	{
		fid1 = fopen("g_all_pop_n.txt","w");
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e ", result_g_pop[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);

		fid1 = fopen("Si_all_pop_n.txt","w");
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e ", Si_grid_all_pop[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);
	}


	// save so pop result
	if (scattering_mechanisms[10] == 1)
	{

		fid1 = fopen("g_all_so_pop_n.txt","w");
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e ", result_g_so_pop[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);

		fid1 = fopen("Si_all_so_pop_n.txt","w");
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e ", Si_grid_all_so_pop[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);
	}


	if(Bfield!=0)
	{
		fid1 = fopen("gH_all_n.txt","w");
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e ", result_gH[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);
	
		fid1 = fopen("hH_all_n.txt","w");
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e ", result_hH[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);
	}
	//*/
	
	/*
	fid1 = fopen("g_original.txt","w");

	for (int i = 0; i < points; i++)
		fprintf(fid1,"%d    %e\n", i+1, g[i]);
	fclose(fid1);

	/* 
	fid1 = fopen("g_th.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, g_th[i]);
	fclose(fid1);

	fid1 = fopen("g_pop.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, g_pop[i]);
	fclose(fid1);


	fid1 = fopen("g_rta.txt","w");
	for (int i = 0; i < points; i++)
	fprintf(fid1,"%d    %e\n", i+1, g_rta[i]);
	fclose(fid1);
	*/
	//cout<<"saved result check here";
	//getchar();
	//getchar();
	//----------------------------------------------------------------------------------
	
}
