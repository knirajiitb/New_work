#include"main.h"

void solve_g(double T)
{

	FILE *fid;	
	char line[1000];	
	
	double S_i_grid_all[points][iterations+1]={0};
	double gplus[points][iterations+1]={0};
	double gminus[points][iterations+1]={0};

//---------------------------------------- solve BTE for g ----------------------------------------------------------
	if(type =="n")
	{	
            for (int i=0;i<points;i++)
            {
		g[i]=0;

		g_rta[i]=0;

		g_old[i]=0;
		g_LO[i]=0;
		g_iv[i]=0;

		g_th[i]=0;
		g_th_old[i]=0;
		g_LO_th[i]=0;


		S_i_grid[i]=0;
		S_iLO_grid[i]=0;
		S_i_th_grid[i]=0;
		S_iLO_th_grid[i]=0;

		beta1[i]=0;
		gH[i]=0;
		hH[i]=0;

		gH_rta[i]=0;
		hH_rta[i]=0;

		gH_LO[i]=0;
		hH_LO[i]=0;

		S_i_grid_g[i]=0;
		S_i_grid_h[i]=0;

		S_iLO_grid_g[i]=0;
		S_iLO_grid_h[i]=0;

	     }
	     	
            for(int counter=0;counter<points;counter++)
            {
		 S_i_grid_all[counter][0] = k_grid[counter];  	
                result_g[counter][0] = k_grid[counter];
                result_g_LO[counter][0] = k_grid[counter];
                result_g_th[counter][0] = k_grid[counter];
            }


            for (int iteration = 0;iteration<iterations;iteration++)
            {
		
		/*
                for (int counter1=0;counter1<points;counter1++)
                {
                    //cout<<"energy_n[counter1] = "<<energy_n[counter1]<<endl;
                    //cout<<"g[counter1] = "<<g[counter1]<<endl;
                    //getchar();
                }
		*/

                    // If POP scattering is included
                    if (scattering_mechanisms[1] == 1)
                    {
		        for (int counter1 = 0;counter1 < points;counter1++)
		        {

                        S_i_grid[counter1] = (N_poph_atT + f_dist[counter1])*lambda_i_minus_grid[counter1]*g[minus_index_pop[counter1]]+
                        (N_poph_atT+1-f_dist[counter1])*lambda_i_plus_grid[counter1]*g[plus_index_pop[counter1]];

                        S_iLO_grid[counter1] = (N_poph_atT+f_dist[counter1])*lambda_i_minus_grid[counter1]*g_LO[minus_index_pop[counter1]]+
                        (N_poph_atT+1-f_dist[counter1])*lambda_i_plus_grid[counter1]*g_LO[plus_index_pop[counter1]];


                        S_i_th_grid[counter1] = (N_poph_atT + f_dist[counter1])*lambda_i_minus_grid[counter1]*g_th[minus_index_pop[counter1]]+
                        (N_poph_atT + 1 - f_dist[counter1])*lambda_i_plus_grid[counter1]*g_th[plus_index_pop[counter1]];


                        S_iLO_th_grid[counter1] = (N_poph_atT + f_dist[counter1])*lambda_i_minus_grid[counter1]*g_LO_th[minus_index_pop[counter1]]+
                        (N_poph_atT+1-f_dist[counter1])*lambda_i_plus_grid[counter1]*g_LO_th[plus_index_pop[counter1]];

                            
			 S_i_grid_all[counter1][iteration+1] = S_i_grid[counter1];  	

                            /*
                            cout<<"counter1 = "<<counter1<<endl;
                            cout<<"S_i_grid[counter1] = "<<S_i_grid[counter1]<<endl;
                            cout<<"S_iLo_grid[counter1] = "<<S_iLO_grid[counter1]<<endl;
                            cout<<"S_i_th_grid[counter1] =    "<<S_i_th_grid[counter1]<<endl;
                            cout<<"S_iLO_th_grid[counter1] =   "<<S_iLO_th_grid[counter1]<<endl;
			     cout<<"N_poph_atT =    "<<N_poph_atT<<endl;
			     cout<<"f_dist[minus_index_pop[counter1]]   "<<f_dist[minus_index_pop[counter1]]<<endl;
			     cout<<"lambda_o_minus_grid[counter1] =    "<<lambda_o_minus_grid[counter1]<<endl;
			     cout<<"lambda_o_plus_grid[counter1]   "<<lambda_o_plus_grid[counter1]<<endl;
			     getchar();
			     //*/	                            
                    }

                }
         

                for (int counter1=0;counter1<points;counter1++)
                {
                    g[counter1] = (S_i_grid[counter1]+electric_driving_force[counter1])/(denom[counter1]);
                    g_th[counter1] = (S_i_th_grid[counter1] + thermal_driving_force[counter1])/(denom[counter1]);
                }

                // If POP scattering is included
                if (scattering_mechanisms[1] == 1)
                {
                    for (int counter1=0;counter1<points;counter1++)
                        g_LO[counter1] = (S_iLO_grid[counter1] + electric_driving_force[counter1])/S_o_grid_total[counter1];
                }

                for (int counter1=0;counter1<points;counter1++)
                {
                    result_g[counter1][iteration+1] = g[counter1];

                    result_g_th[counter1][iteration+1] = g_th[counter1];

                //fprintf('Iteration %d in BTE: at T = %5.2f K. Average change in g = %e \n',iteration,T,sum(g-g_old)/points);

                    g_old[counter1] = g[counter1];
                    g_th_old[counter1] = g_th[counter1];

                    if (iteration==0)
                       g_rta[counter1] = g[counter1] ;
		     
		     /*
		     if(iteration == iterations-1)
		     {	
			     cout<<"counter1 =   "<<counter1<<endl;	
			     cout<<"g[counter1] =   "<<g[counter1]<<endl;
			     cout<<"g_th[counter1] =   "<<g_th[counter1]<<endl;
			     cout<<"g_LO[counter1] =   "<<g_LO[counter1]<<endl;
			     getchar();
		     }
		     //*/
                }
            }   // iteration loop end here
		
		/*
		FILE *fid1;

		fid1 = fopen("electric_driving_force_n.txt","w");
		for (int i = 0; i < points; i++)
			fprintf(fid1,"%e \n", electric_driving_force[i]);	
		fclose(fid1);

		fid1 = fopen("g_all_n.txt","w");


		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e ", result_g[i][j]);
			fprintf(fid1," \n");
		}

		fclose(fid1);
		
		fid1 = fopen("Si_all_n.txt","w");

		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e ", S_i_grid_all[i][j]);
			fprintf(fid1," \n");
		}

		fclose(fid1);
		*/
		
		/*
		for (int counter1=0;counter1<points;counter1++)
		{
			     cout<<"counter1 =   "<<counter1<<endl;	
			     cout<<"g[counter1] =   "<<g[counter1]<<endl;
			     cout<<"g_th[counter1] =   "<<g_th[counter1]<<endl;
			     cout<<"g_LO[counter1] =   "<<g_LO[counter1]<<endl;
			     getchar();
		}
		*/
		
//---------------------------------------- solve BTE for g finished----------------------------------------------------------


//---------------------------------------- Solve BTE for g and h for hall mobility---------------------------------------------


		if(Bfield!=0)
		{
		    for (int iteration = 0;iteration<iterations;iteration++)
		    {

		        double beta1_LO[points];
		     		
	            // If POP scattering is included
	            if (scattering_mechanisms[1] == 1)
	            {
		        for (int counter1 = 0;counter1 < points;counter1++)
		        {

		                S_i_grid_g[counter1] = (N_poph_atT + f_dist[counter1])*lambda_i_minus_grid[counter1]*gH[minus_index_pop[counter1]] +
		                (N_poph_atT+1-f_dist[counter1])*lambda_i_plus_grid[counter1]*gH[plus_index_pop[counter1]];

		                S_iLO_grid_g[counter1] = (N_poph_atT+f_dist[counter1])*lambda_i_minus_grid[counter1]*gH_LO[minus_index_pop[counter1]]+
		                (N_poph_atT+1-f_dist[counter1])*lambda_i_plus_grid[counter1]*gH_LO[plus_index_pop[counter1]];

		                S_i_grid_h[counter1] = (N_poph_atT + f_dist[counter1])*lambda_i_minus_grid[counter1]*hH[minus_index_pop[counter1]]+
		                (N_poph_atT+1-f_dist[counter1])*lambda_i_plus_grid[counter1]*hH[plus_index_pop[counter1]];

		                S_iLO_grid_h[counter1] = (N_poph_atT+f_dist[counter1])*lambda_i_minus_grid[counter1]*hH_LO[minus_index_pop[counter1]]+
		                (N_poph_atT+1-f_dist[counter1])*lambda_i_plus_grid[counter1]*hH_LO[plus_index_pop[counter1]];


		                    /*		
		                    cout<<"counter1 = "<<counter1<<endl;
		                    cout<<"S_i_grid_g[counter1] = "<<S_i_grid_g[counter1]<<endl;
		                    cout<<"S_iLo_grid_g[counter1] = "<<S_iLO_grid_g[counter1]<<endl;
		                    cout<<"S_i_grid_h[counter1] = "<<S_i_grid_h[counter1]<<endl;
		                    cout<<"S_iLo_grid_h[counter1] = "<<S_iLO_grid_h[counter1]<<endl;
				     cout<<"N_poph_atT =    "<<N_poph_atT<<endl;
				     cout<<"f_dist[minus_index_pop[counter1]]   "<<f_dist[minus_index_pop[counter1]]<<endl;
				     cout<<"lambda_o_minus_grid[counter1] =    "<<lambda_o_minus_grid[counter1]<<endl;
				     cout<<"lambda_o_plus_grid[counter1]   "<<lambda_o_plus_grid[counter1]<<endl;
				     getchar();
				     //*/	                            


		            }
		        }
		        


		        for (int counter1=0;counter1<points;counter1++)
		        {

		        	beta1[counter1] = e*(v_n[counter1]*0.01)*Bfield/((h_bar*e)*(k_grid[counter1]*pow(10,9)) * 				(denom[counter1]));
		        	
		        	beta1_LO[counter1] = e*(v_n[counter1]*0.01)*Bfield/((h_bar*e)*(k_grid[counter1]*pow(10,9)) * 				(S_o_grid_total[counter1]));

		            	gH[counter1] = (S_i_grid_g[counter1] + electric_driving_force[counter1] + beta1[counter1] * 			S_i_grid_h[counter1])/((denom[counter1])*(1 + beta1[counter1]*beta1[counter1]));
		            	
		            	hH[counter1] = (S_i_grid_h[counter1] - beta1[counter1] * electric_driving_force[counter1] -                    	 			beta1[counter1] * S_i_grid_g[counter1])/((denom[counter1]) * (1 + 				beta1[counter1]*beta1[counter1]));

		            	gH_LO[counter1] = (S_iLO_grid_g[counter1] + electric_driving_force[counter1] + beta1_LO[counter1] * 			S_iLO_grid_h[counter1])/((S_o_grid_total[counter1] )*(1 + beta1_LO[counter1]*beta1_LO[counter1]));
		            	
		            	hH_LO[counter1] = (S_iLO_grid_h[counter1] - beta1_LO[counter1] * electric_driving_force[counter1] 				- beta1_LO[counter1] * S_iLO_grid_g[counter1])/((S_o_grid_total[counter1]) * (1 + 				       beta1_LO[counter1]*beta1_LO[counter1]));
		        }
		        
		        if (iteration==0)
			{
				for (int counter1=0;counter1<points;counter1++)
				{
					gH_rta[counter1] = gH[counter1] ;
					hH_rta[counter1] = hH[counter1] ;
				}	        	
			}

		    }
		    	
	 	}  // if condition !Bfield terminated

		//----------------------- Solve BTE for g and h for hall mobility finished-------------

		//------------- save data ---------------------------------------------------	    
		    /*		
		    FILE *fid1;
		    fid1 = fopen("g_original.txt","w");
		    for (int i = 0; i < points; i++)
			fprintf(fid1,"%d    %e\n", i+1, g[i]);
		fclose(fid1);
			 
			 /* 
		    fid1 = fopen("g_th.txt","w");
		    for (int i = 0; i < points; i++)
			fprintf(fid1,"%d    %e\n", i+1, g_th[i]);
		fclose(fid1);

		    fid1 = fopen("g_LO.txt","w");
		    for (int i = 0; i < points; i++)
			fprintf(fid1,"%d    %e\n", i+1, g_LO[i]);
		fclose(fid1);

		    fid1 = fopen("g_iv.txt","w");
		    for (int i = 0; i < points; i++)
			fprintf(fid1,"%d    %e\n", i+1, g_iv[i]);
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
	else    // for p type
	{
		
		for (int i=0;i<points;i++)
		{
			g[i]=0;

			g_rta[i]=0;

			g_old[i]=0;
			g_LO[i]=0;
			g_iv[i]=0;

			g_th[i]=0;
			g_th_old[i]=0;
			g_LO_th[i]=0;
						
			S_i_grid[i]=0;
			S_iLO_grid[i]=0;
			S_i_th_grid[i]=0;
			S_iLO_th_grid[i]=0;

			/*		
			// for magnetic field
			beta1[i]=0;
			gH[i]=0;
			hH[i]=0;

			gH_rta[i]=0;
			hH_rta[i]=0;

			gH_LO[i]=0;
			hH_LO[i]=0;

			S_i_grid_g[i]=0;
			S_i_grid_h[i]=0;

			S_iLO_grid_g[i]=0;
			S_iLO_grid_h[i]=0;
			*/	
		}
	     	
		for(int counter=0;counter<points;counter++)
		{	
	 		gplus[counter][0] = k_grid[counter]; 
	 		gminus[counter][0] = k_grid[counter];  	
			S_i_grid_all[counter][0] = k_grid[counter]; 
			result_g[counter][0] = k_grid[counter];
			result_g_LO[counter][0] = k_grid[counter];
			result_g_th[counter][0] = k_grid[counter];
		}

		
		for (int iteration = 0;iteration<iterations;iteration++)
		{
			
			//cout<<"iteration = "<<iteration<<endl<<endl;
			//getchar();
			//getchar();
			//getchar();
			/*
			for (int counter1=0;counter1<points;counter1++)
			{
			//cout<<"energy_p[counter1] = "<<energy_p[counter1]<<endl;
			//cout<<"g[counter1] = "<<g[counter1]<<endl;
			//getchar();
			}
			*/

			// If POP scattering is included
			if (scattering_mechanisms[1] == 1)
			{
				for (int counter1 = 0;counter1 < points;counter1++)
				{        
					gplus[counter1][iteration+1] = g[plus_index_pop[counter1]];
					gminus[counter1][iteration+1] = g[minus_index_pop[counter1]];
					
					/*
					if(iteration!=0 && (counter1 == 0 || counter1 == 4 || counter1 == 100 || counter1 == 200
					|| counter1 == 300   || counter1 == 400))
					{
						
					cout<<"counter1   =  "<<counter1<<endl;
					cout<<"g[plus_index_pop[counter1]]  =   "<<g[plus_index_pop[counter1]]<<endl;
					cout<<"g[minus_index_pop[counter1]]  =   "<<g[minus_index_pop[counter1]]<<endl;

					cout<<"lambda_plus = "<<lambda_i_plus_grid[counter1]<<endl;
					cout<<"lambda_minus = "<<lambda_i_minus_grid[counter1]<<endl;

					cout<<"plus_index_pop[counter1] = "<<plus_index_pop[counter1]<<endl;
					cout<<"minus_index_pop[counter1] = "<<minus_index_pop[counter1]<<endl;

					cout<<"In plus = "<<lambda_i_plus_grid[counter1]*g[plus_index_pop[counter1]]<<endl;
					cout<<"In minus = "<<lambda_i_minus_grid[counter1]*g[minus_index_pop[counter1]]<<endl;

					getchar();
					}
					*/
							
				        S_i_grid[counter1] = lambda_i_minus_grid[counter1]*g[minus_index_pop[counter1]] +
				        lambda_i_plus_grid[counter1]*g[plus_index_pop[counter1]];

				        S_iLO_grid[counter1] = lambda_i_minus_grid[counter1]*g_LO[minus_index_pop[counter1]] +
				        lambda_i_plus_grid[counter1]*g_LO[plus_index_pop[counter1]];


				        S_i_th_grid[counter1] = lambda_i_minus_grid[counter1]*g_th[minus_index_pop[counter1]] +
				        lambda_i_plus_grid[counter1]*g_th[plus_index_pop[counter1]];


				        S_iLO_th_grid[counter1] = lambda_i_minus_grid[counter1]*g_LO_th[minus_index_pop[counter1]] +
				        lambda_i_plus_grid[counter1]*g_LO_th[plus_index_pop[counter1]];
					 
					 S_i_grid_all[counter1][iteration+1] = S_i_grid[counter1];  	
				

					/*
					cout<<"counter1 = "<<counter1<<endl;
					cout<<"S_i_grid[counter1] = "<<S_i_grid[counter1]<<endl;
					cout<<"S_iLo_grid[counter1] = "<<S_iLO_grid[counter1]<<endl;
					cout<<"S_i_th_grid[counter1] =    "<<S_i_th_grid[counter1]<<endl;
					cout<<"S_iLO_th_grid[counter1] =   "<<S_iLO_th_grid[counter1]<<endl;
					cout<<"lambda_o_minus_grid[counter1] =    "<<lambda_o_minus_grid[counter1]<<endl;
					cout<<"lambda_o_plus_grid[counter1]   "<<lambda_o_plus_grid[counter1]<<endl;
					getchar();
					//*/	                            
				}

			}


			for (int counter1=0;counter1<points;counter1++)
			{
				g[counter1] = (S_i_grid[counter1]+electric_driving_force[counter1])/(denom[counter1]);
				g_th[counter1] = (S_i_th_grid[counter1] + thermal_driving_force[counter1])/(denom[counter1]);
			}

			// If POP scattering is included
			if (scattering_mechanisms[1] == 1)
			{
				for (int counter1=0;counter1<points;counter1++)
				g_LO[counter1] = (S_iLO_grid[counter1] + electric_driving_force[counter1])/nu_So_p[counter1][0][0];
			}

			for (int counter1=0;counter1<points;counter1++)
			{
				result_g[counter1][iteration+1] = g[counter1];

				result_g_th[counter1][iteration+1] = g_th[counter1];

				//fprintf('Iteration %d in BTE: at T = %5.2f K. Average change in g = %e \n',iteration,T,sum(g-g_old)/points);

				g_old[counter1] = g[counter1];
				g_th_old[counter1] = g_th[counter1];

				if (iteration==0)
					g_rta[counter1] = g[counter1] ;

				/*
				if(iteration == iterations-1)
				{	
				cout<<"counter1 =   "<<counter1<<endl;	
				cout<<"g[counter1] =   "<<g[counter1]<<endl;
				cout<<"g_th[counter1] =   "<<g_th[counter1]<<endl;
				cout<<"g_LO[counter1] =   "<<g_LO[counter1]<<endl;
				getchar();
				}
				//*/
			}
		}   // iteration loop end here

		/*
		FILE *fid1;
		
		fid1 = fopen("g_all_p.txt","w");		
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e\t", result_g[i][j]);
			fprintf(fid1," \n");
		}	
		fclose(fid1);
		/*

		fid1 = fopen("Si_all_p.txt","w");		
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e\t", S_i_grid_all[i][j]);
			fprintf(fid1," \n");
		}
	
		fclose(fid1);

		fid1 = fopen("gplus.txt","w");		
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e\t", gplus[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);
		
		fid1 = fopen("gminus.txt","w");		
		for (int i = 0; i < points; i++)
		{
			fprintf(fid1,"%d \t", i+1);
			for(int j=0;j<iterations+1;j++)
				fprintf(fid1,"%e\t", gminus[i][j]);
			fprintf(fid1," \n");
		}
		fclose(fid1);

		//*/
	
		/*
		for (int counter1=0;counter1<points;counter1++)
		{
			     cout<<"counter1 =   "<<counter1<<endl;	
			     cout<<"g[counter1] =   "<<g[counter1]<<endl;
			     cout<<"g_th[counter1] =   "<<g_th[counter1]<<endl;
			     cout<<"g_LO[counter1] =   "<<g_LO[counter1]<<endl;
			     getchar();
		}
		*/
		
//---------------------------------------- solve BTE for g finished----------------------------------------------------------
	
	}
}
