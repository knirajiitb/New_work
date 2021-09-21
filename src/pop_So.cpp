#include"main.h"

void pop_So(double T, double efef, int ii, int T_loop)
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
	

	for (int i=0;i<points;i++)
	{
		S_o_grid[i] = 1;
		S_o_grid_total[i]=0;
	}	


//------------------------------------------ S_o_grid calculated ----------------------------------------------
	
	double sum=0;
	for (int counter1 = 0;counter1 < points;counter1++)
	    sum = sum +  S_o_grid[counter1];

	double average_dummy = sum/points;
	//cout<<"average dummy = "<<average_dummy<<endl;
	//cout<<"N_poph_atT   =   "<<N_poph_atT<<endl;
	
	for (int counter1 = 0;counter1 < points;counter1++)
	{
		//cout<<"counter = "<<counter1<<endl;
		//cout<<"plus_index   =   "<<plus_index_pop[counter1]<<endl;
		//cout<<"minus_index   =   "<<minus_index_pop[counter1]<<endl;

		if ((lambda_o_minus_grid[counter1]==0) && (lambda_o_plus_grid[counter1]==0))
		    S_o_grid[counter1] = average_dummy;
		    // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
		else
		    S_o_grid[counter1] = (N_poph_atT+1-f_dist[minus_index_pop[counter1]])*lambda_o_minus_grid[counter1] +
		    (N_poph_atT+f_dist[plus_index_pop[counter1]])*lambda_o_plus_grid[counter1];
		/*
		cout<<"(N_poph_atT+1-f_dist[minus_index_pop[counter1]])  = "<<(N_poph_atT+1-f_dist[minus_index_pop[counter1]])<<endl;
		cout<<"(N_poph_atT+f_dist[plus_index_pop[counter1]])   =   "<<(N_poph_atT+f_dist[plus_index_pop[counter1]])<<endl;
		cout<<"lambda_o_plus_grid[counter1]   = "<<lambda_o_plus_grid[counter1]<<endl;
		cout<<"lambda_o_minus_grid[counter1] = "<<lambda_o_minus_grid[counter1]<<endl;
		cout<<"f_dist[plus_index_pop[counter1]]  = "<<f_dist[plus_index_pop[counter1]]<<endl;
		cout<<"f_dist[minus_index_pop[counter1]]  = "<<f_dist[minus_index_pop[counter1]]<<endl;
		cout<<"S_o_grid[counter1]  =    "<<S_o_grid[counter1]<<endl;
		getchar();
		*/
	}

        for (int counter1=0;counter1<points;counter1++)
            S_o_grid_total[counter1] = S_o_grid[counter1];
}


