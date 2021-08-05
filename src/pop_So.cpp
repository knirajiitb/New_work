#include"main.h"

void pop_So(double T, double efef, int ii)
{

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


