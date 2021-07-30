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

	for (int counter1 = 0;counter1 < points;counter1++)
	{
		double k_dum =	k_grid[counter1];
		//f_dist(counter) = 1/(1+exp((energy(counter)-efef)/(k_B*T)))+g(counter);

		double arr[points];
		for (int i=0;i<points;i++)
		arr[i] = abs(k_grid[i] - kminus_grid[counter1]);
		int minus_index =FindMinInd(arr,points);

		for (int i=0;i<points;i++)
		arr[i] = abs(k_grid[i] - kplus_grid[counter1]);
		int plus_index =FindMinInd(arr,points);

		// If POP scattering is included
		if (scattering_mechanisms[1] == 1)
		{

			if ((lambda_o_minus_grid[counter1]==0) && (lambda_o_plus_grid[counter1]==0))
			    S_o_grid[counter1] = average_dummy;
			    // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
			else
			    S_o_grid[counter1] = (N_poph_atT+1-f_dist[minus_index])*lambda_o_minus_grid[counter1] +
			    (N_poph_atT+f_dist[plus_index])*lambda_o_plus_grid[counter1];


		}
	}

        for (int counter1=0;counter1<points;counter1++)
            S_o_grid_total[counter1] = S_o_grid[counter1];
}


