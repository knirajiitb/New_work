#include"main.h"

void save_results()
{

	FILE *fid1, *fid2;

	int count;
	if(variation==0)
		count = count_t;
	else
		count = count_d;


	// --------- saving mobility results ----------------------- 
	if (ispin == 1)
	{
		fid1 = fopen("mobility.dat","w");		
		if(Bfield!=0) // && type == "n")
			fid2 = fopen("mobility_hall.dat","w");				
	}
	else
	{	
		if(ispin==2 && kk==0)
		{
			fid1 = fopen("mobility_up_spin.dat","w");
			if(Bfield!=0) // && type == "n")
				fid2 = fopen("mobility_up_spin_hall.dat","w");				
		}
		
		if(ispin==2 && kk==1)		
		{
			fid1 = fopen("mobility_down_spin.dat","w");

			if(Bfield!=0) // && type == "n")
				fid2 = fopen("mobility_down_spin_hall.dat","w");				
		}
	}			

	if (variation==0)   // temperature variation
	{
		fprintf(fid1,"#Temperature(K)");
		if(Bfield!=0) // && type == "n")
			fprintf(fid2,"#Temperature(K)");		
	}
	else
	{
		if(geometry==1)
		{
			fprintf(fid1,"#Doping(cm^-3)");

			if(Bfield!=0) // && type == "n")
				fprintf(fid2,"#Doping(cm^-3)");						
		}
		else if(geometry==2)
		{
			fprintf(fid1,"#Doping(cm^-2)");
			if(Bfield!=0) // && type == "n")
				fprintf(fid2,"#Doping(cm^-2)");
		}
	}
		
	
	if(geometry==1 && type=="n")
	{
		fprintf(fid1,"Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");

		for (int i = 0; i <count ; i++)
		fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e    %e\n",
		calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1], calc_mobility_npop[i][1], calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_dis[i][1], calc_mobility_to[i][1],
		calc_mobility_alloy[i][1], calc_mobility_iv[i][1], calc_mobility_neutral[i][1]);


		if(Bfield!=0)
		{				
			for (int i = 0; i <count ; i++)
			    fprintf(fid2,"%e     %e          %e          %e     %e    %e    %e    %e    %e    %e    %e    %e    %e\n",
				    calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1], calc_mobility_hall_npop[i][1], 
				    calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_dis[i][1], calc_mobility_hall_to[i][1],
				    calc_mobility_hall_alloy[i][1], calc_mobility_hall_iv[i][1], calc_mobility_hall_neutral[i][1]);
			fclose(fid2);
		}
	}
	else if (geometry==1 && type=="p")
	{
		fprintf(fid1,"  Mobility(cm^2/V-s)    Mobility_rta     Mobility_remote_ii   Mobility_po  Mobility_npop   Mobility_de\n");

		for (int i = 0; i <count ; i++)
		fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    \n",
		calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], 
		calc_mobility_po[i][1], calc_mobility_npop[i][1], calc_mobility_de[i][1]);

		if(Bfield!=0)
		{				
			fprintf(fid2," Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii   Mobility_po   Mobility_npop	Mobility_de    \n");

			for (int i = 0; i <count ; i++)
			fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    \n",
			calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], calc_mobility_hall_ii[i][1], 
			calc_mobility_hall_po[i][1], calc_mobility_hall_npop[i][1], calc_mobility_hall_de[i][1]);

			fclose(fid2);
		}
	}   // end of geometry ==1 part if 
	else  // geometry==2
	{
		fprintf(fid1," Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     mobility_so_pop \n");

		for (int i = 0; i <count ; i++)
		{
		
		fprintf(fid1,"%e     %e          %e          %e     %e    %e    %e    %e    %e  \n",
		calc_mobility[i][0], calc_mobility[i][1], calc_mobility_rta[i][1], calc_mobility_ii[i][1], calc_mobility_po[i][1],
		calc_mobility_npop[i][1], calc_mobility_de[i][1], calc_mobility_pe[i][1], calc_mobility_so_pop[i][1]);
		
		}

		if(Bfield!=0)
		{				
			for (int i = 0; i <count ; i++)
				fprintf(fid2,"%e     %e          %e          %e     %e    %e    %e    %e    %e  \n",
				calc_mobility_hall[i][0], calc_mobility_hall[i][1], calc_mobility_hall_rta[i][1], 
				calc_mobility_hall_ii[i][1], calc_mobility_hall_po[i][1], calc_mobility_hall_npop[i][1],
				calc_mobility_hall_de[i][1], calc_mobility_hall_pe[i][1], calc_mobility_hall_so_po[i][1]);
			fclose(fid2);
		}

	}  // end of geometry ==2 part 
	
	fclose(fid1);
	// --------- saving mobility results completed -----------------------
	


	// --------- saving conductivity results ----------------------- 
	if (ispin == 1)
	{
		fid1 = fopen("conductivity.dat","w");		
		if(Bfield!=0) // && type == "n")
			fid2 = fopen("conductivity_hall.dat","w");				
	}
	else
	{	
		if (ispin == 2 && kk==0)
		{
			fid1 = fopen("conductivity_up_spin.dat","w");
			if(Bfield!=0) // && type == "n")
				fid2 = fopen("conductivity_up_spin_hall.dat","w");
		}
		else
		{
			fid1 = fopen("conductivity_down_spin.dat","w");
			if(Bfield!=0) // && type == "n")
				fid2 = fopen("conductivity_down_spin_hall.dat","w");
		}
	}			

	if (variation==0)   // temperature variation
	{
		fprintf(fid1,"#Temperature(K)");
		if(Bfield!=0) // && type == "n")
			fprintf(fid2,"#Temperature(K)");		
	}
	else
	{
		if(geometry==1)
		{
			fprintf(fid1,"#Doping(cm^-3)");
			if(Bfield!=0) // && type == "n")
				fprintf(fid2,"#Doping(cm^-3)");							
		}
		else if(geometry==2)
		{
			fprintf(fid1,"#Doping(cm^-2)");
			if(Bfield!=0) // && type == "n")
				fprintf(fid2,"#Doping(cm^-2)");		
		}
	}
		
	fprintf(fid1,"      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA) \n");
	for (int i = 0; i <count; i++)
	    fprintf(fid1," %e         %e              %e \n", calc_sigma[i][0], calc_sigma[i][1], calc_sigma_rta[i][1]);
	fclose(fid1);

	if(Bfield!=0) // && type == "n")
	{
		fprintf(fid2,"      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA) \n");
		for (int i = 0; i <count; i++)
		    fprintf(fid2," %e         %e              %e \n", 
		    calc_sigma_hall[i][0], calc_sigma_hall[i][1], calc_sigma_hall_rta[i][1]);
		fclose(fid2);
	}	
	
	// --------- saving thermopower results ----------------------- 
	if(geometry==1 && type == "n")
	{
		if (ispin == 1)
			fid1 = fopen("thermopower.dat","w");		
		else
		{	if (ispin == 2 && kk==0)
				fid1 = fopen("thermopower_up_spin.dat","w");
			else
				fid1 = fopen("thermopower_down_spin.dat","w");
		}			

		if (variation==0)   // temperature variation
			fprintf(fid1,"#Temperature(K)");
		else
			fprintf(fid1,"#Doping(cm^-3)");

		fprintf(fid1,"#	Thermopower(uV/K) \n");

		for (int i = 0; i <count ; i++)
		    fprintf(fid1," %e        %e\n", calc_thermopower[i][0], calc_thermopower[i][1]);
		fclose(fid1);	
	}
	
	//-------------saving hall factor ----------------------------
	if(Bfield!=0)
	{
	
		if (ispin == 1)
			fid1 = fopen("hall_factor.dat","w");
		if (ispin == 2 && kk == 0)
			fid1 = fopen("hall_factor_up_spin.dat","w");
		if (ispin == 2 && kk == 1)
			fid1 = fopen("hall_factor_down_spin.dat","w");

		if (variation==0)   // temperature variation
			fprintf(fid1,"#Temperature(K)");
		else
		{
			if(geometry==1)
				fprintf(fid1,"#Doping(cm^-3)");
			else if(geometry==2)
				fprintf(fid1,"#Doping(cm^-2)");
		}

		fprintf(fid1,"  	hall_factor	hall_factor_rta \n");

		for (int i = 0; i <count ; i++)
		{
			fprintf(fid1," %e         %e              %e \n",
			hall_factor[i][0], hall_factor[i][1], hall_factor_rta[i][1]);
		}
		fclose(fid1);
	}
	//*/
}	

