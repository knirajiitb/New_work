#include"main.h"

void generate_output_files()
{

    if (variation==0)   // temperature variation
    {
    	if (ispin == 1)
    	{
		FILE *fid1;
		fid1 = fopen("mobility.dat","w");
		if(geometry==1)
		{
		if(type=="n")
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		else
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii      Mobility_po      Mobility_npop	 Mobility_de    \n");
		}
		else  // geometry==2
		{
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     mobility_so_pop \n");
		}
		fclose(fid1);

		fid1 = fopen("conductivity.dat","w");
		fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);

		if(geometry==1)
		{
		fid1 = fopen("thermopower.dat","w");
		fprintf(fid1,"# Temperature(K)     Thermopower(uV/K)\n");
		fclose(fid1);
		}		
		if(Bfield!=0)
		{
			fid1 = fopen("mobility_hall.dat","w");
			if(geometry==1)
			{
			if(type=="n")
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_remote_ii      Mobility_po       Mobility_npop	Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			
			else   // type==p
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_remote_ii      Mobility_po       Mobility_npop	Mobility_de  \n");
			}
			else  // geometry==2
			{
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_remote_ii      Mobility_po       Mobility_npop	Mobility_de     Mobility_pe     Mobility_so_pop \n");
			}
			fclose(fid1);

			fid1 = fopen("conductivity_hall.dat","w");
			fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor.dat","w");
			fprintf(fid1,"# Temperature(K)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);
		}
	}
	
    	if (ispin == 2 )
    	{
		FILE *fid1;
		fid1 = fopen("mobility_up_spin.dat","w");
		if(geometry==1)
		{
		if(type=="n")
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po       Mobility_npop	Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral\n");
		
		else  // type==p
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po       Mobility_npop	Mobility_de   \n");
		}
		else  // geometry==2
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     mobility_so_pop \n");
				
		fclose(fid1);

		fid1 = fopen("conductivity_up_spin.dat","w");
		fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);

		if(geometry==1)
		{
		fid1 = fopen("thermopower_up_spin.dat","w");
		fprintf(fid1,"# Temperature(K)     Thermopower(uV/K)\n");
		fclose(fid1);
		}
		
		fid1 = fopen("mobility_down_spin.dat","w");
		if(geometry==1)
		{
		if(type=="n")
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	 Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		else //type=="p"
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	 Mobility_de \n");
		}
		else  // geometry==2
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     mobility_so_pop \n");
		
		fclose(fid1);

		fid1 = fopen("conductivity_down_spin.dat","w");
		fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);

		if(geometry==1)
		{
		fid1 = fopen("thermopower_down_spin.dat","w");
		fprintf(fid1,"# Temperature(K)     Thermopower(uV/K)\n");
		fclose(fid1);
		}		
		
		if(Bfield!=0)
		{
			fid1 = fopen("mobility_hall_up_spin.dat","w");
			if(geometry==1)
			{
			if(type=="n")
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			else
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de  \n");
			}
			else  //geometry==2
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     mobility_so_pop \n");
			
			fclose(fid1);

			fid1 = fopen("conductivity_hall_up_spin.dat","w");
			fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor_up_spin.dat","w");
			fprintf(fid1,"# Temperature(K)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);


			fid1 = fopen("mobility_hall_down_spin.dat","w");
			if(geometry==1)
			{
			if(type=="n")
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	 Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			else  // type==p
			fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	 Mobility_de   \n");
			}
			else  // geometry==2
		fprintf(fid1,"#Temperature(K)  Mobility(cm^2/V-s)    Mobility_rta          Mobility_remote_ii      Mobility_po      Mobility_npop	 Mobility_de     Mobility_pe     mobility_so_pop \n");
			
			fclose(fid1);

			fid1 = fopen("conductivity_hall_down_spin.dat","w");
			fprintf(fid1,"# Temperature(K)      Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor_down_spin.dat","w");
			fprintf(fid1,"# Temperature(K)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);
		}		
	}

    }
    else          // Doping variation
    {
	if (ispin == 1)
	{
		FILE *fid1;

		fid1 = fopen("mobility.dat","w");
		if(geometry==1)
		{
		if(type=="n")
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		else  // type p
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de    \n");
		}
		else  // geometry==2
		
		fclose(fid1);

		fid1 = fopen("conductivity.dat","w");
		fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);
		
		if(geometry==1)
		{
		fid1 = fopen("thermopower.dat","w");
		fprintf(fid1,"# Doping(cm^-3)      Thermopower(uV/K)\n");
		fclose(fid1);
		}
		if(Bfield!=0)
		{

			fid1 = fopen("mobility_hall.dat","w");
			if(geometry==1)
			{
			if(type=="n")
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			else
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de    \n");
			}
			else // geometry==2
			fprintf(fid1,"#Doping(cm^-3) Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_remote_ii      Mobility_po       Mobility_npop	Mobility_de     Mobility_pe     Mobility_so_pop \n");
			
			fclose(fid1);

			fid1 = fopen("conductivity_hall.dat","w");
			fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor.dat","w");
			fprintf(fid1,"# Doping(cm^-3)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);
		}
	}


	if (ispin == 2)
	{
		FILE *fid1;

		fid1 = fopen("mobility_up_spin.dat","w");
		if(geometry==1)
		{
		if(type=="n")
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
		else
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de     \n");
		}
		else // geometry==2
		fprintf(fid1,"#Doping(cm^-3) Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_remote_ii      Mobility_po       Mobility_npop	Mobility_de     Mobility_pe     Mobility_so_pop \n");
		
		fclose(fid1);

		fid1 = fopen("conductivity_up_spin.dat","w");
		fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);
		
		if(geometry==1)
		{
		fid1 = fopen("thermopower_up_spin.dat","w");
		fprintf(fid1,"# Doping(cm^-3)      Thermopower(uV/K)\n");
		fclose(fid1);
		}
		fid1 = fopen("mobility_down_spin.dat","w");
		if(geometry==1)
		{
		if(type=="n")
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral  \n");
		else
		fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	Mobility_de      \n");
		fclose(fid1);
		}
		else  // geometry==2
		fprintf(fid1,"#Doping(cm^-3) Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_remote_ii      Mobility_po       Mobility_npop	Mobility_de     Mobility_pe     Mobility_so_pop \n");
		
		fid1 = fopen("conductivity_down_spin.dat","w");
		fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
		fclose(fid1);
		
		if(geometry==1)
		{
		fid1 = fopen("thermopower_down_spin.dat","w");
		fprintf(fid1,"# Doping(cm^-3)      Thermopower(uV/K)\n");
		fclose(fid1);
		}
		if(Bfield!=0)
		{
				
			fid1 = fopen("mobility_hall_up_spin.dat","w");
			if(geometry==1)
			{
			if(type=="n")
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	   Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral  \n");
			else
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	   Mobility_de     \n");
			}
			else // geometry==2
			fprintf(fid1,"#Doping(cm^-3) Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_remote_ii      Mobility_po       Mobility_npop	Mobility_de     Mobility_pe     Mobility_so_pop \n");
			
			fclose(fid1);

			fid1 = fopen("conductivity_hall_up_spin.dat","w");
			fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);
			
			fid1 = fopen("hall_factor_up_spin.dat","w");
			fprintf(fid1,"# Doping(cm^-3)      Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);

			
			fid1 = fopen("mobility_hall_down_spin.dat","w");
			if(geometry==1)
			{
			if(type=="n")
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	   Mobility_de     Mobility_pe     Mobility_dis    Mobility_to     Mobility_alloy  Mobility_iv   Mobility_neutral \n");
			else
			fprintf(fid1,"#Doping(cm^-3)   Mobility(cm^2/V-s)    Mobility_rta          Mobility_ii      Mobility_po     Mobility_npop	   Mobility_de     \n");
			}
			else  // geometry==2
			fprintf(fid1,"#Doping(cm^-3) Mobility(cm^2/V-s)    Mobility_hall_rta          Mobility_remote_ii      Mobility_po       Mobility_npop	Mobility_de     Mobility_pe     Mobility_so_pop \n");
			
			fclose(fid1);

			fid1 = fopen("conductivity_hall_down_spin.dat","w");
			fprintf(fid1,"# Doping(cm^-3)       Conductivity(S/cm)(Rode)  Conductivity(S/cm)(RTA)\n");
			fclose(fid1);

			fid1 = fopen("hall_factor_down_spin.dat","w");
			fprintf(fid1,"# Doping(cm^-3)     Hall_factor(Rode)  Hall_factor(RTA)\n");
			fclose(fid1);
		}
	}

    }

}
