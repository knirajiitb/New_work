#include"main.h"

int flag[50]={0};

double Efield_time[limit1];
double omega_s, initial, freq[30];
double J_time[limit1]={0}, sigma_time[limit1]={0}, mobility_time[limit1]={0};
int time_variation;
int time_limit, len_freq;
double g_time[limit2]={0},g_time_old[limit2]={0};
double omegam[30];
double J_freqr[30]={0}, sigma_freqr[30]={0}, Efield_freqr[30]={0}, J_freqi[30]={0}, sigma_freqi[30]={0}, Efield_freqi[30]={0};
double mobility_freqr2[30]={0}, mobility_freqi2[30]={0}, sigma_freqr2[30]={0}, sigma_freqi2[30]={0};
double mobility_drude_freqr2[30]={0}, mobility_drude_freqi2[30]={0}, sigma_drude_freqr2[30]={0}, sigma_drude_freqi2[30]={0};
int freq_variation=0;
double tau;

double nu_deformation_p[limit2][2][2]={0}, nu_ionizedimpurity_p[limit2][2][2]={0}, nu_el_p[limit2][2][2]={0};
double nu_npop_p[limit2][2][2]={0}, nu_So_p[limit2][2][2]={0};


double nu_deformation[limit2]={0}, nu_piezoelectric[limit2]={0}, nu_ionizedimpurity[limit2]={0}, nu_dislocation[limit2]={0}, nu_alloy[limit2]={0};
double nu_neutralimpurity[limit2]={0}, nu_npop[limit2][limit5]={0}, nu_npop_total[limit2]={0}, nu_iv[limit2][limit4]={0}, nu_iv_total[limit2]={0}, nu_el[limit2]={0};


string  type="n";

double n0, Nd1,Na1, efef_n, efef_p, N_ii;
int cc=-1, count_d, count_t, VASP=1;
double mobility_ii, mobility_po, mobility_to, mobility_de, mobility_pe, mobility_dis;
double mobility_alloy, mobility_iv, mobility_neutral, mobility_npop, mobility_avg, mobility, mobility_rta;
double mobility_hall_ii, mobility_hall_po, mobility_hall_to, mobility_hall_de;
double mobility_hall_pe, mobility_hall_dis, mobility_hall_alloy, mobility_hall_iv;
double mobility_hall_neutral, mobility_hall_npop, mobility_hall_avg, mobility_hall, mobility_hall_rta, hall_factor1, hall_factor_rta1;
double sigma_hall_rta, sigma_hall, thermopower, sigma, sigma_rta;

double denom[limit2];
int plus_index_pop[limit2], minus_index_pop[limit2];
double g[limit2], g_rta[limit2], g_old[limit2], g_LO[limit2], g_iv[limit2], g_th[limit2], g_th_old[limit2], g_LO_th[limit2];
double S_o_grid[limit2]={0}, S_o_grid_total[limit2]={0}, S_i_grid[limit2], S_iLO_grid[limit2], S_i_th_grid[limit2], S_iLO_th_grid[limit2];
double result_g[limit2][15+1], result_g_LO[limit2][15+1], result_g_th[limit2][15+1], result_f[limit2][15+1];

double N_poph_atT, df0dz_integral_n, N_e[limit4], beta_constant; 

double k_min, k_trans, k_step_fine, k_step;
int points, points1, points2;
double df0dk_grid[limit2], f0x1_f0[limit2], electric_driving_force[limit2], thermal_driving_force[limit2], f_dist[limit2];

double kplus_grid_pop[limit2], kminus_grid_pop[limit2], betaplus_grid[limit2], betaminus_grid[limit2];

double  Aminus_grid[limit2], Aplus_grid[limit2], lambda_i_plus_grid[limit2], lambda_o_plus_grid[limit2];
double lambda_i_minus_grid[limit2], lambda_o_minus_grid[limit2], lambda_e_plus_grid[limit2][limit4], lambda_e_minus_grid[limit2][limit4];
double lambda_e_plus_grid_npop[limit2][limit5], lambda_e_minus_grid_npop[limit2][limit5], N_npop[limit5];
int npop_number;

double Ed;
int a11[2],b11[2];
double h_bar,CBM,VBM;
int count1,count2;
int count_orbital, count_orbital_p;

double beta1[limit2], gH[limit2], hH[limit2], gH_rta[limit2], hH_rta[limit2];
double gH_LO[limit2], hH_LO[limit2], S_i_grid_g[limit2], S_i_grid_h[limit2];
double S_iLO_grid_g[limit2], S_iLO_grid_h[limit2], S_o_gridH[limit2]={0}, S_o_grid_totalH[limit2]={0};

double mobility_all[10]={0} , calc_mobility[30][2] = {0}, calc_mobility_rta[30][2] = {0};
double calc_thermopower[30][2] = {0}, calc_sigma[30][2] = {0}, calc_sigma_rta[30][2] = {0};

double calc_mobility_pe[30][1] = {0}, calc_mobility_de[30][1] = {0}, calc_mobility_dis[30][1] = {0}, calc_mobility_ii[30][1] = {0};
double calc_mobility_po[30][1] = {0}, calc_mobility_to[30][1] = {0}, calc_mobility_alloy[30][1] = {0}, calc_mobility_iv[30][1] = {0};
double calc_mobility_neutral[30][1] = {0}, calc_mobility_npop[30][1] = {0};

double mobility_hall_all[10]={0}, calc_mobility_hall[30][2] = {0}, calc_mobility_hall_rta[30][2] = {0};
double calc_sigma_hall[30][2] = {0}, calc_sigma_hall_rta[30][2] = {0}, hall_factor[30][2] = {0}, hall_factor_rta[30][2] = {0};

double calc_mobility_hall_pe[30][1] = {0}, calc_mobility_hall_de[30][1] = {0}, calc_mobility_hall_dis[30][1] = {0};
double calc_mobility_hall_ii[30][1] = {0}, calc_mobility_hall_iv[30][1] = {0}, calc_mobility_hall_neutral[30][1] = {0}, calc_mobility_hall_npop[30][1] = {0};
double calc_mobility_hall_po[30][1] = {0}, calc_mobility_hall_to[30][1] = {0}, calc_mobility_hall_alloy[30][1] = {0}; 
double kcbm[3],kvbm[3];

double k_grid[limit2]={0}, v_n[limit2]={0}, v_p[limit2]={0};	
double energy_n[limit2]={0}, energy_p[limit2]={0}, a_n[limit2]={0}, c_n[limit2]={0};

double a_p[limit2]={0}, c_p[limit2]={0}, Ds_n[limit2]={0}, Ds_p[limit2]={0};

double B_ii = 0, D_ii = 0, A_ii = 0;

double N_im_de, Nd_plus,N_im_modified;

int kk, ispin=1;

double T_array[30],epsilon_s[30],epsilon_inf[30],Bgap[30],P_piezo[30],C_piezo_h14[30],n_array[30],Nd[30],Na[30],N_im[30];

double fermi;
double we[limit4],De[limit4];
int nfv[limit4];
double we_npop[limit5],De_npop[limit5];


int degree1;
double fraction[4];
int length_fraction;

int De_ionization,N_cb, N_vb, iterations=10, variation, scattering_mechanisms[10], iv_number, fitting_1, fitting_2, fitting_3;

double rho=0, k_max, N_dis, omega_LO, omega_TO, E_deformation, C_long, C_trans, c_bar, C_11, C_12, C_44,
C_piezo_c14, P_piezo_h14, Uall, V0, xx, m ,m_h, T_trans ;

double Bfield=0;

int free_e=0;

//using namespace std;
 

void read_input_file()
{
	
    //char s[180];
    cout.precision(6);                        //set precision

        cout.setf(ios::scientific);
        cout<<"Data from input.dat file"<<endl;
        cout<< " -------------- "<<endl<<endl;
	
        //-------------------------------------------------------	
	int len_T,len_nn=0, len_na=0, len_nd=0;
	
        C_11=0;
        C_12=0;
        C_44=0;
      //-----------------------------------------------------
        for (int i=0;i<30;i++)
        {
		T_array[i]=0; epsilon_s[i]=0; epsilon_inf[i]=0; Bgap[i]=0; P_piezo[i] = 0;
		C_piezo_h14[i]=0; n_array[i]=0;  Nd[i]=0; Na[i]=0; N_im[i]=0;
        }

        for (int i=0;i<limit4;i++)
        {
	        we[i]=0; nfv[i]=0; De[i]=0;
        }

        De_ionization = 0;
        double Ed = 0.0092; // in eV   ionization energy  For ag 1 paper

        if (De_ionization==0)
	        Ed = 0; // in eV   ionization energy
        
        string str;
        string ss;
        
        ifstream in("input.dat");
        int count=0;

        while(!in.eof())
        {
		in>> str;
		
		// to check comments(if a line is not part of input then start tat line with #)
		if(str=="#")
		{
		  getline(in,ss);
		  cout<< "--COMMENT--" <<ss <<endl<<endl;
		}
		
		if(str=="TYPE")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>ss;
		  if(ss=="N")
		  	type="n";
		  else if(ss=="P")
		  	type="p";
		  else
		  {
		  	cout<<"Proper type is not given. Exit from program"<<endl;
		  	exit(EXIT_FAILURE);	
		  }	
		}
		
		if(str=="INPUT")
		{
			getline(in,ss);
			stringstream tmp(ss);


			tmp>>ss;
			if(ss=="VASP")
				VASP = 1;
			else if(ss=="TABLEFORM")
				VASP = 2;
			else
			{
				cout<<"Proper input are not given for VASP-INPUT. Exit from program"<<endl;
				exit(EXIT_FAILURE);	
			}			  
			//cout<< "VASP-INPUT  " <<VASP<<endl;
		}

		if(str=="TEMPERATURE")
		{
		  flag[0] = 1;  	
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0; 
		  while (tmp >> T_array[i]) 
		  { count++; i++; } 
		  
		  len_T = count;                    
		  // To display temperature array
		  cout<<"Temp array"<<endl;
		  for(int i=0;i<len_T;i++)
		  {
		    cout << T_array[i]<<"K    ";
		  }		
		  cout<<endl<<endl;
		  //cout<< "TEMPERATURE " <<ss<<endl;
		}
		
		if(str=="DONOR-DOPING")
		{
		  flag[1] = 1;	
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0; 
		  while (tmp >> Nd[i]) 
		  { count++; i++; } 
		  
		  len_nd = count;                    
		  // To display donor doping array
		  cout<<"donor doping array "<<endl;
		  for(int i=0;i<len_nd;i++)
		  {
		    cout << Nd[i]<<" per cm^3    ";
		  }		
		  cout<<endl<<endl;
		  //cout<< "DONOR DOPING CONCENTRATION " <<ss<<endl;
		}
		
		if(str=="ACCEPTOR-DOPING")
		{
		  flag[1] = 1;	
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0; 
		  while (tmp >> Na[i]) 
		  { count++; i++; } 
		  
		  len_na = count;                    
		  // To display acceptor doping array 
		  cout<<"acceptor doping array "<<endl;
		  for(int i=0;i<len_na;i++)
		  {
		    cout << Na[i]<<" per cm^3    ";
		  }		
		  cout<<endl<<endl;
		  //cout<< "ACCEPTOR DOPING CONCENTRATION " <<ss<<endl;
		}
		
		if(str=="NEUTRAL-IMPURITY")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  int count = 0;
		  int i=0; 
		  while (tmp >> N_im[i]) 
		  { count++; i++; } 
		  
		  len_nn = count;                    
		  // To display neutral impurity array 
		  cout<<"neutral impurity array "<<endl;
		  for(int i=0;i<len_nn;i++)
		  {
		    cout << N_im[i]<<" per cm^3    ";
		  }		
		  cout<<endl<<endl;
		  //cout<< "check saved"<< N_im[2]<<endl;
		  //cout<< "NEUTRAL IMPURITY CONCENTRATION " <<ss<<endl;
		}
		
		if(str=="DIELECTRIC-CONST-LF")
		{
		  flag[2] = 1;	
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> epsilon_s[0]; 
			  
		  // To display low frequency dielectric constant
		  cout<<"low frequency dielectric constant "<<endl;
		  cout <<epsilon_s[0];

		  cout<<endl<<endl;
		  //cout<< "LOW FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
		}
		
		if(str=="DIELECTRIC-CONST-HF")
		{
		  flag[3] = 1;	
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp >> epsilon_inf[0]; 
		  
		  // To display high frequency dielectric constant
		  cout<<"high frequency dielectric constant "<<endl;
		  cout<<epsilon_inf[0];
		  
		  cout<<endl<<endl;
		 // cout<< "HIGH FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
		}
		
		if(str=="BANDGAP")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>Bgap[0];

		  // To display band-gap
		  cout<<"BAND-GAP "<<endl;
		  cout<<Bgap[0];
				
		  cout<<endl<<endl;
		  
		  //cout<< "ENERGY BANDGAP " <<ss<<endl;
		}
		
		if(str=="N-VALENCE-BANDVALLEYS")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>N_vb;
		  cout<< "NUMBER OF VALENCE-BAND VALLIES IN BZ " <<N_vb<<endl;
		}
		
		if(str=="N-CONDUCTION-BANDVALLEYS")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>N_cb;
		  cout<< "NUMBER OF CONDUCTION-BAND VALLIES IN BZ " <<N_cb<<endl;
		}
		
		if(str=="DENSITYOFSEMICONDUCTOR")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> rho;
		  cout<< "DENSITY OF SEMICONDUCTOR " <<rho<< " gm/cm3 "<<endl;
    		  rho = rho*1000; 
		  //converted from g/cm^3 to kg/m^3
		}

		if(str=="DISLOCATIONDENSITY")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> N_dis;
		  cout<< "DISLOCATION DENSITY " <<N_dis<< " /CM^2 "<<endl;
		}

		if(str=="LPOPFREQUENCY")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>omega_LO;
		  cout<< "LONGITUDINAL POP FREQUENCY " <<omega_LO<< " THz " <<endl;
	          omega_LO = omega_LO*2*pi*1e12;
		}

		if(str=="ADP")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>E_deformation;		  
		  cout<< "ACOUSTIC DEFORMATION POTENTIAL " <<E_deformation<< " eV "<<endl;
		}

		if(str=="PZCOEFFICIENT")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>P_piezo[0];		  
		  cout<< "PIEZOELECTRIC COEFFICIENTS " <<P_piezo[0]<<endl;
		}

		if(str=="CL")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> C_long;
		  cout<< "ELASTIC CONSTANT FOR LONGITUDINAL MODE " <<C_long<<" dyne/cm^2 "<<endl;
		}

		if(str=="CT")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>C_trans;
		  cout<< "ELASTIC CONSTANT FOR TRANSVERSE MODE " <<C_trans<<" dyne/cm^2 "<<endl;
		}

		if(str=="ALLOYPOTENTIAL")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> Uall;
		  cout<< "ALLOY POTENTIAL " <<Uall<< "eV"<<endl;
		}

		if(str=="VOLUME-PRIMITIVECELL")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> V0;
		  cout<< "VOLUME OF PRIMITIVE CELL OF ALLOY  " <<V0<< " (nm)^3 "<<endl;
		}

		if(str=="FRACTIONOFATOM")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>xx;
		  cout<< "FRACTION OF ATOM FOR ALLOY " <<xx<<endl;
		}

		if(str=="INTRAVALLEYPHONONFREQUENCY")
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> we[i]) 
			{ count++; i++; } 
		
			if(iv_number!=0)
			{
				if(iv_number!=count)
				{
					cout<<"Error same number of intra valley constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				iv_number = count;

			cout<< "PHONON FREQUENCY FOR INTRAVALLEY SCATTERING "<<endl;
			for(int i=0;i<iv_number;i++)
			{
				cout<<"frequency["<<i<<"] = "<<we[i]<<"THz  "<<endl;
				we[i] = we[i]*2*pi*1e12;
			}
			cout<<endl<<endl;

		}

		if(str=="COUPLINGCONSTANT-INTRAVALLEYSCATTERING")
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> De[i]) 
			{ count++; i++; } 
		
			if(iv_number!=0)
			{
				if(iv_number!=count)
				{
					cout<<"Error same number of intra valley constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				iv_number = count;

			cout<< "COUPLING CONSTANT FOR INTRAVALLEY SCATTERING  "<<endl;
			for(int i=0;i<iv_number;i++)
			{
				cout<<"Coupling constant["<<i<<"] = "<<De[i]<<"e8 eV/cm  "<<endl;
			}
			cout<<endl;


		}

		if(str=="NUMBEROFVALLEY")
		{
			getline(in,ss);
			stringstream tmp(ss);
			int count = 0; 
			int i=0;
			while (tmp >> nfv[i]) 
			{  count++;   i++;   }

			if(iv_number!=0)
			{
				if(iv_number!=count)
				{
					cout<<"Error same number of intra valley constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				iv_number = count;

			for(int i=0;i<iv_number;i++)
			{
				cout<<"number of final valleys["<<i<<"] = "<<nfv[i]<<endl;
			}
			
		}

		if(str=="DOS")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>>free_e;
		  cout<< "-DENSITY OF STATES " <<free_e<<endl;
		  // free_e = 0, means DOSCAR is used, otherwise free electron desnity is used 
                 if (free_e == 0 )
			cout<<"free_e = false, DOSCAR is used for density of states"<<endl;
		  else
			cout<<"free_e = true, free electron density is used for calculation"<<endl;
		  
		}

		if(str=="SCATTERING-MECHANISM-CONTROL")
		{
		  flag[4] = 1;	
		  getline(in,ss);
		  stringstream smc(ss);
		  int a1[9];
		  for(int i=0;i<9;i++)
		    smc >> a1[i];

		  scattering_mechanisms[0] = a1[0];       // Ionized imourity 
		  scattering_mechanisms[1] = a1[1];     	// Polar Optical phonon scattering due to longitudinal phonon
		  scattering_mechanisms[2] = a1[8];     // npop phonon
		  scattering_mechanisms[3] = a1[2];	// Acoustic deformation scattering
		  scattering_mechanisms[4] = a1[3];	// Piezoelectric scattering          
		  
		  scattering_mechanisms[6] = a1[4];	// Dislocation scattering 
		  scattering_mechanisms[7] = a1[5];	// Alloy scattering
		  scattering_mechanisms[8] = a1[6];	// Intra-valley scattering
		  scattering_mechanisms[9] = a1[7];	// Neutral impurity scattering
		  cout<< "-SCATTERING MECHANISM CONTROL " <<ss<<endl;


		  //cout<< "alloy" << scattering_mechanisms[7]<<endl;
		}

		if(str=="NUMBER-OF-ITERATIONS")
		{
		  getline(in,ss);
		  stringstream tmp(ss);
		  tmp>> iterations;
		  cout<< "NUMBER OF ITERARTIONS " <<iterations<<endl;
		}

		if(str=="FITTING-DOP")
		{
		  getline(in,ss);
		  cout<< "DEGREE OF POLYNOMIAL USED FOR FITTING " <<ss<<endl;
		}
		
		if(str=="K-SEGMENT")
		{
			getline(in,ss);
			stringstream tmp(ss);
			fraction[0] = 0;
			fraction[1] = 0;
			fraction[2] = 0;
			fraction[3] = 0; 

			int count = 0;
			int i=0; 
			while (tmp >> fraction[i]) 
			{ count++; i++; } 

			cout<< "FRACTION array AT WHICH K-POINT IS DIVIDED "<<endl;

			for(int i=0;i<count;i++)
			{
			  cout << fraction[i]<<"    "<<endl;
			}

			length_fraction = count;
			//cout<<"Length of fraction array given by user is  "<<length_fraction<<endl;
			if (count==1 && fraction[0]==0)
			{
				//cout<<"Code will automatically calculate points for division of wavevector for minimum discontinuity"<<endl;
			}	
			else if ((count) > 4 || (count) < 3)
			{
				cout<<"Three or four points should be given for division of k segment. Exit from program "<<endl;
				exit (EXIT_FAILURE);		
			}		
			cout<<endl;
		}
	  
		if(str=="MAGNETIC-FIELD")
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>> Bfield;
			cout<< "MAGNETIC-FIELD " <<Bfield<<"  Tesla "<<endl;
		}
		
		
		if(str=="PHONON-FREQUENCY-NPOP")
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> we_npop[i]) 
			{ count++; i++; } 
		
			if(npop_number!=0)
			{
				if(npop_number!=count)
				{
					cout<<"Error same number of npop constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				npop_number = count;

			cout<< "PHONON-FREQUENCY for NPOP  "<<endl;
			for(int i=0;i<npop_number;i++)
			{
				cout<<"frequency["<<i<<"] = "<<we_npop[i]<<"THz  "<<endl;
				we_npop[i] = we_npop[i]*2*pi*1e12;
			}
			cout<<endl;
		}

		
		if(str=="COUPLING-CONSTANT-NPOP")
		{
			getline(in,ss);
			stringstream tmp(ss);
			int count = 0;
			int i=0; 
			while (tmp >> De_npop[i]) 
			{ count++; i++; } 
		
			if(npop_number!=0)
			{
				if(npop_number!=count)
				{
					cout<<"Error same number of npop constants are required"<<endl;
					exit (EXIT_FAILURE);			
				}
			}
			else
				npop_number = count;
				
			cout<< "COUPLING-CONSTANT-NPOP  "<<endl;
			for(int i=0;i<npop_number;i++)
				cout<<"Coupling constant["<<i<<"] = "<<De_npop[i]<<"e8 eV/cm  "<<endl;			
		}

				
		if(str=="TIME-VARIATION")
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>time_variation;
			cout<< " Time_variation  =   " <<time_variation<<endl;		  
		}

		if(str=="OMEGA-S")
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>omega_s;
			cout<< "OMEGA-S " <<omega_s<<endl;
		}
		
		if(str=="TIMESTEPS")
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>time_limit;
			cout<< "STEPS  for time variations   =  " <<time_limit<<endl;
		}

		if(str=="FREQ-VARIATION")
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>freq_variation;
			cout<< " Freq_variation  =   " <<freq_variation<<endl;		  
		}

		if(str=="OMEGA")
		{
			getline(in,ss);
			stringstream tmp(ss);

			int count = 0;
			int i=0; 
			while (tmp >> freq[i]) 
			{ count++; i++; } 
			
			len_freq = count;

			// To display omega array
			cout<<"Frequency array"<<endl;
			for(int i=0;i<len_freq;i++)
			{
			  cout << freq[i]<<" Hz   ";
			}		
			cout<<endl<<endl;
							
		}

		if(str=="INITIAL")
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>initial;
			cout<< " INITIAL time  =   " <<initial<<endl;		  
		}
		
		if(str=="TAU")
		{
			getline(in,ss);
			stringstream tmp(ss);
			tmp>>tau;
			cout<< " TAU =   " <<tau<<endl;		  
		}

        } // end of while loop 
	
	if(flag[0]==0)
	{
		cout<<"Temeprature is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	}
	else if(flag[1]==0)
	{
		cout<<"Doping is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	
	}
	else if(flag[2]==0)
	{
		cout<<"Low frequency dielectric constant is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	}
	else if(flag[3]==0)
	{
		cout<<"High frequency dielectric constant is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	}
	else if(flag[4]==0)
	{
		cout<<"Scattering Mechanism control is not given as input. Exit from program"<<endl;
		exit(EXIT_FAILURE);
	}
	
	
	if(len_nn==0 && len_na ==0)
	{
		len_nn = len_nd;
		len_na = len_nd;
	}
	else if(len_nn==0 && len_nd ==0)
	{
		len_nn = len_na;
		len_nd = len_na;
	}
	else if(len_nn==0)
	{
		len_nn = len_nd;
	}
	
	// To check array size of donor, neutral impurity and acceptor concentration
	if(len_nn!=len_na || len_nn!=len_nd)
	{
		cout<<"Error Donor array length, acceptor array length and neutral impurity conc. array length must be same. Exit from program "<<endl;
		exit (EXIT_FAILURE);
	}

	// to fix ntype or ptype
	if(type =="n")
	{
		cout<<"Net doping array"<<endl;
		for(int j=0;j<len_nd;j++)
		{
			
			n_array[j] = Nd[j] - Na[j];
			cout<<n_array[j]<<"   "<<"cm^(-3)";

		}
		cout<<endl;		
	}
	else
	{
		cout<<"Net doping array"<<endl;
		for(int j=0;j<len_nd;j++)
		{
			
			n_array[j] = Na[j] - Nd[j];
			cout<<n_array[j]<<"   "<<"cm^(-3)";

		}
		cout<<endl;		
	}		
        
	// to fix value of variation
	if (len_nd==1 & len_T==1)
		variation = 0;    // Any variation can be taken temperature variation
	else if (len_nd==1)
		variation = 0;   // temp variation
	else if (len_T==1)
		variation=1;   // Doping variation
	else
	{	
		cout<<"Error, out of doping and temperature array one should be one length exit from the program"<<endl;
		exit (EXIT_FAILURE);
	}	


	count_d = len_nd;
	count_t = len_T;


	N_cb = 1;
	N_vb = 1;


	// copy bandgap dielectric constant in all elements of array	
	if (len_T > 1)
 	{
		for (int i=1;i<len_T;i++)
		{
			epsilon_s[i]=epsilon_s[0]; 
			epsilon_inf[i]=epsilon_inf[0]; 
			Bgap[i]=Bgap[0]; 
		}
	}

	if(VASP!=1)
	{
		if(rho==0)
		{
			cout<<"Density value is required.  "<<endl;
			cout<<"Exit from program"<<endl;
			exit(EXIT_FAILURE);
		}
	}	

	for(int i=0;i<len_T;i++)
	C_piezo_h14[i] = C_piezo_c14/(epsilon_s[i]*epsilon_0);   // unit - N/C

	if (C_trans == 0)
	{
		C_trans =  (C_11 - C_12 + 3 * C_44)/5;  // in dyne/cm2   Equation 99 from rode book
		printf("\n C_trans =   %e dyne/cm^2 \n", C_trans);
	}
	
	if (C_long == 0)
	{
		C_long =  (3*C_11 + 2 * C_12 + 4 * C_44)/5;         // in dyne/cm2   Equation 100 from rode book
		printf("\n C_long  =   %e dyne/cm^2 \n", C_long);
	}
	
	
	c_bar = (1.0/3.0)*C_long + (2.0/3)*C_trans;   // in dyne/cm^2
	printf("\n c_bar   =   %e dyne/cm^2 \n", c_bar);

	if (P_piezo[0] == 0)
	{
		for(int i=0;i<len_T;i++)
		{
		    P_piezo[i] = (pow(C_piezo_h14[i],2)*epsilon_0*epsilon_s[i]*(12/C_long+16/C_trans)/35*1e1 );
		    P_piezo[i] = pow(P_piezo[i],0.5);
		// P_piezo is unitless and 1e1 is to convert 1/(dyn/cm2) to 1/(N/m2)
		}
		printf("\n P_piezo =   %e \n",P_piezo[0]);
		
	}
	else
	{
		if (len_T>1)		
		for(int i=1;i<len_T;i++)
		{
		    P_piezo[i] = P_piezo[0];
		// P_piezo is unitless and 1e1 is to convert 1/(dyn/cm2) to 1/(N/m2)
		}
		
	}

	k_max = 6;
	T_trans = 40;
	
	//cout<<"End of Reading input.dat file"<<endl;

	//printf("\n P_piezo =   %e \n",P_piezo);

	if (scattering_mechanisms[1] == 0)
	{
		if (iterations != 1 )
		{
		    iterations = 1;
		    cout<<"Since polar optical phonon scattering is not used, so iterations is set to 1"<<endl;
		}
		T_trans = 0;
	}

	if (De_ionization == 1)
	{	// neutral impurity scattering
		if (scattering_mechanisms[9] == 0)
		{
		    scattering_mechanisms[9] = 1;
		    cout<<"De-ionization is to be included so neutral impurity scattering is to included in simulation";
		}
	}

	h_bar =  h_planck / (2 * pi);    // Planck constant divided by 2pi[eV.s]

	fitting_1 = 0;  // FOR BAND
	fitting_2 = 0;  // FOR PROCAR
	fitting_3 = 1;  // FOR DOSCAR
  
}
