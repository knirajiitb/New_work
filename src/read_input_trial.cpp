/*#include<iostream>
#include <bits/stdc++.h>
#include<fstream>
#include<cstring>
#include<cstdio>*/
#include"main.h"

double n0, Nd1,Na1, efef_n, efef_p, N_ii;
int cc=-1, count_d, count_t;
double mobility_ii, mobility_po, mobility_to, mobility_npop, mobility_de, mobility_pe, mobility_dis;
double mobility_alloy, mobility_iv, mobility_neutral, mobility_avg, mobility, mobility_rta;
double mobility_hall_ii, mobility_hall_po, mobility_hall_to, mobility_hall_npop, mobility_hall_de;
double mobility_hall_pe, mobility_hall_dis, mobility_hall_alloy, mobility_hall_iv;
double mobility_hall_neutral, mobility_hall_avg, mobility_hall, mobility_hall_rta, hall_factor1, hall_factor_rta1;
double sigma_hall_rta, sigma_hall, thermopower, sigma, sigma_rta;

double g[2000], g_rta[2000], g_old[2000], g_LO[2000], g_iv[2000], g_th[2000], g_th_old[2000], g_LO_th[2000];
double S_o_grid[2000], S_o_grid_total[2000], S_i_grid[2000], S_iLO_grid[2000], S_i_th_grid[2000], S_iLO_th_grid[2000];
double result_g[2000][15+1], result_g_LO[2000][15+1], result_g_th[2000][15+1], result_f[2000][15+1];

double N_poph_atT, df0dz_integral_n, N_e[10], beta_constant; 

double k_min, k_trans, k_step_fine, k_step;
int points, points1, points2;
double df0dk_grid[2000], f0x1_f0[2000], electric_driving_force[2000], thermal_driving_force[2000], f_dist[2000];

double kplus_grid[2000], kminus_grid[2000], betaplus_grid[2000], betaminus_grid[2000];

double  Aminus_grid[2000], Aplus_grid[2000], lambda_i_plus_grid[2000], lambda_o_plus_grid[2000];
double lambda_i_minus_grid[2000], lambda_o_minus_grid[2000], lambda_e_plus_grid[2000][5], lambda_e_minus_grid[2000][5];

double nu_deformation[2000], nu_piezoelectric[2000], nu_ionizedimpurity[2000], nu_dislocation[2000], nu_alloy[2000];
double nu_neutralimpurity[2000], nu_iv[2000][5], nu_iv_total[2000], nu_el[2000];

double Ed;
int a11[2],b11[2];
double h_bar,CBM,VBM;
int count1,count2;
int count_orbital, count_orbital_p;

double beta1[2000], gH[2000], hH[2000], gH_rta[2000], hH_rta[2000];
double gH_LO[2000], hH_LO[2000], S_i_grid_g[2000], S_i_grid_h[2000];
double S_iLO_grid_g[2000], S_iLO_grid_h[2000], S_o_gridH[2000], S_o_grid_totalH[2000];

double mobility_all[10]={0} , calc_mobility[30][2] = {0}, calc_mobility_rta[30][2] = {0};
double calc_thermopower[30][2] = {0}, calc_sigma[30][2] = {0}, calc_sigma_rta[30][2] = {0};

double calc_mobility_pe[30][1] = {0}, calc_mobility_de[30][1] = {0}, calc_mobility_dis[30][1] = {0}, calc_mobility_ii[30][1] = {0};
double calc_mobility_po[30][1] = {0}, calc_mobility_to[30][1] = {0}, calc_mobility_alloy[30][1] = {0}, calc_mobility_iv[30][1] = {0};
double calc_mobility_neutral[30][1] = {0};

double mobility_hall_all[10]={0}, calc_mobility_hall[30][2] = {0}, calc_mobility_hall_rta[30][2] = {0};
double calc_sigma_hall[30][2] = {0}, calc_sigma_hall_rta[30][2] = {0}, hall_factor[30][2] = {0}, hall_factor_rta[30][2] = {0};

double calc_mobility_hall_pe[30][1] = {0}, calc_mobility_hall_de[30][1] = {0}, calc_mobility_hall_dis[30][1] = {0};
double calc_mobility_hall_ii[30][1] = {0}, calc_mobility_hall_iv[30][1] = {0}, calc_mobility_hall_neutral[30][1] = {0};
double calc_mobility_hall_po[30][1] = {0}, calc_mobility_hall_to[30][1] = {0}, calc_mobility_hall_alloy[30][1] = {0}; 
double kcbm[3],kvbm[3];

double k_grid[2000]={0}, v_n[2000]={0}, v_p[2000]={0};	
double energy_n[2000]={0}, energy_p[2000]={0}, a_n[2000]={0}, c_n[2000]={0};

double a_p[2000]={0}, c_p[2000]={0}, Ds_n[2000]={0}, Ds_p[2000]={0};

double B_ii = 0, D_ii = 0, A_ii = 0;

double N_im_de, Nd_plus,N_im_modified;

int kk, ispin;

double T_array[30],epsilon_s[30],epsilon_inf[30],Bgap[30],P_piezo[30],C_piezo_h14[30],n_array[30],Nd[30],Na[30],N_im[30];

double we[10],De[10];
int nfv[10];

int degree1;
double fraction[4];
int length_fraction;

int De_ionization,N_cb, N_vb, iterations, variation, scattering_mechanisms[10], iv_number, fitting_1, fitting_2, fitting_3;

double rho, k_max, N_dis, omega_LO, omega_TO, E_deformation_n, C_long, C_trans, c_bar, C_11, C_12, C_44,
C_piezo_c14, P_piezo_h14, Uall, V0, xx, m ,m_h, T_trans ;

double Bfield;

int free_e,len_T,len_n;


//using namespace std;
 

int main()
{
    //char s[180];
    cout.precision(6);                        //set precision
        cout.setf(ios::scientific);
        cout<<"Data from input.dat file"<<endl;
        cout<< " -------------- "<<endl<<endl;

        //-------------------------------------------------------	
        C_11=0;
        C_12=0;
        C_44=0;
      //-----------------------------------------------------
        for (int i=0;i<30;i++)
        {
        T_array[i]=0; epsilon_s[i]=0; epsilon_inf[i]=0; Bgap[i]=0; P_piezo[i] = 0;
        C_piezo_h14[i]=0; n_array[i]=0;  Nd[i]=0; Na[i]=0; N_im[i]=0;
        }

        for (int i=0;i<10;i++)
        {
        we[i]=0; nfv[i]=0; De[i]=0;
        }

        De_ionization = 0;
        double Ed = 0.0092; // in eV   ionization energy  For ag 1 paper

        if (De_ionization==0)
        Ed = 0; // in eV   ionization energy
        //---------
        int type_value; //if type_value is 1 then n-type,0 then p-type
        string str;
        string ss;
        
        ifstream in("ammcrinput_TRIAL.dat");
        int count=0;

        while(!in.eof()){
          in>> str;
          // to check comments(if a line is not part of input then start tat line with #)
        if(str=="#"){
          getline(in,ss);
          cout<< "--COMMENT--" <<ss <<endl<<endl;
        }
        
        if(str=="TYPE"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> type_value;
          cout<< "TYPE_VALUE " <<ss<<endl;
          //cout<< "saved " <<type_value<<endl;
        }
        if(str=="TEMPERATURE"){
          getline(in,ss);
          stringstream tmp(ss);
          int count = 0,count1;int i=0; 
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
        
        if(str=="DONOR-DOPING"){
          getline(in,ss);
          stringstream tmp(ss);
          int count = 0,count1;int i=0; 
          while (tmp >> Nd[i]) 
          { count++; i++; } 
          
          len_n = count;                    
          // To display donor doping array
          cout<<"donor doping array "<<endl;
          for(int i=0;i<len_n;i++)
          {
            cout << Nd[i]<<" per cm^3    ";
          }		
          cout<<endl<<endl;
          //cout<< "DONOR DOPING CONCENTRATION " <<ss<<endl;
        }
        if(str=="ACCEPTOR-DOPING"){
          getline(in,ss);
          stringstream tmp(ss);
          int count = 0,count1;int i=0; 
          while (tmp >> Na[i]) 
          { count++; i++; } 
          
          len_n = count;                    
          // To display acceptor doping array 
          cout<<"acceptor doping array "<<endl;
          for(int i=0;i<len_n;i++)
          {
            cout << Na[i]<<" per cm^3    ";
          }		
          cout<<endl<<endl;
          //cout<< "ACCEPTOR DOPING CONCENTRATION " <<ss<<endl;
        }
        
        if(str=="NEUTRAL-IMPURITY"){
          getline(in,ss);
          stringstream tmp(ss);
          int count = 0,count1;int i=0; 
          while (tmp >> N_im[i]) 
          { count++; i++; } 
          
          len_n = count;                    
          // To display neutral impurity array 
          cout<<"neutral impurity array "<<endl;
          for(int i=0;i<len_n;i++)
          {
            cout << N_im[i]<<" per cm^3    ";
          }		
          cout<<endl<<endl;
          //cout<< "check saved"<< N_im[2]<<endl;
          //cout<< "NEUTRAL IMPURITY CONCENTRATION " <<ss<<endl;
        }
        if(str=="DIELECTRIC-CONST-LF"){
          getline(in,ss);
          stringstream tmp(ss);
          int count = 0,count1;int i=0; 
          while (tmp >> epsilon_s[i]) 
          { count++; i++; } 
          
          len_T = count;                    
          // To display low frequency dielectric constant
          cout<<"low frequency dielectric constant "<<endl;
          for(int i=0;i<len_T;i++)
          {
            cout << epsilon_s[i]<<"   ";
          }		
          cout<<endl<<endl;
          //cout<< "LOW FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
        }
        if(str=="DIELECTRIC-CONST-HF"){
          getline(in,ss);
          stringstream tmp(ss);
          int count = 0,count1;int i=0; 
          while (tmp >> epsilon_inf[i]) 
          { count++; i++; } 
          
          len_T = count;                    
          // To display high frequency dielectric constant
          cout<<"high frequency dielectric constant "<<endl;
          for(int i=0;i<len_T;i++)
          {
            cout << epsilon_inf[i]<<"   ";
          }		
          cout<<endl<<endl;
         // cout<< "HIGH FREQUENCY DIELECTRIC CONSTANT " <<ss<<endl;
        }
        if(str=="BANDGAP"){
          getline(in,ss);
          stringstream tmp(ss);
          int count = 0,count1;int i=0; 
          while (tmp >> Bgap[i]) 
          { count++; i++; } 
          
          len_T = count;                    
          // To display band-gap
          cout<<"BAND-GAP "<<endl;
          for(int i=0;i<len_T;i++)
          {
            cout << Bgap[i]<<"   ";
          }		
          cout<<endl<<endl;
          //cout<< "ENERGY BANDGAP " <<ss<<endl;
        }
        if(str=="N-VALENCE-BANDVALLEYS"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> N_vb;
          cout<< "NUMBER OF VALENCE-BAND VALLIES IN BZ " <<N_vb<<endl;
        }
        if(str=="N-CONDUCTION-BANDVALLEYS"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> N_cb;
          cout<< "NUMBER OF CONDUCTION-BAND VALLIES IN BZ " <<N_cb<<endl;
        }
        if(str=="DENSITYOFSEMICONDUCTOR"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> rho;
          cout<< "DENSITY OF SEMICONDUCTOR " <<rho<< " gm/cm3 "<<endl;
        }
        if(str=="DISLOCATIONDENSITY"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> N_dis;
          cout<< "DISLOCATION DENSITY " <<N_dis<< " /CM^2 "<<endl;
        }
        if(str=="LPOPFREQUENCY"){
          getline(in,ss);
          cout<< "LONGITUDINAL POP FREQUENCY " <<ss<< " THz " <<endl;
        }
        if(str=="ADP"){
          getline(in,ss);
          cout<< "ACOUSTIC DEFORMATION POTENTIAL " <<ss<< " eV "<<endl;
        }
        if(str=="PZCOEFFICIENT"){
          getline(in,ss);
          cout<< "PIEZOELECTRIC COEFFICIENTS " <<ss<<endl;
        }
        if(str=="CL"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> C_long;
          cout<< "ELASTIC CONSTANT FOR LONGITUDINAL MODE " <<C_long<<" dyne/cm^2 "<<endl;
        }
        if(str=="CT"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> C_trans;
          cout<< "ELASTIC CONSTANT FOR TRANSVERSE MODE " <<C_trans<<" dyne/cm^2 "<<endl;
        }
        if(str=="ALLOYPOTENTIAL"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> Uall;
          cout<< "ALLOY POTENTIAL " <<Uall<< "eV"<<endl;
        }
        if(str=="VOLUME-PRIMITIVECELL"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> V0;
          cout<< "VOLUME OF PRIMITIVE CELL OF ALLOY  " <<V0<< " (nm)^3 "<<endl;
        }
        if(str=="FRACTIONOFATOM"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> xx;
          cout<< "FRACTION OF ATOM FOR ALLOY " <<xx<<endl;
        }
        if(str=="INTRAVALLEYPHONONFREQUENCY"){
          getline(in,ss);
          cout<< "PHONON FREQUENCY FOR INTRAVALLEY SCATTERING " <<ss<<endl;
        }
        if(str=="COUPLINGCONSTANT-INTRAVALLEYSCATTERING"){
          getline(in,ss);
          cout<< "COUPLING CONSTANT FOR INTRAVALLEY SCATTERING  " <<ss<<endl;
        }
        if(str=="NUMBEROFVALLEY"){
          getline(in,ss);
          stringstream tmp(ss);
          int count3 = 0; int i=0;
          while (tmp >> nfv[i]) 
          {  count3++;   i++;   }

          iv_number = count3; 
          for(int i=0;i<iv_number;i++)
          {
            cout<<"number of final valleys["<<i<<"] = "<<nfv[i]<<endl;
          }
          //cout<< "NUMBER OF FINAL VALLEY FOR SCATTERING " <<ss<<endl;
        }
        if(str=="DOS"){
          getline(in,ss);
          cout<< "-DENSITY OF STATES " <<ss<<endl;
        }
        if(str=="SCATTERING-MECHANISM-CONTROL"){
          getline(in,ss);
          stringstream smc(ss);
          int a1[9];
          for(int i=0;i<9;i++)
            smc >> a1[i];

          scattering_mechanisms[0] = a1[0];       // Ionized imourity 
          scattering_mechanisms[1] = a1[1];     	// Polar Optical phonon scattering due to longitudinal phonon
          scattering_mechanisms[3] = a1[2];	// Acoustic deformation scattering
          scattering_mechanisms[4] = a1[3];	// Piezoelectric scattering          
          scattering_mechanisms[5] = a1[4];	// Dislocation scattering 
          scattering_mechanisms[6] = a1[5];	// Alloy scattering
          scattering_mechanisms[7] = a1[6];	// Intra-valley scattering
          scattering_mechanisms[8] = a1[7];	// Neutral impurity scattering
          scattering_mechanisms[9] = a1[8];     // TO npop phonon
          cout<< "-SCATTERING MECHANISM CONTROL " <<ss<<endl;
          //cout<< "alloy" << scattering_mechanisms[7]<<endl;
        }
        if(str=="NUMBER-OF-ITERATIONS"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> iterations;
          cout<< "NUMBER OF ITERARTIONS " <<iterations<<endl;
        }
        if(str=="FITTING-DOP"){
          getline(in,ss);
          cout<< "DEGREE OF POLYNOMIAL USED FOR FITTING " <<ss<<endl;
        }
        if(str=="K-SEGMENT"){
          getline(in,ss);
          cout<< "FRACTION AT WHICH K-POINT IS DIVIDED " <<ss<<endl;
        }
        if(str=="MAGNETIC-FIELD"){
          getline(in,ss);
          stringstream tmp(ss);
          tmp>> Bfield;
          cout<< "MAGNETIC-FIELD " <<Bfield<<endl;
        }
        if(str=="OMEGA-S"){
          getline(in,ss);
          cout<< "OMEGA-S " <<ss<<endl;
        }
        if(str=="STEPS"){
          getline(in,ss);
          cout<< "STEPS  " <<ss<<endl;
        }
        if(str=="PHONON-FREQUENCY-NPOP"){
          getline(in,ss);
          cout<< "PHONON-FREQUENCY-NPOP  " <<ss<<endl;
        }
        if(str=="COUPLING-CONSTANT-NPOP"){
          getline(in,ss);
          cout<< "COUPLING-CONSTANT-NPOP  " <<ss<<endl;
        }
        if(str=="VASP-INPUT"){
          getline(in,ss);
          cout<< "VASP-INPUT  " <<ss<<endl;
        }


        
          
        }
        

        
        
}
