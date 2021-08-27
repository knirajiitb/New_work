
#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include <cstring>
#include<iostream>
using namespace std;

#include<fstream>

#include<stdio.h>
#include <math.h>
#include <time.h>
#include<stdlib.h>
#include <vector>
#include <sys/stat.h>   // this for FileExists function
#include<bits/stdc++.h>
#include <algorithm>


//#include <opencv/cv.h>

//#include <Python.h>
//#include "pyhelper.hpp"



//using namespace std;
//using namespace cv;

//------------ - *** Some Global constants ***----------------------------------------------------
#define m_e 9.10938291e-31   // Electron mass[kg]
#define e 1.60217657e-19     // Electron charge[C]
#define pi 3.14159265359      // Pi number
#define h_planck 4.135667516e-15    // Planck constant[eV.s]
#define k_B 8.6173324e-5       // Boltzmann constant[eV / K]
#define epsilon_0 8.854187817e-12    // Absolute value of dielectric constant in vacuum[C ^ 2 / m ^ 2N]
//#define h_bar

#define k_min0 0.001
#define k_trans0 0.1
#define k_step_fine0  0.001
#define k_step0 0.01
#define E 1    // Electric field unit cgs V/cm
#define dTdz 10
#define LORBIT 11
#define limit1 6000  // size for time variation
#define limit2 2000  // size for kpoints after fitting
#define limit3 20000  // size of kpoints for reading
#define limit4 5    // for intervalley scattering no. limit 
#define limit5 5    // for npop scattering no. limit 
//------------ - ------------------------------------------------------------------------------

//-----------function for Valence band --------
void nu_ii_p_funct(int T_loop);
void nu_So_p_funct(double T, int T_loop,double omega_LO);  
void nu_npop_p_funct(double T,int T_loop,int npop_number);
void nu_de_p_funct(int T_loop);
//---------------------------

void read_OUTCAR();
void read_input_file();
void copyright(); 
void generate_output_files();

void find_cbm_vbm(int , int);
void calculate_mobility(double T, int ii);

void make_band(int typee);
int* analytical_fitting(double band[limit2][2],int counta, int aa);
double * polyfit1(double x[],double y[], int n, int N);   // n-degree   N - length of x[] and y[] array
void polyval(double p[], double x[], int nn, int mm);
double * linspace(double a, double b, int numbers);
double rsquare(double data[], double fitted[], int length);
int decompose_wave();
int decompose_wave_p();
void n_DOS();
int FindMinInd(double arr[],int length);

double conduction_dispersion(double k,double coefficients[5][7],double kindex[4],int aa[2]);
double admixture_value(double k,int coloumn);
double admixture_value_p(double k,int coloumn);
double DOS_value(double energy, int a);

void find_fermi(double n, double T, int ii);
double f0(double E1, double e_f, double T);
double beta(double T,int T_loop);
double N_poph(double omega, double T);
double dedk(double k, double coeff[5][7], double kindex[4], int aa[2]);

double kplus(int counter, double omega, int points, double energy[]);
double kminus(int counter, double omega, int points, double energy[]);

double Aplus(int counter, double omega, int points, double energy[]);
double Aminus(int counter, double omega, int points, double energy[]);

double betaplus(int counter, double omega, double epsilon_s, double epsilon_inf , int points);
double betaminus(int counter, double omega, double epsilon_s, double epsilon_inf, int points);

double lambda_i_plus(int counter, double omega,double A_plus,double epsilon_s, double epsilon_inf, int points);
double lambda_i_minus(int counter,double omega,double A_minus,double epsilon_s, double epsilon_inf,int points);

double lambda_o_plus(int counter,double omega,double A_plus, double epsilon_s,double epsilon_inf,int points);
double lambda_o_minus(int counter,double omega,double A_minus, double epsilon_s,double epsilon_inf,int points);

double lambda_e_minus(int counter,double omega,double rho,double De,int nfv,int points);
double lambda_e_plus(int counter,double omega,double rho,double De,int nfv,int points);

double nu_ii(double k, int counter, double beta_constant, double v, double epsilon_s);
double nu_de(double k,int counter,double T,double v);                    
double nu_pe(double k,int counter,double T,double P_piezo,double epsilon_s,double v);
double nu_dis(double k, int counter, double T, double beta_constant, double epsilon_s, double v);
double nu_alloy1(double k,double v);
double nu_im(double k_dum, int counter, double epsilon_s, double N_im, double v);

void generate_required_data(double T);

double df0dk(double k,double T,double e_f, double coefficients[5][7],double kindex[], int aa[]);
double df0dz(double k, double e_f, double T, double df0dz_integral, double coefficients[5][7], double kindex[], int aa[]);

double mu_elastic(double e_f,double T,double coefficients[5][7],double kindex[],double nu_elastic[],double g[],int points, int aa[]);
double mu_elasticH(double e_f,double T,double coefficients[5][7],double kindex[],
        double nu_elastic[], int points, int aa[]);

double mu_overall(double e_f,double T,double coefficients[5][7],double kindex[],
        double g[],double nu_el[],int points,int aa[]);
double mu_overallH(double e_f,double T,double coefficients[5][7],double kindex[],
        double g[], double h[], double nu_el[],int points,int aa[]);

void save_scattering_rate();
void save_perturbation();
double mu_po(double e_f,double T,double coefficients[5][7],double kindex[], double g_LO[],double g[],double nu_el[],int points, int aa[]);
double mu_poH(double e_f,double T,double coefficients[5][7],double kindex[], double g_LO[],double h_LO[],double nu_el[],int points, int aa[]);
void components_BTE(double T, int T_loop, double efefn, double efefp, int ii);

double J(double T,double m,double g_th[],  int points);
double DOS_value1(double energy, int a);
int read_ispin();
void fitting_band();
void solve_g(double T);
void save_results();
void initialize_array();
void pop_So(double T, double efef, int ii);

/*
void find_cbm(int spinorbit);
void bubbleSort(double arr[], int n);
void swap(double *xp, double *yp);
void analytical_fitting(int type);
int FindMinInd(double *array, int size);
void polyFit(double *x, double *y, double *z, int n, int l1, int l2);
void polyval(double *p, double *x, double *y, int p1, int x1);
double rsquare(double *data, double *fitted, int n);
double *linspace(double initial_value, double final_value, int samples);
void write_files(char type);
double conduction_dispersion(double k);
double dedk(double k);
double find_fermi();
*/
//------------ - *** Some Global variables ***----------------------------------------------------
extern double h_bar;
//--------------------------------------------------------------------------------------------
extern int cc, VASP;
extern double lm[10][10], volume1, ion_mass1[5];
extern int ion_numbers1[5], spin_orbit_coupling;
extern int npop_number;

extern double CBM, kcbm[3];
extern double kvbm[3], VBM;
//extern double conduction_band[10000][2];
//------------------------------------ make_band ----------------------------------------
extern double cond_band[limit2][2], val_band[limit2][2],poly[4000];
extern int count1,count2,countx;    //count1 for conduction band count2 for valence band
extern int count_orbital, count_orbital_p;
extern int count_DOS_n,count_DOS_p;
extern double E_F,n_e,n_h;
//extern double short_range_band_num[1000][3];
//extern int length_dd, maximum_degree, points ;
//extern double coefficients[10][10], kindex[10];
extern int number;
extern int fitting_1;
extern void read_procar();
//--------------------------------------- analytical_fitting----------------------------------
extern double coefficients[5][7],kindex[4],coefficients_cond[5][7],kindex_cond[4],coefficients_val[5][7],kindex_val[4];
extern int a11[2],b11[2];

//------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------

double mu_overall_time(double e_f,double T,double coefficients[5][7],double kindex[], double g[],double nu_el[],int points,int aa[], int j);
void conductivity_time(double T, int j);
void conductivity_freq();
void conductivity_with_freq(double T);


extern double Efield_time[limit1];  // unit cgs V/cm
extern double omega_s, freq[30], initial, tau;    // unit 1/s

extern double mobility_freqr2[30], mobility_freqi2[30], sigma_freqr2[30], sigma_freqi2[30];
extern double mobility_drude_freqr2[30], mobility_drude_freqi2[30], sigma_drude_freqr2[30], sigma_drude_freqi2[30];

extern int len_freq;
extern double J_time[limit1], sigma_time[limit1], mobility_time[limit1];
extern double J_freqr[30], sigma_freqr[30], Efield_freqr[30], J_freqi[30], sigma_freqi[30], Efield_freqi[30];
extern int time_variation;
extern int time_limit;
extern double g_time[limit2],g_time_old[limit2];
extern int freq_variation;

//-----------------------------------------------------------------------------------


extern int degree1,length_fraction;
extern double fraction[4];

//------------------------------------- inside read_input----------------------------------------------



extern double T_array[30],epsilon_s[30],epsilon_inf[30],Bgap[30],P_piezo[30],C_piezo_h14[30],n_array[30],Nd[30],Na[30],N_im[30];
extern double we[limit4],De[limit4];
extern double we_npop[limit5],De_npop[limit5];
extern int nfv[limit4];
extern int variation;
extern double kcbm[3],kvbm[3];
extern int kk, count_d ,count_t;
extern double k_min, k_trans, k_step_fine, k_step;
extern int points, points1, points2;


//extern double k_min0,k_trans0,k_step_fine0,k_step0;

extern int De_ionization,N_cb, N_vb, iterations, scattering_mechanisms[10], iv_number, fitting_1, fitting_2, fitting_3;

extern double Ed, c_lattice, rho, k_max, N_dis, omega_LO, omega_TO, E_deformation, C_long, C_trans, c_bar, C_11, C_12, C_44,
C_piezo_c14, P_piezo_h14, Uall, V0, xx, m ,m_h, T_trans ;

extern double Bfield;
extern double kcbm[3],kvbm[3];

extern double beta1[limit2], gH[limit2], hH[limit2], gH_rta[limit2], hH_rta[limit2];
extern double gH_LO[limit2], hH_LO[limit2], S_i_grid_g[limit2], S_i_grid_h[limit2];
extern double S_iLO_grid_g[limit2], S_iLO_grid_h[limit2], S_o_gridH[limit2], S_o_grid_totalH[limit2];
            
extern string  type;
extern int free_e, flag[50];
extern int len_T,len_n;
//----------------------------find_cbm_vbm-------------------------------------------------------------------------
extern double ecbm,evbm,kcbm1[3],kvbm1[3];
extern double kcbm[3],kvbm[3],VBM,CBM;
extern int NKPTS,NBVAL,NBTOT;
extern double orbital_decomposedd[1000][4], orbital_decomposedd_p[1000][4];
extern int ispin;
extern double DOS_n[limit2][2],DOS_p[limit2][2];
extern double fermi;

extern int ispin;

extern double k_grid[limit2], v_n[limit2], v_p[limit2], energy_n[limit2], energy_p[limit2];
extern double a_n[limit2], c_n[limit2], a_p[limit2], c_p[limit2], Ds_n[limit2], Ds_p[limit2];
extern double electric_driving_force[limit2], thermal_driving_force[limit2];

extern double df0dk_grid[limit2], f0x1_f0[limit2], electric_driving_force[limit2], thermal_driving_force[limit2], f_dist[limit2];

extern double kplus_grid_pop[limit2], kminus_grid_pop[limit2];

extern double betaplus_grid[limit2], betaminus_grid[limit2];

extern double  Aminus_grid[limit2], Aplus_grid[limit2], lambda_i_plus_grid[limit2], lambda_o_plus_grid[limit2];
extern double lambda_i_minus_grid[limit2], lambda_o_minus_grid[limit2];
extern double lambda_e_plus_grid[limit2][limit4], lambda_e_minus_grid[limit2][limit4];
extern double B_ii, D_ii, A_ii;

extern double lambda_e_plus_grid_npop[limit2][limit5], lambda_e_minus_grid_npop[limit2][limit5], N_npop[limit5];
extern double we_npop[limit5],De_npop[limit5];

extern double N_poph_atT, N_e[limit4], df0dz_integral_n, beta_constant;


extern double nu_deformation[limit2], nu_piezoelectric[limit2], nu_ionizedimpurity[limit2], nu_dislocation[limit2], nu_alloy[limit2];
extern double nu_neutralimpurity[limit2], nu_iv[limit2][limit4], nu_iv_total[limit2], nu_el[limit2], nu_npop[limit2][limit5];
extern double nu_npop_total[limit2];
extern double denom[limit2];

extern double nu_deformation_p[limit2][2][2], nu_ionizedimpurity_p[limit2][2][2], nu_el_p[limit2][2][2];
extern double nu_npop_p[limit2][2][2], nu_So_p[limit2][2][2];

extern double beta1[limit2], gH[limit2], hH[limit2], gH_rta[limit2], hH_rta[limit2];
extern double gH_LO[limit2], hH_LO[limit2], S_i_grid_g[limit2], S_i_grid_h[limit2];
extern double S_iLO_grid_g[limit2], S_iLO_grid_h[limit2], S_o_gridH[limit2], S_o_grid_totalH[limit2];

extern int plus_index_pop[limit2], minus_index_pop[limit2]; 

extern double g[limit2], g_rta[limit2], g_old[limit2], g_LO[limit2], g_iv[limit2], g_th[limit2], g_th_old[limit2], g_LO_th[limit2];
extern double S_o_grid[limit2], S_o_grid_total[limit2], S_i_grid[limit2], S_iLO_grid[limit2], S_i_th_grid[limit2], S_iLO_th_grid[limit2];
extern double result_g[limit2][15+1], result_g_LO[limit2][15+1], result_g_th[limit2][15+1], result_f[limit2][15+1];


extern double mobility_ii, mobility_po, mobility_to, mobility_npop, mobility_de, mobility_pe, mobility_dis;
extern double mobility_alloy, mobility_iv, mobility_neutral, mobility_npop, mobility_avg, mobility, mobility_rta;
extern double mobility_hall_ii, mobility_hall_po, mobility_hall_to, mobility_hall_npop, mobility_hall_de;
extern double mobility_hall_pe, mobility_hall_dis, mobility_hall_alloy, mobility_hall_iv;
extern double mobility_hall_neutral, mobility_hall_npop, mobility_hall_avg, mobility_hall, mobility_hall_rta, hall_factor1, hall_factor_rta1;
extern double sigma_hall_rta, sigma_hall, thermopower, sigma, sigma_rta;

extern double mobility_all[10] , calc_mobility[30][2] , calc_mobility_rta[30][2] ;
extern double calc_thermopower[30][2] , calc_sigma[30][2] , calc_sigma_rta[30][2] ;

extern double calc_mobility_pe[30][1] , calc_mobility_de[30][1] , calc_mobility_dis[30][1] , calc_mobility_ii[30][1] ;
extern double calc_mobility_po[30][1] , calc_mobility_to[30][1] , calc_mobility_alloy[30][1] , calc_mobility_iv[30][1] ;
extern double calc_mobility_neutral[30][1], calc_mobility_npop[30][1] ;

extern double mobility_hall_all[10], calc_mobility_hall[30][2] , calc_mobility_hall_rta[30][2];
extern double calc_sigma_hall[30][2] , calc_sigma_hall_rta[30][2] , hall_factor[30][2] , hall_factor_rta[30][2] ;

extern double calc_mobility_hall_pe[30][1] , calc_mobility_hall_de[30][1] , calc_mobility_hall_dis[30][1] ;
extern double calc_mobility_hall_ii[30][1] , calc_mobility_hall_iv[30][1] , calc_mobility_hall_neutral[30][1] ;
extern double calc_mobility_hall_npop[30][1];
extern double calc_mobility_hall_po[30][1], calc_mobility_hall_to[30][1], calc_mobility_hall_alloy[30][1];
extern double n0,Nd1,Na1, efef_n, efef_p, N_ii;
//-------------------------------------------------------------------------------------------------------------
#endif // MAIN_H_INCLUDED
