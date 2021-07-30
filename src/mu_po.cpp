#include"main.h"

double mu_po(double e_f,double T,double coefficients[5][7],double kindex[], double g_LO[],double g[],double nu_el[],int points, int aa[])
// It gives the effct of polar optical phonon on mobility in units of (cm^2/V.s)
// According to Equation (46) in Semiconductors and Semimetals volume1 10 (Rode's chapter),
// but with calculated group velocity from band structure and DOS both calculated from DFT:
{
    double integral_numerator = 0;
    double integral_denominator = 0;
    int factor = 100;

    double ddf,df,dv,k_step,de;
    //{
    if (T < T_trans)
    {
        for (int counter = 0;counter<points-1;counter++)
        {
            ddf = (df0dk(k_grid[counter+1],T,e_f,coefficients,kindex,aa)-df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa))/factor;
            for (int i = 0;i<=factor-1;i++)
                g[counter] = (-1)*e*E/(h_bar*nu_el[counter])*(df0dk(k_grid[counter],T,e_f,coefficients,kindex,aa)+i*ddf)*6.241509324e11;
                    //The last number is the conversion from convensional units to cancel out to be
                    // unitless (as in g)

        }
    }
    //}

    if (free_e ==1)
    {
        for (int counter = 0;counter<points-1;counter++)
        {
            dv = (v_n[counter+1]-v_n[counter])/factor;
            df = (f0(energy_n[counter+1],e_f,T)-f0(energy_n[counter],e_f,T))/factor;
            k_step = (k_grid[counter+1]-k_grid[counter])/factor;

            for (int i = 0;i<=factor-1;i++)
            {
                integral_numerator = integral_numerator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(v_n[counter]+i*dv)*g_LO[counter]/E;
                    // =1/E*int[g_LO(En)*DOS(En)*v(En)*dEn]
                    integral_denominator = integral_denominator+k_step*pow(((k_grid[counter]+i*k_step)/pi),2)*(f0(energy_n[counter],e_f,T)+i*df);
                // =int[f0(En)*DOS(En)*dEn]
            }
        }
    }
    else
    {
        for (int counter = 0;counter<=points-2;counter++)
        {
            de = (energy_n[counter+1]-energy_n[counter]);
            integral_numerator = integral_numerator + de*Ds_n[counter]/volume1*v_n[counter]*g_LO[counter]/E;
            // =1/E*int[g_LO(En)*DOS(En)*v(En)*dEn]

            integral_denominator = integral_denominator + de*Ds_n[counter]/volume1*f0(energy_n[counter],e_f,T);
                // =int[f0(En)*DOS(En)*dEn]
        }
    }

    double mobility_po = 1/3.0*integral_numerator/integral_denominator;
    // According to equation (46) in Rode's book; units of (cm^2/V.s)
    // work out automatically from other conventional units of group velocity (cm$
    return mobility_po;

}

