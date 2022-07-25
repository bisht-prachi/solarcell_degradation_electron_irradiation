#include<iostream>
#include<fstream>
#include<cstdio>
#include<cmath>

/******
        Degraded IV characteristic under 1 MeV electron irradiation using C++
        by Prachi Bisht

        (see readme file for the description)

        Units: kg, cm, seconds
*******/


const double pi = 3.14159265358979323846;  /* pi */

const double A = 32.0;                  //exposed solar cell area, cm square
const double q = 1.6e-19;               //fundamental charge
const double kb = 1.38e-23;             //Boltzmann constant
const double h = 6.6e-34;               //Planck's constant
const double c = 3e8;                   //speed of light in vacuum
const double eps = 8.85e-12;            //vacuum permittivity
const double me = 9.1e-31;              //mass of electron
const double T = 300.0;                 //room temprature
const double Vt = 0.025875;             //thermal energy at room temp: (kT/q) eV

/*******
        Following parameters are specific to an n/p GaAs solar cell,
        taken from the text "Semiconductors and Semimetals, Vol. 11" by H J Hovel

*******/

const double ni = 2.37e6;               //intrinsic charge density
const double Na = 1e17;                 //acceptor doping concentration in p-type GaAs
const double Nd = 2e18;                 //donor doping concentration in n-type GaAs

const double Dp = 6.12;                 //diffusion coefficient for hole
const double tau_p0 = 0.79e-8;          //diffsuion time for hole
const double Sp = 0;//1e4;              //surface charge velocity (0 for ideal case)

const double Dn = 116;                  //diffusion coefficient for hole
const double tau_n0 = 20e-9;            //diffusion time for hole
const double Sn = 0;//1e4;              //surface charge velocity (0 for ideal case)

const double xj = 0.1e-4;
const double H = 3e-4;
const double R = 0.04;


int main()
{
	double Vd = Vt*log(Nd*Na/(ni*ni));                  //built in potential
	double W = sqrt((2*eps*(Na+Nd)*Vd)/(q*Na*Nd));      //junction width


	std::ofstream fout, gout;
	std::ifstream fin, gin;


	int m = 66;
	double alphas[m] = {0.0};                           //array to store absorption coefficient alpha at lambda
	double Fs[m] = {0.0};                               //array to store photon flux in a wavelength interval

/*******
        store absorption coefficient

*******/

	fin.open("alpha_gaas.dat", std::ios::in);           //file contains alpha vs lambda for GaAs
	long double t = 0.0;
	int j = 0; int k = 0;

	while(true)
	{
		fin>>t;
		if (j %2 == 0)
			++j;
		else
		{
			if (k*10 + 300 < 500) alphas[k] = 0.0;
			else alphas[k] = 1.0*t;
			++j;
			++k;
		}

		if (fin.eof() ) break;

	}
	fin.close();

/*******
        store photon flux

*******/
	fin.open("am0_e.dat", std::ios::in);                //file contains photon flux vs lambda
	t = 0.0;                                            //(for an interval of 10 nm for GaAs
	j = 0; k = 0;

	double lambda = 0.0;
	double w_at_lambda = 0.0;
	while(true)
	{
		fin>>t;
		if (j %2 == 0)
		{
			lambda = t;
			++j;
		}

		else
		{
			w_at_lambda = 10.0*t;			            //rememember to multiply with 10
			Fs[k] = 5*lambda*w_at_lambda*1e11;          //as the data file contains no. of photons
			++j;                                        //in / cm^2 s / 10nm
			++k;
		}

		if (k == m ) break;
	}
	fin.close();

/*******
        defining current terms: photocurrent from p side, photocurrent from n side,
                              : photocurrent from depletion region,
                              : diffusion current from p side,
                              : diffusion current from n side,
                              : recombination current from depletion region
*******/

	double Jp, Jn, Jdr, Jdiffp, Jdiffn, Jrg;


    fout.open("iv_gaas_1e15_1MeV.dat", std::ios::out);
    double phi = 1e15;                                              //input fluence, try different values

    //for electrons
	double KLn = 1.84e-08;                                          //diffusion length damage coefficient (from Hovel)
	double Ln0 = sqrt(Dn*tau_n0);                                   //original diffusion length
	double temp = (1/(Ln0*Ln0)) + KLn*phi;
	double Ln = sqrt(1/temp);                                       //updated diffusion length
	double tau_n = (Ln*Ln)/Dn;                                      //updated diffusion time

	//for holes
	double KLp = 5.57e-07;                                          //from Hovel
	double Lp0 = sqrt(Dp*tau_p0);
	double temp1 = (1/(Lp0*Lp0)) + KLp*phi;
	double Lp = sqrt(1/temp1);
	double tau_p = (Lp*Lp)/Dp;


/*******
        light current: exact solutions of the drift-diffusion ODE given in Hovel
*******/


    double Jl = 0.0;
    for (int i = 0; i<m; ++i)
    {
        double F = Fs[i];
        double alpha = alphas[i];
        Jp = (q*F*(1 - R)*alpha*Lp/(alpha*alpha*Lp*Lp - 1))*( -1.0*alpha*Lp*exp(-1.0*alpha*xj)
            + ( (((Sp*Lp)/Dp + alpha*Lp) - exp(-1.0*alpha*xj)*( ((Sp*Lp)/Dp)*cosh(xj/Lp) + sinh(xj/Lp)  ))/( ((Sp*Lp)/Dp)*sinh(xj/Lp) + cosh(xj/Lp) ) )  );
        Jn = (q*F*(1 - R)*alpha*Ln/(alpha*alpha*Ln*Ln - 1))*exp(-alpha*(xj+W))*( alpha*Ln
            - (( ((Sn*Ln)/Dn)*(cosh(H-xj-W) - exp(-1.0*alpha*(H-xj-W))) + sinh((H-xj-W)/Ln) + alpha*Ln*exp(-alpha*(H-xj-W)) )/( ((Sn*Ln)/Dn)*sinh((H-xj-W)/Ln) + cosh((H-xj-W)/Ln) )) );
        Jdr = q*F*(1 - R)*exp(-1.0*alpha*xj)*(1 - exp(-1.0*alpha*W));

        double Jph = (Jp + Jn + Jdr);
        Jl += Jph;
    }


/*******
        dark current: diffusion (exact solutions of the drift-diffusion ODE given in Hovel)
*******/

    Jdiffp = (q*Dp*ni*ni/(Lp*Nd))*((Sp*Lp*cosh(xj/Lp) + Dp*sinh(xj/Lp))/(Sp*Lp*sinh(xj/Lp) + Dp*cosh(xj/Lp)));
    Jdiffn = (q*Dn*ni*ni/(Ln*Na))*((Sn*Ln*cosh((H-xj-W)/Lp) + Dp*sinh((H-xj-W)/Lp))/(Sp*Lp*sinh((H-xj-W)/Lp) + Dp*cosh((H-xj-W)/Lp)));


/*******
        dark current: recombination (see Sah, Noyce, Shockley (SNS) theory 1957)
*******/

    double f_of_b = pi*0.5;// - arctan(b/(sqrt(b*b + 1)))*(1/(sqrt(b*b + 1)));

    double pmax = 0;
    double Jdiff = 0.0;
    double Jdark = 0.0;

    for (int i = 0; i<210; ++i)
    {
        double Vj = 0.01*i - 1.0;

        Jrg = ((kb*T*ni*W*2*f_of_b)/(sqrt(tau_p*tau_n)))*(sinh((q*Vj)/(2*kb*T))/(Vd-Vj));
        Jdiff = (Jdiffp + Jdiffn)*(exp(Vj/Vt) - 1);

        Jdark = Jdiff + Jrg;
        double J = -1.0*(Jdark - Jl);

        fout<<Vj<<"\t"<<J<<"\n";
        if (Vj == 0.0)
            std::cout<<"Isc (A)"<<" "<<J*A<<"\n";

        if (Vj*J>pmax)
            pmax = Vj*J;
            //std::cout<<"Pmax"<<" "<<pmax<<"\n";



        if (J<0)
        {
            std::cout<<"Voc (V)"<<" "<<Vj<<"\n"; break;
        }


    }
    std::cout<<"Pmax (W)"<<" "<<pmax*A<<"\n";
    std::cout<<"IV characteristics for fluence"<<" "<<phi<<" "<<"generated in the folder"<<"\n";

    fout.close();


}
