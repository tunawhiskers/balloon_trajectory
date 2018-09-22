#include <cstdio>
#include <vector>
#include <cmath>

#ifndef ATM_H
#define ATM_H

using namespace std;

struct atm_ref {
	double h;         //elevation (m)
	double T;         //Temperature K
	double p;         //pressure Pa
	double rho;       // density (kg/m^3)
	double g;         //acceleration due to gravity 
};

class Std_atm {
	private:
		vector <atm_ref> atm_ar;
		//General properties of Air

		//Transport Properties Params
		static const double visc0;  //Pa s
		static const double cond0;   //W/(m^2 K)
		static const double T_ref;   //Reference Temp (above transport properties @ this temp)
		static const double C;         //constants for sutherlands law

	public:
		static const double Rsp_air; //J/(kg K) (constant)
		static const double Cp_air0;  //J/(kg/K) (approx constant)

		Std_atm();
		~Std_atm();
		//Transport property functions
		double get_cond(double T) const;
		double get_td(double T, double rho) const;
		double get_visc(double T) const;
		double get_kin_visc(double T, double rho) const;
		double get_Pr(double T) const;

		//Standard Atm Interpolations
		double get_rho(double h) const;
		double get_T(double h) const;
		double get_g(double h) const;
		double get_P(double h) const;

		//Read in standard atm
		void read_atm();
};
#endif 
