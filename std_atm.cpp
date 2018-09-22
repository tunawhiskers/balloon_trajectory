/*
* Code for finding the state (density, Pressure, Temperature)
* and transport (viscosity and thermal conductivity) properties
* using US 1976 Standard Atm and Sutherlands Law
*/

#include "std_atm.h"
const double Std_atm::Cp_air0 = 1003.8;
const double Std_atm::Rsp_air = 287.058;
const double Std_atm::visc0 = 1.725E-5; 
const double Std_atm::cond0 = 2.428E-2;
const double Std_atm::T_ref = 275;    
const double Std_atm::C = 110.4;

/*
* @param T Temperature (K)
* @param rho Density (kg/m^3)
* @return Sutherland's Law thermal diffusivity (W/(m*K))
*/
double Std_atm::get_cond(double T) const {
   return 0.0241*pow((T/273.15),0.9);
}

/*
* @param T Temperature (K)
* @param rho Density (kg/m^3)
* @return Sutherland's Law thermal diffusivity (m^2/s)
*/
double Std_atm::get_td(double T, double rho) const {
   return get_cond(T)/(rho*Cp_air0);
}

/*
* @param T Temperature (K)
* @param rho Density (kg/m^3)
* @return Sutherland's Law viscosity (Pa*s)
*/
double Std_atm::get_visc(double T) const {
   return 1.458E-6*pow(T,1.5)/(T+110.4);
}


/*
* @param T Temperature (K)
* @param rho Density (kg/m^3)
* @return Sutherland's Law kinematic viscosity (m^2/s)
*/
double Std_atm::get_kin_visc(double T, double rho) const {
   return get_visc(T)/rho;
}

/*
* @param T Temperature (K)
* @return Prantl Number
*/
double Std_atm::get_Pr(double T) const{
	return get_visc(T)*Cp_air0/get_cond(T);
}


/*
* Interpolate between two nearest points in standard atm
* @param h Elevation (m)
* @return Interpolated densiity (kg/m^3)
*/
double Std_atm::get_rho(double h) const{
	unsigned int i = 0;
	for(i = 0; i < atm_ar.size(); i++) {
		if(atm_ar[i].h>h) break; 
   	}
   	double rho = (atm_ar[i].rho - atm_ar[i-1].rho)/(atm_ar[i].h - atm_ar[i-1].h);
  	rho *= (h - atm_ar[i-1].h);
  	rho += atm_ar[i-1].rho;
	return rho;
}

/*
* Interpolate between two nearest points in standard atm
* @param h Elevation (m)
* @return Interpolated Temperature (K)
*/
double Std_atm::get_T(double h) const{
	unsigned int i = 0;
	for(i = 0; i < atm_ar.size(); i++) {
		if(atm_ar[i].h>h) break; 
   	}
   	double T = (atm_ar[i].T - atm_ar[i-1].T)/(atm_ar[i].h - atm_ar[i-1].h);
  	T *= (h - atm_ar[i-1].h);
  	T += atm_ar[i-1].T;
	return T;
}

/*
* Interpolate between two nearest points in standard atm
* @param h Elevation (m)
* @return Interpolated g (m/s^2)
*/
double Std_atm::get_g(double h) const {
	unsigned int i = 0;
	for(i = 0; i < atm_ar.size(); i++) {
		if(atm_ar[i].h>h) break; 
   	}
   	double g = (atm_ar[i].g - atm_ar[i-1].g)/(atm_ar[i].h - atm_ar[i-1].h);
  	g *= (h - atm_ar[i-1].h);
  	g += atm_ar[i-1].g;
	return g;
}

/*
* Interpolate between two nearest points in standard atm
* @param h Elevation (m)
* @return Interpolated pressure (Pa)
*/
double Std_atm::get_P(double h) const{
	unsigned int i = 0;
	for(i = 0; i < atm_ar.size(); i++) {
		if(atm_ar[i].h>h) break; 
   	}
   	double p = (atm_ar[i].p - atm_ar[i-1].p)/(atm_ar[i].h - atm_ar[i-1].h);
  	p *= (h - atm_ar[i-1].h);
  	p += atm_ar[i-1].p;
	return p;
}

/*
* Reads in us1976_std_atm.dat file
*/
void Std_atm::read_atm() {
	char line[BUFSIZ];
	FILE * in = fopen("us1976_std_atm.dat", "r");
	while(fgets(line, BUFSIZ, in)) {
    	atm_ref a;
		if(line[0] == '#') continue;
		sscanf(line, "%lg %lg %lg %lg %lg", &a.h, &a.T, &a.g, &a.p, &a.rho);
		a.T += 273.15;    //convert to K
		a.p *= 1E4;       //convert to Pa
		a.rho /= 10;      //convert to kg/m^3
		atm_ar.push_back(a);	
	}
	fclose(in);
}

Std_atm::Std_atm() {
	read_atm();
}

Std_atm::~Std_atm() {
}
