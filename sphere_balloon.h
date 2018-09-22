#include <cstdlib>
#include <iostream>
#include <fstream>
#include "std_atm.h"
#include "rad_flux.h"
#include "consts.h"

#ifndef SPHERE_BALLOON_H
#define SPHERE_BALLOON_H

//constants

struct drag {
   double Re;
   double Cd;
};

class Sphere_balloon {
	private:

		//balloon properties
      double rhoEnv;       //envelope density (kg/m^3)
		double cpEnv;        //envelop specifiv heat (J/(kgK))
		double envThickness; // (m)
      double vmCoeff;      //virtual mass
		double diameter;     // (m)
      double radAbs;      	//adsorbed radiation fraction
		double radRef;      	//reflected radation fraction
		double radTrans;    	//transitted radiation fraction
      double massPayload; 	// (kg)
      
      //calculated properted		
		double emissEnv;     //envelope emissivity
		double totAbs;      	//total effective absorbtivity
		double surfArea;      
		double projArea;    	//projected area
		double tm;   			//thermal mass (J/K)
		double massEnv; 		// (kg)
   
      vector <drag> cd_arr;

		Std_atm * atm;       // atmosphere object
		Rad_flux * rad;      // radiation flux object
      
	public:
		double get_therm_mass() const;
		double get_virt_mass() const;
      double get_vol() const;
		double get_mass() const;
      double get_Cd(double v, double el) const;
      double get_cs_area() const;

		double get_q_rad(double lat, double el, double h) const;
		double solve_T_i(double T_s, double el) const;
		double get_sum_q_int(double T_s, double T_i, double el) const;
		double get_q_int(double T_s, double T_i, double el) const;
      double get_Nu_int(double Ra) const;
		
      double solve_T_surf(double q_rad, double el, double v) const;
		double get_sum_q_surf(double q_rad, double T_s, double el, double v) const;
		double get_q_ext(double T_s, double el, double v) const;
		double get_Nu_ext(double Ra, double Re, double Pr) const;
		
      void read_balloon_coeffs();
      Sphere_balloon(Std_atm * _atm, Rad_flux * _rad);
		~Sphere_balloon();
};
	
#endif
