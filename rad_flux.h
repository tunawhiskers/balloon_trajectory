#include <cmath>
#include <cstdlib>
#include "consts.h"
#include "std_atm.h"
#ifndef RAD_FLUX_H
#define RAD_FLUX_H

class Rad_flux {
	private:
		static const double I0;		//Solar Intensity in Space (W/m^2)
		static const double ECC;   //Earths's orbital eccentricity
		static const double P0;    //pressure at sealevel (Pa)
      static const double cloudAlbedo; //cloud Albedo
		double albedoGround;  		//Ground Albedo 
		double emissGround;  		//Emissivity of Ground
		double tGround;     			//Temperature of Ground
      double cloudFrac;          //Fraction of cloud coverage
      double cloudElev;          //Cloud elevation (m)
		int doy;
		double decl;
		Std_atm * atm;
		void readCoeffs();
	public:
		double get_SI0() const;
		double get_declination() const;
		double get_zenith(double lat, double h_ang) const;
		double get_h0(double lat) const;
		double get_air_mass(double zen, double el) const;
		double get_trans_atm(double zen, double el) const;
		double get_direct_SI(double zen, double el) const;
		double get_diffuse_SI(double zen, double el) const;
		double get_reflected_SI(double zen, double el) const;
		double get_earth_IR(double press) const;
		double get_sky_IR(double el) const;
      void set_doy(int _doy); 
		Rad_flux(int doy, Std_atm * _atm);
		~Rad_flux();
};

#endif
