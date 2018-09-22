/*
*	Code for finding the radiation fluxes (solar and IR)
*
*
*	References:
*	(1) Numerical Prediction of the Performance of High Altitude Balloons
*		F. Kreith and J. F. Kreider
*		NCAR TECHNICAL NOTE, 1974
*
*	(2)	Balloon Ascent: 3-D Simulation Tool for the Ascent and Float
*		of High-Altitude Balloons
*		R. E. Farley 
*		AIAA Paper # 2005-7142, 2005
*
*	(3) Performance simulation of high altitude scientific balloons
*		Q. Dai, X. Fang, X. Li, L. Tian
*		Advances in Space Research, 49, 2012
**/

#include "rad_flux.h" 
const double Rad_flux::I0 = 1358;
const double Rad_flux::ECC = 0.016708;   
const double Rad_flux::P0 = 101325;       
const double Rad_flux::cloudAlbedo = 0.6;  

/*
* Expression from (2) Eq. 21
*
* @return The incident solar radiation above Earths atm (W/m^2)
*/
double Rad_flux::get_SI0() const {
	double ta = 2.*PI*doy/365.;
	double ec = pow(((1.+ECC)/(1.-ECC)),2) -1.;
	return I0*(1.+0.5*ec*cos(ta));
}

/*
* Expression from http://en.wikipedia.org/wiki/Position_of_the_Sun
*
* @return Approximate solar declination (rad) 
*/
double Rad_flux::get_declination() const {
	return (-23.44/RTD)*cos(2.*PI*(doy+10)/365);
}

/*
* Expression from (1) Eq. 35
*
* @param lat Lattitute (rad)
* @param h Solar hour angle (rad)
* @return The approximate solar zenith angle (rad)
*/
double Rad_flux::get_zenith(double lat, double h_ang) const {
	return acos(sin(lat)*sin(decl)+cos(lat)*cos(decl)*cos(h_ang));
} 

/*
* @param lat Lattitude (rad)
* @return The approximate solar hour angle for sunrise
*/
double Rad_flux::get_h0(double lat) const {
	return acos(-tan(lat)*tan(get_declination()));
}

/*
* Expression from (1) Eq. 38
*
* @param zen The solar zenith angle (rad)
* @param el Elevation
* @return The approximate air mass (unitless)
*/
double Rad_flux::get_air_mass(double zen, double el) const {
	double press = atm->get_P(el);
	return (press/P0)*(sqrt(1229 + pow((614*cos(zen)),2))-614*cos(zen));
}

/*
* Expression from (1) Eq. 37
*
* @param zen The solar zenith angle (rad)
* @param el Elevation (m)
* @return The atmospheric trasmittance (unitless)
*/
double Rad_flux::get_trans_atm(double zen, double el) const {
	if(fabs(zen) > PI/2.) return 0;
	double am = get_air_mass(zen, el);
	return 0.5*(exp(-0.65*am) + exp(-0.095*am));
}

/*
* @param zen The solar zenith angle (rad)
* @param el Elevation (m)
* @return The intensity of the direct solar radiation (W/m^2)
*/
double Rad_flux::get_direct_SI(double zen, double el) const {
	double SI0 = get_SI0();
	double trans = get_trans_atm(zen, el);
	return trans*SI0;
}

/*
* Expression from (3) Eq. (6)
* Note: This expression doesn't originate from (3) -- where does it come from?
* Note: Small type in denominator of original expression ?
*
* @param zen The solar zenith angle (rad)
* @param el Elevation (m)
* @return The intensity of the diffuse solar radiation from the sky (W/m^2)
*/
double Rad_flux::get_diffuse_SI(double zen, double el) const {
	if(zen > PI/2.) return 0;
	double SI0 = get_SI0();
	double trans = get_trans_atm(zen, el);
   if(el < cloudElev) {
      return (1-cloudFrac)*0.5*SI0*sin(PI/2.-zen)*(1.-trans)/(1-1.4*log(trans));
   } else {
      return 0.5*SI0*sin(PI/2.-zen)*(1.-trans)/(1-1.4*log(trans));
   }
}

/*
* Expression from (2) Eq. 26
*
* @param zen The solar zenith angle (rad)
* @param el Elevation (m)
* @return The intensity solar radiation reflected by the Earth (W/m^2)
*/
double Rad_flux::get_reflected_SI(double zen, double el) const {
	if(zen > PI/2.) return 0;
	double incident_SI = get_SI0();
   double tau_atm = get_trans_atm(zen,el);
   double albedo;
   if(el < cloudElev) {
      albedo = (1.-cloudFrac)*albedoGround;
   } else {
      albedo = (1.-cloudFrac)*(1-cloudFrac)*albedoGround + cloudAlbedo*cloudFrac;
   }
	return albedo*tau_atm*incident_SI*sin(PI/2.-zen);
}

/*
* Expression from (2) Eq. 24
* Note: Original expression has typo (0.95 should be 0.095, as below)
*
* @param el Elevation (m)
* @return the intensity of IR radiation emitted from earth (W/m^2)
*/
double Rad_flux::get_earth_IR(double el) const {
	double press = atm->get_P(el);
	double IR_trans = 1.716-0.5*(exp(-0.65*press/P0) +exp(-0.095*press/P0));
   double tEarth;
   if(el < cloudElev) {
      tEarth = tGround;
   } else {
      tEarth = tGround*(1.-cloudFrac) + atm->get_T(cloudElev)*cloudFrac;
   }
	return IR_trans*emissGround*SB_CONST*pow(tEarth,4);
}

/*
* Approximate linear fit from for lat 30-40 deg (1) Fig. 17
*
* @param h Altitude (m)
* @return The intensity of IR radiation emitted from the sky (W/m^2)
*/
double Rad_flux::get_sky_IR(double el) const {
	return fmax(-0.03*el+300.,50.0);
}

/*
* Set the day of the year
* @param _doy Day of year (day)
*/
void Rad_flux::set_doy(int _doy) {
   doy = _doy;
}


/*
* Read coefficients in earth.dat
*/
void Rad_flux::readCoeffs() {
	FILE * in = fopen("earth.dat", "r");
	if(in == NULL) {
		printf("could not find earth.dat\n");
		exit(0);
	}
	char line[BUFSIZ];
	fgets(line, BUFSIZ, in);
	sscanf(line, "%lg", &albedoGround);
	fgets(line, BUFSIZ, in);
	sscanf(line, "%lg", &emissGround);
	fgets(line, BUFSIZ, in);
	sscanf(line, "%lg", &tGround);
   fgets(line, BUFSIZ, in);
	sscanf(line, "%lg", &cloudFrac);
   fgets(line, BUFSIZ, in);
	sscanf(line, "%lg", &cloudElev);
}


Rad_flux::Rad_flux(int _doy, Std_atm * _atm) {
	readCoeffs();
   doy = _doy;
	decl = get_declination();
	atm = _atm;
}

Rad_flux::~Rad_flux() {
}
