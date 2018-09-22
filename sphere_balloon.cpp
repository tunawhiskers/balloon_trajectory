#include "sphere_balloon.h"

/*
*@param Ra Raleigh number
*@return External Nusselt Number according to (1)
*/
double Sphere_balloon::get_Nu_ext(double Ra, double Re, double Pr) const {
	double Nu_n;
	if(Ra < 1.5E8) {
		Nu_n = 2.0 + 0.6*pow(Ra,0.25);
	} else {
		Nu_n = 0.1*pow(Ra, 0.34);
	}
	double Nu_f;
	if(Re < 5E4) {
		Nu_f = 2 + 0.47*sqrt(Re)*pow(Pr, (1./3.));
	} else {
		Nu_f = (0.0262*pow(Re, 0.8) - 615.)*pow(Pr, (1./3.));
	}
	return fmax(Nu_f, Nu_n);
}

/*
*@param T_s Surface temperature of sphere (K)
*@param el Elevation (m)
*@return Power transferred from sphere to surrounding atmosphere due to convection(W)
*/
double Sphere_balloon::get_q_ext(double T_s, double el, double v) const {
	double T_atm = atm->get_T(el);
   double p_atm = atm->get_P(el);
	double rho_atm = atm->get_rho(el);
	double Pr_atm = atm->get_Pr(T_atm);

	double T_avg = 0.5*(T_atm + T_s);
   double rho_avg = p_atm/(atm->Rsp_air*T_avg);
   double Pr_avg = atm->get_Pr(T_avg);
   double exp_coeff = 1./T_avg;
   double kin_visc = atm->get_visc(T_avg)/rho_avg;
   double Ra = Pr_avg*atm->get_g(el)*fabs(T_s-T_atm)*pow(diameter,3)*exp_coeff/(kin_visc*kin_visc);
	double Re = rho_atm*v*diameter/atm->get_visc(T_atm);
	double Nu = get_Nu_ext(Ra, Re, Pr_atm);
	double k = atm->get_cond(T_avg);
	double h = Nu*k/diameter;
	return h*surfArea*(T_s-T_atm);
}

/*
*@param q_rad Power input from external radiation (W)
*@param T_S Surface temperature (K)
*@param el Elevation (m)
*@return The sum of power input to the balloon surface (W)
*/
double Sphere_balloon::get_sum_q_surf(double q_rad, double T_s, double el, double v) const {
	double q_ce = -get_q_ext(T_s, el, v);
	double q_re = -emissEnv*SB_CONST*pow(T_s,4)*surfArea;
	return q_rad + q_ce + q_re;
}

/*
*@param q_rad Power input from external radiation (W)
*@param diameter Balloon diameter (m)
*@param Elevation (m)
*@param Balloon velocity (m/s)
*@return Surface temperature (K), solved to steady state using NR method
*/
double Sphere_balloon::solve_T_surf(double q_rad, double el, double v) const {
	double dT = 1; //initial guesses
	double T_s = atm->get_T(el)+10; //T_s > T_atm
	for(int i = 0; i < 10; i++) { //10 iterations should be enough
		double q2 = get_sum_q_surf(q_rad, T_s+dT, el, v);
		double q1 = get_sum_q_surf(q_rad, T_s, el, v);
		double dqdT = (q2-q1)/dT;
		dT = q1/dqdT;
		T_s -= dT;
		if(fabs(dT) < 1E-10) break;
	}
	return T_s;
}

/*
*@param Ra Internal balloon Raleigh Number
*@return Internal Nusselt number according to (1)
*/
double Sphere_balloon::get_Nu_int(double Ra) const {
	if(Ra < 1.35E8) {
		return 2.5*(2+0.6*pow(Ra,0.25));
	} else {
		return 0.325*pow(Ra, 0.333);
	}
}

/*
*@param T_s Surface temperature (K)
*@param T_i Internal temperature (K)
*@param el Elevation (m)
*@return The power input from the interior to the surface of the balloon due to convection (W)
*/
double Sphere_balloon::get_q_int(double T_s, double T_i, double el) const {
	double p_atm = atm->get_P(el);
	
	double T_avg = 0.5*(T_s+T_i);
	double rho_avg = p_atm/(atm->Rsp_air*T_avg);
	double Pr = atm->get_Pr(T_avg);
	double exp_coeff = 1./T_avg;
	double kin_visc = atm->get_visc(T_avg)/rho_avg;
	double Ra = Pr*atm->get_g(el)*fabs(T_i-T_s)*pow(diameter,3)*exp_coeff/(kin_visc*kin_visc);
	double Nu = get_Nu_int(Ra);
	double k = atm->get_cond(T_avg);
	double h = Nu*k/diameter;
	return h*surfArea*(T_i-T_s);
}

/*
* @param T_s Balloon surface temperature (K)
* @param T_i Internal balloon temperature (K)
* @param el Elevation (m)
* @return The sum of a power inputs to the balloon interior (W)
*/
double Sphere_balloon::get_sum_q_int(double T_s, double T_i, double el) const {
	double q_ci = -get_q_int(T_s, T_i, el);
	//should there even be IR transfer between internal air & balloon surf?
	//i dont think so
   //double q_ri = 0*E_int*SB_CONST*(pow(T_s,4)-pow(T_i,4))*surface_area;
	return q_ci;// + q_ri;
}

/*
* @param T_s Balloon surface temperature (K)
* @param el Elevation (m)
* @return Internal temperature of balloon using NR method to steady state
* This method is pretty useless as T_i will just go to T_s unless their is
* another heat source inside--but if you add some internal heat source (panels)
* it could come in handy
*/
double Sphere_balloon::solve_T_i(double T_s, double el) const {
	double T_i = T_s + 10;
	double dT = 1;
	for(int i = 0; i < 10; i++) {
		double q1 = get_sum_q_int(T_s, T_i, el);
		double q2 = get_sum_q_int(T_s, T_i+dT, el);
		double dqdT = (q2 - q1)/dT;
		dT = q1/dqdT;
		T_i -= dT;
		if(fabs(dT) < 1E-11) break;
	}
	return T_i;
}

/*
* @param lat Balloon Lattitude  (rad)
* @param el Balloon Elevation (m)
* @param h Solar Hour Angle (rad)
* @return The total power input to the surface of the balloon from radiation (W)
*/
double Sphere_balloon::get_q_rad(double lat, double el, double h) const {
		double hca = asin(RE/(RE+el)); //half cone angle
        double vf = 0.5*(1. - cos(hca));

        double zen = rad->get_zenith(lat, h);
        double direct_I = rad->get_direct_SI(zen, el);
        double power_direct = direct_I*totAbs*projArea;

        double diffuse_I = rad->get_diffuse_SI(zen, el);
        double power_diffuse = diffuse_I*totAbs*(1.-vf)*surfArea;

        double reflected_I = rad->get_reflected_SI(zen, el);
        double power_reflected = reflected_I*totAbs*vf*surfArea;

        double earth_IR = rad->get_earth_IR(el);
        double power_earth_IR = earth_IR*totAbs*vf*surfArea;

        double sky_IR = rad->get_sky_IR(el);
        double power_sky_IR = sky_IR*totAbs*(1.-vf)*surfArea;

        return power_direct+power_diffuse+power_reflected+power_earth_IR+power_sky_IR;
}

/*
* @return Balloon volume (m^3);
*/
double Sphere_balloon::get_vol() const {
	return (4./3.)*PI*pow(0.5*diameter,3);
}

double Sphere_balloon::get_cs_area() const{
   return PI*diameter*diameter/4.;
}

/*
* @return Balloon skin thermal mass (J/K)
*/
double Sphere_balloon::get_therm_mass() const {
	return massEnv*cpEnv;
}


/*
* @return Total balloon mass (kg)
*/
double Sphere_balloon::get_mass() const {
    return massEnv + massPayload;
}

/*
* @return The drag coefficient
*/
double Sphere_balloon::get_Cd(double v, double el) const {
   return 0.8;
   double Re = atm->get_rho(el)*fabs(v)*diameter/atm->get_visc(atm->get_T(el));
   unsigned int i = 0;
   for(i = 0; i < cd_arr.size(); i++) {
      if(cd_arr[i].Re>Re) break;
   }
   if(i==0) return cd_arr[0].Cd;
   double Cd = (cd_arr[i].Cd - cd_arr[i-1].Cd)/(cd_arr[i].Re - cd_arr[i-1].Re);
   Cd *= (Re - cd_arr[i-1].Re);
   Cd += cd_arr[i-1].Cd;
   return Cd;
}

/*
* @return The virtual mass coefficient
*/
double Sphere_balloon::get_virt_mass() const {
   return vmCoeff;
}

void Sphere_balloon::read_balloon_coeffs() {
   FILE * in = fopen("balloon.dat", "r");
   char line [BUFSIZ];
	if(in == NULL) {
      printf("missing file: balloon.dat\n");
      exit(0);
   } else {
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &diameter);
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &massPayload);
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &rhoEnv);
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &envThickness);
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &cpEnv);
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &vmCoeff);
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &radAbs);
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &radTrans);
      fgets(line, BUFSIZ, in);
		sscanf(line, "%lg", &radRef);
	}
   fclose(in);

   in = fopen("drag_sphere.dat", "r");
   if(in == NULL) {
      printf("missing file: drag_sphere.dat\n");
      exit(0);
   } else {
      while(fgets(line, BUFSIZ, in)) {
         drag ref;
         sscanf(line, "%lg %lg", &ref.Re, &ref.Cd);
         cd_arr.push_back(ref);
      }
   }
}


Sphere_balloon::Sphere_balloon(Std_atm * _atm, Rad_flux * _rad) {
	read_balloon_coeffs();
   surfArea = PI*diameter*diameter;
	projArea = 0.25*PI*diameter*diameter;
	massEnv = surfArea*rhoEnv*envThickness;
	atm = _atm;
	rad = _rad;
   radRef = radRef + radRef*radRef + radRef*radRef*radRef;
	emissEnv = radAbs;
	totAbs = radAbs + radAbs*radTrans + radAbs*radTrans*radRef;
}

Sphere_balloon::~Sphere_balloon() {
}
