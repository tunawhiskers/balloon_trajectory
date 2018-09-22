#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "rad_flux.h"
#include "std_atm.h"
#include "sphere_balloon.h"
#include "consts.h"	

double lat = 35.08333/RTD;//Lattitude of ABQ,NM
double min_el = 132.6;
int doy = 180;

Std_atm atm;
Rad_flux rad(doy, &atm);
Sphere_balloon bal(&atm, &rad);

//current balloon state
struct State {
    double v; 	//verticle velocity
    double el;	//elevation (m)
    double Ti;	//internal temperature (K)
    double Ts;	//surface temperature (K)
};


//derivative of balloon state wrt time
struct stateDeriv {
    double dv;
    double del;
    double dTi;
    double dTs;
};

/*
* @param state Current balloon state
* @return The verticle acceleration of balloon (m/s^2)
*/	
double get_a(const State & state) {
   double rho_int = atm.get_P(state.el)/(atm.Rsp_air*state.Ti);
   double rho_atm = atm.get_rho(state.el);
   double F_b = (rho_atm - rho_int)*bal.get_vol()*atm.get_g(state.el);
   double F_d = bal.get_Cd(state.v, state.el)*(0.5*rho_atm*fabs(state.v)*state.v)*bal.get_cs_area();
   double vm = bal.get_mass() + rho_atm*bal.get_vol() + bal.get_virt_mass()*rho_atm*bal.get_vol();
   return (F_b  - F_d - bal.get_mass()*atm.get_g(state.el))/vm;
}

/*
/ @param state The current balloon state
* @param h The current hour angle
* @return Rate of change of the surface temperature wrt time (K/s)
*/
double get_dTs(const State & state, const double h) {
   double q_rad  = bal.get_q_rad(lat, state.el, h);
	double q_surf = bal.get_sum_q_surf(q_rad, state.Ts, state.el, state.v);
	double q_int  = bal.get_sum_q_int(state.Ts, state.Ti, state.el);
   return (q_surf-q_int)/bal.get_therm_mass();
}

/*
* @param state Current balloon state
* @return Rate of change of internal temperature (K/s)
*/
double get_dTi(const State & state) {
		double q_int  = bal.get_sum_q_int(state.Ts, state.Ti, state.el);
		double tm_air = atm.get_rho(state.el)*bal.get_vol()*atm.Cp_air0;
      return q_int/tm_air;
}

/*
* @param init The initial state of the balloon
* @param dt Time step (s)
* @param h Current hour angle (rad)
* @param d Balloon state derivative
* @return The balloon state derivative dt in the future, given
*		an initial deriviative d
*/
stateDeriv eval(const State & init, const double dt, const double h,  const stateDeriv & d) {
   State s;
   s.el = init.el + d.del*dt;
   s.v  = init.v  + d.dv*dt;
   if ((s.el < min_el) & (s.v < 0)) {
      s.v  = 0;
      s.el = min_el;
   }
   s.Ti = init.Ti + d.dTi*dt;
   s.Ts = init.Ts + d.dTs*dt;
      
    
   stateDeriv out;
   out.del = s.v;
   out.dv = get_a(s);
   out.dTs = get_dTs(s,h + dt*((PI/6.)/(3600.))); //dt is converted to hour angle
   out.dTi = get_dTi(s);
   return out;
}

/*
* Use RK4 method to integrate a balloon state to dt in the future
* @param state current balloon state
* @param h current hour angle (rad)
* @param dt time step (s)
*/
void integrate(State & state, const double h, const double dt) {
    stateDeriv zero;
    zero.del = 0;
    zero.dv = 0;
    zero.dTs = 0;
    zero.dTi = 0;
    stateDeriv a = eval(state, 0.0, h, zero);
    stateDeriv b = eval(state, 0.5*dt, h, a);
    stateDeriv c = eval(state, 0.5*dt, h, b);
    stateDeriv d = eval(state, dt, h, c);

    double del_dt = (1./6.)*(a.del + 2*b.del + 2*c.del + d.del);
    double dv_dt  = (1./6.)*(a.dv  + 2*b.dv  + 2*c.dv  + d.dv);
    double dTi_dt = (1./6.)*(a.dTi + 2*b.dTi + 2*c.dTi + d.dTi);
    double dTs_dt = (1./6.)*(a.dTs + 2*b.dTs + 2*c.dTs + d.dTs);

    state.el += del_dt*dt;
    state.v += dv_dt*dt;
    
   if ((state.el < min_el) & (state.v < 0)) { //balloon does not descend below el = min_el
      state.v  = 0;
      state.el = min_el;
    }

    state.Ti += dTi_dt*dt;
    state.Ts += dTs_dt*dt;
}


int main(int nargs, char ** args) {
   double h0 = rad.get_h0(lat); //hour angle of dawn/sunset
   State cur;
   cur.el = min_el;
   cur.v = 0;
   cur.Ts = atm.get_T(cur.el);
   cur.Ti = atm.get_T(cur.el);

	double dt = 1.0;
	double dh = 15*dt/(RTD*3600);
	int i = 0;
   //printf("time(hr), elevation(m), internal T(K), internal-ambient T (K), velocity (m/s)\n");

   for(double h = -1.1*h0; h < 2.1*h0; h+= dh) {
	   double t = (((h+h0)*RTD/15.)*3600);
		integrate(cur,h,dt);
      if(i%100==0) {
		   printf("%f %f %f %f %f\n", t/3600, cur.el, cur.Ti, cur.Ti-atm.get_T(cur.el), cur.v);
      }
      if( (h > h0) & (cur.el == min_el)) break;
      i++;
	}
   return 0;
}

