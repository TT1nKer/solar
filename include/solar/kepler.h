#pragma once

#include "solar/body.h"

namespace solar {

// Solve Kepler's equation M = E - e*sin(E) for eccentric anomaly E
double solve_kepler(double M, double e, double tol = 1e-12, int max_iter = 50);

// Mean anomaly -> true anomaly (via eccentric anomaly)
double mean_to_true_anomaly(double M, double e);

// True anomaly -> eccentric anomaly
double true_to_eccentric_anomaly(double nu, double e);

// Eccentric anomaly -> mean anomaly
double eccentric_to_mean_anomaly(double E, double e);

// Convert orbital elements + true anomaly to cartesian state
// mu: gravitational parameter of central body
State elements_to_state(const OrbitalElements& elem, double mu, double true_anomaly);

// Convert orbital elements at a given Julian date to cartesian state
// Propagates mean anomaly from epoch to jd, then converts
State elements_to_state_at(const OrbitalElements& elem, double mu, double jd);

// Convert cartesian state to orbital elements
OrbitalElements state_to_elements(const State& s, double mu, double epoch = 0.0);

// Orbital period from semi-major axis and mu
double orbital_period(double a, double mu);

// Mean motion (rad/s) from semi-major axis and mu
double mean_motion(double a, double mu);

// Vis-viva: orbital speed at distance r for orbit with semi-major axis a
double vis_viva_speed(double r, double a, double mu);

} // namespace solar
