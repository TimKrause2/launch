#ifndef ORBIT_H
#define ORBIT_H

#include "caams.hpp"

struct OrbitalParams
{
    double a; // semi-major axis length in meters
    double e; // eccentricity
    double i; // inclination
    double Omega; // longitude of ascending node
    double omega; // argument of periapsis
    double nu; // true anomoly
    // all angles in radians
};

//
// Modified Equinoctial Orbital Elements
//

struct ModEquOrbitalParams
{
    double p; // semiparameter
    double f; // eccentricity vector x-component
    double g; // eccentricity vector y-component
    double h;
    double k;
    double L; // true longitude
    void FromOrbitalParams(OrbitalParams &c);
    void ToOrbitalParams(OrbitalParams &c);
    Eigen::VectorXd ToVector(void);
    void FromVector(Eigen::VectorXd v);
};

void orbital_plane_elements(
        Eigen::Vector3d const &h,
        Eigen::Vector3d const &e,
        OrbitalParams &params);

void orbital_elements(
        Eigen::Vector3d r, // location of secondary relative to primary in meters
        Eigen::Vector3d v, // velocity of secondary relative to primary in meters/second
        double mass_pri, // mass of the primary in kilograms
        double mass_sec, // mass of the secondary in kilograms
        OrbitalParams &params);

Eigen::Vector3d eccentricity_vector(
    Eigen::Vector3d r, // location of satellite to primary body
    Eigen::Vector3d v, // velocity of satellite relative to primary
    double mass_pri,
    double mass_sec);

void orbit_state(
    OrbitalParams &params,
    double mass_pri,
    double mass_sec,
    Eigen::Vector3d &r,
    Eigen::Vector3d &v);

void orbit_initialize(
        OrbitalParams &params,
        double mass_pri, // mass of the primary in kilograms
        double mass_sec, // mass of the secondary in kilograms
        Eigen::Vector3d &pos_pri, // position of the primary in meters
        Eigen::Vector3d &pos_sec, // position of the secondary in meters
        Eigen::Vector3d &vel_pri, // velocity of the primary in meters/second
        Eigen::Vector3d &vel_sec); // velocity of the secondary in meters/second

void polar_orbit_initialize(
        OrbitalParams &params,
        double mass_pri,
        double mass_sec,
        double &theta,
        double &theta_dot,
        double &r,
        double &r_dot);

void angular_limits(
    double theta1,
    double theta2,
    double &theta_init,
    double &theta_final);

double orbit_time(
        OrbitalParams &params,
        double mass_pri,
        double mass_sec,
        double nu1,
        double nu2);

void orbit_integrate(
    OrbitalParams &params0,
    OrbitalParams &params1,
    double mass_pri,
    double mass_sec,
    double delta_t);

Eigen::Matrix3d orbit_basis(OrbitalParams &params);

void orbit_draw( // draws the orbits of the primary and secondary bodies
        OrbitalParams &params,
        double mass_pri, // mass of the primary in kilograms
        double mass_sec); // mass of the secondary in kilograms

double semi_major(double T, double m1, double m2);

double period(double a, double m1, double m2);

double v_circular(double r, double m1, double m2);

Eigen::Vector3d direction_from_coord(
    double latitude, double longitude);

void coord_from_direction(
    Eigen::Vector3d dir,
    double &latitude,
    double &longitude);

Eigen::Matrix3d vehicle_basis(
    Eigen::Vector3d &r, // location of vehicle relative to orbited body
    Eigen::Vector3d &v);// velocity of vehicle relative to orbited body

double semi_major_from_periapsis_and_apoapsis(
    double r_periapsis, double r_apoapsis);

double eccentricity_from_periapsis_and_apoapsis(
    double r_periapsis, double r_apoapsis);

double nu_from_r(
        OrbitalParams params,
        Eigen::Vector3d r);

Eigen::Vector3d n_r_earth_frame(double latitude, double longitude);

Eigen::Vector3d orbit_n_plane(OrbitalParams const &params);

Eigen::Vector3d orbit_ascending_node(OrbitalParams const &params);

bool orbit_ring_intersect(
        OrbitalParams const &params,
        Eigen::Vector3d const &rs,
        Eigen::Vector3d &ri1,
        Eigen::Vector3d &ri2);

void ballistic_orbit_params(
        OrbitalParams &params,
        Eigen::Vector3d const &r_start,
        Eigen::Vector3d const &r_end,
        double radius_apogee);

Eigen::Matrix3d orbit_basis_90(
        Eigen::Matrix3d A_basis);


#endif
