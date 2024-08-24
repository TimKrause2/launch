#ifndef SATELLITE_H
#define SATELLITE_H

#include "gsim.h"
#include "event.h"
#include "maneuver.h"
#include "orbit.h"

#define THRUST (9.819609*2.0)


class Satellite : public Body
{
public:
    Body* primary_body;
    ManeuverQueue maneuverQueue;
    EventQueue eventQueue;
    Eigen::Vector3d rp; // attachment point for the rocket engine
public:
    Satellite(Satellite &s, Body *primary_body);
    Satellite(
        double mass,
        double radius,
        Eigen::Vector3d &position,
        Eigen::Vector3d &velocity,
        Eigen::Vector4d &p,
        Eigen::Vector4d &pdot,
        Eigen::Matrix3d &Jp,
        char* tex_image,
        Body* primary_body,
        Eigen::Vector3d rp);
    void Prepare(void);
    void Update(double dt);
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
    double TimeStep(void);

    TimeSpec GetCurrentTime(void);

    Eigen::Vector3d RelativeVelocity(void);
    Eigen::Vector3d RelativePosition(void);
    double RadialVelocity(void);
    double Altitude(void);

    void NewOrbit(
            OrbitalParams &params);
    void OrbitalElements(
            OrbitalParams &params);
    Eigen::Vector3d ScheduleLaunchManeuver(
            double latitude,
            double longitude,
            double heading);
    void ScheduleCircularizeAtRadius(double radius_target);
    void ScheduleCircularizeAtPeriapsis(void);
    void ScheduleAdjustAtPariapsis(double apoapsis);
    void ScheduleCircularizeAtApoapsis(void);
    void ScheduleAdjustAtApoapsis(double periapsis);
    void ScheduleInclinationManeuver(
            double nu_axis,
            double theta);
    void ScheduleDirectedAccelerationManeuver(
        Eigen::MatrixXd control,
        Eigen::MatrixXd time,
        double thrust,
        double T_coast);
    void ScheduleConstantAccelerationManeuver(
        double latitude, double longitude,
        Eigen::RowVectorXd &time,
        double thrust,
        double T_coast);
    void ScheduleThrustVectorLaunch(void);
    void ScheduleDelayedLaunch(
            double latitude, double longitude, // location on the primary body
            OrbitalParams const &params, // orbital parameters of the target orbit
            double thrust);
    bool ScheduleBallisticLaunch(
            double lat_vehicle, double long_vehicle,
            double lat_target, double long_target,
            double T1, double T2,
            double thrust);
    void ScheduleICBMLaunchTest(
            Eigen::MatrixXd &controls,
            Eigen::MatrixXd &time,
            double thrust,
            double radius_target,
            double &delta_t,
            double &T_reentry);

    bool OptimizeICBMLaunch(
            Eigen::Vector3d const &r_launch,
            Eigen::Vector3d const &v_launch,
            Eigen::Vector3d const &r_aim,
            double h_apgogee,
            double T1,
            double thrust,
            Eigen::MatrixXd &controls,
            Eigen::MatrixXd &time,
            Eigen::Vector3d &r_strike,
            double &T_reentry);

    bool ScheduleICBMLaunch(
            double lat_vehicle, double long_vehicle,
            double hgt_vehicle,
            double lat_target, double long_target,
            double hgt_target);
};

#endif
