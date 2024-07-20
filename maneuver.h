#ifndef MANEUVER_H
#define MANEUVER_H

#include <list>
#include "caams.hpp"
#include "orbit.h"
#include "sequencer.h"

class Satellite;

class Maneuver
{
public:
    double t;
    double T_duration;
    Maneuver();
    virtual ~Maneuver();
    virtual void Update(double dt);
    double Duration(void);
    virtual void ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,   // force in global coordinate space
            Eigen::Vector3d &torque)=0;
};

class CoastManeuver : public Maneuver
{
public:
    CoastManeuver(double T_duration);
    ~CoastManeuver();
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
};

class LaunchManeuver : public Maneuver
{
public:
    Satellite* satellite;
    Sequencer sequencer;
    Eigen::Vector3d h_norm;
public:
    LaunchManeuver(
            Satellite* satellite,
            double latitude,
            double longitude,
            double heading);
    ~LaunchManeuver();
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
    Eigen::Vector3d get_h_norm(void);
};

class ConstantAccelerationManeuver : public Maneuver
{
    Eigen::Vector3d force;
public:
    ConstantAccelerationManeuver(
            Eigen::Vector3d force,
            double T_duration);
    ~ConstantAccelerationManeuver();
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
};

class DirectedAccelerationManeuver : public Maneuver
{
public:
    Satellite *satellite;
    DirCoordSequencer sequencer;
    double thrust;
    DirectedAccelerationManeuver(
        Satellite *satellite,
        Eigen::MatrixXd control,
        Eigen::RowVectorXd time,
        double thrust);
    ~DirectedAccelerationManeuver(){}
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
};

class DirectedAcceleration3dManeuver : public Maneuver
{
public:
    Satellite *satellite;
    DirCoord3dSequencer sequencer;
    double thrust;
    DirectedAcceleration3dManeuver(
        Satellite *satellite,
        Eigen::MatrixXd control,
        Eigen::RowVectorXd time,
        double thrust);
    ~DirectedAcceleration3dManeuver(){}
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
};

class ThrustVectorManeuver : public Maneuver
{
public:
    Satellite *satellite;
    ThrustVectorSequencer sequencer;
    double thrust;
    Eigen::Vector3d rp;
    ThrustVectorManeuver(
        Satellite *satellite,
        Eigen::MatrixXd control,
        Eigen::RowVectorXd time,
        double thrust,
        Eigen::Vector3d rp);
    ~ThrustVectorManeuver(){}
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
};

class SurfaceManeuver : public Maneuver
{
public:
    double latitude;
    double longitude;
    Satellite *satellite;
    SurfaceManeuver(
            Satellite *satellite,
            double latitude,
            double longitude,
            double T_duration);
    ~SurfaceManeuver(){}
    void Update(double dt);
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
};

class ElevatorManeuver : public Maneuver
{
    Satellite *satellite;
    double thrust;
public:
    ElevatorManeuver(
            Satellite *satellite,
            double thrust,
            double T_duration);
    ~ElevatorManeuver(){}
    void ForceAndTorque(
                double dt,
                Eigen::Vector3d &force,   // force in global coordinate space
                Eigen::Vector3d &torque);
};

enum ReentryTestState
{
    RTS_INIT,
    RTS_COARSE,
    RTS_FINE
};

class ReentryTestManeuver : public Maneuver
{
    ReentryTestState state;
    double radius_target;
    Satellite *satellite;
    double &delta_t;
    double &T_reentry;
public:
    ReentryTestManeuver(
            Satellite *satellite,
            double radius_target,
            double &delta_t,
            double &T_reentry);
    ~ReentryTestManeuver(){}
    void Update(double dt);
    void ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,   // force in global coordinate space
            Eigen::Vector3d &torque);
};

class ManeuverQueue
{
    std::list<Maneuver*> maneuvers;
public:
    void AddManeuver(Maneuver* maneuver);
    void DeleteCurrent(void);
    Maneuver* CurrentManeuver(void);
    void Update(double dt);
    bool Empty(void);
};

#endif
