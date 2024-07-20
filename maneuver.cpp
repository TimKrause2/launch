#include "maneuver.h"
#include "orbit.h"
#include "gsim.h"
#include "satellite.h"
#include <cmath>
#include <iomanip>
#include <iostream>

Maneuver::Maneuver(void)
{
    t=0.0;
}

Maneuver::~Maneuver(){}

void Maneuver::Update(double dt)
{
    t+=dt;
}

double Maneuver::Duration(void)
{
    return T_duration;
}

CoastManeuver::CoastManeuver(double T_duration):
    Maneuver()
{
    Maneuver::T_duration = T_duration;
}

CoastManeuver::~CoastManeuver()
{}

void CoastManeuver::ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,
            Eigen::Vector3d &torque)
{
    force = Eigen::Vector3d::Zero();
    torque = Eigen::Vector3d::Zero();
}

LaunchManeuver::LaunchManeuver(
    Satellite *satellite,
    double latitude,
    double longitude,
    double heading):
    Maneuver(),
    satellite(satellite)
{
    Eigen::Vector3d y_axis;
    y_axis << 0.0, 1.0, 0.0;
    Eigen::Vector3d z_axis;
    z_axis << 0.0, 0.0, 1.0;
    Eigen::Matrix3d rot_lat = caams::AAA(latitude,-y_axis);
    Eigen::Matrix3d rot_long = caams::AAA(longitude,z_axis);
    Eigen::Matrix3d rot_pos = rot_long*rot_lat;
    Eigen::Vector3d n_r = rot_pos.col(0); // x-axis
    Eigen::Vector3d n_theta = rot_pos.col(2)*cos(heading)+
                              rot_pos.col(1)*sin(heading);
    h_norm = n_r.cross(n_theta);
    satellite->m_position = n_r*(satellite->primary_body->m_radius+satellite->m_radius) + satellite->primary_body->m_position;
    satellite->m_velocity = satellite->primary_body->m_velocity;

    T_duration = sequencer.TotalTime();
}

LaunchManeuver::~LaunchManeuver()
{}

void LaunchManeuver::ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,
            Eigen::Vector3d &torque)
{
    Eigen::Vector3d r_sat_earth = satellite->m_rk_position - satellite->primary_body->m_rk_position;
    Eigen::Vector3d n_r = r_sat_earth.normalized();
    Eigen::Vector3d n_theta = h_norm.cross(n_r);

    double beta = sequencer.GetBeta(t+dt);

    Eigen::Vector3d n_acc = n_r*cos(beta) + n_theta*sin(beta);

    force = n_acc*THRUST*satellite->m_mass;
    torque = Eigen::Vector3d::Zero();
}

Eigen::Vector3d LaunchManeuver::get_h_norm(void){
    return h_norm;
}

ConstantAccelerationManeuver::ConstantAccelerationManeuver(Eigen::Vector3d force, double T_duration)
    :force(force)
{
    Maneuver::T_duration = T_duration;
}

ConstantAccelerationManeuver::~ConstantAccelerationManeuver(){}

void ConstantAccelerationManeuver::ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,
            Eigen::Vector3d &torque){
    force = ConstantAccelerationManeuver::force;
    torque = Eigen::Vector3d::Zero();
}

DirectedAccelerationManeuver::DirectedAccelerationManeuver(
    Satellite *satellite,
    Eigen::MatrixXd control,
    Eigen::RowVectorXd time,
    double thrust)
    : satellite(satellite),
    sequencer(control, time),
    thrust(thrust)
{
    T_duration = sequencer.TotalTime();
}

void DirectedAccelerationManeuver::ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,
            Eigen::Vector3d &torque)
{
    Eigen::Vector3d r = satellite->m_rk_position - satellite->primary_body->m_rk_position;
    Eigen::Vector3d v = satellite->m_rk_velocity - satellite->primary_body->m_rk_velocity;
    Eigen::Matrix3d A_vehicle = vehicle_basis(r,v);
    Eigen::Vector3d n_acc = sequencer.GetDir(t + dt);
    force = A_vehicle*n_acc*thrust*satellite->m_mass;
    torque = Eigen::Vector3d::Zero();
}


DirectedAcceleration3dManeuver::DirectedAcceleration3dManeuver(
    Satellite *satellite,
    Eigen::MatrixXd control,
    Eigen::RowVectorXd time,
    double thrust)
    : satellite(satellite),
    sequencer(control, time),
    thrust(thrust)
{
    T_duration = sequencer.TotalTime();
}

void DirectedAcceleration3dManeuver::ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,
            Eigen::Vector3d &torque)
{
    Eigen::Vector3d n_acc = sequencer.GetDir(t + dt);
    force = n_acc*thrust*satellite->m_mass;
    torque = Eigen::Vector3d::Zero();
}


ThrustVectorManeuver::ThrustVectorManeuver(
    Satellite *satellite,
    Eigen::MatrixXd control,
    Eigen::RowVectorXd time,
    double thrust,
    Eigen::Vector3d rp)
    : satellite(satellite),
    sequencer(control, time),
    thrust(thrust),
    rp(rp)
{
    T_duration = sequencer.TotalTime();
}

void ThrustVectorManeuver::ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,
            Eigen::Vector3d &torque)
{
    Eigen::Vector3d Fb = sequencer.GetDir(t + dt) * thrust * satellite->m_mass;
    torque = rp.cross(Fb);
    Eigen::Matrix3d A = caams::Ap(satellite->rk_p);
    force = A*Fb;
}

SurfaceManeuver::SurfaceManeuver(
        Satellite *satellite,
        double latitude,
        double longitude,
        double T_duration):
    satellite(satellite),
    latitude(latitude),
    longitude(longitude)
{
    Maneuver::T_duration = T_duration;
}

void SurfaceManeuver::Update(double dt)
{
    Maneuver::Update(dt);

    // place the satellite on the surface of the primaray body
    Eigen::Vector3d n_r_primary = n_r_earth_frame(latitude, longitude);
    Eigen::Vector3d r_primary = n_r_primary *
            (satellite->primary_body->m_radius +
             satellite->m_radius);
    Eigen::Matrix3d A = caams::Ap(satellite->primary_body->p);
    Eigen::Vector3d r_global = A*r_primary;
    satellite->m_position = satellite->primary_body->m_position + r_global;

    // calculate the velocity
    Eigen::Vector3d omega = caams::omega_p_dot(satellite->primary_body->p,
                                               satellite->primary_body->pdot);
    Eigen::Vector3d v_global = omega.cross(r_global);
    satellite->m_velocity = satellite->primary_body->m_velocity + v_global;
}

void SurfaceManeuver::ForceAndTorque(
        double dt,
        Eigen::Vector3d &force,   // force in global coordinate space
        Eigen::Vector3d &torque)
{
    // balalnce the force of gravity from the primary body
    double mu = G_gravity * satellite->primary_body->m_mass;
    Eigen::Vector3d r = satellite->m_rk_position - satellite->primary_body->m_rk_position;
    double a_gravity = mu/r.squaredNorm();
    force = r.normalized()*a_gravity;

    // add the force due to being constrained to the surface
    Eigen::Vector3d omega = caams::omega_p_dot(
                satellite->primary_body->rk_p,
                satellite->primary_body->rk_pdot);
    Eigen::Vector3d v = omega.cross(r);
    force += omega.cross(v);
    force *= satellite->m_mass;

    torque = Eigen::Vector3d::Zero();
}

ElevatorManeuver::ElevatorManeuver(
        Satellite *satellite,
        double thrust,
        double T_duration):
    satellite(satellite),
    thrust(thrust)
{
    Maneuver::T_duration = T_duration;
}

void ElevatorManeuver::ForceAndTorque(
        double dt,
        Eigen::Vector3d &force,   // force in global coordinate space
        Eigen::Vector3d &torque)
{
    Eigen::Vector3d r = satellite->m_rk_position - satellite->primary_body->m_rk_position;
    force = r.normalized()*thrust*satellite->m_mass;

    torque = Eigen::Vector3d::Zero();
}

ReentryTestManeuver::ReentryTestManeuver(
        Satellite *satellite,
        double radius_target,
        double &delta_t,
        double &T_reentry)
    :
      satellite(satellite),
      radius_target(radius_target),
      delta_t(delta_t),
      T_reentry(T_reentry)
{
    state = RTS_INIT;
}

void ReentryTestManeuver::Update(double dt)
{
    Maneuver::Update(dt);
    if(state==RTS_INIT){
        delta_t = 1.0;
        state = RTS_COARSE;
    }else if(state==RTS_COARSE){
        double dr = satellite->RadialVelocity();
        if(dr<0.0){
            double delta_h = -dr*delta_t;
            double r_altitude = satellite->RelativePosition().norm();
            if((r_altitude-radius_target)<(delta_h*2.0)){
                state = RTS_FINE;
                delta_t = 0.01;
            }
        }
    }else if(state==RTS_FINE){
        double dr = satellite->RadialVelocity();
        double delta_h = -dr*delta_t;
        double r_altitude = satellite->RelativePosition().norm();
        double delta_a = r_altitude - radius_target;
        if(delta_a<delta_h){
            double delta_t = delta_a/-dr;
            T_reentry = t + delta_t;
            satellite->maneuverQueue.DeleteCurrent();
        }
    }
}

void ReentryTestManeuver::ForceAndTorque(
        double dt,
        Eigen::Vector3d &force,   // force in global coordinate space
        Eigen::Vector3d &torque)
{
    force = Eigen::Vector3d::Zero();
    torque = Eigen::Vector3d::Zero();
}


void ManeuverQueue::AddManeuver(Maneuver *maneuver)
{
    maneuvers.push_back(maneuver);
}

void ManeuverQueue::DeleteCurrent()
{
    if(maneuvers.empty()) return;
    delete maneuvers.front();
    maneuvers.pop_front();
}

Maneuver* ManeuverQueue::CurrentManeuver()
{
    if(maneuvers.empty()){
        return NULL;
    }else{
        return maneuvers.front();
    }
}

void ManeuverQueue::Update(double dt)
{
    if(maneuvers.empty()) return;
    maneuvers.front()->Update(dt);
}

bool ManeuverQueue::Empty()
{
    return maneuvers.empty();
}
