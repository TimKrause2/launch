#include "satellite.h"
#include "two_phase.h"
#include "thrust_vector.h"
#include "delayed_launch.h"
#include "ballistic.h"
#include "spheroid.h"
#include "icbm_simple.h"
#include <cmath>

#define RAD_PER_DEG (M_PI/180.0)

Satellite::Satellite(Satellite &s, Body *primary_body)
    :
      Body(
          s.m_mass,
          s.m_radius,
          s.m_position,
          s.m_velocity,
          s.p,
          s.pdot,
          s.Jp,
          NULL),
      primary_body(primary_body),
      rp(s.rp)
{}

Satellite::Satellite(double mass,
        double radius,
        Eigen::Vector3d &position,
        Eigen::Vector3d &velocity,
        Eigen::Vector4d &p,
        Eigen::Vector4d &pdot,
        Eigen::Matrix3d &Jp,
        char* tex_image,
        Body* primary_body,
        Eigen::Vector3d rp):
    Body(
        mass,
        radius,
        position,
        velocity,
        p,
        pdot,
        Jp,
        tex_image),
    primary_body(primary_body),
    rp(rp)
{}

void Satellite::Prepare(void)
{
    eventQueue.Check();
}

void Satellite::Update(double dt)
{
    maneuverQueue.Update(dt);
    eventQueue.Update(dt);
}

void Satellite::ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,   // force in global coordinate space
            Eigen::Vector3d &torque)
{
    Maneuver* maneuver = maneuverQueue.CurrentManeuver();
    if(maneuver){
        maneuver->ForceAndTorque(dt, force, torque);
    }else{
        force = Eigen::Vector3d::Zero();
        torque = Eigen::Vector3d::Zero();
    }
}

double Satellite::TimeStep()
{
    if(eventQueue.Empty()){
        return CurvatureTimeStep();
    }else{
        return fmin(CurvatureTimeStep(),eventQueue.NextEvent());
    }
}

TimeSpec Satellite::GetCurrentTime(void)
{
    return eventQueue.AbsoluteTime(0.0);
}

Eigen::Vector3d Satellite::RelativeVelocity(void)
{
    return m_velocity - primary_body->m_velocity;
}

Eigen::Vector3d Satellite::RelativePosition(void)
{
    return m_position - primary_body->m_position;
}

double Satellite::RadialVelocity(void)
{
    Eigen::Vector3d r = RelativePosition();
    Eigen::Vector3d v = RelativeVelocity();
    return r.dot(v)/r.norm();
}

double Satellite::Altitude(void)
{
    Eigen::Vector3d r = RelativePosition();
    return r.norm() - primary_body->m_radius;
}

void Satellite::NewOrbit(OrbitalParams &params)
{
    if(!maneuverQueue.Empty()) return;
    orbit_initialize(
                params,
                primary_body->m_mass,
                m_mass,
                primary_body->m_position,m_position,
                primary_body->m_velocity,m_velocity);

}


Eigen::Vector3d Satellite::ScheduleLaunchManeuver(double latitude, double longitude, double heading)
{
    LaunchManeuver* maneuver = new LaunchManeuver(
                this,
                latitude,
                longitude,
                heading);

    double T_maneuver = maneuver->Duration();

    TimeSpec ts_sat_end = eventQueue.AbsoluteTime(T_maneuver);
    Event* event = new SatelliteEvent(this,ts_sat_end);

    int Npoints = maneuver->sequencer.GetNpoints();
    for(int i = 0;i < Npoints; i++){
        TimeSpec ts(maneuver->sequencer.GetTime(i));
        Event* ievent = new IntegratorEvent(ts);
        eventQueue.AddEvent(ievent);
    }

    maneuverQueue.AddManeuver(maneuver);
    eventQueue.AddEvent(event);

    return maneuver->get_h_norm();
}

void Satellite::OrbitalElements(OrbitalParams &params)
{
    Eigen::Vector3d r = m_position - primary_body->m_position;
    Eigen::Vector3d v = m_velocity - primary_body->m_velocity;
    orbital_elements(r,v,primary_body->m_mass,m_mass,params);
}

void Satellite::ScheduleCircularizeAtRadius(double radius_target)
{
    if(!maneuverQueue.Empty()) return;
    OrbitalParams params_now;
    OrbitalParams params_target;
    OrbitalParams params_new;
    OrbitalElements(params_now);
    double r_periapsis = params_now.a*(1 - params_now.e);
    double r_apoapsis = params_now.a*(1 + params_now.e);
    if(radius_target>r_apoapsis)
        return;
    if(radius_target<r_periapsis)
        return;
    params_target = params_now;
    params_target.nu = acos((params_now.a*(1-params_now.e*params_now.e)-radius_target)/radius_target/params_now.e);
    double nu = params_target.nu/RAD_PER_DEG;
    params_new = params_now;
    params_new.e = 1e-9;
    params_new.a = radius_target/(1 - params_new.e);
    params_new.omega = params_target.omega + params_target.nu;
    params_new.nu = 0.0;

    Eigen::Vector3d r_target;
    Eigen::Vector3d v_target;
    orbit_state(
        params_target,
        primary_body->m_mass, m_mass,
        r_target, v_target);

    Eigen::Vector3d r_new;
    Eigen::Vector3d v_new;
    orbit_state(
        params_new,
        primary_body->m_mass, m_mass,
        r_new, v_new);

    Eigen::Vector3d d_r = r_new - r_target;
    if(d_r.norm()>1e-3){
        std::cout << "ScheduleCircularizeAtRadius: orbits out of allignment." << std::endl;
        return;
    }

    Eigen::Vector3d delta_v = v_new - v_target;
    Eigen::Matrix3d A = vehicle_basis(r_target, v_target);
    Eigen::Vector3d delta_v_vehicle = A.transpose()*delta_v;
    double latitude;
    double longitude;
    coord_from_direction(
        delta_v_vehicle, latitude, longitude);

    double T_coast;
    if(params_now.nu == params_target.nu)
    {
        T_coast = 0.0;
    }else{
        T_coast = orbit_time(
            params_now,
            primary_body->m_mass,
            m_mass,
            params_now.nu,
            params_target.nu);
    }
    double T_period = period(params_now.a, m_mass, primary_body->m_mass);
    if(T_coast < T_period/128){
        T_coast += T_period;
    }
    int Nphase1_nodes = (int)(T_coast*256 / T_period);
    if(Nphase1_nodes < 2){
        Nphase1_nodes = 2;
    }

    double T_maneuver = delta_v.norm()/THRUST;

    Eigen::MatrixXd phase1_time = linspace(0.0, T_coast, Nphase1_nodes);
    Eigen::MatrixXd phase2_time = linspace(0.0, T_maneuver, 10);
    Eigen::MatrixXd phase2_control(2,10);
    phase2_control.row(0) = ones(1,10)*latitude; // latitude
    phase2_control.row(1) = ones(1,10)*longitude; // longitude

    if(T_maneuver > 0.5)
    {
        ModEquOrbitalParams meparams0;
        ModEquOrbitalParams meparams1;
        meparams0.FromOrbitalParams(params_now);
        meparams1.FromOrbitalParams(params_new);

        double mu = G_gravity*(primary_body->m_mass + m_mass);

        two_phase_optimize_constant_longitude(
            meparams0,
            meparams1,
            phase1_time,
            phase2_time,
            phase2_control,
            mu,
            THRUST);

        T_coast = phase1_time(phase1_time.cols()-1);
    }

    ScheduleDirectedAccelerationManeuver(
        phase2_control,
        phase2_time,
        THRUST,
        T_coast);
}





void Satellite::ScheduleCircularizeAtPeriapsis(void)
{
    if(!maneuverQueue.Empty()) return;
    OrbitalParams params_now;
    OrbitalElements(params_now);
    OrbitalParams params_peri;
    OrbitalParams params_opt;
    OrbitalParams params_circ;
    params_peri = params_now;
    params_circ = params_now;
    params_peri.nu = 0.0;

    Eigen::Vector3d r_peri;
    Eigen::Vector3d v_peri;
    orbit_state(
        params_peri,
        primary_body->m_mass, m_mass,
        r_peri, v_peri);

    double mag_r_peri = r_peri.norm();

    // match the periapsis of a near circular orbit to the
    // original periapsis

    params_circ.e = 1.0e-9;
    params_circ.a = mag_r_peri/(1 - params_circ.e);
    params_circ.nu = 0.0;

    Eigen::Vector3d r_circ;
    Eigen::Vector3d v_circ;
    orbit_state(
        params_circ,
        primary_body->m_mass, m_mass,
        r_circ, v_circ);

    Eigen::Vector3d dr = r_circ - r_peri;
    if(dr.norm() > 1e-3){
        std::cout << "ScheduleCircularizeAtPeriapsis: orbits don't intersect." << std::endl;
        return;
    }

    Eigen::Vector3d delta_v = v_circ - v_peri;
    Eigen::Matrix3d A = vehicle_basis( r_peri, v_peri);
    Eigen::Vector3d delta_v_vehicle = A.transpose()*delta_v;
    double latitude;
    double longitude;
    coord_from_direction(delta_v_vehicle, latitude, longitude);

    double T_maneuver = delta_v.norm()/THRUST;
    double T_coast = 3.0*T_maneuver;
    double T_target = orbit_time(params_now,
                                 primary_body->m_mass,
                                 m_mass,
                                 params_now.nu,
                                 params_peri.nu);;
    double T_opt = T_target - T_coast;
    if(T_opt<0.0){
        return;
    }

    orbit_integrate(params_peri,params_opt,
                    primary_body->m_mass, m_mass, -T_coast);

    Eigen::MatrixXd phase1_time = linspace(0.0, T_coast, 10);
    Eigen::MatrixXd phase2_time = linspace(0.0, T_maneuver, 10);
    Eigen::MatrixXd phase2_control(2,10);
    phase2_control.row(0) = ones(1,10)*latitude; // latitude
    phase2_control.row(1) = ones(1,10)*longitude; // longitude

    ModEquOrbitalParams meparams0;
    ModEquOrbitalParams meparams1;
    meparams0.FromOrbitalParams(params_opt);
    meparams1.FromOrbitalParams(params_circ);

    double mu = G_gravity*(primary_body->m_mass + m_mass);

    two_phase_optimize_constant_longitude(
        meparams0,
        meparams1,
        phase1_time,
        phase2_time,
        phase2_control,
        mu,
        THRUST);

    T_coast = T_opt + phase1_time(phase1_time.cols()-1);

    ScheduleDirectedAccelerationManeuver(
        phase2_control,
        phase2_time,
        THRUST,
        T_coast);
}

void Satellite::ScheduleAdjustAtPariapsis(double apoapsis)
{
    if(!maneuverQueue.Empty()) return;
    OrbitalParams params_now;
    OrbitalParams params_peri;
    OrbitalParams params_new;
    OrbitalElements(params_now);
    double mag_r_peri = params_now.a*(1.0-params_now.e);
    if(apoapsis < mag_r_peri )
        return;
    params_peri = params_now;
    params_peri.nu = 0.0;
    Eigen::Vector3d r_peri;
    Eigen::Vector3d v_peri;
    orbit_state(
        params_peri,
        primary_body->m_mass, m_mass,
        r_peri, v_peri);

    params_new = params_now;
    params_new.a = semi_major_from_periapsis_and_apoapsis(
        mag_r_peri, apoapsis);
    params_new.e = eccentricity_from_periapsis_and_apoapsis(
        mag_r_peri, apoapsis);
    params_new.nu = 0.0;

    Eigen::Vector3d r_new;
    Eigen::Vector3d v_new;
    orbit_state(
        params_new,
        primary_body->m_mass, m_mass,
        r_new, v_new);

    Eigen::Vector3d delta_v = v_new - v_peri;
    Eigen::Matrix3d A = vehicle_basis(r_peri, v_peri);
    Eigen::Vector3d delta_v_vehicle = A.transpose()*delta_v;
    double latitude;
    double longitude;
    coord_from_direction(
        delta_v_vehicle, latitude, longitude);

    double T_coast;
    if(params_now.nu == params_peri.nu)
    {
        T_coast = 0.0;
    }else{
        T_coast = orbit_time(
            params_now,
            primary_body->m_mass,
            m_mass,
            params_now.nu,
            params_peri.nu);
    }
    double T_period = period(params_now.a, m_mass, primary_body->m_mass);
    if(T_coast < T_period/128){
        T_coast += T_period;
    }
    int Nphase1_nodes = (int)(T_coast*256 / T_period);
    if(Nphase1_nodes < 2){
        Nphase1_nodes = 2;
    }
    Nphase1_nodes = 10;

    double T_maneuver = delta_v.norm()/THRUST;

    Eigen::MatrixXd phase1_time = linspace(0.0, T_coast, Nphase1_nodes);
    Eigen::MatrixXd phase2_time = linspace(0.0, T_maneuver, 10);
    Eigen::MatrixXd phase2_control(2,10);
    phase2_control.row(0) = ones(1,10)*latitude; // latitude
    phase2_control.row(1) = ones(1,10)*longitude; // longitude

    ModEquOrbitalParams meparams0;
    ModEquOrbitalParams meparams1;
    meparams0.FromOrbitalParams(params_now);
    meparams1.FromOrbitalParams(params_new);

    double mu = G_gravity*(primary_body->m_mass + m_mass);

    two_phase_optimize_constant_longitude(
        meparams0,
        meparams1,
        phase1_time,
        phase2_time,
        phase2_control,
        mu,
        THRUST);

    T_coast = phase1_time(phase1_time.cols()-1);

    ScheduleDirectedAccelerationManeuver(
        phase2_control,
        phase2_time,
        THRUST,
        T_coast);
}

void Satellite::ScheduleCircularizeAtApoapsis(void)
{
    if(!maneuverQueue.Empty()) return;
    OrbitalParams params_now;
    OrbitalElements(params_now);
    OrbitalParams params_apo;
    OrbitalParams params_circ;
    params_apo = params_now;
    params_circ = params_now;
    params_apo.nu = M_PI;

    Eigen::Vector3d r_apo;
    Eigen::Vector3d v_apo;
    orbit_state(
        params_apo,
        primary_body->m_mass, m_mass,
        r_apo, v_apo);

    double mag_r_apo = r_apo.norm();

    // match the periapsis of a near circular orbit to the
    // original periapsis

    params_circ.e = 1.0e-9;
    params_circ.a = mag_r_apo/(1 + params_circ.e);
    params_circ.nu = M_PI;

    Eigen::Vector3d r_circ;
    Eigen::Vector3d v_circ;
    orbit_state(
        params_circ,
        primary_body->m_mass, m_mass,
        r_circ, v_circ);

    Eigen::Vector3d dr = r_circ - r_apo;
    if(dr.norm() > 1e-3){
        std::cout << "ScheduleCircularizeAtPeriapsis: orbits don't intersect." << std::endl;
        return;
    }

    Eigen::Vector3d delta_v = v_circ - v_apo;
    Eigen::Matrix3d A = vehicle_basis( r_apo, v_apo);
    Eigen::Vector3d delta_v_vehicle = A.transpose()*delta_v;
    double latitude;
    double longitude;
    coord_from_direction(delta_v_vehicle, latitude, longitude);


    double T_maneuver = delta_v.norm()/THRUST;
    double T_coast;
    if(params_now.nu == params_apo.nu)
    {
        T_coast = 0.0;
    }else{
        T_coast = orbit_time(
            params_now,
            primary_body->m_mass,
            m_mass,
            params_now.nu,
            params_apo.nu);
    }
    double T_period = period(params_now.a, m_mass, primary_body->m_mass);
    if(T_coast < T_period/128){
        T_coast += T_period;
    }
    int Nphase1_nodes = (int)(T_coast*256 / T_period);
    if(Nphase1_nodes < 2){
        Nphase1_nodes = 2;
    }

    Eigen::MatrixXd phase1_time = linspace(0.0, T_coast, Nphase1_nodes);
    Eigen::MatrixXd phase2_time = linspace(0.0, T_maneuver, 10);
    Eigen::MatrixXd phase2_control(2,10);
    phase2_control.row(0) = ones(1,10)*latitude; // latitude
    phase2_control.row(1) = ones(1,10)*longitude; // longitude

    ModEquOrbitalParams meparams0;
    ModEquOrbitalParams meparams1;
    meparams0.FromOrbitalParams(params_now);
    meparams1.FromOrbitalParams(params_circ);

    double mu = G_gravity*(primary_body->m_mass + m_mass);

    two_phase_optimize_constant_longitude(
        meparams0,
        meparams1,
        phase1_time,
        phase2_time,
        phase2_control,
        mu,
        THRUST);

    T_coast = phase1_time(phase1_time.cols()-1);

    ScheduleDirectedAccelerationManeuver(
        phase2_control,
        phase2_time,
        THRUST,
        T_coast);
}

void Satellite::ScheduleAdjustAtApoapsis(double periapsis)
{
    if(!maneuverQueue.Empty()) return;
    OrbitalParams params_now;
    OrbitalParams params_apo;
    OrbitalParams params_new;
    OrbitalElements(params_now);
    double mag_r_apo = params_now.a*(1.0 + params_now.e);
    if(periapsis > mag_r_apo )
        return;
    params_apo = params_now;
    params_apo.nu = M_PI;
    Eigen::Vector3d r_apo;
    Eigen::Vector3d v_apo;
    orbit_state(
        params_apo,
        primary_body->m_mass, m_mass,
        r_apo, v_apo);

    params_new = params_now;
    params_new.a = semi_major_from_periapsis_and_apoapsis(
        periapsis, mag_r_apo);
    params_new.e = eccentricity_from_periapsis_and_apoapsis(
        periapsis, mag_r_apo);
    params_new.nu = M_PI;

    Eigen::Vector3d r_new;
    Eigen::Vector3d v_new;
    orbit_state(
        params_new,
        primary_body->m_mass, m_mass,
        r_new, v_new);

    Eigen::Vector3d delta_v = v_new - v_apo;
    Eigen::Matrix3d A = vehicle_basis(r_apo, v_apo);
    Eigen::Vector3d delta_v_vehicle = A.transpose()*delta_v;
    double latitude;
    double longitude;
    coord_from_direction(
        delta_v_vehicle, latitude, longitude);

    double T_coast;
    if(params_now.nu == params_apo.nu)
    {
        T_coast = 0.0;
    }else{
        T_coast = orbit_time(
            params_now,
            primary_body->m_mass,
            m_mass,
            params_now.nu,
            params_apo.nu);
    }
    double T_period = period(params_now.a, m_mass, primary_body->m_mass);
    if(T_coast < T_period/128){
        T_coast += T_period;
    }
    int Nphase1_nodes = (int)(T_coast*256 / T_period);
    if(Nphase1_nodes < 2){
        Nphase1_nodes = 2;
    }

    double T_maneuver = delta_v.norm()/THRUST;

    Eigen::MatrixXd phase1_time = linspace(0.0, T_coast, Nphase1_nodes);
    Eigen::MatrixXd phase2_time = linspace(0.0, T_maneuver, 10);
    Eigen::MatrixXd phase2_control(2,10);
    phase2_control.row(0) = ones(1,10)*latitude; // latitude
    phase2_control.row(1) = ones(1,10)*longitude; // longitude

    ModEquOrbitalParams meparams0;
    ModEquOrbitalParams meparams1;
    meparams0.FromOrbitalParams(params_now);
    meparams1.FromOrbitalParams(params_new);

    double mu = G_gravity*(primary_body->m_mass + m_mass);

    two_phase_optimize_constant_longitude(
        meparams0,
        meparams1,
        phase1_time,
        phase2_time,
        phase2_control,
        mu,
        THRUST);

    T_coast = phase1_time(phase1_time.cols()-1);

    ScheduleDirectedAccelerationManeuver(
        phase2_control,
        phase2_time,
        THRUST,
        T_coast);
}

void Satellite::ScheduleInclinationManeuver(
        double nu_axis,
        double theta)
{
    if(!maneuverQueue.Empty()) return;

    OrbitalParams params_now;
    OrbitalParams params_axis;
    OrbitalParams params_new;

    // get the current orbital elements
    OrbitalElements(params_now);

    params_axis = params_now;
    params_axis.nu = nu_axis;

    Eigen::Vector3d r_axis;
    Eigen::Vector3d v_axis;

    orbit_state(
        params_axis,
        primary_body->m_mass, m_mass,
        r_axis, v_axis);

    Eigen::Vector3d r_new = r_axis;
    Eigen::Vector3d v_new;

    // rotate about the r_new vector
    Eigen::Matrix A_rot = caams::AAA(theta, r_new);

    v_new = A_rot*v_axis;

    orbital_elements(
        r_new, v_new,
        primary_body->m_mass, m_mass,
        params_new);

    // calculate the angle between the axis velocity vector
    // and the new vector
    Eigen::Vector3d n_v_new = v_new.normalized();
    Eigen::Vector3d n_v_axis = v_axis.normalized();
    double theta_v = acos(n_v_new.dot(n_v_axis));

    // calculate the rate of rotation of the velocity vector
    double omega_v = THRUST/v_new.norm();

    // calculate the acceleration duration
    double T_maneuver = theta_v/omega_v;

    double T_coast;
    if(params_now.nu == params_axis.nu)
    {
        T_coast = 0.0;
    }else{
        T_coast = orbit_time(
            params_now,
            primary_body->m_mass,
            m_mass,
            params_now.nu,
            params_axis.nu);
    }
    double T_period = period(params_now.a, m_mass, primary_body->m_mass);
    if(T_coast < T_period/128){
        T_coast += T_period;
    }
    int Nphase1_nodes = (int)(T_coast*256 / T_period);
    if(Nphase1_nodes < 2){
        Nphase1_nodes = 2;
    }

    double latitude = 0.0;
    double longitude = (theta>0.0)?M_PI/2.0:-M_PI/2.0;

    Eigen::MatrixXd phase1_time = linspace(0.0, T_coast, Nphase1_nodes);
    Eigen::MatrixXd phase2_time = linspace(0.0, T_maneuver, 10);
    Eigen::MatrixXd phase2_control(2,10);
    phase2_control.row(0) = ones(1,10)*latitude; // latitude
    phase2_control.row(1) = ones(1,10)*longitude; // longitude

    ModEquOrbitalParams meparams0;
    ModEquOrbitalParams meparams1;
    meparams0.FromOrbitalParams(params_now);
    meparams1.FromOrbitalParams(params_new);

    double mu = G_gravity*(primary_body->m_mass + m_mass);

    two_phase_optimize(
        meparams0,
        meparams1,
        phase1_time,
        phase2_time,
        phase2_control,
        mu,
        THRUST);

    T_coast = phase1_time(phase1_time.cols()-1);

    ScheduleDirectedAccelerationManeuver(
        phase2_control,
        phase2_time,
        THRUST,
        T_coast);
}

void Satellite::ScheduleDirectedAccelerationManeuver(
    Eigen::MatrixXd control,
    Eigen::MatrixXd time,
    double thrust,
    double T_coast)
{
    // make a directed acceleration maneuver
    DirectedAccelerationManeuver *acc_maneuver =
        new DirectedAccelerationManeuver(
            this, control, time, thrust);

    double T_acc_maneuver = acc_maneuver->Duration();

    if(T_coast!=0.0){
        CoastManeuver *coast_maneuver =
            new CoastManeuver(T_coast);
        TimeSpec ts_coast_end = eventQueue.AbsoluteTime(T_coast);
        Event* coast_event = new SatelliteEvent(this, ts_coast_end);
        maneuverQueue.AddManeuver(coast_maneuver);
        eventQueue.AddEvent(coast_event);
    }

    TimeSpec ts_acc_end = eventQueue.AbsoluteTime(T_acc_maneuver+T_coast);
    Event *acc_event = new SatelliteEvent(this, ts_acc_end);

    int Npoints = acc_maneuver->sequencer.GetNpoints();
    for(int i = 0;i < Npoints; i++){
        TimeSpec ts = eventQueue.AbsoluteTime(acc_maneuver->sequencer.GetTime(i) + T_coast);
        Event* ievent = new IntegratorEvent(ts);
        eventQueue.AddEvent(ievent);
    }

    maneuverQueue.AddManeuver(acc_maneuver);
    eventQueue.AddEvent(acc_event);

}

void Satellite::ScheduleConstantAccelerationManeuver(
    double latitude, double longitude,
    Eigen::RowVectorXd &time,
    double thrust,
    double T_coast)
{
    Eigen::Vector3d n_acc = direction_from_coord(
        latitude, longitude);

    Eigen::Vector3d acc = n_acc*thrust;

    double T_maneuver = time(time.cols()-1);

    if(T_coast!=0.0){
        CoastManeuver *coast_maneuver =
            new CoastManeuver(T_coast);
        TimeSpec ts_coast_end = eventQueue.AbsoluteTime(T_coast);
        Event* coast_event = new SatelliteEvent(this, ts_coast_end);
        maneuverQueue.AddManeuver(coast_maneuver);
        eventQueue.AddEvent(coast_event);
    }

    TimeSpec ts_acc_end = eventQueue.AbsoluteTime(T_maneuver+T_coast);
    Event *acc_event = new SatelliteEvent(this, ts_acc_end);

    int Npoints = time.cols();
    for(int i = 0;i < Npoints; i++){
        TimeSpec ts = eventQueue.AbsoluteTime(time(i) + T_coast);
        Event* ievent = new IntegratorEvent(ts);
        eventQueue.AddEvent(ievent);
    }

    ConstantAccelerationManeuver *acc_maneuver =
        new ConstantAccelerationManeuver(acc, T_maneuver);

    maneuverQueue.AddManeuver(acc_maneuver);
    eventQueue.AddEvent(acc_event);
}

void Satellite::ScheduleThrustVectorLaunch(void)
{
    if(!maneuverQueue.Empty()) return;

    double latitude = 2.0 * RAD_PER_DEG;
    double longitude = 0.0;//-78.5 * RAD_PER_DEG;

    Eigen::Vector3d y_axis(0,1,0);
    Eigen::Vector3d z_axis(0,0,1);
    Eigen::Matrix3d rot_lat = caams::AAA(latitude,-y_axis);
    Eigen::Matrix3d rot_long = caams::AAA(longitude,z_axis);
    Eigen::Matrix3d A_basis = rot_long*rot_lat;
    Eigen::Vector3d n_r_earth = A_basis.col(0);
    Eigen::Vector3d r_earth = n_r_earth*(primary_body->m_radius + m_radius);
    Eigen::Matrix3d A_sat_earth;
    A_sat_earth.col(0) = A_basis.col(1);
    A_sat_earth.col(1) = A_basis.col(2);
    A_sat_earth.col(2) = A_basis.col(0);
    Eigen::Matrix3d A_earth = caams::Ap(primary_body->p);
    Eigen::Matrix3d A_sat_world = A_earth*A_sat_earth;
    Eigen::Vector4d p_sat = caams::pA(A_sat_world);
    Eigen::Vector3d omega_earth = caams::omega_p_dot(primary_body->p, primary_body->pdot);
    Eigen::Vector4d pdot_sat = caams::p_dot_omega(p_sat, omega_earth);
    Eigen::Vector3d r_world = A_earth*r_earth;
    Eigen::Vector3d v_world = caams::SS(omega_earth)*r_world;
    p = p_sat;
    pdot = pdot_sat;
    m_position = r_world + primary_body->m_position;
    m_velocity = v_world + primary_body->m_velocity;

    OrbitalParams params;
    double r_peri = primary_body->m_radius + 500e3;
    Eigen::Vector3d n_peri = A_earth*n_r_earth;
    params.e = 1e-6;
    params.a = r_peri/(1-params.e);
    params.i = latitude;
    params.Omega = 270.0 * RAD_PER_DEG;
    params.omega = 90 * RAD_PER_DEG;//atan2(n_peri(1),n_peri(0));
    params.nu = 12.0 * RAD_PER_DEG;

    ThrustVectorData tvdata;
    tvdata.mass = m_mass;
    tvdata.thrust = THRUST;
    tvdata.mu = G_gravity * primary_body->m_mass;
    tvdata.sp = rp;
    tvdata.Jp = Jp;

    Eigen::MatrixXd controls;
    Eigen::MatrixXd time;

    thrust_vector_simple(
                r_world,
                v_world,
                params,
                &tvdata,
                controls,
                time);
}

void Satellite::ScheduleDelayedLaunch(
        double latitude, double longitude, // location on the primary body
        OrbitalParams const &params, // orbital parameters of the target orbit
        double thrust)
{
    if(!maneuverQueue.Empty()) return;

    Eigen::Vector3d n_rp = n_r_earth_frame(latitude, longitude);
    Eigen::Vector3d rp = n_rp*(m_radius + primary_body->m_radius);
    Eigen::Matrix3d Ap = caams::Ap(primary_body->p);
    Eigen::Vector3d rs = Ap*rp;
    Eigen::Vector3d ri1,ri2;
    if(!orbit_ring_intersect(params, rs, ri1, ri2)){
        return;
    }
    double long_rs_0 = atan2(rs(1), rs(0));
    double long_ri1_0 = atan2(ri1(1), ri1(0));
    double long_init;
    double long_final;
    angular_limits(long_rs_0, long_ri1_0, long_init, long_final);
    double dlong = long_final - long_init;
    Eigen::Vector3d omega = caams::omega_p_dot(primary_body->p, primary_body->pdot);
    double mag_omega = omega.norm();
    double T_intersect = dlong/mag_omega;

    DelayedLaunchData *dldata = new DelayedLaunchData;
    dldata->mass = m_mass;
    dldata->thrust = thrust;
    dldata->mu = G_gravity * primary_body->m_mass;
    dldata->mag_omega = mag_omega;
    dldata->omega = omega;
    dldata->T1 = 1200;
    dldata->T2 = 10;
    dldata->T3 = 600;

    Spheroid *spheroid = dynamic_cast<Spheroid*>(primary_body);
    if(spheroid){
        dldata->ring_mass = spheroid->ring_mass;
        dldata->ring_radius = spheroid->ring_radius;
    }else{
        dldata->ring_mass = 0.0;
        dldata->ring_radius = 1e3;
    }
    dldata->G = G_gravity;

    Eigen::Matrix3d Arot = caams::AAA(-mag_omega*dldata->T1, omega);
    dldata->r1i = Arot*ri1;
    dldata->v1i = omega.cross(dldata->r1i);

    double nu_intersect = nu_from_r(params, ri1);
    dldata->params_final = params;
    dldata->params_final.nu = nu_intersect + 12*RAD_PER_DEG;

    MatrixXd controls;
    MatrixXd time;

    double T_opt = T_intersect - dldata->T1;

    delayed_launch(dldata, controls, time);

    double T_surface = T_opt + dldata->T1;

    SurfaceManeuver *surf_maneuver = new SurfaceManeuver(
                this, latitude, longitude, T_surface);
    maneuverQueue.AddManeuver(surf_maneuver);

    ElevatorManeuver *elev_maneuver = new ElevatorManeuver(
                this, thrust, dldata->T2);
    maneuverQueue.AddManeuver(elev_maneuver);

    DirectedAcceleration3dManeuver *acc_maneuver =
            new DirectedAcceleration3dManeuver(
                this, controls, time, thrust);
    maneuverQueue.AddManeuver(acc_maneuver);

    TimeSpec ts_surf_end = eventQueue.AbsoluteTime(T_surface);
    Event *surf_event = new SatelliteEvent(this, ts_surf_end);
    eventQueue.AddEvent(surf_event);

    TimeSpec ts_elev_end = eventQueue.AbsoluteTime(T_surface + dldata->T2);
    Event *elev_event = new SatelliteEvent(this, ts_elev_end);
    eventQueue.AddEvent(elev_event);

    int Npoints = time.cols();
    for(int i = 0;i < Npoints; i++){
        TimeSpec ts = eventQueue.AbsoluteTime(time(i) + T_surface + dldata->T2);
        Event* ievent = new IntegratorEvent(ts);
        eventQueue.AddEvent(ievent);
    }
    TimeSpec ts_acc_end = eventQueue.AbsoluteTime(T_surface + dldata->T2 + time(Npoints-1));
    Event *acc_event = new SatelliteEvent(this, ts_acc_end);
    eventQueue.AddEvent(acc_event);

    delete dldata;
}

bool Satellite::ScheduleBallisticLaunch(
        double lat_vehicle, double long_vehicle,
        double lat_target, double long_target,
        double T1, double T2,
        double thrust)
{
    if(!maneuverQueue.Empty()) return true;

    Eigen::Vector3d n_r_vehicle = n_r_earth_frame(lat_vehicle, long_vehicle);
    Eigen::Vector3d n_r_target = n_r_earth_frame(lat_target, long_target);
    double mag_r = m_radius + primary_body->m_radius;
    Eigen::Vector3d r_p_vehicle = n_r_vehicle*mag_r;
    Eigen::Vector3d r_p_target = n_r_target*mag_r;
    Eigen::Matrix3d A = caams::Ap(primary_body->p);
    Eigen::Vector3d r_vehicle0 = A*r_p_vehicle;
    Eigen::Vector3d r_target0 = A*r_p_target;

    Eigen::Vector3d omega = caams::omega_p_dot(primary_body->p, primary_body->pdot);
    double mag_omega = omega.norm();

    Eigen::Vector3d v_vehicle0 = omega.cross(r_vehicle0);
    Eigen::Vector3d v_target0 = omega.cross(r_target0);

    BallisticData bdata;
    bdata.mass = m_mass;
    bdata.thrust = thrust;
    bdata.mu = G_gravity*primary_body->m_mass;
    bdata.mag_omega = mag_omega;
    bdata.omega = omega;
    bdata.T1 = T1;
    bdata.T2 = T2;
    bdata.r_vehicle0 = r_vehicle0;
    bdata.v_vehicle0 = v_vehicle0;
    bdata.r_target0 = r_target0;
    bdata.v_target0 = v_target0;

    Spheroid *spheroid = dynamic_cast<Spheroid*>(primary_body);
    if(spheroid){
        bdata.ring_mass = spheroid->ring_mass;
        bdata.ring_radius = spheroid->ring_radius;
    }else{
        bdata.ring_mass = 0.0;
        bdata.ring_radius = 1e3;
    }
    bdata.G = G_gravity;

    MatrixXd controls;
    MatrixXd time;

    if(!ballistic_launch(&bdata, controls, time))
        return false;

    DirectedAcceleration3dManeuver *acc_maneuver =
            new DirectedAcceleration3dManeuver(
                this, controls, time, thrust);
    maneuverQueue.AddManeuver(acc_maneuver);

    CoastManeuver *coast_maneuver = new CoastManeuver(bdata.T2);
    maneuverQueue.AddManeuver(coast_maneuver);

    int Npoints = time.cols();
    for(int i = 0;i < Npoints; i++){
        TimeSpec ts = eventQueue.AbsoluteTime(time(i));
        Event* ievent = new IntegratorEvent(ts);
        eventQueue.AddEvent(ievent);
    }
    TimeSpec ts_acc_end = eventQueue.AbsoluteTime(time(Npoints-1));
    Event *acc_event = new SatelliteEvent(this, ts_acc_end);
    eventQueue.AddEvent(acc_event);

    TimeSpec ts_coast_end = eventQueue.AbsoluteTime(bdata.T1 + bdata.T2);
    Event *coast_event = new SatelliteEvent(this, ts_coast_end);
    eventQueue.AddEvent(coast_event);

    m_position = primary_body->m_position + r_vehicle0;
    m_velocity = primary_body->m_velocity + v_vehicle0;

    return true;
}

void Satellite::ScheduleICBMLaunchTest(Eigen::MatrixXd &controls,
                                       Eigen::MatrixXd &time,
                                       double thrust,
                                       double radius_target,
                                       double &delta_t,
                                       double &T_reentry)
{
    DirectedAcceleration3dManeuver *acc_maneuver =
            new DirectedAcceleration3dManeuver(
                this, controls, time, thrust);
    maneuverQueue.AddManeuver(acc_maneuver);

    ReentryTestManeuver *reentry_maneuver =
            new ReentryTestManeuver(
                this, radius_target, delta_t, T_reentry);
    maneuverQueue.AddManeuver(reentry_maneuver);

    double T_acc_maneuver = acc_maneuver->Duration();
    int Npoints = acc_maneuver->sequencer.GetNpoints();
    for(int i = 0;i < Npoints; i++){
        TimeSpec ts = eventQueue.AbsoluteTime(acc_maneuver->sequencer.GetTime(i));
        Event* ievent = new IntegratorEvent(ts);
        eventQueue.AddEvent(ievent);
    }

    TimeSpec ts_acc_end = eventQueue.AbsoluteTime(T_acc_maneuver);
    Event *acc_event = new SatelliteEvent(this, ts_acc_end);
    eventQueue.AddEvent(acc_event);
}

bool Satellite::OptimizeICBMLaunch(
        Eigen::Vector3d const &r_launch,
        Eigen::Vector3d const &v_launch,
        Eigen::Vector3d const &r_aim,
        double h_apogee,
        double T1,
        double thrust,
        Eigen::MatrixXd &controls,
        Eigen::MatrixXd &time,
        Eigen::Vector3d &r_strike,
        double &T_reentry)
{
    ICBMSimpleData bdata;
    OrbitalParams params;
    // Calculate the orbit for the trajectory
    ballistic_orbit_params(
                params,
                r_launch,
                r_aim,
                primary_body->m_radius + h_apogee);

    Eigen::Matrix3d A_basis = orbit_basis(params);
    Eigen::Matrix3d A_basis_90 = orbit_basis_90(A_basis);

    // calculate the estimate of the end of the boost phase
    Eigen::Vector3d r_boost;
    Eigen::Vector3d v_boost;
    params.nu += 3.0*RAD_PER_DEG;
    orbit_state(
                params,
                primary_body->m_mass,
                m_mass,
                r_boost,
                v_boost);

    bdata.mass = m_mass;
    bdata.thrust = thrust;
    bdata.mu = G_gravity*primary_body->m_mass;
    bdata.T1 = T1;
    bdata.r_vehicle0 = r_launch;
    bdata.v_vehicle0 = v_launch;
    bdata.r_boost = r_boost;
    bdata.v_boost = v_boost;
    bdata.params = params;
    bdata.params.i = M_PI/2.0;
    bdata.params.Omega = M_PI/2.0;
    bdata.params.omega = M_PI/2.0;
    bdata.A_g2a = A_basis_90.transpose();

    bool is_spheroid;
    Spheroid *spheroid = dynamic_cast<Spheroid*>(primary_body);
    Body *body;
    if(spheroid){
        is_spheroid = true;
        bdata.ring_mass = spheroid->ring_mass;
        bdata.ring_radius = spheroid->ring_radius;
    }else{
        is_spheroid = false;
        body = primary_body;
        bdata.ring_mass = 0.0;
        bdata.ring_radius = 1e3;
    }
    bdata.G = G_gravity;

    if(!icbm_simple_launch(&bdata, controls, time))
        return false;

    System system;
    Satellite *satellite;
    Body *earth;

    if(is_spheroid){
        earth = new Spheroid(*spheroid);
    }else{
        earth = new Body(*body);
    }
    satellite = new Satellite(*this,earth);
    system.AddBody(satellite);
    system.AddBody(earth);

    double delta_t = 1.0;

    satellite->ScheduleICBMLaunchTest(
                controls,
                time,
                thrust,
                r_aim.norm(),
                delta_t,
                T_reentry);

    satellite->m_position = earth->m_position + r_launch;
    satellite->m_velocity = earth->m_velocity + v_launch;

    while(!satellite->maneuverQueue.Empty()){
        system.rkIntegrate(delta_t);
    }

    r_strike = satellite->RelativePosition();

    delete earth;
    delete satellite;

    return true;
}

bool Satellite::ScheduleICBMLaunch(
        double lat_vehicle, double long_vehicle,
        double hgt_vehicle,
        double lat_target, double long_target,
        double hgt_target)
{
    if(!maneuverQueue.Empty()) return true;

    Eigen::Vector3d n_r_vehicle = n_r_earth_frame(lat_vehicle, long_vehicle);
    Eigen::Vector3d n_r_target = n_r_earth_frame(lat_target, long_target);
    double mag_r_vehicle = primary_body->m_radius + hgt_vehicle;
    double mag_r_target = primary_body->m_radius + hgt_target;
    Eigen::Vector3d r_p_vehicle = n_r_vehicle*mag_r_vehicle;
    Eigen::Vector3d r_p_target = n_r_target*mag_r_target;
    Eigen::Matrix3d A = caams::Ap(primary_body->p);
    Eigen::Vector3d r_vehicle0 = A*r_p_vehicle;
    Eigen::Vector3d r_target0 = A*r_p_target;

    Eigen::Vector3d omega = caams::omega_p_dot(primary_body->p, primary_body->pdot);
    double mag_omega = omega.norm();

    Eigen::Vector3d v_vehicle0 = omega.cross(r_vehicle0);

    MatrixXd controls;
    MatrixXd time;
    double T_reentry;
    Eigen::Vector3d r_strike;
    double h_apogee = 1500.0e3;

    Eigen::Vector3d r_aim = r_target0;

    do{

        if(!OptimizeICBMLaunch(
                    r_vehicle0,
                    v_vehicle0,
                    r_aim,
                    h_apogee,
                    450, // T1
                    THRUST,
                    controls,
                    time,
                    r_strike,
                    T_reentry)){
            return false;
        }

        // find the location of the target after rotation of the earth
        double T_total = time(time.cols()-1) + T_reentry;
        Eigen::Matrix3d A_rot = caams::AAA(mag_omega*T_total, omega);
        Eigen::Vector3d r_target1 = A_rot*r_target0;

        // calculate the angle between the strike point and the target
        Eigen::Vector3d n_target = r_target1.normalized();
        Eigen::Vector3d n_strike = r_strike.normalized();
        double phi = acos(n_target.dot(n_strike));
        std::cout << "phi:" << phi*180/M_PI << std::endl;
        std::cout << "s:" << phi*primary_body->m_radius << std::endl;
        if(phi < (100.0/primary_body->m_radius)){
            break;
        }

        // rotation axis from n_strike to n_target
        Eigen::Vector3d axis = n_strike.cross(n_target);

        // rotation matrix to correct the aim point
        Eigen::Matrix3d A_corr = caams::AAA(phi, axis);
        r_aim = A_corr*r_aim;

    }while(1);

    DirectedAcceleration3dManeuver *acc_maneuver =
            new DirectedAcceleration3dManeuver(
                this, controls, time, THRUST);
    maneuverQueue.AddManeuver(acc_maneuver);
    double T_acc_maneuver = acc_maneuver->Duration();

    CoastManeuver *coast_maneuver = new CoastManeuver(T_reentry);
    maneuverQueue.AddManeuver(coast_maneuver);

    int Npoints = time.cols();
    for(int i = 0;i < Npoints; i++){
        TimeSpec ts = eventQueue.AbsoluteTime(time(i));
        Event* ievent = new IntegratorEvent(ts);
        eventQueue.AddEvent(ievent);
    }
    TimeSpec ts_acc_end = eventQueue.AbsoluteTime(time(Npoints-1));
    Event *acc_event = new SatelliteEvent(this, ts_acc_end);
    eventQueue.AddEvent(acc_event);

    TimeSpec ts_coast_end = eventQueue.AbsoluteTime(T_acc_maneuver + T_reentry);
    Event *coast_event = new SatelliteEvent(this, ts_coast_end);
    eventQueue.AddEvent(coast_event);

    m_position = primary_body->m_position + r_vehicle0;
    m_velocity = primary_body->m_velocity + v_vehicle0;

    return true;
}



