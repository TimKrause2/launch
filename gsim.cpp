#include <math.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "gsim.h"
#include "uv_sphere.h"

glTexture *Body::cursor_tex;
int Body::n_bodies = 0;

Body::Body(Body &b)
    :
      Body(b.m_mass,
           b.m_radius,
           b.m_position,
           b.m_velocity,
           b.p,
           b.pdot,
           b.Jp,
           NULL)
{}

Body::Body(double mass,
           double radius,
           Eigen::Vector3d &position,
           Eigen::Vector3d &velocity,
           Eigen::Vector4d &p,
           Eigen::Vector4d &pdot,
           Eigen::Matrix3d &Jp,
           char* tex_file):
    m_mass(mass),
    m_radius(radius),
    m_position(position),
    m_velocity(velocity),
    p(p),
    pdot(pdot),
    Jp(Jp),
    enabled(true)
{
    if(tex_file){
        texture = texture_load(tex_file);
    }else{
        texture = NULL;
    }
    InitializeCursor();
}

Body::~Body(void)
{
    DeleteCursor();
}

void Body::InitializeCursor(void)
{
    if(!n_bodies){
        cursor_tex = rgba_texture_load("cursor.tiff");
    }
    n_bodies++;
}

void Body::DeleteCursor(void)
{
    if(n_bodies){
        n_bodies--;
        if(!n_bodies){
            texture_free(cursor_tex);
        }
    }
}

void Body::PushState(void)
{
    states.push_front(State(m_position,m_velocity,p,pdot));
}

void Body::PopState(void)
{
    if(states.empty())return;
    states.pop_front();
}

void Body::RestoreState(void)
{
    if(states.empty())
        return;
    *this = states.front();
}

bool Body::Enabled(void)
{
    return enabled;
}

void Body::Enable(void)
{
    enabled = true;
}

void Body::Disable(void)
{
    enabled = false;
}

void Body::draw(glm::mat4 proj, glm::mat4 view)
{
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    int width = vp[2];
    int height = vp[3];

    // determine pixel size of diameter
    glm::vec4 rp(m_position(0), m_position(1), m_position(2), 1.0f);
    glm::vec4 rp_cam = view*rp;
    if(rp_cam.z>=0.0f) return;
    glm::vec4 rp_test(m_radius*2, 0.0f, rp_cam.z, 1.0f);
    glm::vec4 rp_pers = proj*rp_test;
    float pixels = (float)width/2.0f*rp_pers.x/rp_pers.w;

    if(pixels<4.0f || !texture){
        // plot the cursor
        // find the pixel location of the body
        rp_pers = proj*rp_cam;
        float p_x = rp_pers.x / rp_pers.w;
        if((p_x<-1.0f) || (p_x>1.0f)) return;
        float p_y = rp_pers.y / rp_pers.w;
        if((p_y<-1.0f) || (p_y>1.0f)) return;
        float xc = (float)width*(p_x+1.0f)/2.0f;
        float yc = (float)height*(p_y+1.0f)/2.0f;
        texture_sprite(cursor_tex, xc, yc, 0.4f);
    }

    if(!texture)
        return;
    glm::vec3 scale((float)m_radius,(float)m_radius,(float)m_radius);
    glm::mat4 A_scale = glm::scale(glm::mat4(1.0),scale);

    Eigen::Matrix3d A = caams::Ap(p);
    glm::dmat3 A_rot_dbl = E2GLM(A);
    glm::mat3 A_rot_flt(A_rot_dbl);
    glm::mat4 A_rotation(A_rot_flt);

    //glMultMatrixd( glm::value_ptr(A_body) );
    //glm::vec3 translation((float)m_position.data[0],(float)m_position.data[1],(float)m_position.data[2]);
    glm::vec3 translation = E2GLM(m_position);
    glm::mat4 A_translation = glm::translate(glm::mat4(1.0),translation);

    glm::mat4 model = A_translation*A_rotation*A_scale;

    //glMultMatrixd( glm::value_ptr( A_body ) );
    glm::mat4 mvp = proj*view*model;


    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    glEnable(GL_DEPTH_TEST);
    uv_sphere_draw(mvp, texture->tex_id);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

}


Eigen::Vector4d p_ddot_solve(
        Eigen::Vector4d &p,
        Eigen::Vector4d &p_dot,
        Eigen::Matrix3d &Jp,
        Eigen::Vector3d &np)
{
    Eigen::Matrix4d H = 4.0*caams::L(p_dot).transpose()*Jp*caams::L(p);
    Eigen::Matrix4d R;
    R.topLeftCorner(3,4) = caams::L(p)*H;
    R.row(3) = p_dot.transpose();

    Eigen::Vector4d rhs;
    rhs << np, 0.0;

    Eigen::Vector4d y = rhs - R*p_dot;

    Eigen::Matrix4d A;
    A.topLeftCorner(3,4) = 2.0*Jp*caams::L(p);
    A.row(3) = p.transpose();

    Eigen::Vector4d x = A.inverse()*y;

    return x;
}

// void Body::update_rotation( double dt )
// {
//     Eigen::Matrix4d k_p_dot;
//     Eigen::Matrix4d k_p_ddot;
//     Eigen::Vector4d p_norm;
//     k_p_dot.col(0) = pdot;
//     k_p_ddot.col(0) = p_ddot_solve(p,k_p_dot.col(0),Jp);

//     p_norm = p + (dt/2.0)*k_p_dot.col(0);
//     p_norm.normalize();
//     k_p_dot.col(1) = pdot + (dt/2.0)*k_p_ddot.col(0);
//     k_p_ddot.col(1) = p_ddot_solve(p_norm, k_p_dot.col(1), Jp);

//     p_norm = p + (dt/2.0)*k_p_dot.col(1);
//     p_norm.normalize();
//     k_p_dot.col(2) = pdot + (dt/2.0)*k_p_ddot.col(1);
//     k_p_ddot.col(2) = p_ddot_solve(p_norm, k_p_dot.col(2), Jp);

//     p_norm = p + dt*k_p_dot.col(2);
//     p_norm.normalize();
//     k_p_dot.col(3) = pdot + dt*k_p_ddot.col(2);
//     k_p_ddot.col(3) = p_ddot_solve(p_norm, k_p_dot.col(3), Jp);

//     Eigen::Vector4d c;
//     c << 1.0, 2.0, 2.0, 1.0;
//     c *= dt/6.0;

//     p = p + k_p_dot*c;
//     p.normalize();
//     pdot = pdot + k_p_ddot*c;
// }

void Body::Prepare()
{

}

void Body::Update(double dt)
{

}

void Body:: ForceAndTorque(
        double dt,
        const Eigen::VectorXd &y,
        Eigen::Vector3d &force,   // force in global coordinate space
        Eigen::Vector3d &torque)
{
    force = Eigen::Vector3d::Zero();
    torque = Eigen::Vector3d::Zero();
}

Eigen::Vector3d Body::AtmosphericDrag(
        Eigen::Vector3d &r,
        Eigen::Vector3d &v)
{
    return Eigen::Vector3d::Zero();
}

double Body::TimeStep(void)
{
    return 1e9;
}

double Body::CurvatureTimeStep()
{
    //caams::matrix v_x_a(caams::SS(m_velocity)*m_rk_acceleration);
    //double mag_v = caams::norm(m_velocity);
    //double omega = caams::norm(v_x_a)/mag_v/mag_v;
    Eigen::Vector3d v_x_a = m_velocity.cross(m_rk_acceleration);
    double mag_v2 = m_velocity.squaredNorm();
    double omega = v_x_a.norm()/mag_v2;
    double dt;
    if(omega==0.0){
        dt = 1e9;
    }else{
        dt = 2*M_PI/1024/omega;
    }
    return dt;
}

Eigen::Vector3d Body::rk_acceleration(
        Eigen::Vector3d const &rk_position,
        Eigen::Vector3d const &rk_velocity)
{
    Eigen::Vector3d r = rk_position - m_rk_position;
    double mag_r2 = r.squaredNorm();
    double mag_a = G_gravity*m_mass/mag_r2;
    return -r.normalized()*mag_a;
}

void BodyPair::calculate_forces(Eigen::VectorXd const &y)
{
    // vector from body1 to body2
    Eigen::Vector3d r1 = y.segment<3>(body1->y_offset + R_OFFSET);
    Eigen::Vector3d r2 = y.segment<3>(body2->y_offset + R_OFFSET);
    Eigen::Vector3d r_21 = r2 - r1;
    Eigen::Vector3d r_12 = -r_21;
    Eigen::Vector3d v1 = y.segment<3>(body1->y_offset + V_OFFSET);
    Eigen::Vector3d v2 = y.segment<3>(body2->y_offset + V_OFFSET);
    Eigen::Vector3d v_21 = v2 - v1;
    Eigen::Vector3d v_12 = -v_21;
    double mag_r2 = r_21.squaredNorm();
    double F = G_gravity*body1->m_mass*body2->m_mass/mag_r2;
    Eigen::Vector3d F_12 = r_21/sqrt(mag_r2)*F;
    body1->m_rk_force += body2->AtmosphericDrag(r_12, v_12) + F_12;
    body2->m_rk_force += body1->AtmosphericDrag(r_21, v_21) - F_12;
}

System::System(void)
{
    n_proc_threads = std::thread::hardware_concurrency();
    m_threads.resize(n_proc_threads);
    current_y_offset = 0;
}

void System::AddBody(Body* p_body)
{
	m_bodies.push_back(p_body);
    p_body->y_offset = current_y_offset;
    current_y_offset += BODY_OFFSET;
}

void System::rkPrepare()
{
    std::list<Body*>::iterator l_it_b1;
    for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        if(!(*l_it_b1)->Enabled()) continue;
        (*l_it_b1)->Prepare();
    }
}

void System::rkUpdate(double dt)
{
    std::list<Body*>::iterator l_it_b1;
    for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        if(!(*l_it_b1)->Enabled()) continue;
        (*l_it_b1)->Update(dt);
    }
}

double System::rkTimeStep(void)
{
    double r = 1e9;
    std::list<Body*>::iterator l_it_b1;
    for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        r = fmin(r, (*l_it_b1)->TimeStep());
    }
    return r;
}

// void System::rkAccelerations( double dt )
// {
//     std::list<Body*>::iterator l_it_b1;
//     std::list<Body*>::iterator l_it_b2;
//     int n_bodies=0;
// 	// set all accelerations to zero
// 	for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
//         if(!(*l_it_b1)->Enabled()) continue;
//         Eigen::Vector3d force;
//         Eigen::Vector3d torque;
//         (*l_it_b1)->ForceAndTorque(dt, force, torque);
//         (*l_it_b1)->m_rk_force = force;
//         (*l_it_b1)->rk_pddot = p_ddot_solve(
//                     (*l_it_b1)->rk_p,
//                     (*l_it_b1)->rk_pdot,
//                     (*l_it_b1)->Jp,
//                     torque);
//         n_bodies++;
// 	}

//     // clear the space for the pairs list
//     int n_pairs = n_bodies*(n_bodies-1)/2;
//     m_body_pairs.clear();

//     // create the body pair list
//     for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
//         if(!(*l_it_b1)->Enabled()) continue;
//         l_it_b2 = l_it_b1;
//         l_it_b2++;
//         for( ; l_it_b2!=m_bodies.end(); l_it_b2++){
//             if(!(*l_it_b2)->Enabled()) continue;
//             m_body_pairs.emplace_back(*l_it_b1, *l_it_b2);
//         }
//     }

//     // execute the pair list
//     std::list<BodyPair>::iterator bp_it1;
//     std::list<BodyPair>::iterator bp_it2;
//     if(n_pairs<n_proc_threads){
//         bp_it1 = m_body_pairs.begin();
//         bp_it2 = bp_it1;
//         for(int i=0;i<n_pairs;i++)
//             bp_it2++;
//         rkForces(bp_it1, bp_it2);
//     }else{
//         std::list<std::thread>::iterator th_it;
//         th_it = m_threads.begin();
//         bp_it1 = m_body_pairs.begin();
//         int i_bp_last = 0;
//         for(int i=0;i<n_proc_threads;i++,th_it++){
//             int i_bp_next = (i+1)*n_pairs/n_proc_threads;
//             int n = i_bp_next - i_bp_last;
//             bp_it2 = bp_it1;
//             for(int j=0;j<n;j++)
//                 bp_it2++;
//             *th_it = std::thread(&System::rkForces, this, bp_it1, bp_it2);
//             bp_it1 = bp_it2;
//             i_bp_last = i_bp_next;
//         }
//         th_it = m_threads.begin();
//         for(int i=0;i<n_proc_threads;i++, th_it++)
//             (*th_it).join();
//     }


//     // calculate accelerations from the forces
// 	for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
//         (*l_it_b1)->m_rk_acceleration = (*l_it_b1)->m_rk_force/(*l_it_b1)->m_mass;
// 	}
// }

void System::rkInterBodyForces(
        std::list<BodyPair>::iterator first_pair,
        std::list<BodyPair>::iterator last_pair,
        Eigen::VectorXd const &y)
{
    std::list<BodyPair>::iterator it;
    for(it = first_pair; it != last_pair; it++){
        (*it).calculate_forces(y);
    }
}

Eigen::Vector4d System::ortho_p_dot(
        Eigen::Vector4d const &p,
        Eigen::Vector4d const &pdot)
{
    double sigma = p.dot(pdot);
    return pdot - sigma*p;
}

void System::constrainRotations(Eigen::VectorXd &y)
{
    for(int y_offset=0;y_offset<current_y_offset;y_offset+=BODY_OFFSET)
    {
        y.segment<4>(y_offset+P_OFFSET).normalize();
        y.segment<4>(y_offset+PDOT_OFFSET) =
            ortho_p_dot(y.segment<4>(y_offset+P_OFFSET),
                        y.segment<4>(y_offset+PDOT_OFFSET));
    }
}

Eigen::VectorXd System::dy_func(double t, Eigen::VectorXd y)
{
    Eigen::VectorXd dy(current_y_offset);

    // initialize the body generated forces and torques
    // caclulate p_ddot for each body
    // while travesring the list generate the body pair list
    // for gravity and atmospheric drag
    std::list<Body*>::iterator l_it_b1;
    std::list<Body*>::iterator l_it_b2;
    int n_bodies=0;
    // set all accelerations to zero
    for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        Eigen::Vector3d force;
        Eigen::Vector3d torque;
        (*l_it_b1)->ForceAndTorque(t, y, force, torque);
        (*l_it_b1)->m_rk_force = force;
        Eigen::Vector4d p = y.segment<4>((*l_it_b1)->y_offset + P_OFFSET);
        Eigen::Vector4d pdot = y.segment<4>((*l_it_b1)->y_offset + PDOT_OFFSET);
        dy.segment<4>((*l_it_b1)->y_offset + P_OFFSET) =
                pdot;
        dy.segment<4>((*l_it_b1)->y_offset + PDOT_OFFSET) =
                p_ddot_solve(p, pdot, (*l_it_b1)->Jp, torque);
        n_bodies++;
    }

    // clear the space for the pairs list
    int n_pairs = n_bodies*(n_bodies-1)/2;
    m_body_pairs.clear();

    // create the body pair list
    for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        l_it_b2 = l_it_b1;
        l_it_b2++;
        for( ; l_it_b2!=m_bodies.end(); l_it_b2++){
            m_body_pairs.emplace_back(*l_it_b1, *l_it_b2);
        }
    }

    // execute the pair list
    std::list<BodyPair>::iterator bp_it1;
    std::list<BodyPair>::iterator bp_it2;
    if(n_pairs<n_proc_threads){
        bp_it1 = m_body_pairs.begin();
        bp_it2 = bp_it1;
        for(int i=0;i<n_pairs;i++)
            bp_it2++;
        rkInterBodyForces(bp_it1, bp_it2, y);
    }else{
        std::list<std::thread>::iterator th_it;
        th_it = m_threads.begin();
        bp_it1 = m_body_pairs.begin();
        int i_bp_last = 0;
        for(int i=0;i<n_proc_threads;i++,th_it++){
            int i_bp_next = (i+1)*n_pairs/n_proc_threads;
            int n = i_bp_next - i_bp_last;
            bp_it2 = bp_it1;
            for(int j=0;j<n;j++)
                bp_it2++;
            *th_it = std::thread(&System::rkInterBodyForces, this, bp_it1, bp_it2, y);
            bp_it1 = bp_it2;
            i_bp_last = i_bp_next;
        }
        th_it = m_threads.begin();
        for(int i=0;i<n_proc_threads;i++, th_it++)
            (*th_it).join();
    }

    // calculate accelerations from the forces
    for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        dy.segment<3>((*l_it_b1)->y_offset + R_OFFSET) =
                y.segment<3>((*l_it_b1)->y_offset + V_OFFSET);
        dy.segment<3>((*l_it_b1)->y_offset + V_OFFSET) =
                (*l_it_b1)->m_rk_force/(*l_it_b1)->m_mass;
    }

    return dy;
}

void System::rkIntegrate( double p_dt_total )
{
    double t=0.0;
    double dt;

    Eigen::VectorXd y0(current_y_offset);
    Eigen::VectorXd y1(current_y_offset);

    // initialize y0
    std::list<Body*>::iterator l_it_b1;
    for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        int y_offset = (*l_it_b1)->y_offset;
        y0.segment<3>(y_offset + R_OFFSET) = (*l_it_b1)->m_position;
        y0.segment<3>(y_offset + V_OFFSET) = (*l_it_b1)->m_velocity;
        y0.segment<4>(y_offset + P_OFFSET) = (*l_it_b1)->p;
        y0.segment<4>(y_offset + PDOT_OFFSET) = (*l_it_b1)->pdot;
    }

    // for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
    //     Eigen::Vector3d force;
    //     Eigen::Vector3d torque;
    //     (*l_it_b1)->ForceAndTorque(t, y0, force, torque);
    // }

    while(t<p_dt_total){
        double dt_remain = p_dt_total - t;
        rkPrepare();
        double dt_max = rkTimeStep();
        dt = fmin(dt_remain, dt_max);
        integrator.integrate(
            y1, y0,
            [this](double t, Eigen::VectorXd y){
            return this->dy_func(t, y);},
            // std::bind(&System::dy_func,
            //           this,
            //           std::placeholders::_1,
            //           std::placeholders::_2),
            0, dt);
        constrainRotations(y1);
        rkUpdate(dt);
        y0 = y1;
        t+=dt;
    }

    // transfer the results back to the objects
    for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        int y_offset = (*l_it_b1)->y_offset;
        (*l_it_b1)->m_position = y1.segment<3>(y_offset + R_OFFSET);
        (*l_it_b1)->m_velocity = y1.segment<3>(y_offset + V_OFFSET);
        (*l_it_b1)->p = y1.segment<4>(y_offset + P_OFFSET);
        (*l_it_b1)->pdot = y1.segment<4>(y_offset + PDOT_OFFSET);
    }
}

void System::extrapolate(double t)
{
    // save the state
    std::list<Body*>::iterator l_it;
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        (*l_it)->ex_position0 = (*l_it)->m_position;
        (*l_it)->ex_velocity0 = (*l_it)->m_velocity;
    }

    // integrate to the time
    rkIntegrate(t);

    // save the last sate and restore the original state
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        (*l_it)->ex_position1 = (*l_it)->m_position;
        (*l_it)->ex_velocity1 = (*l_it)->m_velocity;
        (*l_it)->m_position = (*l_it)->ex_position0;
        (*l_it)->m_velocity = (*l_it)->ex_velocity0;
    }
}

void System::draw(glm::mat4 proj, glm::mat4 view)
{
    std::list<Body*>::iterator l_it;
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
        (*l_it)->draw(proj, view);
    }
}

void System::origin(Eigen::Vector3d r0 )
{
    std::list<Body*>::iterator l_it;
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        (*l_it)->m_position -= r0;
    }

}

void System::PushState(void)
{
    std::list<Body*>::iterator l_it;
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        (*l_it)->PushState();
    }
}

void System::PopState(void)
{
    std::list<Body*>::iterator l_it;
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        (*l_it)->PopState();
    }
}

void System::RestoreState(void)
{
    std::list<Body*>::iterator l_it;
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        (*l_it)->RestoreState();
    }
}

Eigen::Vector3d System::accelerationAtPoint(Eigen::Vector3d p)
{
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    std::list<Body*>::iterator l_it;
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled())
            continue;
        Eigen::Vector3d r_bp = (*l_it)->m_position - p;
        Eigen::Vector3d n_bp = r_bp.normalized();
        double mag_r2 = r_bp.squaredNorm();
        double mag_acc = (*l_it)->m_mass*G_gravity/mag_r2;
        acc += mag_acc*n_bp;
    }
    return acc;
}



