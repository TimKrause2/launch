#include <math.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "gsim.h"
#include "image.h"
#include "uv_sphere.h"
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
    if(!tex_file || tex_file[0]==0){
        texName = 0;
        return;
    }

    // load the texture for the body
    image* tex_image = image_load(tex_file);
    if(!tex_image){
        printf("couldn't load texture file:\"%s\"",tex_file);
        return;
    }

    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

    glTexImage2D(
                GL_TEXTURE_2D, 0, GL_RGB8,
                tex_image->width, tex_image->height,
                0,
                GL_RGB, GL_UNSIGNED_BYTE,
                (const GLvoid*)tex_image->data );
    image_free(tex_image);
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
    // first set up the transformation matrix
    //glMatrixMode(GL_MODELVIEW);
    //glPushMatrix();

    glm::vec3 scale((float)m_radius,(float)m_radius,(float)m_radius);
    glm::mat4 A_scale = glm::scale(glm::mat4(1.0),scale);

    //caams::matrix A(caams::Ap(p));
    Eigen::Matrix3d A = caams::Ap(p);
    //glm::dmat3 A_rot_dbl = glm::make_mat3((~A).data);
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
    uv_sphere_draw(mvp, texName);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    // set up texturing
//    glEnable( GL_TEXTURE_2D );
//    glEnable(GL_DEPTH_TEST);
//    glEnable( GL_CULL_FACE );
//    glFrontFace( GL_CW );
//    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );
//    glBindTexture( GL_TEXTURE_2D, texName );
//    uv_sphere( SPHERE_N_LATITUDE, SPHERE_N_LONGITUDE );
//    glDisable( GL_TEXTURE_2D );
//    glDisable( GL_DEPTH_TEST );

//    glPopMatrix();
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
        Eigen::Vector3d &force,   // force in global coordinate space
        Eigen::Vector3d &torque)
{
    force = Eigen::Vector3d::Zero();
    torque = Eigen::Vector3d::Zero();
}

double Body::TimeStep(void)
{
    return CurvatureTimeStep();
}

double Body::CurvatureTimeStep()
{
    //caams::matrix v_x_a(caams::SS(m_velocity)*m_rk_acceleration);
    //double mag_v = caams::norm(m_velocity);
    //double omega = caams::norm(v_x_a)/mag_v/mag_v;
    Eigen::Vector3d v_x_a = m_velocity.cross(m_rk_acceleration);
    double mag_v2 = m_velocity.squaredNorm();
    double omega = v_x_a.norm()/mag_v2;
    if(omega==0.0){
        return 1e9;
    }else{
        return 2*M_PI/1024/omega;
    }
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

void System::AddBody(Body* p_body)
{
	m_bodies.push_back(p_body);
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

void System::rkAccelerations( double dt )
{
    std::list<Body*>::iterator l_it_b1;
    std::list<Body*>::iterator l_it_b2;
	// set all accelerations to zero
	for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        if(!(*l_it_b1)->Enabled()) continue;
        Eigen::Vector3d force;
        Eigen::Vector3d torque;
        (*l_it_b1)->ForceAndTorque(dt, force, torque);
        (*l_it_b1)->m_rk_acceleration = force/(*l_it_b1)->m_mass;
        (*l_it_b1)->rk_pddot = p_ddot_solve(
                    (*l_it_b1)->rk_p,
                    (*l_it_b1)->rk_pdot,
                    (*l_it_b1)->Jp,
                    torque);
	}
	// calculate accelerations for all pairs of bodies
	for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
        if(!(*l_it_b1)->Enabled()) continue;
        l_it_b2 = l_it_b1;
		l_it_b2++;
        for( ; l_it_b2!=m_bodies.end(); l_it_b2++){
            if(!(*l_it_b2)->Enabled()) continue;

            (*l_it_b1)->m_rk_acceleration += (*l_it_b2)->rk_acceleration(
                        (*l_it_b1)->m_rk_position,
                        (*l_it_b1)->m_rk_velocity);
            (*l_it_b2)->m_rk_acceleration += (*l_it_b1)->rk_acceleration(
                        (*l_it_b2)->m_rk_position,
                        (*l_it_b2)->m_rk_velocity);
		}
	}
}

void System::ortho_p_dot(Eigen::Vector4d const &p, Eigen::Vector4d &pdot)
{
    double sigma = p.dot(pdot);
    pdot -= sigma*p;
}

void System::rkPhase1Positions( void )
{
    std::list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
        (*l_it)->m_rk_position = (*l_it)->m_position;
        (*l_it)->m_rk_velocity = (*l_it)->m_velocity;
        (*l_it)->rk_p = (*l_it)->p;
        (*l_it)->rk_pdot = (*l_it)->pdot;
	}
}

void System::rkPhase2Positions( double p_dt )
{
    std::list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
        (*l_it)->m_rk_position = (*l_it)->m_position + (*l_it)->m_kr.col(0)*(p_dt*0.5);
        (*l_it)->m_rk_velocity = (*l_it)->m_velocity + (*l_it)->m_kv.col(0)*(p_dt*0.5);
        (*l_it)->rk_p = (*l_it)->p + (*l_it)->k_p_dot.col(0)*(p_dt*0.5);
        (*l_it)->rk_p.normalize();
        (*l_it)->rk_pdot = (*l_it)->pdot + (*l_it)->k_p_ddot.col(0)*(p_dt*0.5);
        ortho_p_dot((*l_it)->rk_p, (*l_it)->rk_pdot);
	}
}

void System::rkPhase3Positions( double p_dt )
{
    std::list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
        (*l_it)->m_rk_position = (*l_it)->m_position + (*l_it)->m_kr.col(1)*(p_dt*0.5);
        (*l_it)->m_rk_velocity = (*l_it)->m_velocity + (*l_it)->m_kv.col(1)*(p_dt*0.5);
        (*l_it)->rk_p = (*l_it)->p + (*l_it)->k_p_dot.col(1)*(p_dt*0.5);
        (*l_it)->rk_p.normalize();
        (*l_it)->rk_pdot = (*l_it)->pdot + (*l_it)->k_p_ddot.col(1)*(p_dt*0.5);
        ortho_p_dot((*l_it)->rk_p, (*l_it)->rk_pdot);
    }
}

void System::rkPhase4Positions( double p_dt )
{
    std::list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
        (*l_it)->m_rk_position = (*l_it)->m_position + (*l_it)->m_kr.col(2)*p_dt;
        (*l_it)->m_rk_velocity = (*l_it)->m_velocity + (*l_it)->m_kv.col(2)*p_dt;
        (*l_it)->rk_p = (*l_it)->p + (*l_it)->k_p_dot.col(2)*p_dt;
        (*l_it)->rk_p.normalize();
        (*l_it)->rk_pdot = (*l_it)->pdot + (*l_it)->k_p_ddot.col(2)*p_dt;
        ortho_p_dot((*l_it)->rk_p, (*l_it)->rk_pdot);
    }
}

void System::rkPhase1Integrate( double p_dt )
{
    std::list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
//        (*l_it)->m_kr.sub( (*l_it)->m_velocity, 1, 1);
//        (*l_it)->m_kv.sub( (*l_it)->m_rk_acceleration, 1, 1);
        (*l_it)->m_kr.col(0) = (*l_it)->m_velocity;
        (*l_it)->m_kv.col(0) = (*l_it)->m_rk_acceleration;
        (*l_it)->k_p_dot.col(0) = (*l_it)->rk_pdot;
        (*l_it)->k_p_ddot.col(0) = (*l_it)->rk_pddot;
    }
}

void System::rkPhase2Integrate( double p_dt )
{
    std::list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
//        (*l_it)->m_kr.sub( (*l_it)->m_velocity + (*l_it)->m_kv.sub(3,1,1,1)*(p_dt*0.5), 1, 2 );
//        (*l_it)->m_kv.sub( (*l_it)->m_rk_acceleration, 1, 2);
        (*l_it)->m_kr.col(1) = (*l_it)->m_velocity + (*l_it)->m_kv.col(0)*(p_dt*0.5);
        (*l_it)->m_kv.col(1) = (*l_it)->m_rk_acceleration;
        (*l_it)->k_p_dot.col(1) = (*l_it)->rk_pdot;
        (*l_it)->k_p_ddot.col(1) = (*l_it)->rk_pddot;
    }
}

void System::rkPhase3Integrate( double p_dt )
{
    std::list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
//        (*l_it)->m_kr.sub( (*l_it)->m_velocity + (*l_it)->m_kv.sub(3,1,1,2)*(p_dt*0.5), 1, 3 );
//        (*l_it)->m_kv.sub( (*l_it)->m_rk_acceleration, 1, 3);
        (*l_it)->m_kr.col(2) = (*l_it)->m_velocity + (*l_it)->m_kv.col(1)*(p_dt*0.5);
        (*l_it)->m_kv.col(2) = (*l_it)->m_rk_acceleration;
        (*l_it)->k_p_dot.col(2) = (*l_it)->rk_pdot;
        (*l_it)->k_p_ddot.col(2) = (*l_it)->rk_pddot;
    }
}

void System::rkPhase4Integrate( double p_dt )
{
    std::list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
        (*l_it)->m_kr.col(3) = (*l_it)->m_velocity + (*l_it)->m_kv.col(2)*p_dt;
        (*l_it)->m_kv.col(3) = (*l_it)->m_rk_acceleration;
        (*l_it)->k_p_dot.col(3) = (*l_it)->rk_pdot;
        (*l_it)->k_p_ddot.col(3) = (*l_it)->rk_pddot;
        Eigen::Vector4d c;
        c << 1.0, 2.0, 2.0, 1.0;
        c *= p_dt/6.0;
        (*l_it)->m_velocity += (*l_it)->m_kv*c;
        (*l_it)->m_position += (*l_it)->m_kr*c;
        (*l_it)->p += (*l_it)->k_p_dot*c;
        (*l_it)->p.normalize();
        (*l_it)->pdot += (*l_it)->k_p_ddot*c;
        ortho_p_dot((*l_it)->p, (*l_it)->pdot);
    }
}

double System::rkTimeStep()
{
    double dt_min;
    std::list<Body*>::iterator l_it;
    dt_min = 1e9;
    for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
        if(!(*l_it)->Enabled()) continue;
        double dt_test = (*l_it)->TimeStep();
        if(dt_test<dt_min){
            dt_min = dt_test;
        }
    }
    return dt_min;
}

void System::rkIntegrate( double p_dt_total )
{
    double t=0.0;
    double dt;

    while(t<p_dt_total){
        double dt_remain = p_dt_total - t;
        rkPrepare();
        rkPhase1Positions();
        rkAccelerations(0.0);
        double dt_max = rkTimeStep();
        if(dt_remain>dt_max){
            dt = dt_max;
        }else{
            dt = dt_remain;
        }
        rkPhase1Integrate(dt);

        rkPhase2Positions(dt);
        rkAccelerations(dt*0.5);
        rkPhase2Integrate(dt);

        rkPhase3Positions(dt);
        rkAccelerations(dt*0.5);
        rkPhase3Integrate(dt);

        rkPhase4Positions(dt);
        rkAccelerations(dt);
        rkPhase4Integrate(dt);

        rkUpdate(dt);

        t+=dt;
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



