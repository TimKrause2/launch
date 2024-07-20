#include <math.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "gsimms.h"
#include "orbit.h"
#include "image.h"

Body::Body(double mass, double radius, char *tex_file)
    :mass(mass),
      radius(radius),
      Jp(caams::J_p_sphere(mass,radius))
{
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
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

    glTexImage2D(
                GL_TEXTURE_2D, 0, GL_RGB,
                tex_image->width, tex_image->height,
                0,
                GL_RGB, GL_UNSIGNED_BYTE,
                (const GLvoid*)tex_image->data );
    image_free(tex_image);

}

State::State(Body* body)
    :body(body),
      position(3,1),
      velocity(3,1),
      p(4,1),
      pdot(4,1),
      a_thrust(3,1,caams::zeros),
      r_thrust(3,1,caams::zeros),
      rk_position(3,1),
      rk_acceleration(3,1),
      m_kr(3,4),
      m_kv(3,4)
{

}

State::State(State &o) :
    body(o.body),
    position(o.position),
    velocity(o.velocity),
    p(o.p),
    pdot(o.pdot),
    a_thrust(3,1,caams::zeros),
    r_thrust(3,1,caams::zeros),
    rk_position(3,1),
    rk_acceleration(3,1),
    m_kr(3,4),
    m_kv(3,4)
{

}



glm::dvec3 spherical( double latitude, double longitude )
{
    glm::dvec3 r(cos(latitude)*cos(longitude),
                             cos(latitude)*sin(longitude),
                             sin(latitude));
    return r;
}

void uv_sphere( int n_latitude, int n_longitude )
{
    double latitude = M_PI/2;
    double dlatitude = -M_PI/n_latitude;
    double dlongitude = 2*M_PI/n_longitude;
    double du=1.0/n_longitude;
    double dv=1.0/n_latitude;
    double v=0;
    int lat,lon;
    glBegin(GL_QUADS);
    for(lat=0;lat<n_latitude;lat++){
        double longitude = -M_PI;
        double u=0.0;
        glm::dvec2 uv00( u, v );
        glm::dvec2 uv10( u+du, v );
        glm::dvec2 uv11( u+du, v+dv );
        glm::dvec2 uv01( u, v+dv );
        glm::dvec3 v00=spherical(latitude,longitude);
        glm::dvec3 v10=spherical(latitude,longitude+dlongitude);
        glm::dvec3 v11=spherical(latitude+dlatitude,longitude+dlongitude);
        glm::dvec3 v01=spherical(latitude+dlatitude,longitude);
        for(lon=0;lon<n_longitude;lon++){
            glTexCoord2dv( glm::value_ptr( uv00 ) );
            glVertex3dv( glm::value_ptr( v00 ) );
            glTexCoord2dv( glm::value_ptr( uv10 ) );
            glVertex3dv( glm::value_ptr( v10 ) );
            glTexCoord2dv( glm::value_ptr( uv11 ) );
            glVertex3dv( glm::value_ptr( v11 ) );
            glTexCoord2dv( glm::value_ptr( uv01 ) );
            glVertex3dv( glm::value_ptr( v01 ) );
            u+=du;
            longitude+=dlongitude;
            uv00 = uv10;
            uv01 = uv11;
            uv10 = glm::dvec2( u+du, v );
            uv11 = glm::dvec2( u+du, v+dv );
            v00 = v10;
            v01 = v11;
            v10 = spherical(latitude,longitude+dlongitude);
            v11 = spherical(latitude+dlatitude,longitude+dlongitude);
        }
        v+=dv;
        latitude+=dlatitude;
    }
    glEnd();
}


void State::draw(void)
{
    // first set up the transformation matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    double radius = body->radius;
    glm::dvec3 scale(radius,radius,radius);
    glm::dmat4 A_scale = glm::scale(glm::dmat4(1.0),scale);

    caams::matrix A(caams::Ap(p));
    glm::dmat3 A_glm = glm::make_mat3((~A).data);
    glm::dmat4 A_rotation(A_glm);

    glm::dvec3 translation(position.data[0],position.data[1],position.data[2]);
    glm::dmat4 A_translation = glm::translate(glm::dmat4(1.0),translation);

    glm::dmat4 A_body = A_translation*A_rotation*A_scale;

    glMultMatrixd( glm::value_ptr( A_body ) );


    // set up texturing
    glEnable( GL_TEXTURE_2D );
    glEnable(GL_DEPTH_TEST);
    glEnable( GL_CULL_FACE );
    glFrontFace( GL_CW );
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );
    glBindTexture( GL_TEXTURE_2D, body->texName );
    uv_sphere( SPHERE_N_LATITUDE, SPHERE_N_LONGITUDE );
    glDisable( GL_TEXTURE_2D );
    glDisable( GL_DEPTH_TEST );

    glPopMatrix();
}


caams::matrix p_ddot_solve(
        caams::matrix p,
        caams::matrix p_dot,
        caams::matrix Jp)
{
    caams::matrix H(4.0*~caams::L(p_dot)*Jp*caams::L(p));
    caams::matrix R(4,4);
    R.sub(caams::L(p)*H,1,1);
    R.sub(~p_dot,4,1);

    caams::matrix y(-R*p_dot);

    caams::matrix A(4,4);
    A.sub(2.0*Jp*caams::L(p),1,1);
    A.sub(~p,4,1);

    caams::matrix x(A.inverse()*y);

    return x;
}

void normalize_p(caams::matrix &p){
    p = (1.0/caams::norm(p))*p;
}

void State::IntegrateRotation( double dt )
{
    caams::matrix k_p_dot(4,4);
    caams::matrix k_p_ddot(4,4);
    caams::matrix p_norm(4,1);
    k_p_dot.sub(pdot,1,1);
    k_p_ddot.sub(p_ddot_solve(p,k_p_dot.sub(4,1,1,1),body->Jp),1,1);

    p_norm = p + (dt/2.0)*k_p_dot.sub(4,1,1,1);
    normalize_p(p_norm);
    k_p_dot.sub(pdot+(dt/2.0)*k_p_ddot.sub(4,1,1,1),1,2);
    k_p_ddot.sub(p_ddot_solve(p_norm,k_p_dot.sub(4,1,1,2),body->Jp),1,2);

    p_norm = p + (dt/2.0)*k_p_dot.sub(4,1,1,2);
    normalize_p(p_norm);
    k_p_dot.sub(pdot+(dt/2.0)*k_p_ddot.sub(4,1,1,2),1,3);
    k_p_ddot.sub(p_ddot_solve(p_norm,k_p_dot.sub(4,1,1,3),body->Jp),1,3);

    p_norm = p + dt*k_p_dot.sub(4,1,1,3);
    normalize_p(p_norm);
    k_p_dot.sub(pdot+dt*k_p_ddot.sub(4,1,1,3),1,4);
    k_p_ddot.sub(p_ddot_solve(p_norm,k_p_dot.sub(4,1,1,4),body->Jp),1,4);

    p = p + (dt/6.0)*k_p_dot*caams::matrix(4,1,1.0,2.0,2.0,1.0);
    normalize_p(p);
    pdot = pdot + (dt/6.0)*k_p_ddot*caams::matrix(4,1,1.0,2.0,2.0,1.0);
}

void System::AddBody(Body *body)
{
    bodies.emplace_back(body);
}

State& System::GetState(Body *body)
{
    std::list<State>::iterator it;
    for(it=bodies.begin();it!=bodies.end();it++){
        if((*it).body == body){
            return *it;
        }
    }
    return *it;
}

System System::Duplicate(void)
{
    System s;
    std::list<State>::iterator it;
    for(it=bodies.begin();it!=bodies.end();it++){
        s.bodies.emplace_back(*it);
    }
    return s;
}

void System::draw(void)
{
    std::list<State>::iterator it;
    for(it=bodies.begin();it!=bodies.end();it++){
        (*it).draw();
    }
}

void System::rkIntegrate(double dt)
{
    rkPhase1Positions();
    rkAccelerations();
    rkPhase1Integrate(dt);

    rkPhase2Positions(dt);
    rkAccelerations();
    rkPhase2Integrate(dt);

    rkPhase3Positions(dt);
    rkAccelerations();
    rkPhase3Integrate(dt);

    rkPhase4Positions(dt);
    rkAccelerations();
    rkPhase4Integrate(dt);

    update_rotations(dt);

}

void System::rkAccelerations( void )
{
    std::list<State>::iterator l_it_b1;
    std::list<State>::iterator l_it_b2;
    // set all accelerations to zero
    for( l_it_b1=bodies.begin(); l_it_b1!=bodies.end(); l_it_b1++){
        (*l_it_b1).rk_acceleration = (*l_it_b1).a_thrust;
    }
    // calculate accelerations for all pairs of bodies
    for( l_it_b1=bodies.begin(); l_it_b1!=bodies.end(); l_it_b1++){
        l_it_b2 = l_it_b1;
        l_it_b2++;
        for( ; l_it_b2!=bodies.end(); l_it_b2++){
            // compute the normal vector of body 1 relative to 2
            caams::matrix  l_r12( (*l_it_b1).rk_position - (*l_it_b2).rk_position );
            double l_magr = caams::norm(l_r12);
            caams::matrix l_n12(l_r12 * (1.0/l_magr));
            l_magr*=l_magr;
            // normal vector of body 2 relative to 1
            caams::matrix l_n21( -l_n12 );
            double l_g=GM/l_magr;
            double l_a1 = (*l_it_b2).body->mass * l_g;
            double l_a2 = (*l_it_b1).body->mass * l_g;

            (*l_it_b1).rk_acceleration += l_a1*l_n21;
            (*l_it_b2).rk_acceleration += l_a2*l_n12;
        }
    }
}

void System::rkPhase1Positions( void )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).rk_position = (*l_it).position;
    }
}

void System::rkPhase2Positions( double p_dt )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).rk_position = (*l_it).position + (*l_it).m_kr.sub(3,1,1,1)*(p_dt*0.5);
    }
}

void System::rkPhase3Positions( double p_dt )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).rk_position = (*l_it).position + (*l_it).m_kr.sub(3,1,1,2)*(p_dt*0.5);
    }
}

void System::rkPhase4Positions( double p_dt )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).rk_position = (*l_it).position + (*l_it).m_kr.sub(3,1,1,3)*p_dt;
    }
}

void System::rkPhase1Integrate( double p_dt )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).m_kr.sub( (*l_it).velocity, 1, 1);
        (*l_it).m_kv.sub( (*l_it).rk_acceleration, 1, 1);
    }
}

void System::rkPhase2Integrate( double p_dt )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).m_kr.sub( (*l_it).velocity + (*l_it).m_kv.sub(3,1,1,1)*(p_dt*0.5), 1, 2 );
        (*l_it).m_kv.sub( (*l_it).rk_acceleration, 1, 2);
    }
}

void System::rkPhase3Integrate( double p_dt )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).m_kr.sub( (*l_it).velocity + (*l_it).m_kv.sub(3,1,1,2)*(p_dt*0.5), 1, 3 );
        (*l_it).m_kv.sub( (*l_it).rk_acceleration, 1, 3);
    }
}

void System::rkPhase4Integrate( double p_dt )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).m_kr.sub( (*l_it).velocity + (*l_it).m_kv.sub(3,1,1,3)*p_dt, 1, 4);
        (*l_it).m_kv.sub( (*l_it).rk_acceleration, 1, 4);
        (*l_it).velocity += p_dt*(*l_it).m_kv*caams::matrix(4,1,1.0/6.0,2.0/6.0,2.0/6.0,1.0/6.0);
        (*l_it).position += p_dt*(*l_it).m_kr*caams::matrix(4,1,1.0/6.0,2.0/6.0,2.0/6.0,1.0/6.0);
    }
}

void System::update_rotations( double p_dt )
{
    std::list<State>::iterator l_it;
    for(l_it=bodies.begin(); l_it!=bodies.end(); l_it++){
        (*l_it).IntegrateRotation( p_dt );
    }
}





