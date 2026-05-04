#include "dormand_prince.h"
#include <cmath>

DormandPrince::DormandPrince()
{
    c << 0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0;
    a2 =  1.0/5.0;
    a3 << 3.0/40.0,       9.0/40.0;
    a4 << 44.0/45.0,      -56.0/15.0,      32.0/9.0;
    a5 << 19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0;
    a6 << 9017.0/3168.0,  -355.0/33.0,     46732.0/5247.0, 49.0/176.0,  -5103.0/18656.0;
    a7 << 35.0/384.0,      0.0,            500.0/1113.0,   125.0/192.0, -2187.0/6784.0,    11.0/84.0;
    bstar
       << 5179.0/57600.0,  0.0,            7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0;
}

void DormandPrince::integrate_step(
        Eigen::VectorXd &y1,
        Eigen::VectorXd y0,
        std::function<Eigen::VectorXd(double, Eigen::VectorXd)> dy_func,
        double t,
        double dt,
        double &error)
{
    Eigen::MatrixXd k(y0.rows(), 7);
    k.col(0) = dy_func(t + c(0)*dt, y0);
    k.col(1) = dy_func(t + c(1)*dt, y0 + k.col(0)*a2*dt);
    k.col(2) = dy_func(t + c(2)*dt, y0 + k.leftCols(2)*a3*dt);
    k.col(3) = dy_func(t + c(3)*dt, y0 + k.leftCols(3)*a4*dt);
    k.col(4) = dy_func(t + c(4)*dt, y0 + k.leftCols(4)*a5*dt);
    k.col(5) = dy_func(t + c(5)*dt, y0 + k.leftCols(5)*a6*dt);
    Eigen::VectorXd dy1(y0.rows());
    dy1 = k.leftCols(6)*a7*dt;
    y1 = y0 + dy1;
    k.col(6) = dy_func(t + c(6)*dt, y1);
    Eigen::VectorXd dy1p(y0.rows());
    dy1p = k*bstar*dt;
    error = 0.0;
    for(int i=0;i<y0.rows();i++)
    {
        if(dy1(i)!=0.0){
            double this_error = fabs((dy1(i)-dy1p(i))/dy1(i));
            if(this_error > error){
                error = this_error;
            }
        }
    }
}

void DormandPrince::integrate(
        Eigen::VectorXd &y1,
        const Eigen::VectorXd &y0_first,
        std::function<Eigen::VectorXd(double, Eigen::VectorXd)> dy_func,
        double t,
        double dt)
{
    Eigen::VectorXd y0(y0_first.rows());
    y0 = y0_first;
    double dt_sign = (dt>0.0)?1.0:-1.0;
    double dt_total = fabs(dt);
    double dt_step = dt_total;
    double dt_current = 0.0;
    while(dt_current < dt_total){
        if(dt_current+dt_step > dt_total){
            dt_step = dt_total - dt_current;
        }
        double error;
        integrate_step(y1, y0, dy_func, t + dt_sign*dt_current, dt_sign*dt_step, error);
        if(error>1e-11){
            dt_step *= 0.5;
        }else{
            dt_current += dt_step;
            y0 = y1;
            if(error < 2e-12){
                dt_step *= 2.0;
            }
        }
    }
}
