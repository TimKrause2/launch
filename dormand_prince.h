#include "integrator.h"

class DormandPrince : public Integrator
{
private:
    Eigen::Matrix<double,7,1> c;
    double a2;
    Eigen::Matrix<double,2,1> a3;
    Eigen::Matrix<double,3,1> a4;
    Eigen::Matrix<double,4,1> a5;
    Eigen::Matrix<double,5,1> a6;
    Eigen::Matrix<double,6,1> a7; // b is the same as a7
    Eigen::Matrix<double,7,1> bstar;
    void integrate_step(
            Eigen::VectorXd &y1,
            Eigen::VectorXd y0,
            std::function<Eigen::VectorXd(double, Eigen::VectorXd)> dy_func,
            double t,
            double dt,
            double &error);
public:
    DormandPrince();
    void integrate(
            Eigen::VectorXd &y1,
            const Eigen::VectorXd &y0,
            std::function<Eigen::VectorXd(double, Eigen::VectorXd)> dy_func,
            double t,
            double dt);
};
