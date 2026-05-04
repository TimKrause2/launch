#include <Eigen/Dense>
#include <functional>

class Integrator
{
public:
    virtual void integrate(
            Eigen::VectorXd &y1,
            const Eigen::VectorXd &y0,
            std::function<Eigen::VectorXd(double, Eigen::VectorXd)> dy_func,
            double t,
            double dt) = 0;
};
