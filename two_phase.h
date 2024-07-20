#include "psopt.h"
#include "orbit.h"

void two_phase_optimize(
    ModEquOrbitalParams &meparams0,
    ModEquOrbitalParams &meparams1,
    Eigen::MatrixXd &phase1_time, // initial time vector for phase 1 and result
    Eigen::MatrixXd &phase2_time, // initial time vector for phase 2 and result
    Eigen::MatrixXd &phase2_control, // control guess and result on return
    double mu,
    double thrust);


#define two_phase_optimize_constant_longitude two_phase_optimize

//void two_phase_optimize_constant_longitude(
//    ModEquOrbitalParams &meparams0,
//    ModEquOrbitalParams &meparams1,
//    Eigen::MatrixXd &phase1_time, // initial time vector for phase 1 and result
//    Eigen::MatrixXd &phase2_time, // initial time vector for phase 2 and result
//    Eigen::MatrixXd &phase2_control, // control guess and result on return
//    double mu,
//    double thrust);

