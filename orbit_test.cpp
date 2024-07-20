#include "orbit.h"
#include <iostream>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <ctime>

void timespec_sub(std::timespec const &ts1,
                  std::timespec const &ts0,
                  std::timespec &dts)
{
    dts.tv_sec = ts1.tv_sec - ts0.tv_sec;
    dts.tv_nsec = ts1.tv_nsec - ts0.tv_nsec;
    if(dts.tv_nsec < 0){
        dts.tv_sec--;
        dts.tv_nsec += 1000000000;
    }
}


#define RAD_PER_DEG (M_PI/180)
#define LATITUDE 28.608530289578532
#define LONGITUDE -80.60392779821008
int main(void)
{
    OrbitalParams params;
    params.a = 10.0;
    params.e = 0.0001;
    params.i = 45 * RAD_PER_DEG;
    params.Omega = 270.0 * RAD_PER_DEG;
    params.omega = 0.0;
    params.nu = 0.0;

    Eigen::Vector3d rs = n_r_earth_frame(0 * RAD_PER_DEG, 45 * RAD_PER_DEG);
    Eigen::Vector3d ri1;
    Eigen::Vector3d ri2;

    for(int Omega=0;Omega<360;Omega++)
    {
        std::cout << "Omega:" << Omega << std::endl;
        params.Omega = Omega*RAD_PER_DEG;

        bool iflag = orbit_ring_intersect(params, rs, ri1, ri2);

        if(iflag){
            std::cout << "intersection." << std::endl;
            std::cout << "ri1:\n" << ri1 << "\nri2:\n" << ri2 << std::endl;
            std::cout << "longitude of 1:" << atan2(ri1(1), ri1(0)) / RAD_PER_DEG << std::endl;
            std::cout << "longitude of 2:" << atan2(ri2(1), ri2(0)) / RAD_PER_DEG << std::endl;
            std::cout << "mag of 1:" << ri1.norm() << std::endl;
            std::cout << "mag of 2:" << ri2.norm() << std::endl;
            Eigen::Vector3d n = orbit_n_plane(params);
            std::cout << "ri1 dot N:" << ri1.dot(n) << std::endl;
            std::cout << "ri2 dot N:" << ri2.dot(n) << std::endl;
        }else{
            std::cout << "no intersection." << std::endl;
        }
    }

    std::timespec ts0;
    std::timespec ts1;
    std::timespec dts;
    std::timespec_get(&ts0, TIME_UTC);

    for(int i=0;i<1000;i++){
        double k = (double)i/1000;
        boost::math::ellint_1(k);
    }

    std::timespec_get(&ts1, TIME_UTC);
    timespec_sub(ts1, ts0, dts);
    double delta_t = dts.tv_sec + dts.tv_nsec*1e-9;
    delta_t/=1000;
    std::cout << "ellint_1(k) time :" << delta_t << " seconds" << std::endl;

    return 0;
}
