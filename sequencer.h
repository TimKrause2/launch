#include "caams.hpp"

class Sequencer
{
    double *beta_data;
    double *time_data;
    int Npoints;
    void load_file(const char *filename, double **data, int *Npoints);
    int i0;
    int i1;
public:
    Sequencer(void);
    ~Sequencer();
    double GetBeta(double t);
    double TotalTime(void);
    int GetNpoints(void);
    double GetTime(int index);
};

class DirCoordSequencer
{
    Eigen::MatrixXd control;
    Eigen::RowVectorXd time;
    int i0;
    int i1;
    int Npoints;
public:
    DirCoordSequencer(
        Eigen::MatrixXd control,
        Eigen::RowVectorXd time);
    Eigen::Vector3d GetDir(double t);
    double TotalTime(void);
    int GetNpoints(void);
    double GetTime(int index);
};

class DirCoord3dSequencer
{
    Eigen::MatrixXd control;
    Eigen::RowVectorXd time;
    int i0;
    int i1;
    int Npoints;
public:
    DirCoord3dSequencer(
        Eigen::MatrixXd control,
        Eigen::RowVectorXd time);
    Eigen::Vector3d GetDir(double t);
    double TotalTime(void);
    int GetNpoints(void);
    double GetTime(int index);
};

class ThrustVectorSequencer
{
    Eigen::MatrixXd control;
    Eigen::RowVectorXd time;
    int i0;
    int i1;
    int Npoints;
public:
    ThrustVectorSequencer(
        Eigen::MatrixXd control,
        Eigen::RowVectorXd time);
    Eigen::Vector3d TVector(Eigen::Vector2d theta);
    Eigen::Vector3d GetDir(double t);
    double TotalTime(void);
    int GetNpoints(void);
    double GetTime(int index);
};

