#include "sequencer.h"
#include "orbit.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

void Sequencer::load_file(const char *filename, double **data, int *Npoints)
{
    std::ifstream  fstr;
    double value;

    fstr.open(filename);
    if(!fstr){
        std::cout << "Sequencer: Couldn't open file \"" << filename << "\"" << std::endl;
        exit(EXIT_FAILURE);
    }

    // count the number of data elements
    int N=0;
    while(fstr >> value){
        N++;
    }

    *Npoints = N;

    // allocate the memory for the elements
    *data = new double[N];

    // rewind and read in the elements
    fstr.clear();
    fstr.seekg(0, fstr.beg);

    double *dst = *data;
    while(fstr >> value){
        *(dst++) = value;
    }

    fstr.close();
}

Sequencer::Sequencer(void)
{
    int N;
    load_file("beta.dat", &beta_data, &N);
    Npoints = N;
    load_file("time.dat", &time_data, &N);
    if(N!=Npoints){
        std::cout << "Sequencer: Number of points in data sets don't match." << std::endl;
        exit(EXIT_FAILURE);
    }
    i0=0;
    i1=1;
}

Sequencer::~Sequencer()
{
    delete [] beta_data;
    delete [] time_data;
}


double Sequencer::GetBeta(double t)
{
    // position the current segment to include t
    // move the segment forward if required
    while(time_data[i1]<=t){
        if(i1!=(Npoints-1)){
            i0++;
            i1++;
        }else{
            printf("Sequencer::GetBeta:extrapolating\n");
            printf("delta:%le\n", t - time_data[i1]);
            break;
        }
    }

    // move the segment back if required
    while(time_data[i0]>t){
        if(i0!=0){
            i0--;
            i1--;
        }else{
            break;
        }
    }

    // interpolate beta in the segment
    double alpha = (t-time_data[i0])/(time_data[i1]-time_data[i0]);

    return beta_data[i0]*(1-alpha) + beta_data[i1]*alpha;

}

double Sequencer::TotalTime()
{
    return time_data[Npoints-1];
}

int Sequencer::GetNpoints(void)
{
    return Npoints;
}

double Sequencer::GetTime(int index)
{
    return time_data[index];
}

DirCoordSequencer::DirCoordSequencer(
    Eigen::MatrixXd control,
    Eigen::RowVectorXd time):
    control(control),
    time(time)
{
    if(control.cols()!=time.cols()){
        std::cout << "DirCoordSequencer: number of columns don't match." << std::endl;
        exit(EXIT_FAILURE);
    }
    if(control.rows()!=2){
        std::cout << "DirCoordSequencer: control doesn't have 2 rows." << std::endl;
        exit(EXIT_FAILURE);
    }

    i0=0;
    i1=1;
    Npoints = control.cols();
}

Eigen::Vector3d DirCoordSequencer::GetDir(double t)
{
    // position the current segment to include t
    // move the segment forward if required
    while(time(i1)<=t){
        if(i1!=(Npoints-1)){
            i0++;
            i1++;
        }else{
            printf("DirCoordSequencer::GetDir:extrapolating!!\n");
            printf("delta:%le\n",t - time(i1));
            break;
        }
    }

    // move the segment back if required
    while(time(i0)>t){
        if(i0!=0){
            i0--;
            i1--;
        }else{
            break;
        }
    }

    // interpolate between the two direction vectors
    Eigen::Vector2d control0 = control.col(i0);
    Eigen::Vector2d control1 = control.col(i1);
    Eigen::Vector3d r0 = direction_from_coord(control0(0),control0(1));
    Eigen::Vector3d r1 = direction_from_coord(control1(0),control1(1));
    Eigen::Vector3d axis = r0.cross(r1);
    double mag_axis = axis.squaredNorm();
    if(mag_axis == 0.0){
        return r0;
    }else{
        double theta = acos(r0.dot(r1));
        double alpha = (t - time(i0))/(time(i1)-time(i0));
        return caams::AAA(theta*alpha, axis)*r0;
    }
}

double DirCoordSequencer::TotalTime(void)
{
    return time(Npoints-1);
}

int DirCoordSequencer::GetNpoints(void)
{
    return Npoints;
}

double DirCoordSequencer::GetTime(int index)
{
    return time(index);
}





DirCoord3dSequencer::DirCoord3dSequencer(
    Eigen::MatrixXd control,
    Eigen::RowVectorXd time):
    control(control),
    time(time)
{
    if(control.cols()!=time.cols()){
        std::cout << "DirCoord3dSequencer: number of columns don't match." << std::endl;
        exit(EXIT_FAILURE);
    }
    if(control.rows()!=3){
        std::cout << "DirCoord3dSequencer: control doesn't have 3 rows." << std::endl;
        exit(EXIT_FAILURE);
    }

    i0=0;
    i1=1;
    Npoints = control.cols();
}

Eigen::Vector3d DirCoord3dSequencer::GetDir(double t)
{
    // position the current segment to include t
    // move the segment forward if required
    while(time(i1)<=t){
        if(i1!=(Npoints-1)){
            i0++;
            i1++;
        }else{
            printf("DirCoord3dSequencer::GetDir:extrapolating!!\n");
            printf("delta:%le\n",t - time(i1));
            break;
        }
    }

    // move the segment back if required
    while(time(i0)>t){
        if(i0!=0){
            i0--;
            i1--;
        }else{
            break;
        }
    }

    // interpolate between the two direction vectors
    Eigen::Vector3d r0 = control.col(i0);
    Eigen::Vector3d r1 = control.col(i1);
    Eigen::Vector3d axis = r0.cross(r1);
    double mag_axis = axis.squaredNorm();
    if(mag_axis == 0.0){
        return r0;
    }else{
        double theta = acos(r0.dot(r1));
        double alpha = (t - time(i0))/(time(i1)-time(i0));
        return caams::AAA(theta*alpha, axis)*r0;
    }
}

double DirCoord3dSequencer::TotalTime(void)
{
    return time(Npoints-1);
}

int DirCoord3dSequencer::GetNpoints(void)
{
    return Npoints;
}

double DirCoord3dSequencer::GetTime(int index)
{
    return time(index);
}



ThrustVectorSequencer::ThrustVectorSequencer(
    Eigen::MatrixXd control,
    Eigen::RowVectorXd time):
    control(control),
    time(time)
{
    if(control.cols()!=time.cols()){
        std::cout << "ThrustVectorSequencer: number of columns don't match." << std::endl;
        exit(EXIT_FAILURE);
    }
    if(control.rows()!=2){
        std::cout << "ThrustVectorSequencer: control doesn't have 2 rows." << std::endl;
        exit(EXIT_FAILURE);
    }

    i0=0;
    i1=1;
    Npoints = control.cols();
}

Eigen::Vector3d ThrustVectorSequencer::TVector(Eigen::Vector2d theta)
{
    return Eigen::Vector3d(cos(theta(0))*sin(theta(1)),
                           -sin(theta(0)),
                           cos(theta(0))*cos(theta(1)));
}

Eigen::Vector3d ThrustVectorSequencer::GetDir(double t)
{
    // position the current segment to include t
    // move the segment forward if required
    while(time(i1)<=t){
        if(i1!=(Npoints-1)){
            i0++;
            i1++;
        }else{
            printf("ThrustVectorSequencer::GetDir:extrapolating!!\n");
            printf("delta:%le\n",t - time(i1));
            break;
        }
    }

    // move the segment back if required
    while(time(i0)>t){
        if(i0!=0){
            i0--;
            i1--;
        }else{
            break;
        }
    }

    // interpolate between the two direction vectors
    Eigen::Vector2d control0 = control.col(i0);
    Eigen::Vector2d control1 = control.col(i1);
    Eigen::Vector3d r0 = TVector(control0);
    Eigen::Vector3d r1 = TVector(control1);
    Eigen::Vector3d axis = r0.cross(r1);
    double mag_axis = axis.squaredNorm();
    if(mag_axis == 0.0){
        return r0;
    }else{
        double theta = acos(r0.dot(r1));
        double alpha = (t - time(i0))/(time(i1)-time(i0));
        return caams::AAA(theta*alpha, axis)*r0;
    }
}

double ThrustVectorSequencer::TotalTime(void)
{
    return time(Npoints-1);
}

int ThrustVectorSequencer::GetNpoints(void)
{
    return Npoints;
}

double ThrustVectorSequencer::GetTime(int index)
{
    return time(index);
}

