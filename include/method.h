#pragma once
#include <list>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/format.hpp>

#include "handle_io.h"

using namespace std;
using namespace Eigen;
using boost::format;

VectorXd rotate(VectorXd coord, Vector3d axis, double angle);

VectorXd fitting(VectorXd coord1, VectorXd coord2);

VectorXd calc_displacement(VectorXd coord1, VectorXd coord2);

double calc_rmsd(VectorXd coord);

double calc_rmsd(VectorXd coord1, VectorXd coord2);

double calc_mindist(VectorXd coord1, VectorXd coord2);

double calc_norm(VectorXd vector);

double calc_average(VectorXd force);

MatrixXd gen_distmat(VectorXd coord);

MatrixXd gen_differ(MatrixXd X, MatrixXd Y);

list<size_t> gen_pocket(double cutoff, VectorXd dist2ligand);

void normal_equation(VectorXd & coefficient, MatrixXd X, VectorXd Y);
