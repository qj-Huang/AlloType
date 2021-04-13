#include "method.h"

VectorXd rotate(VectorXd coord, Vector3d axis, double angle)
{
	if (axis.dot(axis) != 1.0)
		axis /= sqrt(axis.dot(axis));

	Matrix3d rotmat = Matrix3d::Zero();
	rotmat(1, 0) = axis(2);
	rotmat(0, 2) = axis(1);
	rotmat(2, 1) = axis(0);
	rotmat(0, 1) = -rotmat(1, 0);
	rotmat(2, 0) = -rotmat(0, 2);
	rotmat(1, 2) = -rotmat(2, 1);

	rotmat *= sin(angle);
	rotmat += (1 - cos(angle)) * (axis * axis.transpose());
	rotmat += cos(angle) * Matrix3d::Identity();

	Map<Matrix3Xd> xyz(coord.data(), 3, coord.size() / 3);
	Vector3d center = xyz.rowwise().mean();
	xyz.colwise() -= center;

	Matrix3Xd rot_xyz = rotmat * xyz;
	rot_xyz.colwise() += center;
	Map<VectorXd> rot_coord(rot_xyz.data(), rot_xyz.size());

	return rot_coord;
}

VectorXd fitting(VectorXd coord1, VectorXd coord2)
{
	if (coord1.size() == coord2.size())
	{
		size_t vlen = coord1.size() / 3;
		Map<Matrix3Xd> xyz1(coord1.data(), 3, vlen);
		Map<Matrix3Xd> xyz2(coord2.data(), 3, vlen);
		Vector3d c1 = xyz1.rowwise().mean();
		Vector3d c2 = xyz2.rowwise().mean();
		xyz1.colwise() -= c1;
		xyz2.colwise() -= c2;

		Matrix3d r_mat = Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 3; ++j)
				for (size_t k = 0; k < vlen; ++k)
					r_mat(i, j) += xyz2(i, k) * xyz1(j, k);

		Matrix3d sym_mat = r_mat.transpose() * r_mat;

		SelfAdjointEigenSolver<Matrix3d> eigensolver(sym_mat);
		Vector3d eigenvalues = eigensolver.eigenvalues();
		Matrix3d eigenvectors = eigensolver.eigenvectors();

		Matrix3d b_mat = Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			b_mat.row(i) = r_mat * eigenvectors.col(i) / sqrt(eigenvalues[i]);

		Matrix3d a_mat = eigenvectors.transpose();

		Matrix3d u_mat = Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 3; ++j)
				for (size_t k = 0; k < 3; ++k)
					u_mat(i, j) += b_mat(k, j) * a_mat(k, i);

		Matrix3Xd fit_xyz2 = u_mat * xyz2;
		fit_xyz2.colwise() += c1;
		Map<VectorXd> fit_coord2(fit_xyz2.data(), fit_xyz2.size());
		return fit_coord2;
	}
	else
		return VectorXd();
}

VectorXd calc_displacement(VectorXd coord1, VectorXd coord2)
{
	if (coord1.size() == coord2.size())
		return coord2 - coord1;
	else
		return VectorXd();
}

double calc_rmsd(VectorXd coord)
{
	size_t vlen = coord.size() / 3;
	Map<Array3Xd> diffxyz(coord.data(), 3, vlen);
	VectorXd sqdist = diffxyz.pow(2).colwise().sum();
	return sqrt(sqdist.sum() / vlen);
}

double calc_rmsd(VectorXd coord1, VectorXd coord2)
{
	if (coord1.size() == coord2.size())
	{
		VectorXd diff = coord1 - coord2;
		size_t vlen = diff.size() / 3;

		Map<Array3Xd> diffxyz(diff.data(), 3, vlen);

		VectorXd sqdist = diffxyz.pow(2).colwise().sum();
		return sqrt(sqdist.sum() / vlen);
	}
	else
		return 0.0;
}

double calc_mindist(VectorXd coord1, VectorXd coord2)
{
	Map<Matrix3Xd> xyz1(coord1.data(), 3, coord1.size() / 3);
	Map<Matrix3Xd> xyz2(coord2.data(), 3, coord2.size() / 3);

	VectorXd dist(xyz1.cols());
	Vector3d onexyz = Vector3d::Zero();
	Array3Xd diffxyz(3, xyz2.cols());
	for (size_t i = 0; i < size_t(xyz1.cols()); ++i)
	{
		onexyz = xyz1.col(i);
		diffxyz = xyz2.colwise() - onexyz;
		dist(i) = diffxyz.pow(2).colwise().sum().minCoeff();
	}
	return dist.minCoeff();
}

double calc_norm(VectorXd vector)
{
	ArrayXd v = vector;
	return sqrt(v.pow(2).sum());
}

double calc_average(VectorXd force)
{
	Map<Array3Xd> forcexyz(force.data(), 3, force.size() / 3);
	VectorXd sqdist = forcexyz.pow(2).colwise().sum().sqrt();
	return sqdist.sum() / sqdist.size();
}

MatrixXd gen_distmat(VectorXd coord)
{
	size_t resn = coord.size() / 3;
	MatrixXd distmat = MatrixXd::Zero(resn, resn);
	for (size_t i = 0; i < resn; ++i)
		for (size_t j = i + 1; j < resn; ++j)
			distmat(j, i) = distmat(i, j) = sqrt(pow(coord(3 * i) - coord(3 * j), 2) + pow(coord(3 * i + 1) - coord(3 * j + 1), 2) + pow(coord(3 * i + 2) - coord(3 * j + 2), 2));
	return distmat;
}

MatrixXd gen_differ(MatrixXd holo, MatrixXd apo)
{
	size_t len = apo.col(1).size();
	MatrixXd diff = MatrixXd::Zero(len, len);
	for (size_t i = 0; i < len; i++)
		for (size_t j = i; j < len; j++)
			if (apo(i, j) != 0)
				diff(i, j) = holo(i, j) - apo(i, j);
	return diff;
}

list<size_t> gen_pocket(double cutoff, VectorXd dist2ligand)
{
	list<size_t> pocket;
	if (cutoff > dist2ligand.minCoeff())
	{
		for (size_t i = 0; i < size_t(dist2ligand.size()); ++i)
			if (dist2ligand(i) < cutoff)
				pocket.push_back(i);
	}
	else
		handle_error(
			format("Given cutoff is too short. Minimum possible cutoff is %1$.2f.") % dist2ligand.minCoeff()
		);
	return pocket;
}

void normal_equation(VectorXd &coeff, MatrixXd X, VectorXd Y)
{
	coeff = (X.transpose() * X).inverse() * X.transpose() * Y;
}
