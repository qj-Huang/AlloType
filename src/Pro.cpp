#include "Pro.h"

Pro::Pro()
{
}

Pro::Pro(string fpath, bool has_ligand_flag, vector<string> exclude, double k, double cutoff)
{
	with_ligand_flag = has_ligand_flag;
	for (vector<string>::iterator it = exclude.begin(); it != exclude.end(); ++it)
		exclres.emplace(*it);
	k_default = k_inter = k_intra = k;
	cutoff_inter = cutoff_intra = cutoff;
	handle_info(format("Spring constant = %1$.1f Kcal/(mol A^2).") % k);
	handle_info(format("Cutoff = %1$.1f A.") % cutoff);
	read(fpath);
	handle_info(format("Successfully loaded protein at path %1%") % fpath);
	gen_coord();
	gen_distmat();
	gen_contact();
	handle_info("Coordinate matrix, distance matrix, contact map have been generated for this protein.");
	if (with_ligand_flag)
	{
		gen_dist2ligand();
		handle_info("Residue distance to ligand has been calculated.");
	}
}

Pro::~Pro()
{
}

bool Pro::has_ligand()
{
	return with_ligand_flag;
}

void Pro::read(string fpath)
{
	string line;
	ifstream pdb(fpath);
	size_t proid = 0, proatomid = 0, ligandatomid = 0;
	string prev_chain = "";
	size_t prev_resid = 0;
	if (pdb.is_open())
	{
		while (getline(pdb, line))
		{
			string record = read_record(line);
			if (record == "ATOM" || record == "HETATM")
			{
				string resname = read_resname(line);
				if (prores.find(resname) != prores.end())
				{
					if (prev_resid == 0 && prev_chain == "")
					{
						prev_resid = read_resid(line);
						prev_chain = read_chain(line);
					}
					if (prev_resid != read_resid(line) || prev_chain != read_chain(line))
						++proid;
					if (read_atomname(line) == "CA")
					{
						ResInfo res = read_res(line);
						pro[proid] = res;
					}
					if (proatoms.find(proid) != proatoms.end())
						proatoms[proid].push_back(read_atom(line));
					else
					{
						vector<AtomInfo> grp = { read_atom(line) };
						proatoms[proid] = grp;
					}
					prev_resid = read_resid(line);
					prev_chain = read_chain(line);
					++proatomid;
				}
				else if (exclres.find(resname) == exclres.end())
				{
					if (ligand.find(resname) != ligand.end())
						ligand[resname].push_back(read_atom(line));
					else
					{
						vector<AtomInfo> grp = { read_atom(line) };
						ligand[resname] = grp;
					}
					++ligandatomid;
				}
				else
				{
					AtomInfo ex = read_atom(line);
					if (excl.find(resname) != excl.end())
						excl[resname].push_back(ex);
					else
					{
						vector<AtomInfo> grp;
						grp.push_back(ex);
						excl[resname] = grp;
					}
				}
			}
		}
		pdb.close();
		line.clear();
		resn = proid + 1;
		proatomn = proatomid;
		ligandatomn = ligandatomid;
	}
	else
		handle_error(format("Unbale to open file %1%.") % fpath);
}

void Pro::gen_contact()
{
	contact_map = MatrixXi::Zero(resn, resn);
	kmat = ArrayXXd::Zero(resn, resn);
	dist_map = MatrixXd::Zero(resn, resn);
	for (size_t i = 0; i < resn; ++i)
	{
		contact_map(i, i) = 1;
		for (size_t j = i + 1; j < resn; ++j)
		{
			double dist_ij = gen_distmat_flag ? distmat(i, j) : distance(i, j);
			if (pro[i].chain == pro[j].chain && dist_ij < cutoff_intra)
			{
				contact_map(i, j) = contact_map(j, i) = 2;
				contact_pairs.push_back(make_pair(i, j));
				contact_pairs.push_back(make_pair(j, i));
				kmat(i, j) = kmat(j, i) = k_intra;
				dist_map(i, j) = distmat(i, j) = dist_ij;
			}
			else if (pro[i].chain != pro[j].chain && dist_ij < cutoff_inter)
			{
				contact_map(i, j) = contact_map(j, i) = 3;
				contact_pairs.push_back(make_pair(i, j));
				contact_pairs.push_back(make_pair(j, i));
				kmat(i, j) = kmat(j, i) = k_inter;
				dist_map(i, j) = distmat(i, j) = dist_ij;
			}
			else
				contact_map(i, j) = contact_map(j, i) = 0;
		}
	}
	gen_contact_flag = true;
}

void Pro::gen_coord()
{
	procoord = VectorXd::Zero(3 * resn);
	for (size_t i = 0; i < resn; i++)
	{
		procoord(3 * i) = pro[i].x;
		procoord(3 * i + 1) = pro[i].y;
		procoord(3 * i + 2) = pro[i].z;
	}
	ligandcoord = VectorXd::Zero(3 * ligandatomn);
	size_t j = 0;
	for (map<string, vector<AtomInfo>>::iterator it = ligand.begin(); it != ligand.end(); ++it)
	{
		for (vector<AtomInfo>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit)
		{
			ligandcoord(3 * j) = iit->x;
			ligandcoord(3 * j + 1) = iit->y;
			ligandcoord(3 * j + 2) = iit->z;
			++j;
		}
	}
	for (map<size_t, vector<AtomInfo>>::iterator it = proatoms.begin(); it != proatoms.end(); ++it)
	{
		VectorXd rescoord = VectorXd::Zero(3 * it->second.size());
		size_t k = 0;
		for (vector<AtomInfo>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit)
		{
			rescoord(3 * k) = iit->x;
			rescoord(3 * k + 1) = iit->y;
			rescoord(3 * k + 2) = iit->z;
			++k;
		}
		rescoords[it->first] = rescoord;
	}
}

MatrixXd Pro::gen_hessian()
{
	MatrixXd hessian = MatrixXd::Zero(3 * resn, 3 * resn);
	if (gen_contact_flag)
	{
		for (vector<pair<size_t, size_t>>::iterator it = contact_pairs.begin(); it != contact_pairs.end(); ++it)
		{
			size_t pi = it->first;
			size_t pj = it->second;
			double diffx = pro[pi].x - pro[pj].x;
			double diffy = pro[pi].y - pro[pj].y;
			double diffz = pro[pi].z - pro[pj].z;
			double d = pow(diffx, 2) + pow(diffy, 2) + pow(diffz, 2);
			double k = k_default;
			if (get_contact(pi, pj) == 2)
				k = k_intra;
			else if (get_contact(pi, pj) == 3)
				k = k_inter;
			k /= d;
			double hxx = -k * pow(diffx, 2);
			double hyy = -k * pow(diffy, 2);
			double hzz = -k * pow(diffz, 2);
			double hxy = -k * diffx * diffy;
			double hxz = -k * diffx * diffz;
			double hyz = -k * diffy * diffz;
			hessian(3 * pi, 3 * pj) = hxx;
			hessian(3 * pi, 3 * pi) -= hxx;
			hessian(3 * pi + 1, 3 * pj + 1) = hyy;
			hessian(3 * pi + 1, 3 * pi + 1) -= hyy;
			hessian(3 * pi + 2, 3 * pj + 2) = hzz;
			hessian(3 * pi + 2, 3 * pi + 2) -= hzz;
			hessian(3 * pi, 3 * pj + 1) = hxy;
			hessian(3 * pi + 1, 3 * pj) = hxy;
			hessian(3 * pi, 3 * pi + 1) -= hxy;
			hessian(3 * pi + 1, 3 * pi) -= hxy;
			hessian(3 * pi, 3 * pj + 2) = hxz;
			hessian(3 * pi + 2, 3 * pj) = hxz;
			hessian(3 * pi, 3 * pi + 2) -= hxz;
			hessian(3 * pi + 2, 3 * pi) -= hxz;
			hessian(3 * pi + 1, 3 * pj + 2) = hyz;
			hessian(3 * pi + 2, 3 * pj + 1) = hyz;
			hessian(3 * pi + 1, 3 * pi + 2) -= hyz;
			hessian(3 * pi + 2, 3 * pi + 1) -= hyz;
		}
	}
	return hessian;
}

MatrixXd Pro::gen_covariance(MatrixXd hessian)
{
	MatrixXd covariance = MatrixXd::Zero(3 * resn, 3 * resn);
	SelfAdjointEigenSolver<MatrixXd> eigensolver(hessian);
	VectorXd eigenvalues = eigensolver.eigenvalues();
	VectorXd zero2inf_eigenvalues = eigenvalues;
	MatrixXd eigenvectors = eigensolver.eigenvectors();
	vector<size_t> zeromodes, nonzeromodes;
	size_t zeromoden = calc_zero_modes(eigenvalues, zero2inf_eigenvalues);
	IOFormat Clean = IOFormat(4, 0, ", ", "\n");
	if (zeromoden == 6)
	{
		MatrixXd U = ArrayXXd(eigenvectors.transpose()).colwise() / ArrayXd(zero2inf_eigenvalues);
		covariance = (kB * Navo * Temp / k_default) * (eigenvectors * U);
	}
	else
		handle_error(format("Hessian matrix has %1% zero modes. Please check it before constructing covariance matrix.") % zeromoden);
	return covariance;
}

void Pro::gen_distmat()
{
	distmat = MatrixXd::Zero(resn, resn);
	for (size_t i = 0; i < resn; i++)
		for (size_t j = i + 1; j < resn; j++)
			distmat(j, i) = distmat(i, j) = distance(i, j);
	gen_distmat_flag = true;
}

void Pro::gen_dist2ligand()
{
	dist2ligand = VectorXd::Zero(resn);

	if (with_ligand_flag)
	{
		Map<Matrix3Xd> ligandxyz(ligandcoord.data(), 3, ligandatomn);
		for (size_t i = 0; i < resn; ++i)
		{
			size_t resatomn = get_resatomn(i);
			VectorXd rescoord = get_rescoord(i);
			Map<Matrix3Xd> resxyz(rescoord.data(), 3, resatomn);
			VectorXd dist = VectorXd::Zero(resatomn);
			Array3Xd diffxyz(3, ligandxyz.cols());
			for (size_t j = 0; j < resatomn; ++j)
			{
				diffxyz = ligandxyz.colwise() - resxyz.col(j);
				dist(j) = sqrt(diffxyz.pow(2).colwise().sum().minCoeff());
			}
			dist2ligand(i) = dist.minCoeff();
		}
	}
}

double Pro::distance(size_t i, size_t j)
{
	return sqrt(pow(pro[i].x - pro[j].x, 2) + pow(pro[i].y - pro[j].y, 2) + pow(pro[i].z - pro[j].z, 2));
}

size_t Pro::calc_zero_modes(VectorXd eigenvalues, VectorXd &zero2inf_eigenvalues)
{
	size_t count = 0;
	zero2inf_eigenvalues = eigenvalues;
	for (size_t i = 0; i < size_t(eigenvalues.size()); i++)
		if (eigenvalues(i) == 0.0 || abs(eigenvalues(i)) < 1e-8)
		{
			zero2inf_eigenvalues(i) = numeric_limits<double>::infinity();
			++count;
		}
	return count;
}

bool Pro::has_res(size_t id)
{
	if (id < resn)
		return true;
	else
		return false;
}

string Pro::get_resname(size_t id)
{
	if (id < resn)
		return pro[id].resname;
	else
		return string();
}

string Pro::get_chain(size_t id)
{
	if (id < resn)
		return pro[id].chain;
	else
		return string();
}

size_t Pro::get_resid(size_t id)
{
	if (id < resn)
		return pro[id].resid;
	else
		return 0;
}

double Pro::get_x(size_t id)
{
	if (id < resn)
		return pro[id].x;
	else
		return 0.0;
}

double Pro::get_y(size_t id)
{
	if (id < resn)
		return pro[id].y;
	else
		return 0.0;
}

double Pro::get_z(size_t id)
{
	if (id < resn)
		return pro[id].z;
	else
		return 0.0;
}

double Pro::get_bfactor(size_t id)
{
	if (id < resn)
		return pro[id].bfactor;
	else
		return 0.0;
}

size_t Pro::get_resn()
{
	return resn;
}

int Pro::get_contact(size_t i, size_t j)
{
	if (i < resn && j < resn)
		return contact_map(i, j);
	else
		return 0;
}

void Pro::show_contact_pairs()
{
	vector<string> buf;
	for (vector<pair<size_t, size_t>>::iterator it = contact_pairs.begin(); it != contact_pairs.end(); ++it)
		buf.push_back(
			(format("(%1%, %2%)") % it->first % it->second).str()
		);
	handle_result("Contact pairs:", buf);
}

MatrixXi Pro::get_contact_map()
{
	return contact_map;
}

VectorXd Pro::get_procoord()
{
	return procoord;
}

VectorXd Pro::get_ligandcoord()
{
	return ligandcoord;
}

VectorXd Pro::get_dist2ligand()
{
	return dist2ligand;
}

VectorXd Pro::get_rescoord(size_t id)
{
	if (id < resn)
		return rescoords[id];
	else
		return VectorXd();
}

size_t Pro::get_resatomn(size_t id)
{
	if (id < resn)
		return rescoords[id].size() / 3;
	else
		return 0;
}

MatrixXd Pro::get_distmat()
{
	return distmat;
}

MatrixXd Pro::get_dist_map()
{
	return dist_map;
}

ArrayXXd Pro::get_kmat()
{
	return kmat;
}

bool Pro::empty()
{
	if (resn == 0)
		return true;
	else
		return false;
}
