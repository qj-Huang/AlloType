#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <limits>
#include <boost/format.hpp>
#include <Eigen/Dense>

#include "handle_io.h"
#include "read_pdb.h"

using namespace std;
using namespace Eigen;
using boost::format;

constexpr double Navo = 6.02214076e23;
constexpr double kB = 1.380649e-23;
constexpr double Temp = 298.15;
constexpr double h = 6.6260693e-34;

class Pro
{
public:
	Pro();
	Pro(string fpath, bool has_ligand_flag, vector<string> exclude = { }, double k = 1.0, double cutoff = 9.0);
	~Pro();

	bool has_ligand();

	bool has_res(size_t id);
	string get_resname(size_t id);
	string get_chain(size_t id);
	size_t get_resid(size_t id);
	double get_x(size_t id);
	double get_y(size_t id);
	double get_z(size_t id);
	double get_bfactor(size_t id);

	VectorXd get_rescoord(size_t id);
	size_t get_resatomn(size_t id);

	size_t get_resn();

	int get_contact(size_t i, size_t j);

	void show_contact_pairs();
	MatrixXi get_contact_map();

	VectorXd get_procoord();
	VectorXd get_ligandcoord();

	VectorXd get_dist2ligand();

	MatrixXd get_distmat();
	ArrayXXd get_kmat();
	MatrixXd get_dist_map();

	MatrixXd gen_hessian();

	MatrixXd gen_covariance(MatrixXd hessian);

	bool empty();

private:
	double distance(size_t i, size_t j);

	size_t calc_zero_modes(VectorXd eigenvalues, VectorXd &zero2inf_eigenvalues);

	void read(string fpath);

	void gen_contact();

	void gen_coord();

	void gen_distmat();

	void gen_dist2ligand();

	map<size_t, ResInfo> pro;
	map<size_t, vector<AtomInfo>> proatoms;
	map<string, vector<AtomInfo>> ligand;
	map<string, vector<AtomInfo>> excl;

	MatrixXi contact_map;
	vector<pair<size_t, size_t>> contact_pairs;
	ArrayXXd kmat;
	MatrixXd distmat;
	VectorXd procoord;
	map<size_t, VectorXd> rescoords;
	VectorXd ligandcoord;
	VectorXd dist2ligand;
	MatrixXd dist_map;

	set<string> prores = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };
	set<string> ligandres = {};
	set<string> exclres = { "HOH" };

	size_t resn = 0;
	size_t proatomn = 0;
	size_t ligandatomn = 0;
	double cutoff_intra = 9.0;
	double cutoff_inter = 9.0;
	double k_intra = 1.0;
	double k_inter = 1.0;
	double k_default = 1.0;

	bool gen_distmat_flag = false;
	bool gen_contact_flag = false;
	bool with_ligand_flag = false;
};
