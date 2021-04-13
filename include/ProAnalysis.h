#pragma once
#include <algorithm>
#include <fstream>
#include <list>
#include <utility>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "handle_io.h"
#include "Pro.h"

using namespace std;
using namespace Eigen;
using boost::format;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using boost::algorithm::trim;

constexpr double PI = 3.1415926535897932;

class ProAnalysis
{
public:
	ProAnalysis();
	ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex);
	ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex, string hessian_path, string covariance_path);
	~ProAnalysis();

	void interactive_pocket(unsigned int mode);
	void interactive();

	MatrixXd get_hessian()
	{
		return hessian;
	}

	Matrix3d get_hessian(size_t i, size_t j);
	double get_hessian_s(size_t si, size_t sj);
	void write_hessian(string writepath)
	{
		write_matrix(hessian, writepath);
	}
	void write_hessian_binary(string writepath)
	{
		write_matrix_binary(hessian, writepath);
	}
	void read_hessian_binary(string fpath)
	{
		read_matrix_binary(hessian, fpath);
	}

	MatrixXd get_covariance()
	{
		return covariance;
	}
	Matrix3d get_covariance(size_t i, size_t j);
	double get_covariance_s(size_t si, size_t sj);
	void write_covariance(string writepath)
	{
		write_matrix(covariance, writepath);
	}
	void write_covariance_binary(string writepath)
	{
		write_matrix_binary(covariance, writepath);
	}
	void read_covariance_binary(string fpath)
	{
		read_matrix_binary(covariance, fpath);
	}

	list<size_t> get_pocketS() {
		return pocketS;
	}
	list<size_t> get_pocketA() {
		return pocketA;
	}
	list<size_t> get_pocketAS() {
		return pocketAS;
	}

	void show_pocketS() {
		show_pocket(pocketS);
	}
	void show_pocketA() {
		show_pocket(pocketA);
	}
	void show_pocketAS() {
		show_pocket(pocketAS);
	}

	void show_pocketS_force() {
		show_pocket_force(has_pocketS_force_flag, pocketS, pocketS_force);
	}
	void show_pocketA_force() {
		show_pocket_force(has_pocketA_force_flag, pocketA, pocketA_force);
	}
	void show_pocketAS_force() {
		show_pocket_force(has_pocketAS_force_flag, pocketAS, pocketAS_force);
	}

	void test_pocketS() {
		test_pocket(has_pocketS_force_flag, ES_info, pocketS, ES_displacement, pocketS_force, S_fitprocoord);
		calc_model_correlation(has_pocketS_force_flag, ES_displacement, pocketS_force);
	}
	void test_pocketA() {
		test_pocket(has_pocketA_force_flag, EA_info, pocketA, EA_displacement, pocketA_force, A_fitprocoord);
		calc_model_correlation(has_pocketA_force_flag, EA_displacement, pocketA_force);
	}
	void test_pocketAS() {
		test_pocket(has_pocketAS_force_flag, EAS_info, pocketAS, EAS_displacement, pocketAS_force, AS_fitprocoord);
		calc_model_correlation(has_pocketAS_force_flag, EAS_displacement, pocketAS_force);
	}

	bool in_pocketS(size_t id) {
		return in_pocket(pocketS, id);
	}
	bool in_pocketA(size_t id) {
		return in_pocket(pocketA, id);
	}
	bool in_pocketAS(size_t id) {
		return in_pocket(pocketAS, id);
	}

	void add_to_pocketS(size_t id) {
		add_to_pocket(pocketS, id);
	}
	void add_to_pocketA(size_t id) {
		add_to_pocket(pocketA, id);
	}
	void add_to_pocketAS(size_t id) {
		add_to_pocket(pocketAS, id);
	}

	void remove_from_pocketS(size_t id) {
		remove_from_pocket(pocketS, id);
	}
	void remove_from_pocketA(size_t id) {
		remove_from_pocket(pocketA, id);
	}
	void remove_from_pocketAS(size_t id) {
		remove_from_pocket(pocketAS, id);
	}

	void gen_pocketS(double cutoff) {
		if (!ProS.empty())
			gen_pocket(ProS.has_ligand(), pocketS, cutoff, S_dist2ligand);
	}
	void gen_pocketA(double cutoff) {
		if (!ProA.empty())
			gen_pocket(ProA.has_ligand(), pocketA, cutoff, A_dist2ligand);
	}
	void gen_pocketAS(double cutoff) {
		if (!ProAS.empty())
			gen_pocket(ProAS.has_ligand(), pocketAS, cutoff, AS_dist2ligand);
	}

	double calc_ES_rmsd() {
		if (ES_info)
			calc_model_rmsd(has_pocketS_force_flag, pocketS_force, S_fitprocoord);
	}
	double calc_EA_rmsd() {
		if (EA_info)
			calc_model_rmsd(has_pocketA_force_flag, pocketA_force, A_fitprocoord);
	}
	double calc_EAS_rmsd() {
		if (EAS_info)
			calc_model_rmsd(has_pocketAS_force_flag, pocketAS_force, AS_fitprocoord);
	}

	void gen_pocketS_force() {
		if (ES_info)
			gen_pocket_force(has_pocketS_force_flag, pocketS_force, pocketS, ES_displacement);
	}
	void gen_pocketA_force() {
		if (EA_info)
			gen_pocket_force(has_pocketA_force_flag, pocketA_force, pocketA, EA_displacement);
	}
	void gen_pocketAS_force() {
		if (EAS_info)
		{
			if (has_pocketS_force_flag)
				gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketS_force, pocketAS, EAS_displacement);
			else if (has_pocketA_force_flag)
				gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketA_force, pocketAS, EAS_displacement);
			else
				gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketAS, EAS_displacement);
		}
	}

	void gen_free_energy();

private:

	void write_matrix(MatrixXd mat, string writepath);

	void write_matrix_binary(MatrixXd mat, string writepath);

	void read_matrix_binary(MatrixXd & mat, string fpath);

	double calc_model_rmsd(bool flag, VectorXd pocket_force, VectorXd refcoord);

	double calc_model_correlation(bool flag, VectorXd displacement, VectorXd pocket_force);

	void show_pocket(list<size_t> pocket);

	void test_pocket(bool flag, bool info, list<size_t> pocket, VectorXd displacement, VectorXd pocket_force, VectorXd refcoord);

	void show_pocket_force(bool flag, list<size_t> pocket, VectorXd pocket_force);

	bool in_pocket(list<size_t> pocket, size_t id);

	void add_to_pocket(list<size_t> & pocket, size_t id);

	void remove_from_pocket(list<size_t> & pocket, size_t id);

	void gen_pocket(bool has_ligand, list<size_t> & pocket, double cutoff, VectorXd dist2ligand);

	void gen_pocket_force(bool & flag, VectorXd & pocket_force, list<size_t> pocket, VectorXd displacement);

	void gen_pocket_force(bool & flag, VectorXd & pocket_force, VectorXd fixed_force, list<size_t> pocket, VectorXd displacement);

	void calc_energy_unknown(bool flag, double &proenergy, double &pocketenergy, double &totenergy, VectorXd pocket_force, ArrayXd &per_residue_ene);

	void print_energy_results();

	Pro ProE;
	Pro ProS;
	Pro ProA;
	Pro ProAS;

	double pocS_ave;
	double pocA_ave;
	double pocAS_ave;

	ArrayXd ps_per_residue_ene;
	ArrayXd pa_per_residue_ene;
	ArrayXd pas_per_residue_ene;

	VectorXd S_fitprocoord;
	VectorXd S_dist2ligand;
	VectorXd pocketS_force;
	double S_proenergy = 0.0;
	double S_pocketenergy = 0.0;
	double S_energy = 0.0;

	VectorXd A_fitprocoord;
	VectorXd A_dist2ligand;
	VectorXd pocketA_force;
	double A_proenergy = 0.0;
	double A_pocketenergy = 0.0;
	double A_energy = 0.0;

	VectorXd AS_fitprocoord;
	VectorXd AS_dist2ligand;
	VectorXd pocketAS_force;
	double AS_proenergy = 0.0;
	double AS_pocketenergy = 0.0;
	double AS_energy = 0.0;

	MatrixXd hessian;
	MatrixXd covariance;
	MatrixXd mprocoord;
	
	double ddG = 0.0;

	bool ES_info = false;
	VectorXd ES_displacement;
	double ES_rmsd;

	bool EA_info = false;
	VectorXd EA_displacement;
	double EA_rmsd;

	bool EAS_info = false;
	VectorXd EAS_displacement;
	double EAS_rmsd;

	list<size_t> pocketS;
	list<size_t> pocketA;
	list<size_t> pocketAS;

	// Status
	bool has_pocketS_force_flag = false;
	bool has_pocketA_force_flag = false;
	bool has_pocketAS_force_flag = false;

	// Matrix formats
	IOFormat CleanFmt = IOFormat(4, 0, ", ", "\n", "[", "]");

	double cutoff_intra = 9.0;
	double k_intra = 1.0;
};
