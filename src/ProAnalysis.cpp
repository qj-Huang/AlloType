#include "ProAnalysis.h"
#include "Pro.h"
#include "method.h"

ProAnalysis::ProAnalysis()
{
}

ProAnalysis::ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex)
{
	ProE = apo;
	ProS = binding;
	ProA = allostery;
	ProAS = complex;
	if (!ProE.empty())
	{
		hessian = ProE.gen_hessian();
		handle_info("Finish constructing Hessian matrix.");
		covariance = ProE.gen_covariance(hessian);
		handle_info("Finish constructing Covariance matrix.");
		if (!ProS.empty())
		{
			if (ProE.get_resn() != ProS.get_resn())
				handle_error("Apo state protein and binding state protein do not have equal residue numbers.");
			S_dist2ligand = ProS.get_dist2ligand();
			S_fitprocoord = fitting(ProE.get_procoord(), ProS.get_procoord());
			handle_info("Fitting process succeed.");
			ES_displacement = calc_displacement(ProE.get_procoord(), S_fitprocoord);
			handle_info("Calculating displacement succeed.");
			ES_rmsd = calc_rmsd(ES_displacement);
			handle_result(format("RMSD from binding state PDB file: %1$.4f A.") % ES_rmsd);
			ES_info = true;
		}
		if (!ProA.empty())
		{
			if (ProE.get_resn() != ProA.get_resn())
				handle_error("Apo state protein and allostery state protein do not have equal residue numbers.");
			A_dist2ligand = ProA.get_dist2ligand();
			A_fitprocoord = fitting(ProE.get_procoord(), ProA.get_procoord());
			handle_info("Fitting process succeed.");
			EA_displacement = calc_displacement(ProE.get_procoord(), A_fitprocoord);
			handle_info("Calculating displacement succeed.");
			EA_rmsd = calc_rmsd(EA_displacement);
			handle_result(format("RMSD from allostery state PDB file: %1$.4f A.") % EA_rmsd);
			EA_info = true;
		}
		if (!ProAS.empty())
		{
			if (ProE.get_resn() != ProAS.get_resn())
				handle_error("Apo state protein and complex state protein do not have equal residue numbers.");
			AS_dist2ligand = ProAS.get_dist2ligand(); 
			AS_fitprocoord = fitting(ProE.get_procoord(), ProAS.get_procoord());
			handle_info("Fitting process succeed.");
			EAS_displacement = calc_displacement(ProE.get_procoord(), AS_fitprocoord);
			handle_info("Calculating displacement succeed.");
			EAS_rmsd = calc_rmsd(EAS_displacement);
			handle_result(format("RMSD from complex state PDB file: %1$.4f A.") % EAS_rmsd);
			EAS_info = true;
		}
	}
}

ProAnalysis::ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex, string hessian_path, string covariance_path)
{
	ProE = apo;
	ProS = binding;
	ProA = allostery;
	ProAS = complex;
	if (!ProE.empty())
	{
		read_hessian_binary(hessian_path);
		handle_info("Finish reading Hessian matrix.");
		read_covariance_binary(covariance_path);
		handle_info("Finish reading Covariance matrix.");
		if (!ProS.empty())
		{
			if (ProE.get_resn() != ProS.get_resn())
				handle_error("Apo state protein and binding state protein do not have equal residue numbers.");
			S_dist2ligand = ProS.get_dist2ligand();
			S_fitprocoord = fitting(ProE.get_procoord(), ProS.get_procoord());
			handle_info("Fitting process succeed.");
			ES_displacement = calc_displacement(ProE.get_procoord(), S_fitprocoord);
			handle_info("Calculating displacement succeed.");
			ES_rmsd = calc_rmsd(ES_displacement); 
			handle_result(format("RMSD from binding state PDB file: %1$.4f A.") % ES_rmsd);
			ES_info = true;
		}
		if (!ProA.empty())
		{
			if (ProE.get_resn() != ProA.get_resn())
				handle_error("Apo state protein and allostery state protein do not have equal residue numbers.");
			A_dist2ligand = ProA.get_dist2ligand();
			A_fitprocoord = fitting(ProE.get_procoord(), ProA.get_procoord());
			handle_info("Fitting process succeed.");
			EA_displacement = calc_displacement(ProE.get_procoord(), A_fitprocoord);
			handle_info("Calculating displacement succeed.");
			EA_rmsd = calc_rmsd(EA_displacement); 
			handle_result(format("RMSD from allostery state PDB file: %1$.4f A.") % EA_rmsd);
			EA_info = true;
		}
		if (!ProAS.empty())
		{
			if (ProE.get_resn() != ProAS.get_resn())
				handle_error("Apo state protein and complex state protein do not have equal residue numbers.");
			AS_dist2ligand = ProAS.get_dist2ligand();
			AS_fitprocoord = fitting(ProE.get_procoord(), ProAS.get_procoord());
			handle_info("Fitting process succeed.");
			EAS_displacement = calc_displacement(ProE.get_procoord(), AS_fitprocoord);
			handle_info("Calculating displacement succeed.");
			EAS_rmsd = calc_rmsd(EAS_displacement); 
			handle_result(format("RMSD from complex state PDB file: %1$.4f A.") % EAS_rmsd);
			EAS_info = true;
		}
	}
}

ProAnalysis::~ProAnalysis()
{
}

void ProAnalysis::interactive_pocket(unsigned int mode)
{
	string buf, cmd, label;
	vector<string> para;
	switch (mode)
	{
	case 0:
		label = "S";
		break;
	case 1:
		label = "A";
		break;
	case 2:
		label = "AS";
		break;
	}

	while (true)
	{
		cout << '[' << mode << ']' << " >>> ";
		getline(cin, buf);
		trim(buf);
		split(para, buf, boost::is_any_of(" "));
		cmd = para[0];
		if (cmd == "back")
			break;
		else if (cmd == "add")
		{
			size_t resid = 0;
			bool convert_flag = true;
			try
			{
				resid = lexical_cast<size_t>(para[1]) - 1;
			}
			catch (bad_lexical_cast)
			{
				convert_flag = false;
				handle_hint("No residue ID given. Please enter again.");
			}

			if (convert_flag)
			{
				switch (mode)
				{
				case 0:
					add_to_pocketS(resid);
					break;
				case 1:
					add_to_pocketA(resid);
					break;
				case 2:
					add_to_pocketAS(resid);
					break;
				}
			}
		}
		else if (cmd == "del")
		{
			size_t resid = 0;
			bool convert_flag = true;
			try
			{
				resid = lexical_cast<size_t>(para[1]) - 1;
			}
			catch (bad_lexical_cast)
			{
				convert_flag = false;
				handle_hint("No residue ID given. Please enter again.");
			}

			if (convert_flag)
			{
				switch (mode)
				{
				case 0:
					remove_from_pocketS(resid);
					break;
				case 1:
					remove_from_pocketA(resid);
					break;
				case 2:
					remove_from_pocketAS(resid);
					break;
				}
			}
		}
		else if (cmd == "gen-pocket")
		{
			double cutoff = 0.0;
			bool convert_flag = true;
			try
			{
				cutoff = lexical_cast<double>(para[1]);
			}
			catch (bad_lexical_cast)
			{
				convert_flag = false;
				handle_hint("Illegal cutoff length given. Please enter again.");
			}

			if (convert_flag)
			{
				switch (mode)
				{
				case 0:
					gen_pocketS(cutoff);
					break;
				case 1:
					gen_pocketA(cutoff);
					break;
				case 2:
					gen_pocketAS(cutoff);
					break;
				}
			}
		}
		else if (cmd == "show")
		{
			switch (mode)
			{
			case 0:
				show_pocketS();
				break;
			case 1:
				show_pocketA();
				break;
			case 2:
				show_pocketAS();
				break;
			}
		}
		else if (cmd == "test")
		{
			switch (mode)
			{
			case 0:
				test_pocketS();
				break;
			case 1:
				test_pocketA();
				break;
			case 2:
				test_pocketAS();
				break;
			}
		}
		else if (cmd == "gen-force")
		{
			switch (mode)
			{
			case 0:
				gen_pocketS_force();
				break;
			case 1:
				gen_pocketA_force();
				break;
			case 2:
				gen_pocketAS_force();
				break;
			}
		}
		else if (cmd == "show-force")
		{
			switch (mode)
			{
			case 0:
				show_pocketS_force();
				break;
			case 1:
				show_pocketA_force();
				break;
			case 2:
				show_pocketAS_force();
				break;
			}
		}
		else
			handle_hint("Unknown command.");
		cmd.clear();
	}
}

void ProAnalysis::interactive()
{
	string buf, cmd;
	vector<string> para;
	while (true)
	{
		cout << ">>> ";
		getline(cin, buf);
		trim(buf);
		split(para, buf, boost::is_any_of(" "));
		cmd = para[0];
		if (cmd == "exit")
			break;
		else if (cmd == "pocketS")
			interactive_pocket(0);
		else if (cmd == "pocketA")
			interactive_pocket(1);
		else if (cmd == "pocketAS")
			interactive_pocket(2);
		else if (cmd == "energy")
			gen_free_energy();
		else
			handle_hint("Unknown command.");
		cmd.clear();
	}
}

Matrix3d ProAnalysis::get_hessian(size_t i, size_t j)
{
	if (i < size_t(hessian.rows() / 3) && j < size_t(hessian.cols() / 3))
		return Matrix3d(hessian.block(3 * i, 3 * j, 3, 3));
	else
		return Matrix3d();
}

double ProAnalysis::get_hessian_s(size_t si, size_t sj)
{
	if (si < size_t(hessian.rows()) && sj < size_t(hessian.cols()))
		return hessian(si, sj);
	else
		return 0.0;
}

Matrix3d ProAnalysis::get_covariance(size_t i, size_t j)
{
	if (i < size_t(covariance.rows() / 3) && j < size_t(covariance.cols() / 3))
		return Matrix3d(covariance.block(3 * i, 3 * j, 3, 3));
	else
		return Matrix3d();
}

double ProAnalysis::get_covariance_s(size_t si, size_t sj)
{
	if (si < size_t(covariance.rows()) && sj < size_t(covariance.cols()))
		return covariance(si, sj);
	else
		return 0.0;
}

void ProAnalysis::gen_free_energy()
{
	if (ES_info && EA_info)
	{
		calc_energy_unknown(has_pocketS_force_flag, S_proenergy, S_pocketenergy, S_energy, pocketS_force, ps_per_residue_ene);
		calc_energy_unknown(has_pocketA_force_flag, A_proenergy, A_pocketenergy, A_energy, pocketA_force, pa_per_residue_ene);
		pocketAS_force = pocketS_force + pocketA_force;
		write_matrix(pocketAS_force, "force.txt");
		has_pocketAS_force_flag = true;
		calc_energy_unknown(has_pocketAS_force_flag, AS_proenergy, AS_pocketenergy, AS_energy, pocketAS_force, pas_per_residue_ene);
		for (list<size_t>::iterator it = pocketS.begin(); it != pocketS.end(); ++it)
		{
			add_to_pocketAS(*it);
		}
		for (list<size_t>::iterator it = pocketA.begin(); it != pocketA.end(); ++it)
		{
			add_to_pocketAS(*it);
		}
		ddG = AS_energy - S_energy - A_energy;
		VectorXd ddG_per_residue = pas_per_residue_ene - pa_per_residue_ene - ps_per_residue_ene;
		write_matrix(ddG_per_residue, "ddG_per_residue.txt");
		pocS_ave = S_pocketenergy / pocketS.size() / kB / Temp / Navo * 4186;
		pocA_ave = A_pocketenergy / pocketA.size() / kB / Temp / Navo * 4186;
		pocAS_ave = AS_pocketenergy / pocketAS.size() / kB / Temp / Navo * 4186;
		print_energy_results();
	}
	else if (EAS_info && ES_info)
	{
		calc_energy_unknown(has_pocketS_force_flag, S_proenergy, S_pocketenergy, S_energy, pocketS_force, ps_per_residue_ene);
		calc_energy_unknown(has_pocketAS_force_flag, AS_proenergy, AS_pocketenergy, AS_energy, pocketAS_force, pas_per_residue_ene);
		write_matrix(pocketAS_force, "force.txt");
		pocketA_force = VectorXd::Zero(pocketAS_force.size());
		for (list<size_t>::iterator it = pocketAS.begin(); it != pocketAS.end(); ++it)
		{
			add_to_pocketA(*it);
		}
		pocketA_force = pocketAS_force - pocketS_force;
		has_pocketA_force_flag = true;
		calc_energy_unknown(has_pocketA_force_flag, A_proenergy, A_pocketenergy, A_energy, pocketA_force, pa_per_residue_ene);
		ddG = AS_energy - S_energy - A_energy;
		ArrayXd ddG_per_residue = pas_per_residue_ene - pa_per_residue_ene - ps_per_residue_ene;
		write_matrix(ddG_per_residue, "ddG_per_residue.txt");
		pocS_ave = S_pocketenergy / pocketS.size() / kB / Temp / Navo * 4186;
		pocA_ave = A_pocketenergy / pocketA.size() / kB / Temp / Navo * 4186;
		pocAS_ave = AS_pocketenergy / pocketAS.size() / kB / Temp / Navo * 4186;
		print_energy_results();
	}
	else if (EAS_info && EA_info)
	{
		calc_energy_unknown(has_pocketA_force_flag, A_proenergy, A_pocketenergy, A_energy, pocketA_force, pa_per_residue_ene);
		calc_energy_unknown(has_pocketAS_force_flag, AS_proenergy, AS_pocketenergy, AS_energy, pocketAS_force, pas_per_residue_ene);
		write_matrix(pocketAS_force, "force.txt");
		pocketS_force = VectorXd::Zero(pocketAS_force.size());
		for (list<size_t>::iterator it = pocketAS.begin(); it != pocketAS.end(); ++it)
		{
			add_to_pocketS(*it);
		}
		pocketS_force = pocketAS_force - pocketA_force;
		has_pocketS_force_flag = true;
		calc_energy_unknown(has_pocketS_force_flag, S_proenergy, S_pocketenergy, S_energy, pocketS_force, ps_per_residue_ene);
		ddG = AS_energy - S_energy - A_energy;
		VectorXd ddG_per_residue = pas_per_residue_ene - pa_per_residue_ene - ps_per_residue_ene;
		write_matrix(ddG_per_residue, "ddG_per_residue.txt");
		pocS_ave = S_pocketenergy / pocketS.size() / kB / Temp / Navo * 4186;
		pocA_ave = A_pocketenergy / pocketA.size() / kB / Temp / Navo * 4186;
		pocAS_ave = AS_pocketenergy / pocketAS.size() / kB / Temp / Navo * 4186;
		print_energy_results();
	}
	else
		handle_error("Lack necessary protein information.");
}

void ProAnalysis::write_matrix(MatrixXd mat, string writepath)
{
	ofstream matf(writepath, ios::out);
	if (matf.is_open())
	{
		matf << mat.format(CleanFmt);
		matf.close();
		handle_info(format("Matrix has been written to %1%.") % writepath);
	}
	else
		handle_warning(format("Can not open file %1%.") % writepath);
}

void ProAnalysis::write_matrix_binary(MatrixXd mat, string writepath)
{
	ofstream matf(writepath, ios::out | ios::binary);
	if (matf.is_open())
	{
		size_t Msize = mat.size();
		size_t Mrown = mat.rows();
		size_t Mcoln = mat.cols();
		double* M = new double[Msize];
		for (size_t i = 0; i < Mrown; ++i)
			for (size_t j = 0; j < Mcoln; ++j)
				M[i * Mcoln + j] = mat(i, j);
		matf.write((char *)&Msize, sizeof(size_t));
		matf.write((char *)&Mrown, sizeof(size_t));
		matf.write((char *)&Mcoln, sizeof(size_t));
		matf.write((char *)&M[0], Msize * sizeof(double));
		matf.close();
		delete[] M;
		handle_info(format("Matrix (binary data) has been written to %1%.") % writepath);
	}
	else
		handle_warning(format("Can not open file %1%.") % writepath);
}

void ProAnalysis::read_matrix_binary(MatrixXd & mat, string fpath)
{
	ifstream matf(fpath, ios::in | ios::binary);
	if (matf.is_open())
	{
		size_t Msize = 0;
		size_t Mrown = 0;
		size_t Mcoln = 0;
		matf.read((char *)&Msize, sizeof(size_t));
		matf.read((char *)&Mrown, sizeof(size_t));
		matf.read((char *)&Mcoln, sizeof(size_t));
		double* M = new double[Msize];
		matf.read((char *)&M[0], Msize * sizeof(double));
		mat = MatrixXd::Zero(Mrown, Mcoln);
		for (size_t i = 0; i < Mrown; ++i)
			for (size_t j = 0; j < Mcoln; ++j)
				mat(i, j) = M[i * Mcoln + j];
		matf.close();
		delete[] M;
	}
	else
		handle_warning(format("Can not open file %1%.") % fpath);
}

double ProAnalysis::calc_model_rmsd(bool flag, VectorXd pocket_force, VectorXd refcoord)
{
	if (flag)
	{
		mprocoord = (covariance / kB / Temp / Navo) * pocket_force  + ProE.get_procoord();
		VectorXd fitmprocoord = fitting(refcoord, mprocoord);
		return calc_rmsd(refcoord, fitmprocoord);
	}
	else
	{
		handle_warning("Can not calculate model RMSD without pocket force generated. Call \"gen_pocket*_force\" function first.");
		return 0.0;
	}
}

double ProAnalysis::calc_model_correlation(bool flag, VectorXd displacement, VectorXd pocket_force)
{
	if (flag)
	{
		mprocoord = (covariance / kB / Temp / Navo) * pocket_force + ProE.get_procoord();
		VectorXd new_displacement = calc_displacement(ProE.get_procoord(), mprocoord);
		if (displacement.size() == new_displacement.size())
		{
			double xBar = calc_average(new_displacement);
			double yBar = calc_average(displacement);
			double varX = 0;
			double varY = 0;
			double SSR = 0;
			double SST = 0;
			size_t vlen = displacement.size() / 3;
			VectorXd X = VectorXd::Zero(vlen);
			VectorXd Y = VectorXd::Zero(vlen);
			VectorXd diffX = VectorXd::Zero(vlen);
			VectorXd diffY = VectorXd::Zero(vlen);
			for (size_t i = 0; i < vlen; ++i)
			{
				X(i) = sqrt(pow(new_displacement(i * 3), 2) + pow(new_displacement(i * 3 + 1), 2) + pow(new_displacement(i * 3 + 2), 2));
				Y(i) = sqrt(pow(displacement(i * 3), 2) + pow(displacement(i * 3 + 1), 2) + pow(displacement(i * 3 + 2), 2));
				diffX(i) = X(i) - xBar;
				diffY(i) = Y(i) - yBar;
				SSR += (diffX(i) * diffY(i));
				varX += pow(diffX(i), 2);
				varY += pow(diffY(i), 2);
			}				
			SST = sqrt(varX * varY);
			return SSR / SST;
		}
		else
			return 0.0;
	}
	else
	{
		handle_warning("Can not calculate model correlation without pocket force generated. Call \"gen_pocket*_force\" function first.");
		return 0.0;
	}
}

void ProAnalysis::show_pocket(list<size_t> pocket)
{
	string buf = "Pocket residues: ";
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		buf += to_string(*it + 1) + " ";
	}
	handle_result(buf);
}

void ProAnalysis::test_pocket(bool flag, bool info, list<size_t> pocket, VectorXd displacement, VectorXd pocket_force, VectorXd refcoord)
{
	if (!pocket.empty() && flag)
	{
		show_pocket(pocket);
		if (info)
		{
			handle_result(format("RMSD between real structure and structure calculated according to current pocket: %1% A.") % calc_model_rmsd(flag, pocket_force, refcoord));
			handle_result(format("Pearson between real structure and structure calculated according to current pocket: %1% A.") % calc_model_correlation(flag, displacement, pocket_force));
		}
		else
			handle_warning("Lack necessary protein information.");
	}
	else
		handle_warning("The binding pocket domain is not specificed.");
}

void ProAnalysis::show_pocket_force(bool flag, list<size_t> pocket, VectorXd pocket_force)
{
	vector<string> buf;
	if (flag)
	{
		for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
		{
			Vector3d resforce = Vector3d::Zero();
			resforce << pocket_force(*it * 3), pocket_force(*it * 3 + 1), pocket_force(*it * 3 + 2);
			buf.push_back(
				(format("RES %1% FORCE %2$.4f") % (*it + 1) % calc_norm(resforce)).str()
			);
		}
		handle_result("Pocket force:", buf);
	}
	else
		handle_warning("The binding pocket domain is not specificed.");
}

bool ProAnalysis::in_pocket(list<size_t> pocket, size_t id)
{
	bool find_element_flag = false;
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
		if (*it == id)
			find_element_flag = true;
	return find_element_flag;
}

void ProAnalysis::add_to_pocket(list<size_t> & pocket, size_t id)
{
	if (id < ProE.get_resn())
	{
		if (in_pocket(pocket, id))
			handle_info("Given residue ID already in pocket.");
		else
		{
			pocket.push_back(id);
			handle_info(format("Add residue %1% to pocket.") % (id + 1));
		}
	}
	else
		handle_warning("Given residue ID out of range.");
}

void ProAnalysis::remove_from_pocket(list<size_t> & pocket, size_t id)
{
	if (id < ProE.get_resn())
	{
		if (in_pocket(pocket, id))
		{
			for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
				if (*it == id)
				{
					pocket.erase(it);
					handle_info(format("Remove residue %1% from pocket.") % (id + 1));
				}
		}
		else
			handle_info("Given residue ID not in pocket. Can not erase it.");
	}
	else
		handle_warning("Given residue ID out of range.");
}

void ProAnalysis::gen_pocket(bool has_ligand, list<size_t> &pocket, double cutoff, VectorXd dist2ligand)
{
	if (has_ligand)
	{
		if (cutoff > dist2ligand.minCoeff())
		{
			if (!pocket.empty())
				pocket.clear();
			for (size_t i = 0; i < size_t(dist2ligand.size()); ++i)
				if (dist2ligand(i) < cutoff)
					pocket.push_back(i);
		}
		else
			handle_warning(format("Given cutoff is too short. Minimum possible cutoff is %1$.2f A.") % dist2ligand.minCoeff());
	}
	else
		handle_warning("Can not find ligand information.");
}

void ProAnalysis::gen_pocket_force(bool & flag, VectorXd & pocket_force, list<size_t> pocket, VectorXd displacement)
{
	pocket_force = VectorXd::Zero(covariance.rows());  
	size_t ndim = pocket.size() * 3;	
	MatrixXd X = MatrixXd::Zero(covariance.rows(), ndim);
	VectorXd Y = displacement;
	VectorXd coeff = VectorXd::Zero(ndim);
	size_t i = 0;
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)  
	{
		X.col(3 * i) = covariance.col(*it * 3);
		X.col(3 * i + 1) = covariance.col(*it * 3 + 1);
		X.col(3 * i + 2) = covariance.col(*it * 3 + 2);	
		++i;
	} 
	X /= (Navo * Temp * kB);
	normal_equation(coeff, X, Y);
	i = 0;
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		pocket_force(*it * 3) = coeff(i * 3);
		pocket_force(*it * 3 + 1) = coeff(i * 3 + 1);
		pocket_force(*it * 3 + 2) = coeff(i * 3 + 2);
		++i;
	}
	flag = true;
}

void ProAnalysis::gen_pocket_force(bool & flag, VectorXd & pocket_force, VectorXd fixed_force, list<size_t> pocket, VectorXd displacement)
{
	VectorXd equiv_displacement = displacement - covariance * fixed_force / kB / Temp / Navo;
	gen_pocket_force(flag, pocket_force, pocket, equiv_displacement);
	pocket_force += fixed_force;
	flag = true;
}

void ProAnalysis::calc_energy_unknown(bool flag, double &proenergy, double &pocketenergy, double &totenergy, VectorXd pocket_force, ArrayXd &per_residue_ene)
{
	if (flag)
	{
		VectorXd procoord = covariance / kB / Temp / Navo * pocket_force + ProE.get_procoord();
		ArrayXXd distdiffmat = gen_distmat(procoord) - ProE.get_distmat();
		proenergy = (distdiffmat.pow(2) * ProE.get_kmat()).sum() / 4;
		VectorXd dis = covariance / kB / Temp / Navo * pocket_force;
		pocketenergy = -pocket_force.transpose() * covariance / kB / Temp / Navo * pocket_force;
		totenergy = proenergy + pocketenergy;
		size_t vlen = pocket_force.size()/3;
		per_residue_ene = (distdiffmat.pow(2) * ProE.get_kmat()).rowwise().sum()/2;	
		for (size_t i = 0; i < vlen; ++i)
		{
			per_residue_ene(i) -= (pocket_force(i * 3) * dis(i * 3)+ pocket_force(i * 3 + 1) * dis(i * 3 + 1 )+ pocket_force(i * 3 + 2) * dis(i * 3 + 2));
		}
	}
	else
		handle_warning("Can not calculate energy without pocket force generated. Call \"gen_pocket*_force\" function first.");
}

void ProAnalysis::print_energy_results()
{
	vector<string> buf;
	buf.push_back((format("Free energy for binding state structure S: %1$.4f Kcal/mol.") % S_energy).str());
	buf.push_back((format("Pro: %1$.4f Kcal/mol.") % S_proenergy).str());
	buf.push_back((format("Pocket: %1$.4f Kcal/mol.") % S_pocketenergy).str());
	buf.push_back((format("Free energy for allostery state structure A: %1$.4f Kcal/mol.") % A_energy).str());
	buf.push_back((format("Pro: %1$.4f Kcal/mol.") % A_proenergy).str());
	buf.push_back((format("Pocket: %1$.4f Kcal/mol.") % A_pocketenergy).str());
	buf.push_back((format("Free energy for complex state structure AS: %1$.4f Kcal/mol.") % AS_energy).str());
	buf.push_back((format("Pro: %1$.4f Kcal/mol.") % AS_proenergy).str());
	buf.push_back((format("Pocket: %1$.4f Kcal/mol.") % AS_pocketenergy).str());
	buf.push_back((format("Change of free energy: : %1$.4f Kcal/mol.") % ddG).str());
	buf.push_back((format("S_ave: : %1$.4f Kcal/mol.") % pocS_ave).str());
	buf.push_back((format("A_ave: : %1$.4f Kcal/mol.") % pocA_ave).str());
	buf.push_back((format("AS_ave: : %1$.4f Kcal/mol.") % pocAS_ave).str());
	handle_result("All predict free energy results: ", buf);
}
