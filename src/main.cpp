#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "handle_io.h"
#include "Pro.h"
#include "ProAnalysis.h"

using namespace std;
using namespace boost::algorithm;
using boost::lexical_cast;
using boost::filesystem::path;

vector<string> get_exclude_res()
{
	string buf = "";
	cout << "Enter residues to be excluded: ";
	getline(cin, buf);
	trim(buf);
	vector<string> excl = {};
	split(excl, buf, boost::is_any_of(" "));
	return excl;
}

double get_spring_constant()
{
	string buf;
	cout << "Enter spring constant: ";
	getline(cin, buf);
	trim(buf);
	return lexical_cast<double>(buf);
}

double get_cutoff()
{
	string buf;
	cout << "Enter cutoff: ";
	getline(cin, buf);
	trim(buf);
	return lexical_cast<double>(buf);
}

int main()
{
	string dataset = "", pro = "";
	cout << "Enter dataset path: ";
	getline(cin, dataset);
	cout << "Enter protein family name: ";
	getline(cin, pro);
	path workdir = path(dataset) / pro;

	if (!boost::filesystem::is_directory(workdir))
		handle_error(boost::format("Can not find directory %1%.") % workdir.string());

	Pro Apo;
	Pro Binding;
	Pro Allostery;
	Pro Complex;

	string aponame = "", bindingname = "", allosteryname = "", complexname = "";
	path pdbf = "";
	vector<string> excl = {};
	double spring_constant = 1.0;
	double cutoff = 9.0;

	cout << "Enter PDB name for apo state:";
	getline(cin, aponame);
	trim(aponame);
	if (aponame.empty())
		handle_error("Apo state must be loaded.");
	else
	{
		if (!ends_with(aponame, ".pdb"))
			aponame += ".pdb";
		pdbf = workdir / aponame;
		excl = get_exclude_res();
		spring_constant = get_spring_constant();
		cutoff = get_cutoff();
		Apo = Pro(pdbf.string(), false, excl, spring_constant, cutoff);
	}

	cout << "Enter PDB name for binding state:";
	getline(cin, bindingname);
	trim(bindingname);
	if (!bindingname.empty())
	{
		if (!boost::algorithm::ends_with(bindingname, ".pdb"))
			bindingname += ".pdb";
		pdbf = workdir / bindingname;
		excl = get_exclude_res();
		spring_constant = get_spring_constant();
		cutoff = get_cutoff();
		Binding = Pro(pdbf.string(), true, excl, spring_constant, cutoff);
	}

	cout << "Enter PDB name for allostery state:";
	getline(cin, allosteryname);
	trim(allosteryname);
	if (!allosteryname.empty())
	{
		if (!boost::algorithm::ends_with(allosteryname, ".pdb"))
			allosteryname += ".pdb";
		pdbf = workdir / allosteryname;
		excl = get_exclude_res();
		spring_constant = get_spring_constant();
		cutoff = get_cutoff();
		Allostery = Pro(pdbf.string(), true, excl, spring_constant, cutoff);
	}

	cout << "Enter PDB name for complex state:";
	getline(cin, complexname);
	trim(complexname);
	if (!complexname.empty())
	{
		if (!boost::algorithm::ends_with(complexname, ".pdb"))
			complexname += ".pdb";
		pdbf = workdir / complexname;
		excl = get_exclude_res();
		spring_constant = get_spring_constant();
		cutoff = get_cutoff();
		Complex = Pro(pdbf.string(), true, excl, spring_constant, cutoff);
	}
	
	ProAnalysis Cycle(Apo, Binding, Allostery, Complex);

	Cycle.write_hessian_binary(aponame + ".hessian");
	Cycle.write_covariance_binary(aponame + ".covariance");
	
	//ProAnalysis Cycle(Apo, Binding, Allostery, Complex, aponame + ".hessian", aponame + ".covariance");

	Cycle.interactive();
   
	return 0;
}