#pragma once
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using boost::lexical_cast;
using boost::algorithm::trim;

struct ResInfo
{
	string resname;
	string chain;
	size_t resid;
	double x;
	double y;
	double z;
	double bfactor;
};
struct AtomInfo
{
	string atomname;
	string resname;
	string chain;
	size_t atomid;
	size_t resid;
	double x;
	double y;
	double z;
	double bfactor;
};

string sslice(size_t begin, size_t size, string in);

size_t read_resid(string line);

string read_record(string line);

string read_resname(string line);

string read_chain(string line);

string read_atomname(string line);

ResInfo read_res(string line);

AtomInfo read_atom(string line);
