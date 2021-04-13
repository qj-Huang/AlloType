#include "read_pdb.h"

string sslice(size_t begin, size_t end, string in)
{
	string out = in.substr(begin, end - begin);
	trim(out);
	return out;
}

size_t read_resid(string line)
{
	return lexical_cast<size_t>(sslice(22, 26, line));
}

string read_record(string line)
{
	return sslice(0, 6, line);
}

string read_resname(string line)
{
	return sslice(17, 20, line);
}

string read_chain(string line)
{
	return sslice(21, 22, line);
}

string read_atomname(string line)
{
	return sslice(12, 16, line);
}

ResInfo read_res(string line)
{
	ResInfo res;
	res.resname = sslice(17, 20, line);
	res.chain = sslice(21, 22, line);
	res.resid = lexical_cast<size_t>(sslice(22, 26, line));
	res.x = lexical_cast<double>(sslice(30, 38, line));
	res.y = lexical_cast<double>(sslice(38, 46, line));
	res.z = lexical_cast<double>(sslice(46, 54, line));
	res.bfactor = lexical_cast<double>(sslice(60, 66, line));
	return res;
}

AtomInfo read_atom(string line)
{
	AtomInfo atom;
	atom.atomname = sslice(12, 16, line);
	atom.resname = sslice(17, 20, line);
	atom.chain = sslice(21, 22, line);
	atom.resid = lexical_cast<size_t>(sslice(22, 26, line));
	atom.x = lexical_cast<double>(sslice(30, 38, line));
	atom.y = lexical_cast<double>(sslice(38, 46, line));
	atom.z = lexical_cast<double>(sslice(46, 54, line));
	atom.bfactor = lexical_cast<double>(sslice(60, 66, line));
	return atom;
}