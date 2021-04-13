#include "handle_io.h"

void handle_info(string info)
{
	cout << "[Info] " << info << endl;
}

void handle_info(boost::format info)
{
	cout << "[Info] " << info << endl;
}

void handle_warning(string info)
{
	cout << "[Warning] " << info << endl;
}

void handle_warning(boost::format info)
{
	cout << "[Warning] " << info << endl;	
}

void handle_error(string info)
{
	cin.get();
	cout << "[Error] " << info << endl;
	cout << "Program interrupted. Press any key to exit...";
	cin.get();
	exit(1);
}

void handle_error(boost::format info)
{
	cin.get();
	cout << "[Error] " << info << endl;
	cout << "Program interrupted. Press any key to exit...";
	cin.get();
	exit(1);
}

void handle_result(string info)
{
	cout << "[Result] " << info << endl;
}

void handle_result(boost::format info)
{
	cout << "[Result] " << info << endl;
}


void handle_result(string info, vector<string> buf)
{
	cout << "[Result] " << info << endl;
	for (vector<string>::iterator it = buf.begin(); it != buf.end(); ++it)
		cout << *it << endl;
}

void handle_hint(string info)
{
	cout << info << endl;
}

void handle_hint(boost::format info)
{
	cout << info << endl;
}
