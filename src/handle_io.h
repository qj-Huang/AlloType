#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <boost/format.hpp>

using namespace std;

void handle_info(string info);
void handle_info(boost::format info);

void handle_warning(string info);
void handle_warning(boost::format info);

void handle_error(string info);
void handle_error(boost::format info);

void handle_result(string info);
void handle_result(boost::format info);
void handle_result(string info, vector<string> buf);

void handle_hint(string info);

void handle_hint(boost::format info);
