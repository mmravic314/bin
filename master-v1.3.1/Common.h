/*
----------------------------------------------------------------------------
This file is part of MASTER.

MASTER is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

MASTER is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with MASTER.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2014 Jianfu Zhou, Gevorg Grigoryan
----------------------------------------------------------------------------
*/

#if !defined(COMMON_H)
#define COMMON_H

#undef DEBUG_BB
#undef DEBUG_CMDLINE
#undef DEBUG_GAPLEN
#undef DEBUG_PDBINFO
#undef DEBUG_PS
#undef DEBUG_QS
#undef DEBUG_TS
#undef DEBUG_ZSCORE
#undef TIME_OUTPUT
#undef TIME_SEARCH

#define IMPOSSIBLE_ANGLE	999.9
#define IMPOSSIBLE_COORD	999999.999

#define NUM_BBA				4
#define NUM_COORDS			3

#define N_IDX				0
#define CA_IDX				1
#define C_IDX				2
#define O_IDX				3

#define NOT_MATCH			0
#define MATCH_L1			1
#define MATCH_L2			2
#define MATCH_L3			3

#define LEN_AA_CODE			3

// Protein Structure
#define SEC_BBCOOR			0
#define SEC_CACOOR			1
#define SEC_NUMRES			2
#define SEC_PDBINFO			3
#define SEC_SEQ				4

#define NUM_SEC_PS			5

// Query Structure
#define SEC_BEFBRK			5
#define SEC_CENRES			6
#define SEC_CRDIHED			7
#define SEC_NUMSEG			8

#define NUM_SEC_QS			9

// Target Structure
#define SEC_DIHEDDISTR		5
#define SEC_DISTDISTR		6

#define NUM_SEC_TS			7

#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstddef> 
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <getopt.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

string basename(const string &, const bool &);
void error(const string &);
bool exist(const char *);
void file2array(string, vector<string> &);
int isBackboneAtom(const char *);
bool isDigit(const string &);
void openFileC (FILE* &, const string &, const char *);
void openFileCPP(fstream &, const string &, const ios_base::openmode &);
string optionUsage(string, string, int, int, int);

template <class T>
string toString(const T & val)
{
	return static_cast<ostringstream*>( &(ostringstream() << val) )->str();
}

char* trim(char *);
char* trimSpace(char *, const size_t &);

#define ASSERT(cond, msg) {	if (!(cond)) { error(msg); }}

#endif
