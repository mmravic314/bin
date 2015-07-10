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
along with MaDCaT.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2014 Jianfu Zhou, Gevorg Grigoryan             
----------------------------------------------------------------------------
*/

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
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

void error(const string & msg)
{
  fprintf(stdout, "Error: %s.\n", msg.c_str());
  exit(-1);
}

#define ASSERT(cond, msg) {if (!(cond)) error(msg);};

void openFileC (FILE* & fp, const string & fname, const char* mode)
{
  fp = fopen(fname.c_str(), mode);
  ASSERT(fp != NULL, "could not open file " + fname);
}

void openFileCPP(fstream & fs, const string & fname, const ios_base::openmode & mode)
{
	fs.open(fname.c_str(), mode);
	ASSERT(fs.is_open(), "could not open file " + fname);
}

template <class T>
string toString (T val)
{
	return static_cast<ostringstream*>( &(ostringstream() << val) )->str();
}

char* trim(char* str) {
  char* nptr; int i;
  for (i = 0; i < strlen(str); i++) {
    if ((str[i] != '\n') && (str[i] != '\t') && (str[i] != ' ')) { break; }
  }
  nptr = str + i;
  if (i == strlen(str)) { return nptr; }

  for (i = strlen(str)-1; i >= 0; i--) {
    if ((str[i] != '\n') && (str[i] != '\t') && (str[i] != ' ')) {
      str[i+1] = '\0';
      break;
    }
  }
  return nptr;
}

void file2array(string _filename, vector<string> & lines) {
  FILE* ifp;
  int maxline = 1000;
  char *line, *tline;
  line = (char*) malloc(sizeof(char)*maxline);

  ifp = fopen(_filename.c_str(), "r");
  if (ifp == NULL) { error("unable to open file " + _filename); }

  while (fgets(line, maxline, ifp) != NULL) {
    if (line[strlen(line)-1] != '\n') { error("lines in file " + _filename + " are over " + toString(maxline) + "  long - increase max line limit and recompile."); }
    tline = trim(line);
    if (strlen(tline) > 0) { lines.push_back(line); }
    if (feof(ifp)) { break; }
  }
  if (ferror(ifp)) { error("an error occurred while reading file " + _filename); }
  free(line);
}

string fileBase(string fn) {
  if (fn.find_last_of(".") == string::npos) return fn;
  else return fn.substr(0, fn.find_last_of("."));
}

vector<string> tokenize(const std::string & _input, const std::string & _delimiter=" ", bool _allowEmtpy=false){
	vector<string> results;

	if (_input == "") {
		if (_allowEmtpy) {
			results.push_back(_input);
		}
		return results;
	}
	
	if (_allowEmtpy) {
		size_t prePos = 0;
		size_t pos  = _input.find(_delimiter);
		unsigned int delimiterSize = _delimiter.size();
		string left = _input, right;

		while (pos != std::string::npos) {
			results.push_back(left.substr(prePos, pos));
			if( pos + delimiterSize <= left.size() ) {
				left = left.substr(pos + delimiterSize, left.size() );
			} else {
				left = "";
			}
			pos  = left.find(_delimiter);
		}

		results.push_back(left);
	} else {
		int start  = _input.find_first_not_of(_delimiter);
		int end    = 0;
		string cur = _input;


		while (start != std::string::npos){
			end    = _input.find_first_of(_delimiter, start);
			results.push_back(_input.substr(start, end-start));
			start  = _input.find_first_not_of(_delimiter, end);
		}
	}
	return results;

}

string getFileName(string fullpath){

  string name = fullpath;

  // Assume linux paths '/' , sorry windows users!
  vector<string> paths = tokenize(name, "/");


  if (paths.size() > 0){
    name = paths[paths.size()-1];
  }

  int cut;
  if (name.find_last_of(".") > 1000 || name.find_last_of(".") <= 0){
	  cut = name.length();
  } else {
	  cut = name.length() - name.find_last_of(".");
  }

  name.erase(name.length() - cut);

  return name;
}

string optionUsage(string opt, string mes, int w, int p1, int p2) {
  // first print the name of the option
  string text(p1, ' ');
  text += opt;
  if (p2 > text.size()) text += string(p2 - text.size(), ' ');

  // next print the description text
  int i = 0, k, L = text.size(), n;
  while (i < mes.size()) {
    k = mes.find_first_of(" ", i);
    if (k == string::npos) k = mes.size();
    n = k - i;
    if ((L + n >= w) && (L > 0)) { text += "\n" + string(p2, ' '); L = p2; }
    text += mes.substr(i, n) + " ";
    L += n + 1;
    i = k+1;
  }
  return text;
}

void extractNumberofResidues(FILE* & ifp, int & numres, const string & psfname)
{
	numres = 0;
	ASSERT(fseek(ifp, sizeof(int) * SEC_NUMRES, SEEK_SET) == 0, "could not seek number of residues in file " + psfname);	
	ASSERT(fread((void*) &numres, sizeof(int), 1, ifp) == 1, "could not read number of residues in file " + psfname);
	ASSERT(numres > 0, "invalid number of residues in file " + psfname);
}

void extractEntirePdbInformation(FILE* & ifp, vector<string> & pdbinfo, const int & numres, const string & psfname)
{
	pdbinfo.clear();
	
	int i, offset, os_pdbinfo, pos, reslen, maxreslen = 10000; // guess as to what the max possible residue string will be
	char* resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
	
	ASSERT(fseek(ifp, sizeof(int) * SEC_PDBINFO, SEEK_SET) == 0, "could not seek offset of PDB information section in file " + psfname);
	ASSERT(fread((void*) &os_pdbinfo, sizeof(int), 1, ifp) == 1, "could not read offset of PDB information section in file " + psfname);
	ASSERT(fseek(ifp, os_pdbinfo, SEEK_SET) == 0, "could not seek PDB information section in file " + psfname);
	ASSERT(fread((void*) &offset, sizeof(int), 1, ifp) == 1, "could not read offset of PDB information in file " + psfname);
	pos = os_pdbinfo + offset;
	ASSERT(fseek(ifp, pos, SEEK_SET) == 0, "could not seek PDB information in file " + psfname);
	for (i = 0; i < numres; i++)
	{
		ASSERT(fread((void*) &reslen, sizeof(int), 1, ifp) == 1, "could not read length of PDB information in file " + psfname);
		if (reslen > maxreslen)
		{
			maxreslen = reslen;
			free((void*) resinfo);
			resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
		}
		ASSERT(fread((void*) resinfo, sizeof(char), reslen, ifp) == reslen, "could not read PDB information in file " + psfname);
		resinfo[reslen] = '\0';
		pdbinfo.push_back(string(resinfo));
	}
	free((void*) resinfo);
}

void usage()
{
	int w = 80, p1 = 3, p2 = p1 + 12; // based on the length of the longest option name
	cout << optionUsage("--pds", "input PDS file.", w, p1, p2) << endl;
	cout << optionUsage("--pdsList", "list of input PDS files, one per line.", w, p1, p2) << endl;
	cout << optionUsage("--pdb", "optional: output PDB file name for input PDS file. By default, takes base name of PDS file and appends \".pdb\".", w, p1, p2) << endl;
	cout << optionUsage("--pdbList", "optional: a file with a list of output PDB file names (one per input PDS file).", w, p1, p2) << endl;
}

class CreateOptions
{
	public:		
		string ext;
		vector<string> pdbfnames;
		vector<string> pdsfnames;

		CreateOptions();
		~CreateOptions() {}
};

CreateOptions::CreateOptions()
{
	ext = ".pdb";
	pdbfnames.clear();
	pdsfnames.clear();
}

void parseCommandLine(int argc, char** argv, CreateOptions & copts)
{
	map<string, bool> spec;

	while (1)
	{
		int oind = 0;
    	static struct option opts[] =
    	{
			{"pds", 1, 0, 1},
	  		{"pdsList", 1, 0, 2},
			{"pdb", 1, 0, 3},
			{"pdbList", 1, 0, 4},
      		{0, 0, 0, 0}
    	};

    	int c = getopt_long (argc, argv, "", opts, &oind);
    	if (c == -1)
    	{
			break;
    	}

    	switch (c) {
			case 1:
				if (0 != copts.pdsfnames.size())
				{
					usage();
					error("input PDS files already exist");
				}
        		copts.pdsfnames.push_back(string(optarg));
		        spec[string(opts[oind].name)] = true;
		        break;

			case 2:
				if (0 != copts.pdsfnames.size())
				{
					usage();
					error("input PDS files already exist");
				}
		        file2array(string(optarg), copts.pdsfnames);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 3:
				if (0 != copts.pdbfnames.size())
				{
					usage();
					error("output PDB files already exist");
				}
        		copts.pdbfnames.push_back(string(optarg));
		        spec[string(opts[oind].name)] = true;
		        break;

			case 4:
				if (0 != copts.pdbfnames.size())
				{
					usage();
					error("output PDB files already exist");
				}
		        file2array(string(optarg), copts.pdbfnames);
		        spec[string(opts[oind].name)] = true;
		        break;
		
      		case '?':
				usage();
				exit(-1);

      		default:
        		printf ("?? getopt returned character code %d ??\n", c);
				usage();
				exit(-1);
    	}
	}

  	// make sure all required options have been specified
  	if ((spec.find("pds") == spec.end()) && (spec.find("pdsList") == spec.end()))
  	{
		usage();
	    error("not all required options specified");
	}

	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
    	{
			printf ("%s ", argv[optind++]);
    	}
		printf ("\n");
		exit(-1);
	}

	if (0 == copts.pdbfnames.size())
	{
	  	for (int i = 0; i < copts.pdsfnames.size(); i++)
	  	{
	  		copts.pdbfnames.push_back(fileBase(copts.pdsfnames[i]) + copts.ext);
	  	}
	}
	ASSERT(copts.pdsfnames.size() == copts.pdbfnames.size(), "numbers of input PDS files and output PDB files not consistent");
}

int main(int argc, char *argv[])
{	
	CreateOptions copts;
	int i, j, numres;
	FILE *ifp = NULL;
	fstream ofs;
	vector<string> pdbinfo;
	
	parseCommandLine(argc, argv, copts);	

	for (i = 0; i < copts.pdsfnames.size(); i++)
	{
		openFileC(ifp, copts.pdsfnames[i], "rb");
		extractNumberofResidues(ifp, numres, copts.pdsfnames[i]);
		extractEntirePdbInformation(ifp, pdbinfo, numres, copts.pdsfnames[i]);
		fclose(ifp);
		openFileCPP(ofs, copts.pdbfnames[i], (ios::out | ios::trunc));
		for (j = 0; j < pdbinfo.size(); j++)
		{
			ofs << pdbinfo[j];
		}		
		ofs.close();
	}

	return 0;
}
