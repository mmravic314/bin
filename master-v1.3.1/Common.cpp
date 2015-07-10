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

#include "Common.h"

string basename(const string & fname, const bool & ext)
{
	if (fname.find_last_of("/") == string::npos)
	{
		if (fname.find_last_of(".") == string::npos)
		{
			return fname;
		}
		else
		{
			if (ext)
			{
				return fname;
			}
			else
			{
				return fname.substr(0, fname.find_last_of("."));
			}
		}
	}
	else
	{
		if (fname.find_last_of(".") == string::npos)
		{
			return fname.substr(fname.find_last_of("/") + 1);
		}
		else
		{
			if (ext)
			{
				return fname.substr(fname.find_last_of("/") + 1);
			}
			else
			{
				return fname.substr(fname.find_last_of("/") + 1, 
				fname.find_last_of(".") - 1 - fname.find_last_of("/"));
			}
		}
	}
}

void error(const string & msg)
{
  fprintf(stdout, "Error: %s.\n", msg.c_str());
  exit(-1);
}

bool exist(const char *path)
{
	struct stat buffer;
	if (stat(path, &buffer) == 0)
	{
		return true;
	}
	return false;
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
  fclose(ifp);
  free(line);
}

int isBackboneAtom(const char* an)
{
	if (strcmp(an, "N") == 0)
	{
		return 0;
	}
	if (strcmp(an, "NT") == 0)
	{
		return 10;
	}
	if (strcmp(an, "CA") == 0)
	{
		return 1;
	}
	if (strcmp(an, "C") == 0)
	{
		return 2;
	}
	if (strcmp(an, "O") == 0)
	{
		return 3;
	}
	if (strcmp(an, "OT1") == 0)
	{
		return 13;
	}
	if (strcmp(an, "OT2") == 0)
	{
		return 23;
	}
	if (strcmp(an, "OXT") == 0)
	{
		return 33;
	}

	return -1;
}

bool isDigit(const string & str)
{
	if (str.length() <= 0)
	{
		return false;
	}	
	int i;
	for (i = 0; i < str.length(); i++)
	{
		if (isdigit(str[i]))
		{
			continue;
		}
		return false;
	}
	return true;
}

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

// places a string end character after the last non-space character and returns a new
// point that points to the first non-space character
/*
char* trimSpace(char* str) {
  char* nptr; int i;
  for (i = 0; i < strlen(str); i++) {
    if (str[i] != ' ') { break; }
  }
  nptr = str + i;
  if (i == strlen(str)) { return nptr; }

  for (i = strlen(str)-1; i >= 0; i--) {
    if (str[i] != ' ') {
      str[i+1] = '\0';
      break;
    }
  }
  return nptr;
}
*/
char* trimSpace(char* str, const size_t & sz)
{
	char* nptr;
	int i;
	for (i = 0; i < sz; i++)
	{
		if (str[i] != ' ')
		{
			break;
		}
	}
	nptr = str + i;

	for (i = sz - 1; i >= 0; i--)
	{
		if (str[i] != ' ')
		{
			str[i + 1] = '\0';
			break;
		}
	}

	return nptr;
}

