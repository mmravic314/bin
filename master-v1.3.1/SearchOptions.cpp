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

#include "SearchOptions.h"

extern void usage();

SearchOptions::SearchOptions()
{
	_bbrmsd = false;
	_dmode = 0;
	_ddzscore = false;
	_phieps = 180.0;
	_psieps = 180.0;
	_rmode = 0;
	_otype = "match";
	_tune = 0.5;
	_topn = 0;
}

void SearchOptions::setDistDevCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_deps)) == 1) && (_deps >= 0.0)))
	{
		usage();
		error("bad greedy distance deviation cutoff value");
	}
	_dmode = 1;
}

void SearchOptions::setGapLen(const string & str)
{
	ASSERT(_gaplen.size() == 0, "gap length constraints have been specified");
	
	string subs, subss;
	size_t cur, next, subc, subn;

	// the first gap length constraint
	_gaplen.push_back(vector<int>());
	// looking for the first delimiter (i.e., semicolon)
	cur = 0;
	next = str.find_first_of(";", cur);
	while (next != string::npos)
	{
		// found a delimiter, indicating the end of the current gap length constraint		
		subs = str.substr(cur, next - cur);
		if (subs.length() > 0)
		{
			// looking for the connector (i.e., minus)
			subc = 0;
			subn = subs.find_first_of("-", subc);
			while (subn != string::npos)
			{
				// found a connector, indicating the end of the current component of the current gap length constraint
				subss = subs.substr(subc, subn - subc);
				if (!isDigit(subss))
				{
					usage();
					error("bad loop length value");
				}
				// found a valid component of the current gap length constraint
				_gaplen[_gaplen.size() - 1].push_back(atoi(subss.c_str()));
				// looking for the next connector, though only one connector is allowed here
				subc = subn + 1;
				subn = subs.find_first_of("-", subc);				
			}
			// found no connector
			if (subc < subs.length())
			{
				// found the last component of the current gap length constraint
				subss = subs.substr(subc);
				if (!isDigit(subss))
				{
					usage();
					error("bad loop length value");
				}
				// found a valid component of the current gap length constraint
				_gaplen[_gaplen.size() - 1].push_back(atoi(subss.c_str()));
			}
			// a connector should be followed by a vaild component
			if ('-' == subs[subs.length() - 1])
			{
				usage();
				error("bad loop length value");
			}
			// at most two valid components are allowed here
			if (_gaplen[_gaplen.size() - 1].size() > 2)
			{
				usage();
				error("bad loop length value");
			}
		}
		// a delimiter should be followed by the next gap length constraint		
		_gaplen.push_back(vector<int>());
		// looking for the next delimiter		
		cur = next + 1;
		next = str.find_first_of(";", cur);		
	}
	// found no delimiter
	if (cur < str.length())
	{
		// found the last gap length constraint
		subs = str.substr(cur);
		// looking for the connector (i.e., minus)
		subc = 0;
		subn = subs.find_first_of("-", subc);
		while (subn != string::npos)
		{
			// found a connector, indicating the end of the current component of the current gap length constraint
			subss = subs.substr(subc, subn - subc);
			if (!isDigit(subss))
			{
				usage();
				error("bad loop length value");
			}
			// found a valid component of the current gap length constraint
			_gaplen[_gaplen.size() - 1].push_back(atoi(subss.c_str()));
			// looking for the next connector, though only one connector is allowed here
			subc = subn + 1;
			subn = subs.find_first_of("-", subc);				
		}
		// found no connector
		if (subc < subs.length())
		{
			// found the last component of the current gap length constraint
			subss = subs.substr(subc);
			if (!isDigit(subss))
			{
				usage();
				error("bad loop length value");
			}
			// found a valid component of the current gap length constraint
			_gaplen[_gaplen.size() - 1].push_back(atoi(subss.c_str()));
		}
		// a connector should be followed by a vaild component
		if ('-' == subs[subs.length() - 1])
		{
			usage();
			error("bad loop length value");
		}
		// at most two valid components are allowed here
		if (_gaplen[_gaplen.size() - 1].size() > 2)
		{
			usage();
			error("bad loop length value");
		}
	}
#if defined(DEBUG_GAPLEN)
	int i, j;
	for (i = 0; i < _gaplen.size(); i++)
	{
		for (j = 0; j < _gaplen[i].size(); j++)
		{
			cout << " " << _gaplen[i][j];
		}
		cout << "\n";
	}
	exit(-1);
#endif
}

void SearchOptions::setPhiDevCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_phieps)) == 1) && (_phieps >= 0.0)))
	{
		usage();
		error("bad Phi angle deviation cutoff value");
	}
}

void SearchOptions::setPsiDevCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_psieps)) == 1) && (_psieps >= 0.0)))
	{
		usage();
		error("bad Psi angle deviation cutoff value");
	}
}

void SearchOptions::setRmsdBoundMode(const char * s)
{
	if (!((sscanf(s, "%d", &(_rmode)) == 1) && ((_rmode == 0) || (_rmode == 1) || (_rmode == 2))))
	{
		usage();
		error("bad CA RMSD bounding mode value");
	}
}

void SearchOptions::setRmsdCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_rthresh)) == 1) && (_rthresh > 0.0)))
	{
		usage();
		error("bad CA RMSD cutoff value");
	}
}

void SearchOptions::setStructOutDir(const string & dir)
{
	_sodname = dir;
	if ((!_sodname.empty()) && (!exist(_sodname.c_str())))
  	{
  		ASSERT(mkdir(_sodname.c_str(), 0755) == 0, "could not make directory " + _sodname);
  	}
}

void SearchOptions::setOutType(const string & type)
{
	_otype = type;
	if ((_otype.compare("full") != 0) && (_otype.compare("match") != 0) && (_otype.compare("wgap") != 0))
	{
		usage();
		error("bad output type value");
	}
}

void SearchOptions::setTopN(const char * s)
{
	if (!((sscanf(s, "%d", &(_topn)) == 1) && (_topn >= 0)))
	{
		usage();
		error("bad topN value");
	}
}

void SearchOptions::setTsFile(const string & fn)
{
	ASSERT(_tsfnames.size() == 0, "targets already exist");
	_tsfnames.push_back(fn);
}

void SearchOptions::setTsFiles(const string & list)
{
	ASSERT(_tsfnames.size() == 0, "targets already exist");
	file2array(list, _tsfnames);
}

void SearchOptions::setTuningParam(const char * s)
{
	if (!((sscanf(s, "%lf", &(_tune)) == 1) && (_tune >= 0.0) && (_tune <= 1.0)))
	{
		usage();
		error("bad tuning parameter value");
	}
}
