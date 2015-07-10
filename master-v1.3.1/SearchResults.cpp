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

#include "Match.h"
#include "SearchResults.h"

SearchResults::SearchResults(const int & len)
{
	_nMatFnd = 0;
	setMatchLimit(len);
}

SearchResults::SearchResults(const string & mifname, const vector<vector<int> > & gaplen, const int & len)
{
	_nMatFnd = 0;
	setMatchLimit(len);
	setMatchList(mifname, gaplen);
}

SearchResults::~SearchResults()
{
	clear();
}

int SearchResults::addTsFile(const string & fn)
{
	if (_tsfmap.find(fn) == _tsfmap.end())
	{
		int i = _tsfmap.size();
		_tsfmap[fn] = i;
		_tsfiles.push_back((string *) &(_tsfmap.find(fn)->first));
	}
	return _tsfmap[fn];
}

void SearchResults::clear()
{
	_L.clear();
	_limit = 0;
	_nMatFnd = 0;
	_pL.clear();	
	_qs.closeProteinStructFile();
	for (int i = 0; i < _seq.size(); i++)
	{
		delete [] _seq[i];
	}
	_seq.clear();
	_tsfiles.clear();
	_tsfmap.clear();
	_tsnr.clear();
}

void SearchResults::setPointerVector()
{
	_pL.assign(_L.size(), NULL);
	set<Match>::iterator it;
	int i = 0;
	for (it = _L.begin(); it != _L.end(); ++it)
	{
		_pL[i++] = (Match *) &(*it);
	}
}

void SearchResults::setTsFiles(const vector < string > & list)
{
	int i;
	_tsfiles.assign(list.size(), NULL);
	for (i = 0; i < list.size(); i++)
	{
		_tsfiles[i] = (string *) &(list[i]);
	}
}

void SearchResults::setTsNumRes(const int & i, const int & nr)
{
	ASSERT(_tsnr[i] <= 0, "each target should set its number of residues only once");
	_tsnr[i] = nr;
}

// set match list
bool SearchResults::full()
{
	if ((_limit > 0) && (_L.size() >= _limit))
	{
		return true;
	}
	return false;
}

void SearchResults::insertMatch(const Match & m)
{
	_L.insert(m);
	_nMatFnd++;
	if ((_limit > 0) && (_L.size() > _limit))
	{
		_L.erase(--_L.end());
	}
}

void SearchResults::setMatchLimit(const int & len)
{
	_limit = len;
}

void SearchResults::setMatchList(const string & mifname, const vector<vector<int> > & gaplen)
{
	ASSERT(_L.empty(), "match list already exists");
	
	FILE* ifp;
	int maxline = 1000;
	char *line = (char*) malloc(maxline * sizeof(char));
	Match m(this);

	openFileC(ifp, mifname, "r");

	while (fgets(line, maxline, ifp) != NULL)
	{
		ASSERT(line[strlen(line) - 1] == '\n', "lines in file " + mifname + " are over " + toString(maxline) + " long. Please increase max line limit and recompile.");
		if (strlen(line) > 0)
		{
			if (0 == m.parseMatch(line, gaplen))
			{
				insertMatch(m);
			}
		}
		if (feof(ifp))
		{
			break;
		}
	}

	ASSERT(0 == ferror(ifp), "could not read file " + mifname);	
	fclose(ifp);
	free(line);
}

double SearchResults::worstRmsd() const
{
	return _L.rbegin()->getRmsd();
}


// vector match list
ostream& operator<<(ostream & os, SearchResults & L)
{
	int i;
	for (i = 0; i < L.size(); i++)
	{
		os << L[i] << "\n";
	}

	return os;
}

void SearchResults::sortByTarget()
{
	sort(_pL.begin(), _pL.end(), cmpMatchByTarget);
}

void SearchResults::sortByRmsd(vector<int>* old)
{
	if (old == NULL) {
		sort(_pL.begin(), _pL.end(), cmpMatchByRmsd);
	} else {
		vector<pair<Match*, int> > toSort(_pL.size());
		int i;
		for (i = 0; i < _pL.size(); i++)
		{
			toSort[i].first = _pL[i];
			toSort[i].second = i + 1;
		}
		sort(toSort.begin(), toSort.end(), cmpMatchPairByRmsd);
		old->assign(toSort.size(), -1);
		for (i = 0; i < toSort.size(); i++)
		{
			_pL[i] = toSort[i].first;
			(*old)[i] = toSort[i].second;
		}
	}
}

