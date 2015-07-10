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

#if !defined(SEARCHRESULTS_H)
#define SEARCHRESULTS_H

#include "Common.h"
#include "QueryStruct.h"

class Match;

class SearchResults
{
	public:
		SearchResults(const int & len = 0);
		SearchResults(const string & mifname, const vector<vector<int> > & gaplen, const int & len = 0);
		~SearchResults();

		int addTsFile(const string &); // NOTE: (low priority) may want to make this function protected and make Match a friend of SearchResults
		void clear();
		unsigned long int getNumMatFnd() const { return _nMatFnd; }
		int getNumTs() const { return _tsfiles.size(); }
		QueryStruct& getQs() { return _qs; }
		int getQsNumSeg() const { return _qs.getNumSeg(); }
		int getQsSegLen(const int & i) const { return (_qs.getResBefBrk(i + 1) - _qs.getResBefBrk(i)); }
		vector<char*>& getSeq() { return _seq; }
		char* getSeq(const int & i) { return _seq[i]; }
		string getTsFile(const int & i) const { return *(_tsfiles[i]); }
		int getTsNumRes(const int & i) const { return _tsnr[i]; }
		void initTsNumRes() { _tsnr.assign(_tsfiles.size(), -1); }
		void setPointerVector();
		void setQs(const string & fn) { _qs.readQueryStruct(fn); }
		void setTsFiles(const vector<string> &);
		void setTsNumRes(const int &, const int &);
		
		// set match list
		bool empty() const { return _L.empty(); }
		bool full();
		int getMatchLimit() const { return _limit; }
		void insertMatch(const Match &);
		void setMatchLimit(const int &);
		void setMatchList(const string &, const vector<vector<int> > &);
		unsigned int size() const { return _L.size(); }
		double worstRmsd() const;

		// vector match list
		Match& operator[](const unsigned int & i) { return *(_pL[i]); }
		friend ostream& operator<<(ostream &, SearchResults &);
		void sortByTarget();
		void sortByRmsd(vector<int>* old = NULL);
	
	protected:
		set<Match> _L;
		int _limit;
		unsigned long int _nMatFnd;
		vector<Match*> _pL;
		QueryStruct _qs;
		vector<char*> _seq;		
		vector<string*> _tsfiles;
		map<string, int> _tsfmap;
		vector<int> _tsnr;
};

#endif
