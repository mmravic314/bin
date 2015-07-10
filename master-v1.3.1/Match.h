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

#if !defined(MATCH_H)
#define MATCH_H

#include "Common.h"
#include "SearchResults.h"

class Match
{
	public:
		Match(const Match &);  // NOTE: (important) Match _must_ have a copy constructor so that insertion into set works properly (a copy needs to be made and Match has pointers)
		Match(SearchResults* p = NULL, const int & idx = -1);
		~Match();

		void delBegRes();
		void delRotation();
		void delTranslation();
		int* getBegRes() { return _bres; }
		int getBegRes(const int & i) const { return _bres[i]; }
		int getEndRes(const int & i) const { return (_bres[i] - 1 + _sr->getQsSegLen(i)); }
		unsigned long int getInsertTime() const { return _inTime; }
		double getMaxDistDev() const { return _maxdd; }
		int getNumSeg() const { return _sr->getQsNumSeg(); }
		double getRmsd() const { return _rmsd; }
		double** getRotation() { return _u; }
		SearchResults* getSearchResults() { return _sr; }
		int getSegLen(const int & i) const { return _sr->getQsSegLen(i); }
		char* getSeq() { return _sr->getSeq(_tsidx); }
		double* getTranslation() { return _t; }
		string getTsFile() const { return _sr->getTsFile(_tsidx); }
		int getTsIdx() const { return _tsidx; }
		int getTsNumRes() const { return _sr->getTsNumRes(_tsidx); }
		void initMatch(SearchResults* p = NULL, const int & idx = -1);
		void newBegRes(const int & num) { _bres = new int [num]; }
		void newRotation();
		void newTranslation() { _t = new double [3]; }
		bool operator< (const Match &) const;
		friend ostream& operator<<(ostream &, const Match &);
		int parseMatch(char*, const vector<vector<int> > &);
		void setBegRes(int* br, const int & num) { memcpy((void *)_bres, (const void *)br, num * sizeof(int)); }
		void setBegRes(const int & i, const int & idx) { _bres[i] = idx; }
		void setInsertTime(const unsigned long int & t) { _inTime = t; }
		void setMaxDistDev(const double & d) { _maxdd = d; } // NOTE: (low priority) remove this function and replace with updateMxDistDev(double d) { if (d > _maxdd) _maxdd = d; }
		void setRmsd(const double & r) { _rmsd = r; }
		void setRotation(double**);
		void setSearchResults(SearchResults* p) { _sr = p; }
		void setTranslation(double* t) { memcpy((void *)_t, (const void *)t, 3 * sizeof(double)); }
		void setTsIdx(const int & idx) { _tsidx = idx; }
		void setTsNumRes(const int & nr) { _sr->setTsNumRes(_tsidx, nr); }

	protected:		
		int* _bres; // beginning indices for each segment
		unsigned long int _inTime; // when this match was inserted (unique for this match)
		double _maxdd; // max inter-segment center-to-center distance
		double _rmsd; // full-BB or CA RMSD
		SearchResults* _sr; // pointer to search results which contain this match
		double* _t; // translation vector
		int _tsidx; // index of the target structure from which this match originates
		double** _u; // rotation matrix
};

bool cmpMatchByTarget(Match* a, Match* b);
bool cmpMatchByRmsd(Match* a, Match* b);
bool cmpMatchPairByRmsd(pair<Match*, int> a, pair<Match*, int> b);

#endif
