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

#if !defined(CENTERRESIDUE_H)
#define CENTERRESIDUE_H

#include "Common.h"

class CenterResidue
{
	public:		
		CenterResidue() {}
		~CenterResidue() {}

		int getBegRes() const { return _bres; }
		double getCurDistDevCut() const { return _cddc; }
		int getCurMatchCandIdx() const { return _cmci; }
		int getEndRes() const { return _eres; }
		vector<int>& getGapLen() { return _gaplen; }
		vector<int>& getMatchCand() { return _mcand; }
		int getMatchCand(const int & i) const { return _mcand[i]; }
		int getMatchCandBegRes() const { return _mcbres; }
		int getMatchCandEndRes() const { return _mceres; }
		double getPhi() const { return _phi; }
		double getPsi() const { return _psi; }
		float getRank() const { return _rank; }
		int getResIdx() const { return _rind; }
		int getSegIdx() const { return _sind; }
		int getSegLen() const { return _len; }
		double getSelfRmsdCut() const { return _srthresh; }
		int getSofarLen() const { return _lsofar; }
		double getSofarRmsd() const { return _rsofar; }
		double getSofarRmsdCut() const { return _mrthresh; }
		int getTsResFlag(const int & i) const { return _tsrflag[i]; }
		double getTsResRmsd(const int & i) const { return _tsrrmsd[i]; }
		void initTsResFlag(const int & num) { _tsrflag.assign(num, NOT_MATCH); }
		void initTsResRmsd(const int & num) { _tsrrmsd.assign(num, DBL_MAX); }
		void setBegRes(const int & idx) { _bres = idx; }
		void setCurDistDevCut(const double & ddcut) { _cddc = ddcut; }
		void setCurMatchCandIdx(const int & idx) { _cmci = idx; }
		void setEndRes(const int & idx) { _eres = idx; }
		void setGapLen(const vector<int> & len) { _gaplen = len; }
		void setMatchCandBegRes(const int & idx) { _mcbres = idx; }
		void setMatchCandEndRes(const int & idx) { _mceres = idx; }
		void setPhi(const double & phi) { _phi = phi; }
		void setPsi(const double & psi) { _psi = psi; }
		void setRank(const float & r) { _rank = r; }
		void setResIdx(const int & idx) { _rind = idx; }
		void setSegIdx(const int & idx) { _sind = idx; }
		void setSegLen(const int & len) { _len = len; }
		void setSelfRmsdCut(const double & rcut) { _srthresh = rcut; }
		void setSofarLen(const int & len) { _lsofar = len; }
		void setSofarRmsd(const double & r) { _rsofar = r; }
		void setSofarRmsdCut(const double & rcut) { _mrthresh = rcut; }
		void setTsResFlag(const int & i, const int & f) { _tsrflag[i] = f; }
		void setTsResRmsd(const int & i, const double & r) { _tsrrmsd[i] = r; }

	protected:
		int _bres; // begin residue of segment
		double _cddc; // current distance deviation cutoff
		int _cmci; // current match candidate index
		int _eres; // end residue of segment
		int _len; // length of segment
		vector<int> _gaplen; // gap length
		int _lsofar; // total length of segments so far
		vector<int> _mcand; // match candidates
		int _mcbres, _mceres; // begin and end residues of segment around match candidate
		double _mrthresh; // rmsd threshold for all the segments so far, assuming the remaining part has perfect match		
		double _phi, _psi; // dihedral angles of central residue
		float _rank; // (# match candidates / segment length)
		int _rind; // central residue index
		double _rsofar; // RMSD of segments so far
		int _sind; // segment index
		double _srthresh; // rmsd threshold for single segment, assuming all the other segments have perfect matches
		vector<int> _tsrflag; // target structure residue flags
							 // 0: not pass some filter
						   	 // 1: pass dihedral angles filter and segment length filter
						   	 // 2: pass 1 and individual RMSD filter
		vector<double> _tsrrmsd; // RMSD of candidate matching segment centering at each residue of target structure
};

bool cmpCenResByRank(const CenterResidue & a, const CenterResidue & b);
#endif
