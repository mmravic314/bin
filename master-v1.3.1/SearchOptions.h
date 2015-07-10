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

#if !defined(SEARCHOPTIONS_H)
#define SEARCHOPTIONS_H

#include "Common.h"

class SearchOptions
{
	public:		
		SearchOptions();
		~SearchOptions() {}

		bool getBBRmsd() const { return _bbrmsd; }
		int getDistDevBoundMode() const { return _dmode; }
		double getDistDevCut() const { return _deps; }
		bool getDistDevZscore() const { return _ddzscore; }
		vector<vector<int> >& getGapLen() { return _gaplen; }
		string getMatchInFile() const { return _mifname; }
		string getMatchOutFile() const { return _mofname; }		
		double getPhiDevCut() const { return _phieps; }
		double getPsiDevCut() const { return _psieps; }
		string getQsFile() const { return _qsfname; }
		int getRmsdBoundMode() const { return _rmode; }
		double getRmsdCut() const { return _rthresh; }
		string getSeqOutFile() const { return _sofname; }
		string getStructOutDir() const { return _sodname; }
		string getOutType() const { return _otype; }
		int getTopN() const { return _topn; }
		vector<string>& getTsFiles() { return _tsfnames; }
		double getTuningParam() const { return _tune; }
		void setBBRmsd(const bool & f) { _bbrmsd = f; }
		void setDistDevCut(const char*);
		void setDistDevZscore(const bool & f) { _ddzscore = f; }
		void setGapLen(const string &);
		void setMatchInFile(const string & fn) { _mifname = fn; }
		void setMatchOutFile(const string & fn) { _mofname = fn; }
		void setPhiDevCut(const char*);
		void setPsiDevCut(const char*);
		void setQsFile(const string & fn) { _qsfname = fn; }
		void setRmsdBoundMode(const char*);
		void setRmsdCut(const char*);
		void setSeqOutFile(const string & fn) { _sofname = fn; }
		void setStructOutDir(const string &);
		void setOutType(const string &);
		void setTopN(const char*);
		void setTsFile(const string &);
		void setTsFiles(const string &);
		void setTuningParam(const char*);

	protected:
		bool _bbrmsd; // flag for full-backbone RMSD
		string _ddofname;
		bool _ddzscore; // flag for distance deviation Z-score
		double _deps; // distance deviation cutoff
		int _dmode; // distance deviation bounding mode
		vector<vector<int> > _gaplen; // user-defined inter-segment gap lengths
		string _mifname; // input match list file
		string _mofname; // output match list file
		double _phieps; // Phi angle deviation cutoff
		double _psieps; // Psi angle deviation cutoff
		string _qsfname; // query structure file		
		int _rmode; // CA RMSD bounding mode
		double _rthresh; // CA RMSD cutoff
		string _sodname; // output structure directory
		string _sofname; // output sequence file
		string _otype; // output structure type
		double _tune; // tuning parameter for CA RMSD bounding
		int _topn; // keep the best this many hits during search
		vector<string> _tsfnames; // target structure files
};

#endif
