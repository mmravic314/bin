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

#if !defined(PROTEINSTRUCT_H)
#define PROTEINSTRUCT_H

#include "Common.h"

class ProteinStruct
{	
	public:
		ProteinStruct() { _ifp = NULL; }
		ProteinStruct(const string &);
		~ProteinStruct() {}

		double calcDistBB(const int &, const int &);
		double calcDistCA(const int &, const int &);
		void closeProteinStructFile();
		bool fullBB();
		vector<vector<double> >& getBBCoords() { return _bbcoor; }
		vector<vector<double> >& getCACoords() { return _cacoor; }
		string getFileName() { return _fname; }
		bool getFullBB(const int & i) const { return _fullbb[i]; }
		int getNumRes() const { return _numres; }
		vector<string>& getPdbInfo() { return _pdbinfo; }	
		void initPdbInfo();
		void readBBCoords();
		void readCACoords();
		void readNumRes();
		void readPdbInfo();
		void readPdbInfo(const vector<int> &);
		void readProteinStruct(const string &);
		void readSeq(char* &);
	
	protected:
			string _fname;		
			FILE* _ifp;
	
			vector<vector<double> > _bbcoor;
			vector<vector<double> > _cacoor;
			vector<bool> _fullbb;
			int _numres;
			vector<string> _pdbinfo;		
};

#endif
