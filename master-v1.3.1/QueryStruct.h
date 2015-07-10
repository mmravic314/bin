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

#if !defined(QUERYSTRUCT_H)
#define QUERYSTRUCT_H

#include "Common.h"
#include "ProteinStruct.h"

class QueryStruct : public ProteinStruct
{
	public:
		QueryStruct() : ProteinStruct() {}
		QueryStruct(const string & qsfname);
		~QueryStruct() {}

		int getCenRes(const int & i) const { return _cenres[i]; }
		vector<int>& getCenRes() { return _cenres; }
		vector<vector<double> >& getCenResDihed() { return _crdihed; }
		int getNumSeg() const { return _numseg; }
		int getResBefBrk(const int & i) const { return _befbrk[i]; }
		vector<int>& getResBefBrk() { return _befbrk; }		

		void readCenRes();
		void readCenResDihed();
		void readNumSeg();
		void readQueryStruct(const string & qsfname);
		void readResBefBrk();

	protected:
		int _numseg;
		vector<int> _befbrk;
		vector<int> _cenres;
		vector<vector<double> > _crdihed;
};

#endif
