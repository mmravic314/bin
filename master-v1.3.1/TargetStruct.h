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

#if !defined(TARGETSTRUCT_H)
#define TARGETSTRUCT_H

#include "Common.h"
#include "ProteinStruct.h"

class TargetStruct : public ProteinStruct
{	
	public:
		TargetStruct() : ProteinStruct() {}
		TargetStruct(const string & tsfname);
		~TargetStruct() {}

		void initDihedDistr();
		void initDistDistr();
		void readDihedDistr(vector<int> & mcand, const double & phi, const double & psi, const double & phieps, const double & psieps);
		void readDistDistr(vector<int> & mcand, const int & ri, const double & d, const double & deps);
		void readTargetStruct(const string & tsfname);
	
	protected:
			int _os_diheddistr, _os_distdistr;
			double _phistep, _psistep, _dcut, _dstep;		
			vector<vector<vector<int> > > _diheddistr;
			vector<vector<vector<int> > > _distdistr;		
};

#endif
