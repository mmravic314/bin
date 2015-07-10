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

#if !defined(CENRESGRP_H)
#define CENRESGRP_H

#include "CenterResidue.h"
#include "Common.h"

class CenResGrp
{
	public:
		CenResGrp() {}
		CenResGrp(const int & num);
		~CenResGrp() {}

		CenterResidue& operator[](const unsigned int & i) { return _cres[i]; }
		void initCenResGrp(const int & num);
		unsigned int size() const { return _cres.size(); }
		void sortByRank() { sort(_cres.begin(), _cres.end(), cmpCenResByRank); }

	protected:
		vector<CenterResidue> _cres;
};

#endif
