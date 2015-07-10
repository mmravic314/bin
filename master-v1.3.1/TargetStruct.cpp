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

#include "TargetStruct.h"

TargetStruct::TargetStruct(const string & tsfname)
{
	readTargetStruct(tsfname);
}

void TargetStruct::initDihedDistr()
{
	ASSERT(fseek(_ifp, sizeof(int) * SEC_DIHEDDISTR, SEEK_SET) == 0, "could not seek offset of dihedral angles distribution section in file " + _fname);
	ASSERT(fread((void*) &_os_diheddistr, sizeof(int), 1, _ifp) == 1, "could not read offset of dihedral angles distribution section in file " + _fname);
	ASSERT(fseek(_ifp, _os_diheddistr, SEEK_SET) == 0, "could not seek dihedral angles distribution section in file " + _fname);
	ASSERT(fread((void*) &_phistep, sizeof(double), 1, _ifp) == 1, "could not read phi step in file " + _fname);
	ASSERT(fread((void*) &_psistep, sizeof(double), 1, _ifp) == 1, "could not read psi step in file " + _fname);
	int phinumbin = int(ceil(360.0 / _phistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	int psinumbin = int(ceil(360.0 / _psistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	_diheddistr.assign(phinumbin, vector<vector<int> >(psinumbin));
}

void TargetStruct::initDistDistr()
{
	ASSERT(fseek(_ifp, sizeof(int) * SEC_DISTDISTR, SEEK_SET) == 0, "could not seek offset of distance distribution section in file " + _fname);
	ASSERT(fread((void*) &_os_distdistr, sizeof(int), 1, _ifp) == 1, "could not read offset of distance distribution section in file " + _fname);
	ASSERT(fseek(_ifp, _os_distdistr, SEEK_SET) == 0, "could not seek distance distribution section in file " + _fname);
	ASSERT(fread((void*) &_dcut, sizeof(double), 1, _ifp) == 1, "could not read distance cutoff in file " + _fname);
	ASSERT(fread((void*) &_dstep, sizeof(double), 1, _ifp) == 1, "could not read distance step in file " + _fname);
	int numbin = int(ceil(_dcut / _dstep));
	_distdistr.assign(_numres, vector<vector<int> >(numbin));
}

void TargetStruct::readDihedDistr(vector<int> & mcand, const double & phi, const double & psi, const double & phieps, const double & psieps)
{
	mcand.clear();

	int i;
	if ((phieps >= 180.0) && (psieps >= 180.0))
	{
		for (i = 0; i < _numres; i++)
		{
			mcand.push_back(i);
		}
		return;
	}

	int numphibin = int(ceil(360.0 / _phistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	int numpsibin = int(ceil(360.0 / _psistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	vector<int> phibi, psibi; // bin indices
	int j;
	if (phi >= IMPOSSIBLE_ANGLE)
	{
		for (i = 0; i < numphibin; i++)
		{
			phibi.push_back(i);
		}
	}
	else
	{
		int bphibin = int(ceil((phi - phieps + 180.0) / _phistep)) - 1;
		int ephibin = int(ceil((phi + phieps + 180.0) / _phistep)) - 1;
		if (ephibin - bphibin >= numphibin - 2)
		{
			for (i = 0; i < numphibin; i++)
			{
				phibi.push_back(i);
			}
		}
		else
		{
			for (i = 0; (bphibin + i) <= ephibin; i++)
			{
				j = bphibin + i;
				while ((j < 0) || (j > (numphibin - 2)))
				{
					if (j < 0)
					{
						j += (numphibin - 1);
					}
					if (j > (numphibin - 2))
					{
						j -= (numphibin - 1);
					}
				}
				phibi.push_back(j);
			}
			sort(phibi.begin(), phibi.end());
			phibi.push_back(numphibin - 1);
		}
	}	
	if (psi >= IMPOSSIBLE_ANGLE)
	{
		for (i = 0; i < numpsibin; i++)
		{
			psibi.push_back(i);
		}
	}
	else
	{
		int bpsibin = int(ceil((psi - psieps + 180.0) / _psistep)) - 1;
		int epsibin = int(ceil((psi + psieps + 180.0) / _psistep)) - 1;
		if (epsibin - bpsibin >= numpsibin - 2)
		{
			for (i = 0; i < numpsibin; i++)
			{
				psibi.push_back(i);
			}
		}
		else
		{
			for (i = 0; (bpsibin + i) <= epsibin; i++)
			{
				j = bpsibin + i;
				while ((j < 0) || (j > (numpsibin - 2)))
				{
					if (j < 0)
					{
						j += (numpsibin - 1);
					}
					if (j > (numpsibin - 2))
					{
						j -= (numpsibin - 1);
					}
				}
				psibi.push_back(j);
			}
			sort(psibi.begin(), psibi.end());
			psibi.push_back(numpsibin - 1);
		}
	}	

	if ((phibi.size() >= numphibin) && (psibi.size() >= numpsibin))
	{
		for (i = 0; i < _numres; i++)
		{
			mcand.push_back(i);
		}
		return;
	}
	
	vector<int> offset;
	int numelem, pos;
	int* mc = (int*) malloc(_numres * sizeof(int));
	for (i = 0; i < phibi.size(); i++)
	{
		offset.assign(psibi.size(), 0);

		for (j = 0; j < psibi.size(); j++)
		{
			if (_diheddistr[phibi[i]][psibi[j]].size() > 0)
			{
				if (_diheddistr[phibi[i]][psibi[j]].size() > 1)
				{
					mcand.insert(mcand.end(), _diheddistr[phibi[i]][psibi[j]].begin() + 1, _diheddistr[phibi[i]][psibi[j]].end());
				}
				continue;
			}
			pos = _os_diheddistr + sizeof(double) * 2 + sizeof(int) * (phibi[i] * numpsibin + psibi[j]);
			if (pos != ftell(_ifp))
			{
				ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek offset of dihedral angles distribution");
			}
			ASSERT(fread((void*) &offset[j], sizeof(int), 1, _ifp) == 1, "could not read offset of dihedral angles distribution");
			if (offset[j] == 0)
			{
				_diheddistr[phibi[i]][psibi[j]].push_back(offset[j]);
			}
		}

		for (j = 0; j < psibi.size(); j++)
		{
			if (offset[j] == 0)
			{
				continue;
			}
			pos = _os_diheddistr + offset[j];
			if (pos != ftell(_ifp))
			{
				ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek dihedral angles distribution");
			}			
			ASSERT(fread((void*) &numelem, sizeof(int), 1, _ifp) == 1, "could not read number of elements in dihedral angles distribution");
			_diheddistr[phibi[i]][psibi[j]].push_back(numelem);
			ASSERT(fread((void*) mc, sizeof(int), numelem, _ifp) == numelem, "could not read dihedral angles distribution");
			_diheddistr[phibi[i]][psibi[j]].insert(_diheddistr[phibi[i]][psibi[j]].end(), mc, mc + numelem);
			mcand.insert(mcand.end(), mc, mc + numelem);
		}
	}	
	free((void*) mc);
}

void TargetStruct::readDistDistr(vector<int> & mcand, const int & ri, const double & d, const double & deps)
{
	mcand.clear();
	if (d - deps > _dcut)
	{
		return;
	}
		
	int bbin = int(ceil(((d - deps < _dstep) ? _dstep : (d - deps)) / _dstep)) - 1;
	int ebin = int(ceil(((d + deps > _dcut) ? _dcut : (d + deps)) / _dstep)) - 1;
	int numbin = int(ceil(_dcut / _dstep));
		
	vector<int> offset(ebin - bbin + 1, 0);	
	int i, pos;
	for (i = 0; (bbin + i) <= ebin; i++)
	{
		if (_distdistr[ri][bbin + i].size() > 0)
		{
			if (_distdistr[ri][bbin + i].size() > 1)
			{
				mcand.insert(mcand.end(), _distdistr[ri][bbin + i].begin() + 1, _distdistr[ri][bbin + i].end());
			}
			continue;
		}

		pos = _os_distdistr + sizeof(double) * 2 + sizeof(int) * (ri * numbin + bbin + i);
		if (pos != ftell(_ifp))
		{
			ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek offset of distance distribution");
		}
		ASSERT(fread((void*) &offset[i], sizeof(int), 1, _ifp) == 1, "could not read offset of distance distribution");
		if (offset[i] == 0)
		{
			_distdistr[ri][bbin + i].push_back(offset[i]);
		}
	}
	
	int numelem;
	int* mc = (int*) malloc(_numres * sizeof(int));
	for (i = 0; (bbin + i) <= ebin; i++)
	{
		if (offset[i] == 0)
		{
			continue;
		}

		pos = _os_distdistr + offset[i];
		if (pos != ftell(_ifp))
		{
			ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek distance distribution");
		}		
		ASSERT(fread((void*) &numelem, sizeof(int), 1, _ifp) == 1, "could not read number of elements in distance distribution");
		_distdistr[ri][bbin + i].push_back(numelem);
		ASSERT(fread((void*) mc, sizeof(int), numelem, _ifp) == numelem, "could not read distance distribution");
		_distdistr[ri][bbin + i].insert(_distdistr[ri][bbin + i].end(), mc, mc + numelem);
		mcand.insert(mcand.end(), mc, mc + numelem);
	}
	free((void*) mc);
}

void TargetStruct::readTargetStruct(const string & tsfname)
{
	readProteinStruct(tsfname);
	initDihedDistr();
	initDistDistr();
}

