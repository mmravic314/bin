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

#include "ProteinStruct.h"

ProteinStruct::ProteinStruct(const string & psfname)
{
	readProteinStruct(psfname);
}

double ProteinStruct::calcDistBB(const int & ri,const int & rj)
{
	if (ri == rj)
	{
		return 0.0;
	}

	int cai = ri * NUM_BBA + CA_IDX;
	int caj = rj * NUM_BBA + CA_IDX;
	return sqrt((_bbcoor[cai][0] - _bbcoor[caj][0]) * (_bbcoor[cai][0] - _bbcoor[caj][0])
		+ (_bbcoor[cai][1] - _bbcoor[caj][1]) * (_bbcoor[cai][1] - _bbcoor[caj][1])
		+ (_bbcoor[cai][2] - _bbcoor[caj][2]) * (_bbcoor[cai][2] - _bbcoor[caj][2]));
}

double ProteinStruct::calcDistCA(const int & ri, const int & rj)
{
	if (ri == rj)
	{
		return 0.0;
	}

	return sqrt((_cacoor[ri][0] - _cacoor[rj][0]) * (_cacoor[ri][0] - _cacoor[rj][0])
		+ (_cacoor[ri][1] - _cacoor[rj][1]) * (_cacoor[ri][1] - _cacoor[rj][1])
		+ (_cacoor[ri][2] - _cacoor[rj][2]) * (_cacoor[ri][2] - _cacoor[rj][2]));
} 

void ProteinStruct::closeProteinStructFile()
{
	if (_ifp != NULL)
	{
		fclose(_ifp);
	}
}

bool ProteinStruct::fullBB()
{
	for (int i = 0; i < _fullbb.size(); i++)
	{
		if (!_fullbb[i])
		{
			return false;
		}
	}
	return true;
}

void ProteinStruct::initPdbInfo()
{
	_pdbinfo.assign(_numres, string());
}

void ProteinStruct::readBBCoords()
{
	int i, j, offset, pos;
	const int coordlen = _numres * NUM_BBA * NUM_COORDS;
	double* xyz = (double*) malloc(coordlen * sizeof(double));
	_bbcoor.assign(_numres * NUM_BBA, vector<double>());
	_fullbb.assign(_numres, true);

	ASSERT(fseek(_ifp, sizeof(int) * SEC_BBCOOR, SEEK_SET) == 0, "could not seek offset of BB coords section in file " + _fname);
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of BB coords section in file " + _fname);
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek BB coords section in file " + _fname);
	ASSERT(fread((void*) xyz, sizeof(double), coordlen, _ifp) == coordlen, "could not read BB coords section in file " + _fname);
	for (i = 0; i < _numres; i++)
	{
		for (j = 0; j < NUM_BBA; j++)
		{
			pos = i * NUM_BBA * NUM_COORDS + j * NUM_COORDS;
			_bbcoor[i * NUM_BBA + j].assign(xyz + pos, xyz + pos + NUM_COORDS);

			if (_fullbb[i] && (IMPOSSIBLE_COORD == _bbcoor[i * NUM_BBA + j][0]))
			{
				_fullbb[i] = false;
			}
		}
	}
	
	free((void*) xyz);

#if defined(DEBUG_BB)
	cout << "BB coords:\n";
	for (i = 0; i < _numres; i++)
	{
		printf("N\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 0][0], _bbcoor[4 * i + 0][1], _bbcoor[4 * i + 0][2]);
		printf("CA\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 1][0], _bbcoor[4 * i + 1][1], _bbcoor[4 * i + 1][2]);
		printf("C\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 2][0], _bbcoor[4 * i + 2][1], _bbcoor[4 * i + 2][2]);
		printf("O\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 3][0], _bbcoor[4 * i + 3][1], _bbcoor[4 * i + 3][2]);
		cout << _fullbb[i] << "\n";
		cout << "\n";
	}
	exit(-1);
#endif

}

void ProteinStruct::readCACoords()
{
	int i, offset, pos;
	const int coordlen = _numres * NUM_COORDS;
	double* xyz = (double*) malloc(coordlen * sizeof(double));
	_cacoor.assign(_numres, vector<double>());
	
	ASSERT(fseek(_ifp, sizeof(int) * SEC_CACOOR, SEEK_SET) == 0, "could not seek offset of CA coords section in file " + _fname);
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of CA coords section in file " + _fname);
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek CA coords section in file " + _fname);
	ASSERT(fread((void*) xyz, sizeof(double), coordlen, _ifp) == coordlen, "could not read CA coords section in file " + _fname);
	for (i = 0; i < _numres; i++)
	{
		pos = i * NUM_COORDS;
		_cacoor[i].assign(xyz + pos, xyz + pos + NUM_COORDS);
	}

	free((void*) xyz);
}

void ProteinStruct::readNumRes()
{
	ASSERT(fseek(_ifp, sizeof(int) * SEC_NUMRES, SEEK_SET) == 0, "could not seek number of residues in file " + _fname);	
	ASSERT(fread((void*) &_numres, sizeof(int), 1, _ifp) == 1, "could not read number of residues in file " + _fname);
}

void ProteinStruct::readPdbInfo()
{
	_pdbinfo.clear();
	
	int i, offset, os_pdbinfo, pos, reslen, maxreslen = 10000; // guess as to what the max possible residue string will be
	char* resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
	
	ASSERT(fseek(_ifp, sizeof(int) * SEC_PDBINFO, SEEK_SET) == 0, "could not seek offset of PDB info section in file " + _fname);
	ASSERT(fread((void*) &os_pdbinfo, sizeof(int), 1, _ifp) == 1, "could not read offset of PDB info section in file " + _fname);
	ASSERT(fseek(_ifp, os_pdbinfo, SEEK_SET) == 0, "could not seek PDB info section in file " + _fname);
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of PDB info in file " + _fname);
	pos = os_pdbinfo + offset;
	ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek PDB info in file " + _fname);
	for (i = 0; i < _numres; i++)
	{
		ASSERT(fread((void*) &reslen, sizeof(int), 1, _ifp) == 1, "could not read length of PDB info in file " + _fname);
		if (reslen > maxreslen)
		{
			maxreslen = reslen;
			free((void*) resinfo);
			resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
		}
		ASSERT(fread((void*) resinfo, sizeof(char), reslen, _ifp) == reslen, "could not read PDB info in file " + _fname);
		resinfo[reslen] = '\0';
		_pdbinfo.push_back(string(resinfo));
	}
	free((void*) resinfo);
}

void ProteinStruct::readPdbInfo(const vector < int > & residx)
{
	int i;
	vector<int> ri;
	// check which residues need to actually be read
	for (i = 0; i < residx.size(); i++)
	{
		if (_pdbinfo[residx[i]].empty())
		{
			ri.push_back(residx[i]);
		}
	}
	if (ri.empty())
	{
		return;
	}

	vector<int> offset(ri.size());
	int os_pdbinfo, pos, reslen, maxreslen = 10000; // guess as to what the max possible residue string will be
	char* resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'

	ASSERT(fseek(_ifp, sizeof(int) * SEC_PDBINFO, SEEK_SET) == 0, "could not seek offset of PDB info section in file " + _fname);
	ASSERT(fread((void*) &os_pdbinfo, sizeof(int), 1, _ifp) == 1, "could not read offset of PDB info section in file " + _fname);

	for (i = 0; i < ri.size(); i++)
	{
		pos = os_pdbinfo + sizeof(int) * ri[i];
		if (pos != ftell(_ifp))
		{
			ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek offset of PDB info in file " + _fname);
		}
		ASSERT(fread((void*) &offset[i], sizeof(int), 1, _ifp) == 1, "could not read offset of PDB info in file " + _fname);
	}
	for (i = 0; i < ri.size(); i++)
	{
		pos = os_pdbinfo + offset[i];
		if (pos != ftell(_ifp))
		{
			ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek PDB info in file " + _fname);
		}
		ASSERT(fread((void*) &reslen, sizeof(int), 1, _ifp) == 1, "could not read length of PDB info in file " + _fname);
		if (reslen > maxreslen)
		{
			maxreslen = reslen;
			free((void*) resinfo);
			resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
		}
		ASSERT(fread((void*) resinfo, sizeof(char), reslen, _ifp) == reslen, "could not read PDB info in file " + _fname);
		resinfo[reslen] = '\0';
		_pdbinfo[ri[i]] = string(resinfo);
	}
	free((void*) resinfo);
}

void ProteinStruct::readProteinStruct(const string & psfname)
{
	_fname = psfname;

	openFileC(_ifp, _fname, "rb");

	readNumRes();

#if defined(DEBUG_PS)
	cout << "file name: " << _fname << "\n";
	cout << "number of residues: " << _numres << "\n";
	exit(-1);
#endif
}

void ProteinStruct::readSeq(char * & seq)
{
	delete [] seq;

	int offset;
	const int seqlen = LEN_AA_CODE * _numres;
	seq = new char [seqlen];

	ASSERT(fseek(_ifp, sizeof(int) * SEC_SEQ, SEEK_SET) == 0, "could not seek offset of sequence section in file " + _fname);
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of sequence section in file " + _fname);
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek sequence section in file " + _fname);
	ASSERT(fread((void*) seq, sizeof(char), seqlen, _ifp) == seqlen, "could not read sequence section in file " + _fname);
}
