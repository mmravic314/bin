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

#include "QueryStruct.h"

QueryStruct::QueryStruct(const string & qsfname)
{
	readQueryStruct(qsfname);
}

void QueryStruct::readCenRes()
{
	int offset;
	ASSERT(fseek(_ifp, sizeof(int) * SEC_CENRES, SEEK_SET) == 0, "could not seek offset of center residues section in file " + _fname);
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of center residues section in file " + _fname);
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek center residues section in file " + _fname);	
	int* cr = (int*) malloc(_numseg * sizeof(int));
	ASSERT(fread((void*) cr, sizeof(int), _numseg, _ifp) == _numseg, "could not read center residues section in file " + _fname);
	_cenres.assign(cr, cr + _numseg);
	free((void*) cr);
}

void QueryStruct::readCenResDihed()
{
	_crdihed.clear();
	int i, offset;
	vector<double> dihed;
	ASSERT(fseek(_ifp, sizeof(int) * SEC_CRDIHED, SEEK_SET) == 0, "could not seek offset of center residue dihedral angles section in file " + _fname);
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of center residue dihedral angles section in file " + _fname);
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek center residue dihedral angles section in file " + _fname);
	double* phipsi = (double*) malloc(_numseg * 2 * sizeof(double));
	ASSERT(fread((void*) phipsi, sizeof(double), (_numseg * 2), _ifp) == (_numseg * 2), "could not read center residue dihedral angles section in file " + _fname);
	for (i = 0; i < (_numseg * 2); i += 2)
	{
		dihed.assign(phipsi + i, phipsi + i + 2);
		_crdihed.push_back(dihed);
	}
	free((void*) phipsi);
}

void QueryStruct::readNumSeg()
{
	ASSERT(fseek(_ifp, sizeof(int) * SEC_NUMSEG, SEEK_SET) == 0, "could not seek number of segments in file " + _fname);	
	ASSERT(fread((void*) &_numseg, sizeof(int), 1, _ifp) == 1, "could not read number of segments in file " + _fname);
}

void QueryStruct::readQueryStruct(const string & qsfname)
{
	readProteinStruct(qsfname);
	readNumSeg();
	readCenRes();
	readCenResDihed();
	readResBefBrk();

#if defined(DEBUG_QS)
	int i;
	cout << "number of segments: " << _numseg << "\n";	
	cout << "center residues: ";
	for (i = 0; i < _cenres.size(); i++)
	{
		cout << " " << _cenres[i];
	}
	cout << "\n";
	cout << "center residue dihedral angles: ";
	for (i = 0; i < _crdihed.size(); i++)
	{
		printf(" (%f,%f)", _crdihed[i][0], _crdihed[i][1]);
	}
	cout << "\n";
	cout << "residues before breaks: ";
	for (i = 0; i < _befbrk.size(); i++)
	{
		cout << " " << _befbrk[i];
	}
	cout << "\n";
	exit(-1);
#endif
}

void QueryStruct::readResBefBrk()
{
	int offset;		
	ASSERT(fseek(_ifp, sizeof(int) * SEC_BEFBRK, SEEK_SET) == 0, "could not seek offset of residues before breaks section in file " + _fname);
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of residues before breaks section in file " + _fname);
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek residues before breaks section in file " + _fname);
	int* rbb = (int*) malloc((_numseg + 1) * sizeof(int));
	ASSERT(fread((void*) rbb, sizeof(int), (_numseg + 1), _ifp) == (_numseg + 1), "could not read indices before breaks section in file " + _fname);	
	_befbrk.assign(rbb, rbb + _numseg + 1);
	free((void*) rbb);
}

