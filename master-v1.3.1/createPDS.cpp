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

#undef DEBUG_BB
#undef DEBUG_CA
#undef DEBUG_CMDLINE
#undef DEBUG_CONTACT
#undef DEBUG_CORCHA
#undef DEBUG_DIST
#undef DEBUG_NR
#undef DEBUG_PARSEPDB
#undef DEBUG_PDBINFO
#undef DEBUG_PS
#undef DEBUG_QS
#undef DEBUG_SEQALI

#define IMPOSSIBLE_ANGLE	999.9
#define IMPOSSIBLE_COORD	999999.999

#define NUM_BBA				4
#define NUM_COORDS			3

#define N_IDX				0
#define CA_IDX				1
#define C_IDX				2
#define O_IDX				3

#define LEN_AA_CODE			3

// Protein Structure
#define SEC_BBCOOR			0
#define SEC_CACOOR			1
#define SEC_NUMRES			2
#define SEC_PDBINFO			3
#define SEC_SEQ				4

#define NUM_SEC_PS			5

// Query Structure
#define SEC_BEFBRK			5
#define SEC_CENRES			6
#define SEC_CRDIHED			7
#define SEC_NUMSEG			8

#define NUM_SEC_QS			9

// Target Structure
#define SEC_DIHEDDISTR		5
#define SEC_DISTDISTR		6

#define NUM_SEC_TS			7

// STL Includes
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <getopt.h>
#include <unistd.h>

// MSL Includes
#include "AtomContainer.h"
#include "AtomSelection.h"
#include "MslTools.h"
#include "System.h"

using namespace MSL;
using namespace std;

class CorrChains;
class CreateOptions;
class ProteinStruct;
class QueryStruct;
class TargetStruct;

bool alignChains(const vector<string> &, const vector<string> &, const double &, AtomPointerVector &, AtomPointerVector &, const int &, const double &, vector<double> &, vector<vector<double> > &, const double &, const vector<double> &, const vector<vector<double> > &);
bool alignSeqNeedlemanWunsch(const vector<string> &, const vector<string> &, const double &, vector<string> &, vector<string> &);
bool areUnitsInContact(AtomPointerVector*, AtomPointerVector*);
int calcNumOvlpChains(Chain*, vector<Chain*> &, map<Chain*, vector<Chain*> > &);
double calcRmsdKabsch(AtomPointerVector &, AtomPointerVector &, int, vector<double> &, vector<vector<double> > &);
void error(const string &);
void file2array(string, vector<string> &);
string fileBase(string);
int findMostOvlpChainSet(const int &, const int &, vector<vector<Chain*> > &, map<Chain*, vector<Chain*> > &, vector<Chain*> &);
map<Chain*, AtomPointerVector> getChainCa(System &);
map<Chain*, vector<double> > getChainCaCentroid(map<Chain*, AtomPointerVector> &);
map<Chain*, vector<Chain*> > getChainContacts(System &);
map<Chain*, vector<string> > getChainSeq(System &);
vector<CorrChains> getCorrChains(vector<Chain*> &, vector<Chain*> &, map<Chain*, vector<double> > &, const vector<double> &, const vector<vector<double> > &);
string getFileName(string);
void gridPoint(Atom* a, CartesianPoint& c, double gs, int* ip=NULL, int* jp=NULL, int* kp=NULL);
bool hasBreak(Residue &, Residue &);
bool isProtein(const vector<string> &, const string &);
void mergeChainGroups(Chain*, Chain*, map<Chain*, int> &, vector<vector<Chain*> > &);
void openFileC (FILE* &, const char*, const char*);
void openFileCPP(fstream &, const string &, const ios_base::openmode &);
string optionUsage(string, string, int, int, int);
string pad(string, int);
void parseCommandLine(int, char**, CreateOptions &);
vector<string> tokenize(const std::string & _input, const std::string & _delimiter=" ", bool _allowEmtpy=false);

template <class T>
string toString (const T & val)
{
	return static_cast<ostringstream*>( &(ostringstream() << val) )->str();
}

char* trim(char*);
void usage();
template <class T> void writeDatum (fstream &, T, const bool &);
void writeString (fstream &, const string &, const bool &);

#define ASSERT(cond, msg) {if (!(cond)) error(msg);};

class CorrChains
{
	public:
		Chain* chnA;
		Chain* chnB;
		double centDist;

		CorrChains() {}
		~CorrChains() {}
};

struct centDistComp
{
	bool operator() (const CorrChains & ccA, const CorrChains & ccB)
	{
		if (ccA.centDist != ccB.centDist)
		{
			return (ccA.centDist < ccB.centDist);
		}
		else
		{
			return ((&ccA) < (&ccB));
		}
	}
};

class CreateOptions
{
	public:
		CreateOptions();
		~CreateOptions() {}

		bool getBinary() const { return _bin; }
		double getContRmsdThresh() const { return _contRmsdThresh; }
		double getDistCut() const { return _dcut; }
		double getDistStep() const { return _dstep; }
		string getFileExt() const { return _ext; }
		vector<string>& getLegalAA() { return _legalAA; }
		bool getNonRed() const { return _nr; }
		bool getNonRedPdb() const { return _nrPdb; }
		double getOptRmsdThresh() const { return _optRmsdThresh; }
		vector<string>& getPdbFiles() { return _pdbfnames; }
		vector<string>& getPdsFiles() { return _pdsfnames; }
		string getPdsType() const { return _type; }
		double getPhiStep() const { return _phistep; }
		vector<string>& getPostPdbFiles() { return _opdbfnames; }
		double getPsiStep() const { return _psistep; }		
		double getSeqIdenThresh() const { return _seqIdenThresh; }
		string getWordTer() const { return _ter; }
		string getWordSep() const { return _sep; }
				
		void setContRmsdThresh(const char*);
		void setDistCut(const char*);
		void setDistStep(const char*);
		void setNaturalAA();
		void setNonRed(const bool & f) { _nr = f; }
		void setNonRedPdb(const bool & f) { _nrPdb = f; }
		void setOptRmsdThresh(const char*);
		void setPdbFile(const string &);
		void setPdbFiles(const string &);
		void setPdsFile(const string & fn, const bool & chk = true);
		void setPdsFiles(const string &);
		void setPdsType(const string &);
		void setPhiStep(const char*);
		void setPostPdbFile(const string &);
		void setPostPdbFiles(const string &);
		void setPsiStep(const char*);	
		void setSeqIdenThresh(const char*);
		void setUnnaturalAA();

	protected:
		bool _bin;
		double _contRmsdThresh; // --gRMSD
		double _dcut;
		double _dstep;
		string _ext;
		vector<string> _legalAA;
		bool _nr;
		bool _nrPdb;
		vector<string> _opdbfnames;
		double _optRmsdThresh; // --lRMSD
		vector<string> _pdbfnames;
		vector<string> _pdsfnames;
		double _phistep;
		double _psistep;
		string _sep;
		double _seqIdenThresh; // --seqID
		string _ter;
		string _type;
};

CreateOptions::CreateOptions()
{
	_bin = true;
	_contRmsdThresh = 2.0; //--gRMSD
	_dcut = 25.0;
	_dstep = 5.0;
	_ext = ".pds";
	_optRmsdThresh = 1.0; //--lRMSD
	_phistep = 10.0;
	_psistep = 10.0;
	_sep = "";
	_seqIdenThresh = 0.9; // --seqID
	_nr = false; // --nr
	_nrPdb = false; // --nrPDB
	_ter = "";

	setNaturalAA();
}

void CreateOptions::setContRmsdThresh(const char* s)
{
	if (!((sscanf(s, "%lf", &(_contRmsdThresh)) == 1) && (_contRmsdThresh >= 0.0)))
	{
		usage();
		error("bad contact RMSD threshold value");
	}
}

void CreateOptions::setDistCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_dcut)) == 1) && (_dcut > 0.0)))
	{
		usage();
		error("bad distance cutoff value");
	}
}

void CreateOptions::setDistStep(const char * s)
{
	if (!((sscanf(s, "%lf", &(_dstep)) == 1) && (_dstep > 0.0)))
	{
		usage();
		error("bad distance step value");
	}
}

void CreateOptions::setNaturalAA()
{
	// legal residue names that are considered "protein" here
	_legalAA.clear();
	_legalAA.push_back("ALA");
	_legalAA.push_back("ARG");
	_legalAA.push_back("ASN");
	_legalAA.push_back("ASP");
	_legalAA.push_back("CYS");
	_legalAA.push_back("GLN");
	_legalAA.push_back("GLU");
	_legalAA.push_back("GLY");
	_legalAA.push_back("HIS"); _legalAA.push_back("HSC"); _legalAA.push_back("HSD"); _legalAA.push_back("HSE"); _legalAA.push_back("HSP");
	_legalAA.push_back("ILE");
	_legalAA.push_back("LEU");
	_legalAA.push_back("LYS");
	_legalAA.push_back("MET"); _legalAA.push_back("MSE");
	_legalAA.push_back("PHE");
  	_legalAA.push_back("PRO");
	_legalAA.push_back("SER");
	_legalAA.push_back("THR");
	_legalAA.push_back("TRP");
	_legalAA.push_back("TYR");
	_legalAA.push_back("VAL");
}

void CreateOptions::setOptRmsdThresh(const char * s)
{
	if (!((sscanf(s, "%lf", &(_optRmsdThresh)) == 1) && (_optRmsdThresh >= 0.0)))
	{
		usage();
		error("bad optimal RMSD threshold value");
	}
}

void CreateOptions::setPdbFile(const string & fn)
{
	if (_pdbfnames.size() != 0)
	{
		usage();
		error("input PDB files already exist");
	}
	_pdbfnames.push_back(fn);
}

void CreateOptions::setPdbFiles(const string & list)
{
	if (_pdbfnames.size() != 0)
	{
		usage();
		error("input PDB files already exist");
	}
	file2array(list, _pdbfnames);
}

void CreateOptions::setPdsFile(const string & fn, const bool & chk)
{
	if (chk)
	{
		if (_pdsfnames.size() != 0)
		{
			usage();
			error("output PDS files already exist");
		}
	}
	_pdsfnames.push_back(fn);
}

void CreateOptions::setPdsFiles(const string & list)
{
	if (_pdsfnames.size() != 0)
	{
		usage();
		error("output PDS files already exist");
	}
	file2array(list, _pdsfnames);
}

void CreateOptions::setPdsType(const string & t)
{
	_type = t;
	if (!((_type == "query") || (_type == "target")))
	{
		usage();
		error("bad output PDS type value");
	}
}

void CreateOptions::setPhiStep(const char * s)
{
	if (!((sscanf(s, "%lf", &(_phistep)) == 1) && (_phistep > 0.0)))
	{
		usage();
		error("bad phi step value");
	}
}

void CreateOptions::setPostPdbFile(const string & fn)
{
	if (_opdbfnames.size() != 0)
	{
		usage();
		error("output post-processed PDB files already exist");
	}
	_opdbfnames.push_back(fn);
}

void CreateOptions::setPostPdbFiles(const string & list)
{
	if (_opdbfnames.size() != 0)
	{
		usage();
		error("output post-processed PDB files already exist");
	}
	file2array(list, _opdbfnames);
}

void CreateOptions::setPsiStep(const char * s)
{
	if (!((sscanf(s, "%lf", &(_psistep)) == 1) && (_psistep > 0.0)))
	{
		usage();
		error("bad psi step value");
	}
}

void CreateOptions::setSeqIdenThresh(const char * s)
{
	if (!((sscanf(s, "%lf", &(_seqIdenThresh)) == 1) && (_seqIdenThresh >= 0.0) && (_seqIdenThresh <= 1.0)))
	{
		usage();
		error("bad sequence identity threshold value");
	}
}

void CreateOptions::setUnnaturalAA()
{
	// unnatural AA's
	_legalAA.push_back("CSO");
	_legalAA.push_back("HIP");
	_legalAA.push_back("PTR");
	_legalAA.push_back("SEC");
	_legalAA.push_back("SEP");
	_legalAA.push_back("TPO");
}

class ProteinStruct
{
	public:
		ProteinStruct() {}
		ProteinStruct(const string &, const vector<string> &);
		ProteinStruct(const string &, const vector<string> &, const bool &, const double &, const double &, const double &, const bool &);
		~ProteinStruct() {}		
	
		void checkBBCoords(FILE* &, const string &, const CreateOptions &);
		void checkCACoords(FILE* &, const string &, const CreateOptions &);
		void checkProteinStructFile(FILE* &, const string &, const CreateOptions &);
		System& getProteinSys() { return _sys; }
		void parsePdb(const string &, const vector<string> &);
		void parsePdb(const string &, const vector<string> &, const bool &, const double &, const double &, const double &, const bool &);
		void setBBCoords();
		void setCACoords();
		void setDihed();
		void setDist();
		void setFileName(const string &);
		void setNumRes();
		void setPdbInfo();
		void setProteinStruct();
		void setSeq();
		void writeBBCoords(fstream &, const CreateOptions &);
		void writeCACoords(fstream &, const CreateOptions &);
//			void writeDihedralAngles(fstream &, const CreateOptions &);
//			void writeDistance(fstream &, const CreateOptions &);
		void writePdbInfo(fstream &, const CreateOptions &);
		void writeProteinStruct(fstream &, const CreateOptions &);
		void writeProteinStructFile(fstream &, const CreateOptions &);
		int writeProteinStructFileHeader(fstream &, const CreateOptions &, const int &);
		void writeSeq(fstream &, const CreateOptions &);

	protected:
		vector<vector<double> > _bbcoor;
		vector<vector<double> > _cacoor;
		vector<pair<double, double> > _dihed;
		map<pair<int, int>, double> _dist;
		string _fname;
		int _numres;		
		vector<vector<string> > _pdbinfo;
		AtomPointerVector _pselca;
		vector<string> _seq;
		System _sys;
};

ProteinStruct::ProteinStruct(const string & pdbfname, const vector<string> & legalAA)
{
	parsePdb(pdbfname, legalAA);
	setProteinStruct();
}

ProteinStruct::ProteinStruct(const string & pdbfname, const vector<string> & legalAA, const bool & nr, const double & seqIdenThresh, const double & optRmsdThresh, const double & contRmsdThresh, const bool & nrPdb)
{
	parsePdb(pdbfname, legalAA, nr, seqIdenThresh, optRmsdThresh, contRmsdThresh, nrPdb);
	setProteinStruct();
}

void ProteinStruct::checkBBCoords(FILE * & ifp, const string & psfname, const CreateOptions & copts)
{
	int i, j, numres, coordlen, offset, pos;
	double *xyz = NULL;
	vector<vector<double> > bbcoor;

	if (copts.getBinary())
	{
		ASSERT(fseek(ifp, sizeof(int) * SEC_NUMRES, SEEK_SET) == 0, "could not seek number of residues in file " + psfname);
		ASSERT(fread((void*) &numres, sizeof(int), 1, ifp) == 1, "could not read number of residues in file " + psfname);
		ASSERT(numres == _numres, "number of residues not consistent in file " + psfname);

		bbcoor.assign(numres * NUM_BBA, vector<double>());
		coordlen = numres * NUM_BBA * NUM_COORDS;
		xyz = (double*) malloc(coordlen * sizeof(double));
		ASSERT(fseek(ifp, sizeof(int) * SEC_BBCOOR, SEEK_SET) == 0, "could not seek offset of BB coords section in file " + psfname);
		ASSERT(fread((void*) &offset, sizeof(int), 1, ifp) == 1, "could not read offset of BB coords section in file " + psfname);
		ASSERT(fseek(ifp, offset, SEEK_SET) == 0, "could not seek BB coords section in file " + psfname);
		ASSERT(fread((void*) xyz, sizeof(double), coordlen, ifp) == coordlen, "could not read BB coords section in file " + psfname);
		for (i = 0; i < numres; i++)
		{
			for (j = 0; j < NUM_BBA; j++)
			{
				pos = i * NUM_BBA * NUM_COORDS + j * NUM_COORDS;
				bbcoor[i * NUM_BBA + j].assign(xyz + pos, xyz + pos + NUM_COORDS);
			}
		}

		ASSERT(bbcoor.size() == _bbcoor.size(), "number of BB coords not consistent in file " + psfname);
		for (i = 0; i < bbcoor.size(); i++)
		{
			ASSERT((_bbcoor[i][0] == bbcoor[i][0]) && (_bbcoor[i][1] == bbcoor[i][1]) && (_bbcoor[i][2] == bbcoor[i][2]), "BB coords section not consistent in file " + psfname);
		}

		free((void*) xyz);
	}
	else
	{
		error("could not check the correctness of BB coords in text file " + psfname);
	}
	
	cout << "BB coords section in file " << psfname << " checked." << endl;
}

void ProteinStruct::checkCACoords(FILE * & ifp, const string & psfname, const CreateOptions & copts)
{
	int i, numres, coordlen, offset, pos;
	double *xyz = NULL;
	vector<vector<double> > cacoor;

	if (copts.getBinary())
	{
		ASSERT(fseek(ifp, sizeof(int) * SEC_NUMRES, SEEK_SET) == 0, "could not seek number of residues in file " + psfname);
		ASSERT(fread((void*) &numres, sizeof(int), 1, ifp) == 1, "could not read number of residues in file " + psfname);
		ASSERT(numres == _numres, "number of residues not consistent in file " + psfname);

		cacoor.assign(numres, vector<double>());
		coordlen = numres * NUM_COORDS;
		xyz = (double*) malloc(coordlen * sizeof(double));
		ASSERT(fseek(ifp, sizeof(int) * SEC_CACOOR, SEEK_SET) == 0, "could not seek offset of CA coords section in file " + psfname);
		ASSERT(fread((void*) &offset, sizeof(int), 1, ifp) == 1, "could not read offset of CA coords section in file " + psfname);
		ASSERT(fseek(ifp, offset, SEEK_SET) == 0, "could not seek CA coords section in file " + psfname);
		ASSERT(fread((void*) xyz, sizeof(double), coordlen, ifp) == coordlen, "could not read CA coords section in file " + psfname);
		for (i = 0; i < numres; i++)
		{
			pos = i * NUM_COORDS;
			cacoor[i].assign(xyz + pos, xyz + pos + NUM_COORDS);
		}

		ASSERT(cacoor.size() == _cacoor.size(), "number of CA coords not consistent in file " + psfname);
		for (i = 0; i < cacoor.size(); i++)
		{
			ASSERT((_cacoor[i][0] == cacoor[i][0]) && (_cacoor[i][1] == cacoor[i][1]) && (_cacoor[i][2] == cacoor[i][2]), "CA coords section not consistent in file " + psfname);
		}

		free((void*) xyz);
	}
	else
	{
		error("could not check the correctness of CA coords in text file " + psfname);
	}
	
	cout << "CA coords section in file " << psfname << " checked." << endl;
}

void ProteinStruct::checkProteinStructFile(FILE* & ifp, const string & psfname, const CreateOptions & copts)
{
	checkBBCoords(ifp, psfname, copts);
	checkCACoords(ifp, psfname, copts);
}

void ProteinStruct::parsePdb(const string & pdbfname, const vector<string> & legalAA)
{
	_sys.reset();
	_pselca.deletePointers();

	setFileName(pdbfname);

	System S;
	int i;
	
	ASSERT(S.readPdb(pdbfname, true), "could not read file " + pdbfname);
	for (i = 0; i < S.positionSize(); i++)
	{
		Residue& r = S.getResidue(i);

		if (!isProtein(legalAA, r.getResidueName()))
		{
			fprintf(stdout, "Warning: residue '%s,%d%s' has been skipped due to unidentified AA type '%s'.\n",
				r.getChainId().c_str(), r.getResidueNumber(), r.getResidueIcode().c_str(), r.getResidueName().c_str());
			continue;
		}

		if (!r.atomExists("CA"))
		{
			fprintf(stdout, "Warning: residue '%s,%d%s' has been skipped due to missing CA atom.\n",
				r.getChainId().c_str(), r.getResidueNumber(), r.getResidueIcode().c_str());
			continue;
		}

		_sys.addAtoms(r.getAtomPointers(), true);
	}

	AtomSelection sel(_sys.getAtomPointers());
	_pselca = sel.select("name CA");
	ASSERT(_pselca.size() == _sys.positionSize(), "numbers of CA atoms and residues not consistent in file " + pdbfname);
	ASSERT(_pselca.size() > 0, "no valid residues in file " + pdbfname);

#if defined(DEBUG_PARSEPDB)
	cout << "file name: " << _fname << endl;
	cout << "full atoms:" << endl;
	cout << _sys.getAtomPointers() << endl;
	cout << "CA atoms:" << endl;
	cout << _pselca << endl;
	exit(-1);
#endif
}

void ProteinStruct::parsePdb(const string & pdbfname,const vector < string > & legalAA,const bool & nr,const double & seqIdenThresh,const double & optRmsdThresh,const double & contRmsdThresh,const bool & nrPdb)
{
	_sys.reset();
	_pselca.deletePointers();

	setFileName(pdbfname);
	
	System S0, S;
	int i;

	ASSERT(S0.readPdb(pdbfname, true), "could not read file " + pdbfname);
	for (i = 0; i < S0.positionSize(); i++)
	{
		Residue& r = S0.getResidue(i);

		if (!isProtein(legalAA, r.getResidueName()))
		{
			fprintf(stdout, "Warning: residue '%s,%d%s' has been skipped due to unidentified AA type '%s'.\n",
				r.getChainId().c_str(), r.getResidueNumber(), r.getResidueIcode().c_str(), r.getResidueName().c_str());
			continue;
		}

		if (!r.atomExists("CA"))
		{
			fprintf(stdout, "Warning: residue '%s,%d%s' has been skipped due to missing CA atom.\n",
				r.getChainId().c_str(), r.getResidueNumber(), r.getResidueIcode().c_str());
			continue;
		}

		S.addAtoms(r.getAtomPointers(), true);
	}

	if (nr)
	{
		cout << "Removing redundant chains...\n";
		
		map<Chain*, AtomPointerVector> chainCa = getChainCa(S);
		map<Chain*, vector<double> > chainCaCentroid = getChainCaCentroid(chainCa);
		map<Chain*, vector<Chain*> > chainContacts = getChainContacts(S);
		map<Chain*, int> chainGrpIdx;
		vector<vector<Chain*> > chainGrps;
		vector<vector<Chain*> > chainGrpsRefined;
		map<Chain*, bool> chainKept;
		map<Chain*, vector<string> > chainSeq = getChainSeq(S);
		map<Chain*, bool> chainRemoved;
		for (i = 0; i < S.chainSize(); i++)
		{
			chainRemoved[&(S.getChain(i))] = false; // True: can be removed safely
													// False: can't be removed
		}
		Chain* chn;
		vector<CorrChains> coChains;
		vector<Chain*> curChainSet, mostOvlpChainSet;
		bool f;
		map<Chain*, vector<Chain*> >::iterator it, itA, itB;
		int j, k, minChainIdx, minGrpIdx, minNumChains, ttlNumChains;
		vector<double> t(3), ct(3);
		vector<vector<double> > u(3, vector<double>(3)), cu(3, vector<double>(3));
		
		for (itA = chainContacts.begin(); itA != chainContacts.end(); ++itA)
		{		
			if (chainGrpIdx.find(itA->first) == chainGrpIdx.end())
			{
				chainGrps.push_back(vector<Chain*>());
				chainGrps[chainGrps.size() - 1].push_back(itA->first);
				chainGrpIdx[itA->first] = chainGrps.size() - 1; 		
			}

			it = itA;
			++it;
			for (itB = it; itB != chainContacts.end(); ++itB)
			{		
				if ((chainGrpIdx.find(itA->first) != chainGrpIdx.end())
					&& (chainGrpIdx.find(itB->first) != chainGrpIdx.end())
					&& (chainGrpIdx[itA->first] == chainGrpIdx[itB->first]))
				{
					continue;
				}
				
				if (alignChains(chainSeq[itA->first], chainSeq[itB->first], seqIdenThresh, chainCa[itA->first], chainCa[itB->first], 1, optRmsdThresh, ct, cu, contRmsdThresh, t, u))
				{
					// two central chains are redundant
					coChains = getCorrChains(itA->second, itB->second, chainCaCentroid, ct, cu);
					f = true;
					for (i = 0; i < coChains.size(); i++)
					{
						if (!alignChains(chainSeq[coChains[i].chnA], chainSeq[coChains[i].chnB], seqIdenThresh, chainCa[coChains[i].chnA], chainCa[coChains[i].chnB], 0, optRmsdThresh, t, u, contRmsdThresh, ct, cu))
						{
							// two corresponding neighbor chains are not redundant
							f = false;
							break;
						}
					}
					if (f)
					{
						// chain A plus its neighbors and chain B plus its neighbors are redundant
						if (itA->second.size() < itB->second.size())
						{
							// chain B plus its neighbors contains chain A plus its neighbors
							// therefore safely remove chain A and its neighbors					
							chainRemoved[itA->first] = true;
						}
						else
						{
							if (itA->second.size() > itB->second.size())
							{
								// chain A plus its neighbors contains chain B plus its neighbors
								// therefore safely remove chain B and its neighbors
								chainRemoved[itB->first] = true;
							}
							else
							{
								// chain A plus its neighbors == chain B plus its neighbors
								if ((chainGrpIdx.find(itA->first) != chainGrpIdx.end())
									&& (chainGrpIdx.find(itB->first) != chainGrpIdx.end()))
								{
									// merge two different chain groups
									mergeChainGroups(itA->first, itB->first, chainGrpIdx, chainGrps);
								}
								else
								{
									if ((chainGrpIdx.find(itA->first) != chainGrpIdx.end())
										&& (chainGrpIdx.find(itB->first) == chainGrpIdx.end()))
									{
										// add chain B to chain A's group
										chainGrps[chainGrpIdx[itA->first]].push_back(itB->first);
										chainGrpIdx[itB->first] = chainGrpIdx[itA->first];
									}
									else
									{
										if ((chainGrpIdx.find(itA->first) == chainGrpIdx.end())
											&& (chainGrpIdx.find(itB->first) != chainGrpIdx.end()))
										{
											// add chain A to chain B's group
											chainGrps[chainGrpIdx[itB->first]].push_back(itA->first);
											chainGrpIdx[itA->first] = chainGrpIdx[itB->first];
										}
										else
										{
											if ((chainGrpIdx.find(itA->first) == chainGrpIdx.end())
												&& (chainGrpIdx.find(itB->first) == chainGrpIdx.end()))
											{
												// assign a new group for chain A and B
												chainGrps.push_back(vector<Chain*>());											
												chainGrps[chainGrps.size() - 1].push_back(itA->first);
												chainGrps[chainGrps.size() - 1].push_back(itB->first);
												chainGrpIdx[itA->first] = chainGrps.size() - 1;
												chainGrpIdx[itB->first] = chainGrps.size() - 1; 										
											}
										}
									}
								}
							}
						}					
					}
				}
				// 1. two central chains are not redundant
				// 2. two corresponding neighbor chains are not redundant
				// 3. chain A plus its neighbors contains chain B plus its neighbors, or vice versa 		
				if (chainGrpIdx.find(itB->first) == chainGrpIdx.end())
				{
					chainGrps.push_back(vector<Chain*>());
					chainGrps[chainGrps.size() - 1].push_back(itB->first);
					chainGrpIdx[itB->first] = chainGrps.size() - 1; 			
				}
			}		
		}

		for (i = 0; i < chainGrps.size(); i++)
		{
			if (chainGrps[i].size() <= 0)
			{
				continue;
			}
			
			f = false;		
			for (j = 0; j < chainGrps[i].size(); j++)
			{
				if (chainRemoved[chainGrps[i][j]])
				{
					f = true;
					break;
				}
			}
			if (f)
			{
				for (j = 0; j < chainGrps[i].size(); j++)
				{
					chainRemoved[chainGrps[i][j]] = true;

#if defined(DEBUG_NR)
					chn = chainGrps[i][j];
					AtomContainer ac;
					ac.addAtoms(chn->getAtomPointers());
					for (k = 0; k < chainContacts[chn].size(); k++)
					{
						ac.addAtoms(chainContacts[chn][k]->getAtomPointers());
					}
					ac.writePdb(getFileName(pdbfname) + "." + toString(i + 1) + "." + toString(j + 1) + ".r.pdb");
#endif				
				}
				chainGrps[i].clear();
			}
		}

		minNumChains = S.chainSize() + 1;
		minGrpIdx = -1;
		for (i = 0; i < chainGrps.size(); i++)
		{
			if (chainGrps[i].size() <= 0)
			{
				continue;
			}

#if defined(DEBUG_NR)
			for (j = 0; j < chainGrps[i].size(); j++)
			{
				chn = chainGrps[i][j];
				AtomContainer ac;
				ac.addAtoms(chn->getAtomPointers());
				for (k = 0; k < chainContacts[chn].size(); k++)
				{
					ac.addAtoms(chainContacts[chn][k]->getAtomPointers());
				}
				ac.writePdb(getFileName(pdbfname) + "." + toString(i + 1) + "." + toString(j + 1) + ".k.pdb");			
			}
#endif
			
			chainGrpsRefined.push_back(chainGrps[i]);
			if (chainGrpsRefined[chainGrpsRefined.size() - 1].size() < minNumChains)
			{
				minNumChains = chainGrpsRefined[chainGrpsRefined.size() - 1].size();
				minGrpIdx = chainGrpsRefined.size() - 1;
			}
		}
		ASSERT((minNumChains > 0) && (minNumChains <= S.chainSize())
			&& (minGrpIdx >= 0) && (minGrpIdx < chainGrpsRefined.size()),
			"could not find group with minimum number of chains");

#if defined(DEBUG_NR)
		printf("group %d contains %d chains, which is minimum.\n", minGrpIdx + 1, minNumChains);
		for (i = 0; i < chainGrpsRefined.size(); i++)
		{
			printf("group %d contains %d chains.\n", i + 1, chainGrpsRefined[i].size());			
		}	
		exit(-1);
#endif

		minNumChains = S.chainSize() + 1;
		minChainIdx = -1;
		for (i = 0; i < chainGrpsRefined[minGrpIdx].size(); i++)
		{
			ttlNumChains = findMostOvlpChainSet(minGrpIdx, i, chainGrpsRefined, chainContacts, curChainSet);
			if (ttlNumChains < minNumChains)
			{
				minNumChains = ttlNumChains;
				minChainIdx = i;
				mostOvlpChainSet = curChainSet;
			}
		}
		ASSERT((minNumChains > 0) && (minNumChains <= S.chainSize())
			&& (minChainIdx >= 0) && (minChainIdx < chainGrpsRefined[minGrpIdx].size()),
			"could not find most overlapping chain set");

		for (i = 0; i < mostOvlpChainSet.size(); i++)
		{
			chn = mostOvlpChainSet[i];
			
			if (nrPdb)
			{			
				AtomContainer ac;
				ac.addAtoms(chn->getAtomPointers());
				for (j = 0; j < chainContacts[chn].size(); j++)
				{
					ac.addAtoms(chainContacts[chn][j]->getAtomPointers());
				}
				ac.writePdb(getFileName(pdbfname) + "." + toString(i + 1) + ".f.pdb");
			}

			chainKept[chn] = true;
			for (j = 0; j < chainContacts[chn].size(); j++)
			{
				chainKept[chainContacts[chn][j]] = true;
			}
		}
		for (i = 0; i < S.chainSize(); i++)
		{
			if (chainKept.find(&(S.getChain(i))) == chainKept.end())
			{
				continue;
			}
			_sys.addAtoms(S.getChain(i).getAtomPointers(), true);
		}
	}
	else
	{
		_sys.addAtoms(S.getAtomPointers(), true);
	}

	AtomSelection sel(_sys.getAtomPointers());
	_pselca = sel.select("name CA");
	ASSERT(_pselca.size() == _sys.positionSize(), "numbers of CA atoms and residues not consistent in file " + pdbfname);
	ASSERT(_pselca.size() > 0, "no valid residues in file " + pdbfname);

#if defined(DEBUG_PARSEPDB)
	cout << "file name: " << _fname << endl;
	cout << "full atoms:" << endl;
	cout << _sys.getAtomPointers() << endl;
	cout << "CA atoms:" << endl;
	cout << _pselca << endl;
	exit(-1);
#endif
}

void ProteinStruct::setBBCoords()
{
	_bbcoor.assign(NUM_BBA * _sys.positionSize(), vector<double>(NUM_COORDS));
	
	for (int i = 0; i < _sys.positionSize(); i++)
	{
		Residue& r = _sys.getResidue(i);

		// N/NT
		if (r.atomExists("N"))
		{
			Atom& a = r.getAtom("N");
			_bbcoor[NUM_BBA * i + N_IDX][0] = a.getX();
			_bbcoor[NUM_BBA * i + N_IDX][1] = a.getY();
			_bbcoor[NUM_BBA * i + N_IDX][2] = a.getZ();
		}
		else
		{
			if (r.atomExists("NT"))
			{
				Atom& a = r.getAtom("NT");
				_bbcoor[NUM_BBA * i + N_IDX][0] = a.getX();
				_bbcoor[NUM_BBA * i + N_IDX][1] = a.getY();
				_bbcoor[NUM_BBA * i + N_IDX][2] = a.getZ();
			}
			else
			{
				_bbcoor[NUM_BBA * i + N_IDX][0] = IMPOSSIBLE_COORD;
				_bbcoor[NUM_BBA * i + N_IDX][1] = IMPOSSIBLE_COORD;
				_bbcoor[NUM_BBA * i + N_IDX][2] = IMPOSSIBLE_COORD;
			}
		}
		// CA
		if (r.atomExists("CA"))
		{
			Atom& a = r.getAtom("CA");
			_bbcoor[NUM_BBA * i + CA_IDX][0] = a.getX();
			_bbcoor[NUM_BBA * i + CA_IDX][1] = a.getY();
			_bbcoor[NUM_BBA * i + CA_IDX][2] = a.getZ();
		}
		else
		{
			error("CA atom should have existed");
		}
		// C
		if (r.atomExists("C"))
		{
			Atom& a = r.getAtom("C");
			_bbcoor[NUM_BBA * i + C_IDX][0] = a.getX();
			_bbcoor[NUM_BBA * i + C_IDX][1] = a.getY();
			_bbcoor[NUM_BBA * i + C_IDX][2] = a.getZ();
		}
		else
		{
			_bbcoor[NUM_BBA * i + C_IDX][0] = IMPOSSIBLE_COORD;
			_bbcoor[NUM_BBA * i + C_IDX][1] = IMPOSSIBLE_COORD;
			_bbcoor[NUM_BBA * i + C_IDX][2] = IMPOSSIBLE_COORD;
		}
		// O/OT1/OT2/OXT
		if (r.atomExists("O"))
		{
			Atom& a = r.getAtom("O");
			_bbcoor[NUM_BBA * i + O_IDX][0] = a.getX();
			_bbcoor[NUM_BBA * i + O_IDX][1] = a.getY();
			_bbcoor[NUM_BBA * i + O_IDX][2] = a.getZ();
		}
		else
		{
			if (r.atomExists("OT1"))
			{
				Atom& a = r.getAtom("OT1");
				_bbcoor[NUM_BBA * i + O_IDX][0] = a.getX();
				_bbcoor[NUM_BBA * i + O_IDX][1] = a.getY();
				_bbcoor[NUM_BBA * i + O_IDX][2] = a.getZ();
			}
			else
			{
				if (r.atomExists("OT2"))
				{
					Atom& a = r.getAtom("OT2");
					_bbcoor[NUM_BBA * i + O_IDX][0] = a.getX();
					_bbcoor[NUM_BBA * i + O_IDX][1] = a.getY();
					_bbcoor[NUM_BBA * i + O_IDX][2] = a.getZ();
				}
				else
				{
					if (r.atomExists("OXT"))
					{
						Atom& a = r.getAtom("OXT");
						_bbcoor[NUM_BBA * i + O_IDX][0] = a.getX();
						_bbcoor[NUM_BBA * i + O_IDX][1] = a.getY();
						_bbcoor[NUM_BBA * i + O_IDX][2] = a.getZ();
					}
					else
					{
						_bbcoor[NUM_BBA * i + O_IDX][0] = IMPOSSIBLE_COORD;
						_bbcoor[NUM_BBA * i + O_IDX][1] = IMPOSSIBLE_COORD;
						_bbcoor[NUM_BBA * i + O_IDX][2] = IMPOSSIBLE_COORD;
					}
				}
			}
		}
	}

#if defined(DEBUG_BB)
	for (int i = 0; i < _sys.positionSize(); i++)
	{
		cout << "\n" << _sys.getResidue(i).getAtomPointers();
		printf("N\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 0][0], _bbcoor[4 * i + 0][1], _bbcoor[4 * i + 0][2]);
		printf("CA\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 1][0], _bbcoor[4 * i + 1][1], _bbcoor[4 * i + 1][2]);
		printf("C\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 2][0], _bbcoor[4 * i + 2][1], _bbcoor[4 * i + 2][2]);
		printf("O\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 3][0], _bbcoor[4 * i + 3][1], _bbcoor[4 * i + 3][2]);
	}
	exit(-1);
#endif
}

void ProteinStruct::setCACoords()
{
	_cacoor.assign(_pselca.size(), vector<double>(NUM_COORDS));
	for (int i = 0; i < _pselca.size(); i++)
	{
		_cacoor[i][0] = _pselca[i]->getX();
		_cacoor[i][1] = _pselca[i]->getY();
		_cacoor[i][2] = _pselca[i]->getZ();
	}

#if defined(DEBUG_CA)
	for (int i = 0; i < _cacoor.size(); i++)
	{
		cout << *(_pselca[i]) << endl;
		printf("%f %f %f\n", _cacoor[i][0], _cacoor[i][1], _cacoor[i][2]);
	}
	exit(-1);
#endif
}

void ProteinStruct::setDihed()
{
	_dihed.clear();

	int i;
	double phi, psi;
	for (i = 0; i < _sys.positionSize(); i++)
	{
		phi = IMPOSSIBLE_ANGLE;
		psi = IMPOSSIBLE_ANGLE;

		if (i > 0)
		{
			Residue& r0 = _sys.getResidue(i - 1);
			Residue& r1 = _sys.getResidue(i);
			if ((!hasBreak(r0, r1)) && (r0.atomExists("C")) && (r1.atomExists("CA")) && (r1.atomExists("C")))
			{
				if (r1.atomExists("N"))
				{
					phi = r0.getAtom("C").dihedral(r1.getAtom("N"), r1.getAtom("CA"), r1.getAtom("C"));
				}
				else
				{
					if (r1.atomExists("NT"))
					{
						phi = r0.getAtom("C").dihedral(r1.getAtom("NT"), r1.getAtom("CA"), r1.getAtom("C"));
					}
				}
			}
		}
		if (i < (_sys.positionSize() - 1))
		{
			Residue& r1 = _sys.getResidue(i);
			Residue& r2 = _sys.getResidue(i + 1);
			if ((!hasBreak(r1, r2)) && (r1.atomExists("CA")) && (r1.atomExists("C")))
			{
				if (r1.atomExists("N"))
				{
					if (r2.atomExists("N"))
					{
						psi = r1.getAtom("N").dihedral(r1.getAtom("CA"), r1.getAtom("C"), r2.getAtom("N"));
					}
					else
					{
						if (r2.atomExists("NT"))
						{
							psi = r1.getAtom("N").dihedral(r1.getAtom("CA"), r1.getAtom("C"), r2.getAtom("NT"));
						}
					}
				}
				else
				{
					if (r1.atomExists("NT"))
					{
						if (r2.atomExists("N"))
						{
							psi = r1.getAtom("NT").dihedral(r1.getAtom("CA"), r1.getAtom("C"), r2.getAtom("N"));
						}
						else
						{
							if (r2.atomExists("NT"))
							{
								psi = r1.getAtom("NT").dihedral(r1.getAtom("CA"), r1.getAtom("C"), r2.getAtom("NT"));
							}
						}
					}
				}
			}
		}		
	
		phi = (phi == -180.0 ? 180.0 : phi);
		psi = (psi == -180.0 ? 180.0 : psi);
		_dihed.push_back(make_pair(phi, psi));
	}
}

void ProteinStruct::setDist()
{
	_dist.clear();

	int i, j;
	for (i = 0; i < _pselca.size(); i++)
	{
		for (j = i + 1; j < _pselca.size(); j++)
		{
			_dist.insert(make_pair(make_pair(i, j), _pselca[i]->distance(*_pselca[j])));
		}
	}

#if defined(DEBUG_DIST)
	for (i = 0 ; i < _pselca.size(); i++)
	{
		for (j = i + 1; j < _pselca.size(); j++)
		{
			printf("%d %d %f\n", i, j, _dist[make_pair(i, j)]);
		}
	}
	exit(-1);
#endif
}

void ProteinStruct::setFileName(const string & pdbfname)
{
	_fname = getFileName(pdbfname);
}

void ProteinStruct::setNumRes()
{
	_numres = _sys.positionSize();
}

void ProteinStruct::setPdbInfo()
{
	_pdbinfo.clear();

	vector<string> resinfo;
	PDBFormat::AtomData atom;
	int i, j, count = 0;
	for (i = 0; i < _sys.positionSize(); i++)
	{
		resinfo.clear();

		Residue& r = _sys.getResidue(i);
		for (j = 0; j < r.atomSize(); j++)
		{
			count++;
			atom.clear();
			atom = PDBFormat::createAtomData(r.getAtom(j));
			atom.D_SERIAL = count;
			resinfo.push_back(PDBFormat::createAtomLine(atom) + string("\n"));
		}

		_pdbinfo.push_back(resinfo);
	}
	ASSERT(_pdbinfo.size() == _sys.positionSize(), "PDB information section missing some residues");
}

void ProteinStruct::setProteinStruct()
{
	setBBCoords();
	setCACoords();	
	setDihed();
	setDist();
	setNumRes();
	setPdbInfo();
	setSeq();	

#if defined(DEBUG_PS)
	int i, j;
	cout << "CA atom coords:" << endl;
	for (i = 0; i < _cacoor.size(); i++)
	{
		printf("%f %f %f\n", _cacoor[i][0], _cacoor[i][1], _cacoor[i][2]);
	}
	cout << "dihedral angles:" << endl;
	for (i = 0; i < _dihed.size(); i++)
	{
		printf("%f %f\n", _dihed[i].first, _dihed[i].second);
	}
	cout << "distance:" << endl;
	for (i = 0 ; i < _numres; i++)
	{
		for (j = i + 1; j < _numres; j++)
		{
			printf("%d %d %f\n", i, j, _dist[make_pair(i, j)]);
		}
	}
	cout << "number of residues: " << _numres << endl;
	cout << "PDB information:" << endl;
	for (i = 0; i < _pdbinfo.size(); i++)
	{
		for (j = 0; j < _pdbinfo[i].size(); j++)
		{
			cout << _pdbinfo[i][j];
		}
	}
	cout << "sequence:" << endl;
	for (i = 0; i < _seq.size(); i++)
	{
		cout << _seq[i] << endl;
	}
	exit(-1);
#endif
}

void ProteinStruct::setSeq()
{
	_seq.clear();
	
	int i;
	for (i = 0; i < _sys.positionSize(); i++)
	{
		_seq.push_back(_sys.getResidue(i).getResidueName());
	}
}

void ProteinStruct::writeBBCoords(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "BBCOORDS\n", copts.getBinary());
	}
	for (i = 0; i < _bbcoor.size(); i++)
	{
		writeDatum(ofs, double(_bbcoor[i][0]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_bbcoor[i][1]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_bbcoor[i][2]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void ProteinStruct::writeCACoords(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "CACOORDS\n", copts.getBinary());
	}
	for (i = 0; i < _cacoor.size(); i++)
	{
		writeDatum(ofs, double(_cacoor[i][0]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_cacoor[i][1]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_cacoor[i][2]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

//	void ProteinStruct::writeDihedralAngles(fstream & ofs,const CreateOptions & copts)
//	{
//		int i;
//		if (!copts.bin)
//		{
//			writeString(ofs, "DIHEDRALANGLES\n", copts.bin);
//		}
//		for (i = 0; i < _dihed.size(); i++)
//		{
//			writeDatum(ofs, double(_dihed[i].first), copts.bin);
//			writeString(ofs, copts.sep, copts.bin);
//			writeDatum(ofs, double(_dihed[i].second), copts.bin);
//			writeString(ofs, copts.ter, copts.bin);
//		}
//		if (!copts.bin)
//		{
//			writeString(ofs, "END\n", copts.bin);
//		}
//	}
//	
//	void ProteinStruct::writeDistance(fstream & ofs,const CreateOptions & copts)
//	{
//		int i, j;
//		if (!copts.bin)
//		{
//			writeString(ofs, "DISTANCE\n", copts.bin);
//		}
//		for (i = 0; i < _numres; i++)
//		{
//			for (j = i + 1; j < _numres; j++)
//			{
//				writeDatum(ofs, int(i), copts.bin);
//				writeString(ofs, copts.sep, copts.bin);
//				writeDatum(ofs, int(j), copts.bin);
//				writeString(ofs, copts.sep, copts.bin);
//				writeDatum(ofs, double(_dist[make_pair(i, j)]), copts.bin);
//				writeString(ofs, copts.ter, copts.bin);
//			}
//		}
//		if (!copts.bin)
//		{
//			writeString(ofs, "END\n", copts.bin);
//		}
//	}

void ProteinStruct::writePdbInfo(fstream & ofs,const CreateOptions & copts)
{
	int i, j, countbytes;
	if (copts.getBinary())
	{
		countbytes = sizeof(int) * _pdbinfo.size();
		for (i = 0; i < _pdbinfo.size(); i++)
		{
			writeDatum(ofs, int(countbytes), copts.getBinary());
			countbytes += sizeof(int); // for number of chars
			for (j = 0; j < _pdbinfo[i].size(); j++)
			{
				countbytes += sizeof(char) * _pdbinfo[i][j].size();
			}
		}
		for (i = 0; i < _pdbinfo.size(); i++)
		{
			countbytes = 0;
			for (j = 0; j < _pdbinfo[i].size(); j++)
			{
				countbytes += sizeof(char) * _pdbinfo[i][j].size();
			}
			writeDatum(ofs, int(countbytes), copts.getBinary());
			for (j = 0; j < _pdbinfo[i].size(); j++)
			{
				writeString(ofs, _pdbinfo[i][j], copts.getBinary());
			}
		}
	}
	else
	{
		writeString(ofs, "PDBINFO\n", copts.getBinary());
		for (i = 0; i < _pdbinfo.size(); i++)
		{
			writeDatum(ofs, int(_pdbinfo[i].size()), copts.getBinary());
			writeString(ofs, copts.getWordTer(), copts.getBinary());
			for (j = 0; j < _pdbinfo[i].size(); j++)
			{
				writeString(ofs, _pdbinfo[i][j], copts.getBinary());
			}
		}
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void ProteinStruct::writeProteinStruct(fstream & ofs,const CreateOptions & copts)
{
	writeBBCoords(ofs, copts);
	writeCACoords(ofs, copts);
//		writeDihedralAngles(ofs, copts);
//		writeDistance(ofs, copts);
	writePdbInfo(ofs, copts);
	writeSeq(ofs, copts);
}

void ProteinStruct::writeProteinStructFile(fstream & ofs,const CreateOptions & copts)
{
	int countbytes = writeProteinStructFileHeader(ofs, copts, NUM_SEC_PS);
	writeProteinStruct(ofs, copts);

	if (copts.getBinary())
	{
		streampos begin, end;
		ofs.seekp(0, ios::beg);
		begin = ofs.tellp();
		ofs.seekp(0, ios::end);
		end = ofs.tellp();
		ASSERT(countbytes == (end - begin), "size of protein structure file not consistent");
	}

	cout << "Protein structure file created." << endl;
}

int ProteinStruct::writeProteinStructFileHeader(fstream & ofs, const CreateOptions & copts, const int & numsec)
{
	if (!copts.getBinary())
	{
		writeDatum(ofs, int(_numres), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
		return -1;
	}

	int i, j, countbytes = sizeof(int) * numsec;

	// offset of BB coords section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * NUM_COORDS * _bbcoor.size();
	// offset of CA coords section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * NUM_COORDS * _cacoor.size();

//		// offset of dihed section
//		writeDatum(ofs, int(countbytes), copts.bin);
//		countbytes += sizeof(double) * 2 * _dihed.size();
//		// offset of dist section
//		writeDatum(ofs, int(countbytes), copts.bin);
//		countbytes += (sizeof(int) * 2 + sizeof(double)) * _dist.size();

	// number of residues
	writeDatum(ofs, int(_numres), copts.getBinary());	
	// offset of PDB info section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(int) * _pdbinfo.size();
	for (i = 0; i < _pdbinfo.size(); i++)
	{
		countbytes += sizeof(int); // for number of chars
		for (j = 0; j < _pdbinfo[i].size(); j++)
		{
			countbytes += sizeof(char) * _pdbinfo[i][j].size();
		}
	}	
	// offset of seq section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(char) * LEN_AA_CODE * _seq.size();

	return countbytes;
}

void ProteinStruct::writeSeq(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "SEQ\n", copts.getBinary());
	}
	for (i = 0; i < _seq.size(); i++)
	{
		writeString(ofs, pad(_seq[i], LEN_AA_CODE), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

class QueryStruct : public ProteinStruct
{
	public:
		QueryStruct() {}
		QueryStruct(const string &, const vector<string> &, const double &);
		~QueryStruct() {}

		void checkQueryStructFile(FILE* &, const string &, const CreateOptions &);
		void setCenRes(const double &);
		void setNumSeg();
		void setQueryStruct(const double &);
		void setResBefBrk();
		void writeCenRes(fstream &, const CreateOptions &);
		void writeCenResDihed(fstream &, const CreateOptions &);
		void writeQueryStruct(fstream &, const CreateOptions &);
		void writeQueryStructFile(fstream &, const CreateOptions &);
		int writeQueryStructFileHeader(fstream &, const CreateOptions &);
		void writeResBefBrk(fstream &, const CreateOptions &);

	protected:
		vector<int> _befbrk;
		vector<int> _cenres;
		int _numseg;
};

QueryStruct::QueryStruct(const string & pdbfname, const vector<string> & legalAA, const double & dcut)
{
	parsePdb(pdbfname, legalAA);
	setQueryStruct(dcut);
}

void QueryStruct::checkQueryStructFile(FILE* & ifp, const string & qsfname,const CreateOptions & copts)
{
	checkProteinStructFile(ifp, qsfname, copts);
}

void QueryStruct::setCenRes(const double & dcut)
{
	_cenres.clear();

	int i, j, k, crind;
	double radius, maxdist, d;
	for (i = 0; i < _numseg; i++)
	{
		radius = DBL_MAX;
		crind = -1;
		for (j = _befbrk[i] + 1; j <= _befbrk[i + 1]; j++)
		{
			maxdist = 0.0;
			for (k = _befbrk[i] + 1; k <= _befbrk[i + 1]; k++)
			{
				d = _pselca[j]->distance(*_pselca[k]);
				if (d > maxdist)
				{
					maxdist = d;
				}						
			}
			if (maxdist < radius)
			{
				radius = maxdist;
				crind = j;
			}
		}
		ASSERT(crind > -1, "could not find a center residue");
		_cenres.push_back(crind);
	}	

	// distance between any pair of center residues should be <= distance cutoff
	// otherwise, not a good query structure
	for (i = 0; i < _cenres.size(); i++)
	{
		for (j = i + 1; j < _cenres.size(); j++)
		{
			ASSERT(_pselca[_cenres[i]]->distance(*_pselca[_cenres[j]]) <= dcut, "distance between certain pair of center residues > distance cutoff, " + _fname + " not a good query");
		}
	}
}

void QueryStruct::setNumSeg()
{
	_numseg = _befbrk.size() - 1;
}

void QueryStruct::setQueryStruct(const double & dcut)
{
	setProteinStruct();
	setResBefBrk();
	setNumSeg();
	setCenRes(dcut);

#if defined(DEBUG_QS)
	cout << "residues before breaks: ";
	for (int i = 0; i < _befbrk.size(); i++)
	{
		cout << _befbrk[i] << " ";
	}
	cout << endl;
	cout << "center residues: ";
	for (int i = 0; i < _cenres.size(); i++)
	{
		cout << _cenres[i] << " ";
	}
	cout << endl;
	cout << "number of segments: " << _numseg << endl;
	exit(-1);
#endif	
}

void QueryStruct::setResBefBrk()
{
	_befbrk.clear();
	
	_befbrk.push_back(-1);
	int i;
	for (i = 0; i < (_sys.positionSize() - 1); i++)
	{
		if (hasBreak(_sys.getResidue(i), _sys.getResidue(i + 1)))
		{
			_befbrk.push_back(i);
		}
	}
	_befbrk.push_back(_sys.positionSize() - 1);
}

void QueryStruct::writeCenRes(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "CENRES\n", copts.getBinary());
	}
	for (i = 0; i < _cenres.size(); i++)
	{
		writeDatum(ofs, int(_cenres[i]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void QueryStruct::writeCenResDihed(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "CENTERRESIDUEDIHEDRALANGLES\n", copts.getBinary());
	}
	for (i = 0; i < _cenres.size(); i++)
	{
		writeDatum(ofs, double(_dihed[_cenres[i]].first), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_dihed[_cenres[i]].second), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void QueryStruct::writeQueryStruct(fstream & ofs, const CreateOptions & copts)
{
	writeProteinStruct(ofs, copts);
	writeResBefBrk(ofs, copts);
	writeCenRes(ofs, copts);
	writeCenResDihed(ofs, copts);
}

void QueryStruct::writeQueryStructFile(fstream & ofs, const CreateOptions & copts)
{
	int countbytes = writeQueryStructFileHeader(ofs, copts);
	writeQueryStruct(ofs, copts);

	if (copts.getBinary())
	{
		streampos begin, end;
		ofs.seekp(0, ios::beg);
		begin = ofs.tellp();
		ofs.seekp(0, ios::end);
		end = ofs.tellp();
		ASSERT(countbytes == (end - begin), "size of query structure file not consistent");
	}

	cout << "Query structure " << _fname << " created." << endl;
}

int QueryStruct::writeQueryStructFileHeader(fstream & ofs, const CreateOptions & copts)
{
	int countbytes = writeProteinStructFileHeader(ofs, copts, NUM_SEC_QS);

	if (!copts.getBinary())
	{
		writeDatum(ofs, int(_numseg), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
		return -1;
	}

	// offset of residues before breaks section		
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(int) * _befbrk.size();
	// offset of center residues section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(int) * _cenres.size();
	// offset of center residue dihedral angles section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * 2 * _cenres.size();
	// number of segments
	writeDatum(ofs, int(_numseg), copts.getBinary());	

	return countbytes;
}

void QueryStruct::writeResBefBrk(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "RESIDUESBEFOREBREAKS\n", copts.getBinary());
	}
	for (i = 0; i < _befbrk.size(); i++)
	{
		writeDatum(ofs, int(_befbrk[i]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

class TargetStruct : public ProteinStruct
{
	public:
		TargetStruct() {}		
		TargetStruct(const string &, const vector < string > &, const bool &, const double &, const double &, const double &, const bool &, const double &, const double &, const double &, const double &);
		~TargetStruct() {}

		void checkTargetStructFile(FILE* &, const string &, const CreateOptions &);
		void setDihedDistr(const double &, const double &);
		void setDistDistr(const double &, const double &);
		void setTargetStruct(const double &, const double &, const double &, const double &);
		void writeDihedDistr(fstream &, const CreateOptions &);
		void writeDistDistr(fstream &, const CreateOptions &);
		void writeTargetStruct(fstream &, const CreateOptions &);
		void writeTargetStructFile(fstream &, const CreateOptions &);
		int writeTargetStructFileHeader(fstream &, const CreateOptions &);

	protected:
		double _dcut;		
		vector<vector<vector<int> > > _diheddistr;
		vector<vector<vector<int> > > _distdistr;
		double _dstep;
		double _phistep;
		double _psistep;
};

TargetStruct::TargetStruct(const string & pdbfname, const vector < string > & legalAA, const bool & nr, const double & seqIdenThresh, const double & optRmsdThresh, const double & contRmsdThresh, const bool & nrPdb, const double & phistep, const double & psistep, const double & dcut, const double & dstep)
{
	parsePdb(pdbfname, legalAA, nr, seqIdenThresh, optRmsdThresh, contRmsdThresh, nrPdb);
	setTargetStruct(phistep, psistep, dcut, dstep);
}

void TargetStruct::checkTargetStructFile(FILE * & ifp,const string & tsfname,const CreateOptions & copts)
{
	checkProteinStructFile(ifp, tsfname, copts);
}

void TargetStruct::setDihedDistr(const double & phistep, const double & psistep)
{
	_diheddistr.clear();
	_phistep = phistep;
	_psistep = psistep;

	int i, j, ri, numphibin, numpsibin;
	double phi, psi;

	// (0, step], (step, 2step], (2step, 3step], ...
	numphibin = int(ceil(360.0 / phistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	_diheddistr.resize(numphibin);
	numpsibin = int(ceil(360.0 / psistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	for (i = 0; i < _diheddistr.size(); i++)
	{
		_diheddistr[i].resize(numpsibin);
	}
	for (ri = 0; ri < _dihed.size(); ri++)
	{
		phi = _dihed[ri].first;
		psi = _dihed[ri].second;
		if (phi >= IMPOSSIBLE_ANGLE)
		{
			i = numphibin - 1;
		}
		else
		{
			i = int(ceil((phi + 180.0) / phistep)) - 1;
		}		
		if (psi >= IMPOSSIBLE_ANGLE)
		{
			j = numpsibin - 1;
		}
		else
		{
			j = int(ceil((psi + 180.0) / psistep)) - 1;
		}		
		_diheddistr[i][j].push_back(ri);
	}
}

void TargetStruct::setDistDistr(const double & dcut, const double & dstep)
{
	_distdistr.clear();
	_dcut = dcut;
	_dstep = dstep;

	int i, j, numbin, bi;
	double d;
			
	_distdistr.resize(_numres);
	// (0, dstep], (dstep, 2dstep], (2dstep, 3dstep], ...
	numbin = int(ceil(dcut / dstep));
	for (i = 0 ; i < _distdistr.size(); i++)
	{
		_distdistr[i].resize(numbin);
	}
	for (i = 0; i < _numres; i++)
	{
		for (j = i + 1; j < _numres; j++)
		{
			d = _dist[make_pair(i, j)];
			if (d > dcut)
			{
				continue;
			}
			if (d <= 0.0)
			{
				bi = 0;
			}
			else
			{
				bi = int(ceil(d / dstep)) - 1;
			}
			_distdistr[i][bi].push_back(j);
			_distdistr[j][bi].push_back(i);
		}
	}
}

void TargetStruct::setTargetStruct(const double & phistep, const double & psistep, const double & dcut, const double & dstep)
{
	setProteinStruct();
	setDihedDistr(phistep, psistep);
	setDistDistr(dcut, dstep);
}

void TargetStruct::writeDihedDistr(fstream & ofs,const CreateOptions & copts)
{
	int i, j, k, countbytes;
	if (copts.getBinary())
	{
		writeDatum(ofs, double(_phistep), copts.getBinary());
		writeDatum(ofs, double(_psistep), copts.getBinary());
		countbytes = sizeof(double) * 2 + sizeof(int) * _diheddistr.size() * _diheddistr[0].size(); // 2 for phi and psi steps
		for (i = 0; i < _diheddistr.size(); i++)
		{
			for (j = 0; j < _diheddistr[i].size(); j++)
			{
				if (_diheddistr[i][j].size() == 0)
				{
					writeDatum(ofs, int(_diheddistr[i][j].size()), copts.getBinary());
				}
				else
				{
					writeDatum(ofs, int(countbytes), copts.getBinary());
					countbytes += sizeof(int) * (1 + _diheddistr[i][j].size()); // 1 for number of elements
				}
			}
		}
		for (i = 0; i < _diheddistr.size(); i++)
		{
			for (j = 0; j < _diheddistr[i].size(); j++)
			{
				if (_diheddistr[i][j].size() > 0)
				{
					writeDatum(ofs, int(_diheddistr[i][j].size()), copts.getBinary());
					for (k = 0; k < _diheddistr[i][j].size(); k++)
					{
						writeDatum(ofs, int(_diheddistr[i][j][k]), copts.getBinary());
					}
				}
			}
		}
	}
	else
	{
		writeString(ofs, "DIHEDRALANGLESDISTRIBUTION\n", copts.getBinary());
		writeDatum(ofs, double(_phistep), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_psistep), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
		for (i = 0; i < _diheddistr.size(); i++)
		{
			for (j = 0; j < _diheddistr[i].size(); j++)
			{
				writeDatum(ofs, int(_diheddistr[i][j].size()), copts.getBinary());
				writeString(ofs, copts.getWordTer(), copts.getBinary());
				for (k = 0; k < _diheddistr[i][j].size(); k++)
				{
					writeDatum(ofs, int(_diheddistr[i][j][k]), copts.getBinary());
					writeString(ofs, copts.getWordTer(), copts.getBinary());
				}
			}
		}
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void TargetStruct::writeDistDistr(fstream & ofs,const CreateOptions & copts)
{
	int i, j, k, countbytes;
	if (copts.getBinary())
	{
		writeDatum(ofs, double(_dcut), copts.getBinary());
		writeDatum(ofs, double(_dstep), copts.getBinary());
		countbytes = sizeof(double) * 2 + sizeof(int) * _distdistr.size() * _distdistr[0].size(); // 2 for distance cutoff and step
		for (i = 0; i < _distdistr.size(); i++)
		{
			for (j = 0; j < _distdistr[i].size(); j++)
			{
				if (_distdistr[i][j].size() == 0)
				{
					writeDatum(ofs, int(_distdistr[i][j].size()), copts.getBinary());
				}
				else
				{
					writeDatum(ofs, int(countbytes), copts.getBinary());
					countbytes += sizeof(int) * (1 + _distdistr[i][j].size()); // 1 for number of elements
				}
			}
		}
		for (i = 0; i < _distdistr.size(); i++)
		{
			for (j = 0; j < _distdistr[i].size(); j++)
			{
				if (_distdistr[i][j].size() > 0)
				{
					writeDatum(ofs, int(_distdistr[i][j].size()), copts.getBinary());
					for (k = 0; k < _distdistr[i][j].size(); k++)
					{
						writeDatum(ofs, int(_distdistr[i][j][k]), copts.getBinary());
					}
				}
			}
		}
	}
	else
	{
		writeString(ofs, "DISTANCEDISTRIBUTION\n", copts.getBinary());
		writeDatum(ofs, double(_dcut), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_dstep), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
		for (i = 0; i < _distdistr.size(); i++)
		{
			for (j = 0; j < _distdistr[i].size(); j++)
			{
				writeDatum(ofs, int(_distdistr[i][j].size()), copts.getBinary());
				writeString(ofs, copts.getWordTer(), copts.getBinary());
				for (k = 0; k < _distdistr[i][j].size(); k++)
				{
					writeDatum(ofs, int(_distdistr[i][j][k]), copts.getBinary());
					writeString(ofs, copts.getWordTer(), copts.getBinary());
				}
			}
		}
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void TargetStruct::writeTargetStruct(fstream & ofs, const CreateOptions & copts)
{
	writeProteinStruct(ofs, copts);
	writeDihedDistr(ofs, copts);
	writeDistDistr(ofs, copts);
}

void TargetStruct::writeTargetStructFile(fstream & ofs, const CreateOptions & copts)
{
	int countbytes = writeTargetStructFileHeader(ofs, copts);
	writeTargetStruct(ofs, copts);

	if (copts.getBinary())
	{
		streampos begin, end;
		ofs.seekp(0, ios::beg);
		begin = ofs.tellp();
		ofs.seekp(0, ios::end);
		end = ofs.tellp();
		ASSERT(countbytes == (end - begin), "size of target structure file not consistent");
	}

	cout << "Target structure " << _fname << " created." << endl;	
}

int TargetStruct::writeTargetStructFileHeader(fstream & ofs, const CreateOptions & copts)
{
	int countbytes = writeProteinStructFileHeader(ofs, copts, NUM_SEC_TS);

	if (!copts.getBinary())
	{
		return -1;
	}

	int i, j;
	// offset of dihedral angles distribution section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * 2 + sizeof(int) * _diheddistr.size() * _diheddistr[0].size(); // 2 for phi and psi steps
	for (i = 0; i < _diheddistr.size(); i++)
	{
		for (j = 0; j < _diheddistr[i].size(); j++)
		{
			if (_diheddistr[i][j].size() > 0)
			{
				countbytes += sizeof(int) * (1 + _diheddistr[i][j].size()); // 1 for number of elements
			}
		}
	}
	// offset of distance distribution section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * 2 + sizeof(int) * _distdistr.size() * _distdistr[0].size(); // 2 for distance cutoff and step
	for (i = 0; i < _distdistr.size(); i++)
	{
		for (j = 0; j < _distdistr[i].size(); j++)
		{
			if (_distdistr[i][j].size() > 0)
			{
				countbytes += sizeof(int) * (1 + _distdistr[i][j].size()); // 1 for number of elements
			}
		}
	}
	
	return countbytes;
}

bool alignChains(const vector<string> & seqA, const vector<string> & seqB, const double & seqIdenThresh, AtomPointerVector & caA, AtomPointerVector & caB, const int & mode, const double & optRmsdThresh, vector<double> & t, vector<vector<double> > & u, const double & contRmsdThresh, const vector<double> & ct, const vector<vector<double> > & cu)
{
	// True: two chains are redundant
	// False: two chains are non-redundant
	
	// mode=1: align central chains, and calculate optimal RMSD as well as rotation and translation matrices
	// mode=0: align corresponding neighbors, and calculate optimal RMSD as well as contact RMSD
	
	vector<string> seqAAligned, seqBAligned;	
	if (alignSeqNeedlemanWunsch(seqA, seqB, seqIdenThresh, seqAAligned, seqBAligned))
	{
		AtomPointerVector caAAligned, caBAligned;
		int i, j = 0, k = 0;
		for (i = 0; i < seqAAligned.size(); i++)
		{
			if ((seqAAligned[i] != "---") && (seqBAligned[i] != "---"))
			{
				caAAligned.push_back(caA[j]);
				caBAligned.push_back(caB[k]);
			}
			if (seqAAligned[i] != "---")
			{
				j++;
			}
			if (seqBAligned[i] != "---")
			{
				k++;
			}
		}
		ASSERT((j == caA.size()) && (k == caB.size()), "number of Ca atoms not consistent");
		if (calcRmsdKabsch(caAAligned, caBAligned, mode, t, u) <= optRmsdThresh)
		{
			if (mode)
			{
				return true;
			}
			double contRmsd = 0.0, xa, ya, za, xb, yb, zb;
			for (i = 0; i < caAAligned.size(); i++)
			{
				xa = ct[0] + cu[0][0] * caAAligned[i]->getX() + cu[0][1] * caAAligned[i]->getY() + cu[0][2] * caAAligned[i]->getZ();
				ya = ct[1] + cu[1][0] * caAAligned[i]->getX() + cu[1][1] * caAAligned[i]->getY() + cu[1][2] * caAAligned[i]->getZ();
				za = ct[2] + cu[2][0] * caAAligned[i]->getX() + cu[2][1] * caAAligned[i]->getY() + cu[2][2] * caAAligned[i]->getZ();
				xb = caBAligned[i]->getX();
				yb = caBAligned[i]->getY();
				zb = caBAligned[i]->getZ();
				contRmsd += ((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb));
			}
			contRmsd = sqrt(contRmsd / double(caAAligned.size()));
			if (contRmsd <= contRmsdThresh)
			{
				return true;
			}
		}
	}

	return false;
}

bool alignSeqNeedlemanWunsch(const vector<string> & seqA, const vector<string> & seqB, const double & seqIdenThresh, vector<string> & seqAAligned, vector<string> & seqBAligned)
{
	// True: two sequences are redundant
	// False: two sequences are not redundant
	seqAAligned.clear();
	seqBAligned.clear();
/*	
	const int PAM250[20][20] = 
	{
		{ 2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0}, // A: ALA
		{-2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2}, // R: ARG
		{ 0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2}, // N: ASN
		{ 0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2}, // D: ASP
		{-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2}, // C: CYS
		{ 0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2}, // Q: GLN
		{ 0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2}, // E: GLU
		{ 1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1}, // G: GLY
		{-1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2}, // H: HIS, HSC, HSD, HSE, HSP
		{-1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4}, // I: ILE
		{-2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2}, // L: LEU
		{-1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2}, // K: LYS
		{-1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2}, // M: MET, MSE
		{-3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1}, // F: PHE
		{ 1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1}, // P: PRO
		{ 1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1}, // S: SER
		{ 1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0}, // T: THR
		{-6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6}, // W: TRP
		{-3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2}, // Y: TYR
		{ 0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4}  // V: VAL
	};
*/
	const int PAM250[21][21] = 
	{
		{ 2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0, -2}, // A: ALA
		{-2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -2}, // R: ARG
		{ 0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2, -2}, // N: ASN
		{ 0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2, -2}, // D: ASP
		{-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -2}, // C: CYS
		{ 0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2, -2}, // Q: GLN
		{ 0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2, -2}, // E: GLU
		{ 1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1, -2}, // G: GLY
		{-1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2, -2}, // H: HIS, HSC, HSD, HSE, HSP
		{-1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2}, // I: ILE
		{-2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -2}, // L: LEU
		{-1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2, -2}, // K: LYS
		{-1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2}, // M: MET, MSE
		{-3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -2}, // F: PHE
		{ 1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -2}, // P: PRO
		{ 1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1, -2}, // S: SER
		{ 1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0, -2}, // T: THR
		{-6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -2}, // W: TRP
		{-3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -2}, // Y: TYR
		{ 0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2}, // V: VAL
		{-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  6}  // Unnatural: CSO, HIP, PTR, SEC, SEP, TPO
	};

	map<string, int> aaToIdx;
	aaToIdx["ALA"] = 0;
	aaToIdx["ARG"] = 1;
	aaToIdx["ASN"] = 2;
	aaToIdx["ASP"] = 3;
	aaToIdx["CYS"] = 4;
	aaToIdx["GLN"] = 5;
	aaToIdx["GLU"] = 6;
	aaToIdx["GLY"] = 7;
	aaToIdx["HIS"] = 8; aaToIdx["HSC"] = 8; aaToIdx["HSD"] = 8; aaToIdx["HSE"] = 8; aaToIdx["HSP"] = 8;
	aaToIdx["ILE"] = 9;
	aaToIdx["LEU"] = 10;
	aaToIdx["LYS"] = 11;
	aaToIdx["MET"] = 12; aaToIdx["MSE"] = 12;
	aaToIdx["PHE"] = 13;
	aaToIdx["PRO"] = 14;
	aaToIdx["SER"] = 15;
	aaToIdx["THR"] = 16;
	aaToIdx["TRP"] = 17;
	aaToIdx["TYR"] = 18;
	aaToIdx["VAL"] = 19;
	// unnatural AA's
	aaToIdx["CSO"] = 20; aaToIdx["HIP"] = 20; aaToIdx["PTR"] = 20; aaToIdx["SEC"] = 20; aaToIdx["SEP"] = 20; aaToIdx["TPO"] = 20;
	
	int i, j, maxelem, maxi, maxj, nIden;
	const int gap = -1;
	double seqIden;
	vector<vector<int> > alignMat(seqA.size() + 1); // row
	for (i = 0; i < alignMat.size(); i++)
	{
		alignMat[i].resize(seqB.size() + 1); // column
	}
	// complete dynamic programming matrix
	alignMat[0][0] = 0;
	for (i = 1; i <= seqA.size(); i++)
	{
		alignMat[i][0] = alignMat[i - 1][0] + gap; // first column
	}
	for (j = 1; j <= seqB.size(); j++)
	{
		alignMat[0][j] = alignMat[0][j - 1] + gap; // first row
	}	
	for (i = 1; i < alignMat.size(); i++)
	{
		for (j = 1; j < alignMat[i].size(); j++)
		{
			alignMat[i][j] = max(alignMat[i - 1][j - 1] + PAM250[aaToIdx[seqA[i - 1]]][aaToIdx[seqB[j - 1]]], max(alignMat[i - 1][j] + gap, alignMat[i][j - 1] + gap));
		}
	}
	maxelem = alignMat[seqA.size()][seqB.size()];
	maxi = seqA.size();
	maxj = seqB.size();
	for (i = 0; i < seqA.size(); i++)
	{
		if (alignMat[i][seqB.size()] > maxelem) // last column
		{
			maxi = i;
			maxj = seqB.size();
			maxelem = alignMat[i][seqB.size()];
		}
	}
	for (j = 0; j < seqB.size(); j++)
	{
		if (alignMat[seqA.size()][j] > maxelem) // last row
		{
			maxi = seqA.size();
			maxj = j;
			maxelem = alignMat[seqA.size()][j];
		}
	}
	ASSERT((maxi == seqA.size()) || (maxj == seqB.size()), "sequence alignment not complete");
	// trace back
	for (i = seqA.size(); i > maxi; i--) // last column
	{
		seqAAligned.push_back(seqA[i - 1]);
		seqBAligned.push_back("---");
	}
	for (j = seqB.size(); j > maxj; j--) // last row
	{
		seqBAligned.push_back(seqB[j - 1]);
		seqAAligned.push_back("---");
	}
	i = maxi;
	j = maxj;
	while ((i > 0) && (j > 0))
	{
		if (alignMat[i][j] == alignMat[i - 1][j - 1] + PAM250[aaToIdx[seqA[i - 1]]][aaToIdx[seqB[j - 1]]])
		{
			seqAAligned.push_back(seqA[i - 1]);
			seqBAligned.push_back(seqB[j - 1]);
			i--;
			j--;
			continue;
		}
		if (alignMat[i][j] == alignMat[i - 1][j] + gap)
		{
			seqAAligned.push_back(seqA[i - 1]);
			seqBAligned.push_back("---");
			i--;
			continue;
		}
		if (alignMat[i][j] == alignMat[i][j - 1] + gap)
		{
			seqAAligned.push_back("---");
			seqBAligned.push_back(seqB[j - 1]);
			j--;
			continue;
		}
	}
	while (i > 0) // first column
	{
		seqAAligned.push_back(seqA[i - 1]);
		seqBAligned.push_back("---");
		i--;
	}
	while (j > 0) // first row
	{
		seqBAligned.push_back(seqB[j - 1]);
		seqAAligned.push_back("---");
		j--;
	}
	reverse(seqAAligned.begin(), seqAAligned.end());
	reverse(seqBAligned.begin(), seqBAligned.end());
	
#if defined(DEBUG_SEQALI)
	cout << "sequence A:" << endl;
	for (i = 0; i < seqA.size(); i++)
	{
		cout << " " << seqA[i];
	}
	cout << endl;
	cout << "sequence B:" << endl;
	for (j = 0; j < seqB.size(); j++)
	{
		cout << " " << seqB[j];
	}
	cout << endl;
	cout << "aligned sequence A:" << endl;
	for (i = 0; i < seqAAligned.size(); i++)
	{
		cout << " " << seqAAligned[i];
	}
	cout << endl;
	cout << "aligned sequence B:" << endl;
	for (j = 0; j < seqBAligned.size(); j++)
	{
		cout << " " << seqBAligned[j];
	}
	cout << endl;
#endif

	ASSERT(seqAAligned.size() == seqBAligned.size(), "length of aligned sequences not consistent");
	
	// calculate sequence identity
	nIden = 0;
	for (i = 0; i < seqAAligned.size(); i++)
	{
		if ((seqAAligned[i] != "---") && (seqBAligned[i] != "---") && (seqAAligned[i] == seqBAligned[i]))
		{
			nIden++;
		}
	}
	seqIden = double(nIden) / double(seqAAligned.size());

#if defined(DEBUG_SEQALI)
	cout << "sequence identity: " << seqIden << endl;
	exit(-1);
#endif
	
	if (seqIden < seqIdenThresh)
	{
		return false;
	}
	return true;
}

// defines whether two crystal units are in contact
// count how many pairs of atoms from two units are within a distance cutoff
bool areUnitsInContact(AtomPointerVector* A, AtomPointerVector* B)
{
	double dcut = 6.0;
  	double dcut2 = dcut * dcut;
  	int i, n = 0, x, y, z, xi, yi, zi, ai;
  	int nthresh = min(20, int(ceil(0.5 * (A->size() + B->size()) / 10)));
	string atomname;

  	// create a grid and place B into it
  	map<int, map<int, map<int, AtomPointerVector> > > grid;
  	CartesianPoint c = A->getGeometricCenter();
  	for (i = 0; i < B->size(); i++)
  	{
		atomname = (*B)[i]->getName();
		if (atomname[0] == 'H')
		{
			continue;
		}
		
		gridPoint((*B)[i], c, dcut, &x, &y, &z);
    	if (grid.find(x) == grid.end())
    	{
			grid[x] = map<int, map<int, AtomPointerVector> >();
    	}
    	if (grid[x].find(y) == grid[x].end())
    	{
			grid[x][y] = map<int, AtomPointerVector> ();
    	}
    	if (grid[x][y].find(z) == grid[x][y].end())
    	{
			grid[x][y][z] = AtomPointerVector ();
    	}
    	grid[x][y][z].push_back((*B)[i]);
  	}

  	// visit all atoms in A
  	for (i = 0; i < A->size(); i++)
  	{
		atomname = (*A)[i]->getName();
		if (atomname[0] == 'H')
		{
			continue;
		}
		
    	// what grid box does this atom fit into?
    	gridPoint((*A)[i], c, dcut, &x, &y, &z);
    	for (xi = x - 1; xi <= x + 1; xi++)
    	{
      		if (grid.find(xi) == grid.end())
      		{
	  			continue;
      		}
      		map<int, map<int, AtomPointerVector> >& gridX = grid[xi];
      		for (yi = y - 1; yi <= y + 1; yi++)
      		{
        		if (gridX.find(yi) == gridX.end())
        		{
					continue;
        		}
        		map<int, AtomPointerVector>& gridXY = gridX[yi];
        		for (zi = z - 1; zi <= z + 1; zi++)
        		{
          			if (gridXY.find(zi) == gridXY.end())
          			{
		  				continue;
          			}
          			AtomPointerVector& gridXYZ = gridXY[zi];
          			for (ai = 0; ai < gridXYZ.size(); ai++)
          			{
            			if ((*A)[i]->distance2(*(gridXYZ[ai])) < dcut2)
            			{
							n++;
            			}
          			}
        		}
      		}
    	}
    	if (n > nthresh)
    	{
			return true;
    	}
  	}
  	return false;
}

int calcNumOvlpChains(Chain* curChain, vector<Chain*> & curChainSet, map<Chain*, vector<Chain*> > & chainContacts)
{
	int i, j, numOvlpChains = 0;
	map<Chain*, bool> chainInc;

	for (i = 0; i < curChainSet.size(); i++)
	{
		chainInc[curChainSet[i]] = true;
		for (j = 0; j < chainContacts[curChainSet[i]].size(); j++)
		{
			chainInc[chainContacts[curChainSet[i]][j]] = true;
		}
	}

	if (chainInc.find(curChain) != chainInc.end())
	{
		numOvlpChains++;
	}
	for (i = 0; i < chainContacts[curChain].size(); i++)
	{
		if (chainInc.find(chainContacts[curChain][i]) != chainInc.end())
		{
			numOvlpChains++;
		}
	}

	return numOvlpChains;
}

/**************************************************************************
  Implemetation of Kabsch algoritm for finding the best rotation matrix
---------------------------------------------------------------------------
  x    - x(i,m) are coordinates of atom m in set x            (input)
  y    - y(i,m) are coordinates of atom m in set y            (input)
  n    - n is number of atom pairs                            (input)
  mode  - 0:calculate rmsd only                               (input)
          1:calculate rmsd,u,t                                (takes longer)
  rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
  u    - u(i,j) is   rotation  matrix for best superposition  (output)
  t    - t(i)   is translation vector for best superposition  (output)
**************************************************************************/
double calcRmsdKabsch(AtomPointerVector &_align, AtomPointerVector &_ref, int mode, vector<double> & t, vector<vector<double> > & u)
{
	int i, j, m, m1, l, k;
	double e0, rms1, d, h, g;
	double cth, sth, sqrth, p, det, sigma;  
	double xc[3], yc[3];
	double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
	double sqrt3=1.73205080756888, tol=0.01;
	int ip[]={0, 1, 3, 1, 2, 4, 3, 4, 5};
	int ip2312[]={1, 2, 0, 1};
	
	int a_failed=0, b_failed=0;
	double epsilon=0.00000001;
	
	int n=_ref.size();
	ASSERT((n >= 0) && (n == _align.size()), "invalid protein length for calculating RMSD");

	if (0 == n)
	{
		return 0.0;
	}
	
	//initializtation
	double rmsd=0;
	rms1=0;
	e0=0;
	for (i=0; i<3; i++) {
		xc[i]=0.0;
		yc[i]=0.0;
		t[i]=0.0;
		for (j=0; j<3; j++) {
			u[i][j]=0.0;
			r[i][j]=0.0;
			a[i][j]=0.0;
			if (i==j) {
				u[i][j]=1.0;
				a[i][j]=1.0;
			}
		}
	} 

	//compute centers for vector sets x, y
	for(i=0; i<n; i++){
		xc[0] += _align[i]->getX();
		xc[1] += _align[i]->getY();
		xc[2] += _align[i]->getZ();
		
		yc[0] += _ref[i]->getX();
		yc[1] += _ref[i]->getY();
		yc[2] += _ref[i]->getZ();
	}
	for(i=0; i<3; i++){
		xc[i] = xc[i]/(double)n;
		yc[i] = yc[i]/(double)n;        
	}
	
	//compute e0 and matrix r
	for (m=0; m<n; m++) {
		e0 += (_align[m]->getX()-xc[0])*(_align[m]->getX()-xc[0]) \
		  +(_ref[m]->getX()-yc[0])*(_ref[m]->getX()-yc[0]);
		e0 += (_align[m]->getY()-xc[1])*(_align[m]->getY()-xc[1]) \
		  +(_ref[m]->getY()-yc[1])*(_ref[m]->getY()-yc[1]);
		e0 += (_align[m]->getZ()-xc[2])*(_align[m]->getZ()-xc[2]) \
		  +(_ref[m]->getZ()-yc[2])*(_ref[m]->getZ()-yc[2]);
		r[0][0] += (_ref[m]->getX() - yc[0])*(_align[m]->getX() - xc[0]);
		r[0][1] += (_ref[m]->getX() - yc[0])*(_align[m]->getY() - xc[1]);
		r[0][2] += (_ref[m]->getX() - yc[0])*(_align[m]->getZ() - xc[2]);
		r[1][0] += (_ref[m]->getY() - yc[1])*(_align[m]->getX() - xc[0]);
		r[1][1] += (_ref[m]->getY() - yc[1])*(_align[m]->getY() - xc[1]);
		r[1][2] += (_ref[m]->getY() - yc[1])*(_align[m]->getZ() - xc[2]);
		r[2][0] += (_ref[m]->getZ() - yc[2])*(_align[m]->getX() - xc[0]);
		r[2][1] += (_ref[m]->getZ() - yc[2])*(_align[m]->getY() - xc[1]);
		r[2][2] += (_ref[m]->getZ() - yc[2])*(_align[m]->getZ() - xc[2]);
	}
	//compute determinat of matrix r
	det = r[0][0] * ( r[1][1]*r[2][2] - r[1][2]*r[2][1] )		\
	- r[0][1] * ( r[1][0]*r[2][2] - r[1][2]*r[2][0] )		\
	+ r[0][2] * ( r[1][0]*r[2][1] - r[1][1]*r[2][0] ); 
	sigma=det;
	
	//compute tras(r)*r
	m = 0;
	for (j=0; j<3; j++) {
		for (i=0; i<=j; i++) {            
			rr[m]=r[0][i]*r[0][j]+r[1][i]*r[1][j]+r[2][i]*r[2][j];
			m++;
		}
	}
	
	double spur=(rr[0]+rr[2]+rr[5]) / 3.0;
	double cof = (((((rr[2]*rr[5] - rr[4]*rr[4]) + rr[0]*rr[5])	\
		  - rr[3]*rr[3]) + rr[0]*rr[2]) - rr[1]*rr[1]) / 3.0;
	det = det*det; 
	
	for (i=0; i<3; i++){
		e[i]=spur;
	}
	
	if (spur>0) {
		d = spur*spur;
		h = d - cof;
		g = (spur*cof - det)/2.0 - spur*h;

		if (h>0) {
			sqrth = sqrt(h);
			d = h*h*h - g*g;
			if(d<0.0) d=0.0;
			d = atan2( sqrt(d), -g ) / 3.0;			
			cth = sqrth * cos(d);
			sth = sqrth*sqrt3*sin(d);
			e[0]= (spur + cth) + cth;
			e[1]= (spur - cth) + sth;            
			e[2]= (spur - cth) - sth;

			if (mode!=0) {//compute a                
				for (l=0; l<3; l=l+2) {
					d = e[l];  
					ss[0] = (d-rr[2]) * (d-rr[5])  - rr[4]*rr[4];
					ss[1] = (d-rr[5]) * rr[1]      + rr[3]*rr[4];
					ss[2] = (d-rr[0]) * (d-rr[5])  - rr[3]*rr[3];
					ss[3] = (d-rr[2]) * rr[3]      + rr[1]*rr[4];
					ss[4] = (d-rr[0]) * rr[4]      + rr[1]*rr[3];                
					ss[5] = (d-rr[0]) * (d-rr[2])  - rr[1]*rr[1]; 

					if (fabs(ss[0])<=epsilon) ss[0]=0.0;
					if (fabs(ss[1])<=epsilon) ss[1]=0.0;
					if (fabs(ss[2])<=epsilon) ss[2]=0.0;
					if (fabs(ss[3])<=epsilon) ss[3]=0.0;
					if (fabs(ss[4])<=epsilon) ss[4]=0.0;
					if (fabs(ss[5])<=epsilon) ss[5]=0.0;

					if (fabs(ss[0]) >= fabs(ss[2])) {
						j=0;                    
						if( fabs(ss[0]) < fabs(ss[5])){
							j = 2;
						}
					} else if ( fabs(ss[2]) >= fabs(ss[5]) ){
						j = 1;
					} else {
						j = 2;
					}

					d = 0.0;
					j = 3 * j;
					for (i=0; i<3; i++) {
						k=ip[i+j];
						a[i][l] = ss[k];
						d = d + ss[k]*ss[k];						
					} 


					//if( d > 0.0 ) d = 1.0 / sqrt(d);
					if (d > epsilon) d = 1.0 / sqrt(d);
					else d=0.0;
					for (i=0; i<3; i++) {
						a[i][l] = a[i][l] * d;
					}               
				}//for l

				d = a[0][0]*a[0][2] + a[1][0]*a[1][2] + a[2][0]*a[2][2];
				if ((e[0] - e[1]) > (e[1] - e[2])) {
					m1=2;
					m=0;
				} else {
					m1=0;
					m=2;                
				}
				p=0;
				for(i=0; i<3; i++){
					a[i][m1] = a[i][m1] - d*a[i][m];
					p = p + a[i][m1]*a[i][m1];
				}
				if (p <= tol) {
					p = 1.0;
					for (i=0; i<3; i++) {
						if (p < fabs(a[i][m])){
							continue;
						}
						p = fabs( a[i][m] );
						j = i;                    
					}
					k = ip2312[j];
					l = ip2312[j+1];
					p = sqrt( a[k][m]*a[k][m] + a[l][m]*a[l][m] ); 
					if (p > tol) {
						a[j][m1] = 0.0;
						a[k][m1] = -a[l][m]/p;
						a[l][m1] =  a[k][m]/p;                                                       
					} else {//goto 40
						a_failed=1;
					}     
				} else {//if p<=tol
					p = 1.0 / sqrt(p);
					for(i=0; i<3; i++){
						a[i][m1] = a[i][m1]*p;
					}                                  
				}//else p<=tol  
				if (a_failed!=1) {
					a[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
					a[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
					a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];       
				}                                   
			}//if(mode!=0)       
		}//h>0

		//compute b anyway
		if (mode!=0 && a_failed!=1) {//a is computed correctly
			//compute b
			for (l=0; l<2; l++) {
				d=0.0;
				for(i=0; i<3; i++){
					b[i][l] = r[i][0]*a[0][l] + r[i][1]*a[1][l] + r[i][2]*a[2][l];
					d = d + b[i][l]*b[i][l];
				}
				//if( d > 0 ) d = 1.0 / sqrt(d);
				if (d > epsilon) d = 1.0 / sqrt(d);
				else d=0.0;
				for (i=0; i<3; i++) {
					b[i][l] = b[i][l]*d;
				}                
			}            
			d = b[0][0]*b[0][1] + b[1][0]*b[1][1] + b[2][0]*b[2][1];
			p=0.0;

			for (i=0; i<3; i++) {
				b[i][1] = b[i][1] - d*b[i][0];
				p += b[i][1]*b[i][1];
			}

			if (p <= tol) {
				p = 1.0;
				for (i=0; i<3; i++) {
					if (p<fabs(b[i][0])) {
						continue;
					}
					p = fabs( b[i][0] );
					j=i;
				}
				k = ip2312[j];
				l = ip2312[j+1];
				p = sqrt( b[k][0]*b[k][0] + b[l][0]*b[l][0] ); 
				if (p > tol) {
					b[j][1] = 0.0;
					b[k][1] = -b[l][0]/p;
					b[l][1] =  b[k][0]/p;        
				} else {
					//goto 40
					b_failed=1;
				}                
			} else {//if( p <= tol )
				p = 1.0 / sqrt(p);
				for(i=0; i<3; i++){
					b[i][1]=b[i][1]*p;
				}
			}            
			if (b_failed!=1){
				b[0][2] = b[1][0]*b[2][1] - b[1][1]*b[2][0];
				b[1][2] = b[2][0]*b[0][1] - b[2][1]*b[0][0];
				b[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0]; 
				//compute u
				for (i=0; i<3; i++){
					for(j=0; j<3; j++){
						u[i][j] = b[i][0]*a[j][0] + b[i][1]*a[j][1]	\
								+ b[i][2]*a[j][2];
					}
				}
			}

			//compute t
			for(i=0; i<3; i++){
				t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1])	\
						- u[i][2]*xc[2];
			}            
		}//if(mode!=0 && a_failed!=1)
	} else {//spur>0, just compute t and errors
		//compute t
		for (i=0; i<3; i++) {
			t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1]) - u[i][2]*xc[2];
		}
	} //else spur>0 

	//compute rmsd
	for(i=0; i<3; i++){
		if( e[i] < 0 ) e[i] = 0;
		e[i] = sqrt( e[i] );           
	}            
	d = e[2];
	if( sigma < 0.0 ){
		d = - d;
	}
	d = (d + e[1]) + e[0];
	rms1 = (e0 - d) - d; 
	if( rms1 < 0.0 ) rms1 = 0.0;  

	rmsd=sqrt(rms1/(double)n);

	return rmsd;
}

void error(const string & msg)
{
	fprintf(stdout, "Error: %s.\n", msg.c_str());
	exit(-1);
}

void file2array(string _filename, vector<string> & lines) {
  FILE* ifp;
  int maxline = 1000;
  char *line, *tline;
  line = (char*) malloc(sizeof(char)*maxline);

  ifp = fopen(_filename.c_str(), "r");
  if (ifp == NULL) { error("unable to open file " + _filename); }

  while (fgets(line, maxline, ifp) != NULL) {
    if (line[strlen(line)-1] != '\n') { error("lines in file " + _filename + " are over " + toString(maxline) + "  long - increase max line limit and recompile."); }
    tline = trim(line);
    if (strlen(tline) > 0) { lines.push_back(line); }
    if (feof(ifp)) { break; }
  }
  if (ferror(ifp)) { error("an error occurred while reading file " + _filename); }
  free(line);
}

string fileBase(string fn) {
  if (fn.find_last_of(".") == string::npos) return fn;
  else return fn.substr(0, fn.find_last_of("."));
}

int findMostOvlpChainSet(const int & minGrpIdx, const int & curChainIdx, vector<vector<Chain*> > & chainGrpsRefined, map<Chain*, vector<Chain*> > & chainContacts, vector<Chain*> & curChainSet)
{
	curChainSet.clear();

	int i, j, maxi, maxj, maxNumOvlpChains, numGrpsLeft = chainGrpsRefined.size(), numOvlpChains, ttlNumChains;
	vector<bool> chosen(chainGrpsRefined.size(), false);
	curChainSet.push_back(chainGrpsRefined[minGrpIdx][curChainIdx]);
	chosen[minGrpIdx] = true;
	numGrpsLeft--;
	ttlNumChains = chainContacts[chainGrpsRefined[minGrpIdx][curChainIdx]].size() + 1; // 1 for central chain

	while (numGrpsLeft > 0)
	{
		maxNumOvlpChains = -1;
		maxi = -1;
		maxj = -1;
		for (i = 0; i < chainGrpsRefined.size(); i++)
		{
			if (chosen[i])
			{
				continue;
			}

			for (j = 0; j < chainGrpsRefined[i].size(); j++)
			{
				numOvlpChains = calcNumOvlpChains(chainGrpsRefined[i][j], curChainSet, chainContacts);
				if (numOvlpChains > maxNumOvlpChains)
				{
					maxNumOvlpChains = numOvlpChains;
					maxi = i;
					maxj = j;
				}
			}
		}

		curChainSet.push_back(chainGrpsRefined[maxi][maxj]);
		chosen[maxi] = true;
		numGrpsLeft--;
		ttlNumChains += chainContacts[chainGrpsRefined[maxi][maxj]].size() + 1; // 1 for central chain
		ttlNumChains -= maxNumOvlpChains; // overlapping chains have been counted twice
	}

	return ttlNumChains;
}

map<Chain*, AtomPointerVector> getChainCa(System & S)
{
	map<Chain*, AtomPointerVector> chainCa;
	int i, j;
	for (i = 0; i < S.chainSize(); i++)
	{
		Chain& chn = S.getChain(i);
		for (j = 0; j < chn.positionSize(); j++)
		{
			Residue& r = chn.getResidue(j);
			Atom& ca = r.getAtom("CA");
		  	chainCa[&chn].push_back(&ca);
		}
	}
	return chainCa;
}

map<Chain*, vector<double> > getChainCaCentroid(map<Chain*, AtomPointerVector> & chainCa)
{
	map<Chain*, vector<double> > chainCaCentroid;
	map<Chain*, AtomPointerVector>::iterator it;
	double x, y, z;
	int i;	
	for (it = chainCa.begin(); it != chainCa.end(); ++it)
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
		for (i = 0; i < it->second.size(); i++)
		{
			x += it->second[i]->getX();
			y += it->second[i]->getY();
			z += it->second[i]->getZ();
		}
		x = x / (it->second.size());
		chainCaCentroid[it->first].push_back(x);
		y = y / (it->second.size());
		chainCaCentroid[it->first].push_back(y);
		z = z / (it->second.size());
		chainCaCentroid[it->first].push_back(z);
	}

	return chainCaCentroid;
}

map<Chain*, vector<Chain*> > getChainContacts(System& S)
{
  	map<Chain*, vector<Chain*> > chainContacts;
	int i, j;
  	for (i = 0; i < S.chainSize(); i++)
  	{
		if (chainContacts.find(&(S.getChain(i))) == chainContacts.end())
		{
			// deal with the situation that there is only one chain
			chainContacts[&(S.getChain(i))] = vector<Chain*>();
		}
    	for (j = i + 1; j < S.chainSize(); j++)
    	{
      		if (areUnitsInContact(&(S.getChain(i).getAtomPointers()), &(S.getChain(j).getAtomPointers())))
      		{
        		chainContacts[&(S.getChain(i))].push_back(&(S.getChain(j)));
        		chainContacts[&(S.getChain(j))].push_back(&(S.getChain(i)));
      		}
    	}
  	}
  	return chainContacts;
}

map<Chain*, vector<string> > getChainSeq(System & S)
{
	map<Chain*, vector<string> > chainSeq;
	int i, j;
	for (i = 0; i < S.chainSize(); i++)
	{
		Chain& chn = S.getChain(i);
		for (j = 0; j < chn.positionSize(); j++)
		{
			chainSeq[&chn].push_back((chn.getResidue(j)).getResidueName());
		}
	}
	return chainSeq;
}

vector<CorrChains> getCorrChains(vector<Chain*> & contactA, vector<Chain*> & contactB, map<Chain*, vector<double> > & chnCaCentroid, const vector<double> & t, const vector<vector<double> > & u)
{
	vector<CorrChains> coChains;
	set<CorrChains, centDistComp> ccEns;
	set<CorrChains, centDistComp>::iterator it;
	double xa, ya, za, xb, yb, zb;
	int i, j;
	CorrChains cc;
	for (i = 0; i < contactA.size(); i++)
	{
		xa = t[0] + u[0][0] * chnCaCentroid[contactA[i]][0] + u[0][1] * chnCaCentroid[contactA[i]][1] + u[0][2] * chnCaCentroid[contactA[i]][2];
		ya = t[1] + u[1][0] * chnCaCentroid[contactA[i]][0] + u[1][1] * chnCaCentroid[contactA[i]][1] + u[1][2] * chnCaCentroid[contactA[i]][2];
		za = t[2] + u[2][0] * chnCaCentroid[contactA[i]][0] + u[2][1] * chnCaCentroid[contactA[i]][1] + u[2][2] * chnCaCentroid[contactA[i]][2];
		for (j = 0; j < contactB.size(); j++)
		{
			xb = chnCaCentroid[contactB[j]][0];
			yb = chnCaCentroid[contactB[j]][1];
			zb = chnCaCentroid[contactB[j]][2];

			cc.chnA = contactA[i];
			cc.chnB = contactB[j];
			cc.centDist = sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb));
			ccEns.insert(cc);
		}
	}
	while (!ccEns.empty())
	{
		it = ccEns.begin();
		coChains.push_back(*it);
		for (it = ccEns.begin(); it != ccEns.end(); ++it)
		{
			if ((it->chnA == coChains[coChains.size() - 1].chnA) || (it->chnB == coChains[coChains.size() - 1].chnB))
			{
				ccEns.erase(it);
			}
		}
	}	

	return coChains;
}

string getFileName(string fullpath){

  string name = fullpath;

  // Assume linux paths '/' , sorry windows users!
  vector<string> paths = tokenize(name, "/");


  if (paths.size() > 0){
    name = paths[paths.size()-1];
  }

  int cut;
  if (name.find_last_of(".") > 1000 || name.find_last_of(".") <= 0){
	  cut = name.length();
  } else {
	  cut = name.length() - name.find_last_of(".");
  }

  name.erase(name.length() - cut);

  return name;
}

void gridPoint(Atom* a, CartesianPoint& c, double gs, int* ip, int* jp, int* kp) {
  CartesianPoint d = a->getCoor() - c;
  int i = int(d.getX() < 0 ? ceil(d.getX()/gs): floor(d.getX()/gs));
  int j = int(d.getY() < 0 ? ceil(d.getY()/gs): floor(d.getY()/gs));
  int k = int(d.getZ() < 0 ? ceil(d.getZ()/gs): floor(d.getZ()/gs));
  if (ip != NULL) *ip = i;
  if (jp != NULL) *jp = j;
  if (kp != NULL) *kp = k;
}

bool hasBreak(Residue & cr, Residue & nr)
{
	if (cr.atomExists("C"))
	{
		if (nr.atomExists("N"))
		{
			if (cr.getAtom("C").distance(nr.getAtom("N")) > 2.0)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		if (nr.atomExists("NT"))
		{
			if (cr.getAtom("C").distance(nr.getAtom("NT")) > 2.0)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	if (!((cr.getChainId() == nr.getChainId()) && ((cr.getResidueNumber() == nr.getResidueNumber() - 1) || (cr.getResidueNumber() == nr.getResidueNumber()))))
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool isProtein(const vector<string> & legalAA, const string & resname)
{
	for(int i = 0; i < legalAA.size(); i++)
	{
		if (0 == legalAA[i].compare(resname))
		{
			return true;
		}
	}
	return false;
}

int main(int argc, char *argv[])
{
	int i;
	fstream ofs;
	FILE* ifp;
	CreateOptions copts;
	parseCommandLine(argc, argv, copts);

	if ("query" == copts.getPdsType())
	{
		for (i = 0; i < copts.getPdbFiles().size(); i++)
		{
			QueryStruct qs(copts.getPdbFiles()[i], copts.getLegalAA(), copts.getDistCut());
			openFileCPP(ofs, copts.getPdsFiles()[i], copts.getBinary() ? (ios::out | ios::trunc | ios::binary) : (ios::out | ios::trunc));
			qs.writeQueryStructFile(ofs, copts);
			ofs.close();
			openFileC(ifp, copts.getPdsFiles()[i].c_str(), copts.getBinary() ? "rb" : "r");
			qs.checkQueryStructFile(ifp, copts.getPdsFiles()[i], copts);
			fclose(ifp);
			if (copts.getPostPdbFiles().size() > 0)
			{
				System& sys = qs.getProteinSys();
				sys.writePdb(copts.getPostPdbFiles()[i]);
			}
		}
	}
	else
	{
		if ("target" == copts.getPdsType())
		{
			for (i = 0; i < copts.getPdbFiles().size(); i++)
			{
				TargetStruct ts(copts.getPdbFiles()[i], copts.getLegalAA(), copts.getNonRed(), copts.getSeqIdenThresh(), copts.getOptRmsdThresh(), copts.getContRmsdThresh(), copts.getNonRedPdb(), copts.getPhiStep(), copts.getPsiStep(), copts.getDistCut(), copts.getDistStep());
				openFileCPP(ofs, copts.getPdsFiles()[i], copts.getBinary() ? (ios::out | ios::trunc | ios::binary) : (ios::out | ios::trunc));
				ts.writeTargetStructFile(ofs, copts);
				ofs.close();
				openFileC(ifp, copts.getPdsFiles()[i].c_str(), copts.getBinary() ? "rb" : "r");
				ts.checkTargetStructFile(ifp, copts.getPdsFiles()[i], copts);
				fclose(ifp);
				if (copts.getPostPdbFiles().size() > 0)
				{
					System& sys = ts.getProteinSys();
					sys.writePdb(copts.getPostPdbFiles()[i]);
				}
			}
		}
		else
		{
			usage();
			error("bad output PDS type value");
		}
	}

	return 0;
}

void mergeChainGroups(Chain* chnA, Chain* chnB, map<Chain*, int> & chainGrpIdx, vector<vector<Chain*> > & chainGrps)
{
	if (chainGrpIdx[chnA] == chainGrpIdx[chnB])
	{
		return;
	}
	ASSERT((chainGrps[chainGrpIdx[chnA]].size() > 0) && (chainGrps[chainGrpIdx[chnB]].size() > 0), "number of chains in a group should be > 0");
	int i, idxA = chainGrpIdx[chnA], idxB = chainGrpIdx[chnB];
	vector<Chain*>::iterator it, beg, end;
	if (chainGrps[idxA].size() < chainGrps[idxB].size())
	{		
		it = chainGrps[idxB].end();
		beg = chainGrps[idxA].begin();
		end = chainGrps[idxA].end();
		chainGrps[idxB].insert(it, beg, end);
		for (i = 0; i < chainGrps[idxA].size(); i++)
		{
			chainGrpIdx[chainGrps[idxA][i]] = idxB;
		}
		chainGrps[idxA].clear();
	}
	else
	{		
		it = chainGrps[idxA].end();
		beg = chainGrps[idxB].begin();
		end = chainGrps[idxB].end();
		chainGrps[idxA].insert(it, beg, end);
		for (i = 0; i < chainGrps[idxB].size(); i++)
		{
			chainGrpIdx[chainGrps[idxB][i]] = idxA;
		}
		chainGrps[idxB].clear();
	}
}

void openFileC (FILE* & fp, const char* fname, const char* mode)
{
	fp = fopen(fname, mode);
	ASSERT(fp != NULL, "could not open file " + string(fname));
}

void openFileCPP(fstream & fs, const string & fname, const ios_base::openmode & mode)
{
	fs.open(fname.c_str(), mode);
	ASSERT(fs.is_open(), "could not open file " + fname);
}

string optionUsage(string opt, string mes, int w, int p1, int p2) {
  // first print the name of the option
  string text(p1, ' ');
  text += opt;
  if (p2 > text.size()) text += string(p2 - text.size(), ' ');

  // next print the description text
  int i = 0, k, L = text.size(), n;
  while (i < mes.size()) {
    k = mes.find_first_of(" ", i);
    if (k == string::npos) k = mes.size();
    n = k - i;
    if ((L + n >= w) && (L > 0)) { text += "\n" + string(p2, ' '); L = p2; }
    text += mes.substr(i, n) + " ";
    L += n + 1;
    i = k+1;
  }
  return text;
}

string pad(string str, int len)
{
	if (str.size() == len)
  	{
		return str;
  	}
  	else
  	{
		if (str.size() < len)
		{
    		for (int i = 0; i < len - str.size(); i++)
    		{
				str += " ";
    		}
    		return str;
  		}
		else
		{
    		return str.substr(0, len);
		}
  	}
}

void parseCommandLine(int argc, char** argv, CreateOptions & copts)
{
	map<string, bool> spec;

	while (1)
	{
		int oind = 0;
    	static struct option opts[] =
    	{
			{"pdb", 1, 0, 1},
			{"pdbList", 1, 0, 2},
			{"pds", 1, 0, 3},
			{"pdsList", 1, 0, 4},
			{"phiStep", 1, 0, 5},
	  		{"psiStep", 1, 0, 6},
	  		{"dCut", 1, 0, 7},
	  		{"dStep", 1, 0, 8},
			{"type", 1, 0, 9},
			{"opdb", 1, 0, 10},
			{"opdbList", 1, 0, 11},
			{"nr", 0, 0, 12},
			{"gRMSD", 1, 0, 13},
			{"lRMSD", 1, 0, 14},
			{"seqID", 1, 0, 15},
			{"nrPDB", 0, 0, 16},
			{"unnaturalAA", 0, 0, 17},
			{0, 0, 0, 0}
    	};

    	int c = getopt_long (argc, argv, "", opts, &oind);
    	if (c == -1)
    	{
			break;
    	}

    	switch (c) {
			case 1:				
				copts.setPdbFile(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 2:
				copts.setPdbFiles(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 3:
				copts.setPdsFile(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 4:
				copts.setPdsFiles(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;
			
      		case 5:
			  	copts.setPhiStep(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

      		case 6:
		        copts.setPsiStep(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

      		case 7:
		        copts.setDistCut(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

      		case 8:
			  	copts.setDistStep(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 9:
				copts.setPdsType(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;
			
			case 10:
				copts.setPostPdbFile(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 11:
				copts.setPostPdbFiles(string(optarg));
				spec[string(opts[oind].name)] = true;
				break;

			case 12:
				copts.setNonRed(true);
				spec[string(opts[oind].name)] = true;
				break;

			case 13:
		        copts.setContRmsdThresh(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 14:
		        copts.setOptRmsdThresh(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 15:
		        copts.setSeqIdenThresh(optarg);
		        spec[string(opts[oind].name)] = true;
		        break;

			case 16:
				copts.setNonRedPdb(true);
				spec[string(opts[oind].name)] = true;
				break;

			case 17:
				copts.setUnnaturalAA();
				spec[string(opts[oind].name)] = true;
				break;
			
      		case '?':
				usage();
				exit(-1);

      		default:
        		printf ("?? getopt returned character code %d ??\n", c);
				usage();
				exit(-1);
    	}
	}

  	// make sure all required options have been specified
  	if (!(((spec.find(string("pdb")) != spec.end()) || (spec.find(string("pdbList")) != spec.end()))
		&& (spec.find(string("type")) != spec.end())))
  	{
		usage();
		error("not all required options specified");
  	}

	if (optind < argc)
	{
		usage();
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
    	{
			printf ("%s ", argv[optind++]);
    	}
		printf ("\n");
		exit(-1);
	}

  	int i;
	if (0 == copts.getPdsFiles().size())
	{
		for (i = 0; i < copts.getPdbFiles().size(); i++)
		{
			copts.setPdsFile(fileBase(copts.getPdbFiles()[i]) + copts.getFileExt(), false);
		}
	}
	ASSERT(copts.getPdsFiles().size() == copts.getPdbFiles().size(), "numbers of input PDB files and output PDS files not consistent");
	if (copts.getPostPdbFiles().size() > 0)
  	{
  		ASSERT(copts.getPostPdbFiles().size() == copts.getPdbFiles().size(), "numbers of input PDB files and output post-processed PDB files not consistent");
		for (i = 0; i < copts.getPostPdbFiles().size(); i++)
		{
			copts.getPostPdbFiles()[i] = fileBase(copts.getPostPdbFiles()[i]) + ".pdb";
		}
	}

#if defined(DEBUG_CMDLINE)
	cout << "output binary file? " << copts.getBinary() << endl;
	cout << "contact RMSD threshold: " << copts.getContRmsdThresh() << endl;
	cout << "distance cutoff: " << copts.getDistCut() << endl;
	cout << "distance bin size: " << copts.getDistStep() << endl;
	cout << "output PDS file extension: " << copts.getFileExt() << endl;
	cout << "optimal RMSD threshold: " << copts.getOptRmsdThresh() << endl;
	cout << "phi bin size: " << copts.getPhiStep() << endl;
	cout << "psi bin size: " << copts.getPsiStep() << endl;
	cout << "separator: " << copts.getWordSep() << endl;
	cout << "sequence identity threshold: " << copts.getSeqIdenThresh() << endl;	
	cout << "non-redundant? " << copts.getNonRed() << endl;
	cout << "output non-redundant neighbor PDBs? " << copts.getNonRedPdb() << endl;
	cout << "terminator: " << copts.getWordTer() << endl;
	cout << "output PDS type: " << copts.getPdsType() << endl;
	exit(-1);
#endif
}

vector<string> tokenize(const std::string & _input, const std::string & _delimiter, bool _allowEmtpy){
	vector<string> results;

	if (_input == "") {
		if (_allowEmtpy) {
			results.push_back(_input);
		}
		return results;
	}
	
	if (_allowEmtpy) {
		size_t prePos = 0;
		size_t pos  = _input.find(_delimiter);
		unsigned int delimiterSize = _delimiter.size();
		string left = _input, right;

		while (pos != std::string::npos) {
			results.push_back(left.substr(prePos, pos));
			if( pos + delimiterSize <= left.size() ) {
				left = left.substr(pos + delimiterSize, left.size() );
			} else {
				left = "";
			}
			pos  = left.find(_delimiter);
		}

		results.push_back(left);
	} else {
		int start  = _input.find_first_not_of(_delimiter);
		int end    = 0;
		string cur = _input;


		while (start != std::string::npos){
			end    = _input.find_first_of(_delimiter, start);
			results.push_back(_input.substr(start, end-start));
			start  = _input.find_first_not_of(_delimiter, end);
		}
	}
	return results;

}

char* trim(char* str) {
  char* nptr; int i;
  for (i = 0; i < strlen(str); i++) {
    if ((str[i] != '\n') && (str[i] != '\t') && (str[i] != ' ')) { break; }
  }
  nptr = str + i;
  if (i == strlen(str)) { return nptr; }

  for (i = strlen(str)-1; i >= 0; i--) {
    if ((str[i] != '\n') && (str[i] != '\t') && (str[i] != ' ')) {
      str[i+1] = '\0';
      break;
    }
  }
  return nptr;
}

void usage()
{
	int w = 80, p1 = 3, p2 = p1 + 16; // the length of the longest option name + 3
	cout << optionUsage("", "Produces PDS file(s) from given PDB file(s). Different kinds of PDS files will be written depending on whether the input structures are meant as queries (i.e., --type query) or targets to be searched against (--type target).", w, p1, p2) << endl;
	cout << optionUsage("--type", "whether input PDB structures are meant as future queries (query) or target structures (target).", w, p1, p2) << endl;
	cout << optionUsage("--pdb", "input PDB file.", w, p1, p2) << endl;
	cout << optionUsage("--pdbList", "list of input PDB files, one per line.", w, p1, p2) << endl;
	cout << optionUsage("--pds", "optional: output PDS file name for input PDB file. By default, takes base name of PDB file and appends \".pds\".", w, p1, p2) << endl;
	cout << optionUsage("--pdsList", "optional: a file with a list of output PDS file names (one per input PDB file).", w, p1, p2) << endl;
	cout << optionUsage("--opdb", "optional: output post-processed PDB file (useful for keeping track of how all PDB weirdnesses got parsed).", w, p1, p2) << endl;
	cout << optionUsage("--opdbList", "optional: a file with a list of file names for post-processed PDBs (one per input PDB file).", w, p1, p2) << endl;
	cout << optionUsage("--nr", "optional: if specified, output a PDS file corresponding to the non-redundant part of the structure (defined in terms of non-redundant chains and interfaces; see documentation), instead of the entire input structure.", w, p1, p2) << endl;
	cout << optionUsage("--unnaturalAA", "optional: admit unnatural AA's in addition to the 20 natural ones: CSO, HIP, PTR, SEC, SEP, TPO. By default, all residues that are not the 20 natural amino acids will be ignored.", w, p1, p2) << endl;
	cout << optionUsage("--dCut", "optional: upper limit for tabulating CA-CA distances (must be identical between queries and targets); default is 25.0.", w, p1, p2) << endl;
	cout << optionUsage("--dStep", "optional: bin size for tabulating CA-CA distances (only matters for target structures); default is 5.0.", w, p1, p2) << endl;
	cout << optionUsage("--phiStep", "optional: bin size for tabulating phi angles (only matters for target structures); default is 10.0.", w, p1, p2) << endl;
	cout << optionUsage("--psiStep", "optional: bin size for tabulating psi angles (only matters for target structures); default is 10.0.", w, p1, p2) << endl;
	cout << optionUsage("--seqID", "optional: chain-to-chain sequence identity threshold for considering two chains to be \"identical\" with --nr; default is 0.9.", w, p1, p2) << endl;
	cout << optionUsage("--lRMSD", "optional: chain-to-chain RMSD threshold for considering two chains to be \"identical\" with --nr; default is 1.0.", w, p1, p2) << endl;
	cout << optionUsage("--gRMSD", "optional: RMSD threshold for neighboring chains upon superimposing the central chain, used with --nr; default is 2.0.", w, p1, p2) << endl;
//	cout << optionUsage("--nrPDB", "optional: if specified, dump PDB file of each unique neighborhood for large structures. Only works when --nr has been given.", w, p1, p2) << endl;
	
//	cout << optionUsage("--binary", "optional: output binary file? 1 for binary (default), 0 for text.", w, p1, p2) << endl;	  
}

template <class T>
void writeDatum (fstream & ofs, T val, const bool & bin)
{
  if (bin)
  {
    ofs.write((char*) &val, sizeof(T));
  }
  else
  {
    ofs << val;
  }
}

void writeString (fstream & ofs, const string & str, const bool & bin)
{
  if (bin)
  {
    ofs.write(str.c_str(), sizeof(char) * str.size());
  }
  else
  {
    ofs << str;
  }
}
