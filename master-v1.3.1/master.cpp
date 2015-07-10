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
// POST V1.1
#include "CenResGrp.h"
#include "CenterResidue.h"
#include "Common.h"
#include "Match.h"
#include "ProteinStruct.h"
#include "QueryStruct.h"
#include "SearchOptions.h"
#include "SearchResults.h"
#include "TargetStruct.h"

void auxSearchByDistDistrBB(SearchResults &, CenResGrp &, QueryStruct &, TargetStruct &, const int &, const double &, const double &, const int &, const double &, const int &, const double &, const double &);
void auxSearchByDistDistrCA(SearchResults &, CenResGrp &, QueryStruct &, TargetStruct &, const int &, const double &, const double &, const int &, const double &, const int &, const double &, const double &);
void beforeOut(SearchResults &, const string &, const bool &, const string &, const string &, const string &, const bool &, vector<int>* &, const vector<vector<int> > &);
double calcMaxDistDev(const vector<vector<double> > &, const vector<int> &, const vector<vector<double> > &, const vector<int> &);
double calcRmsdKabsch(const vector<vector<double> > &, const vector<int> &, const vector<vector<double> > &, const vector<int> &, const int &, double*, double**);
double distDevCut(const int &, const double &, const int &, const double &, const int &, const double &);
void extractSeq(string &, const Match &, char*, const string &, const vector<vector<int> > &);
void outDistDevZscore(SearchResults &, const bool &, const double &);
void outMatch(SearchResults &, const string &);
void outSeq(SearchResults &, const string &, const string &, const vector<vector<int> > &);
void parseCmdLine(int, char**, SearchOptions &);
void renameStruct(SearchResults &, vector<int>* &, const string &, const string &);
double rmsdCut(const int &, const double &, const int &, const double &, const int &, const int &, const double &);
void searchByDistDistr(SearchResults &, const bool &, const double &, const double &, const int &, const double &, const int &, const double &, const double &, const vector<vector<int> > &);
void transformStructByRMSD(vector<string> &, const vector<string> &, const vector<int> &, double*, double**);
void usage();

void auxSearchByDistDistrBB(SearchResults & mlist, CenResGrp & cres, QueryStruct & qs, TargetStruct & ts, const int & tsidx, const double & phieps, const double & psieps, const int & dmode, const double & deps, const int & rmode, const double & rthresh, const double & tune)
{
	int b, c, cri, e, i, j, k, gl, qri, tri;
	double qd, td, rmsd;
	vector<int> qbbai, tbbai;
	bool chose, f;
	Match m(&mlist, tsidx);
	double rcut, *t, **u;
	t = new double [3];
	u = new double* [3];
	for (i = 0; i < 3; ++i)
	{
		u[i] = new double [3];
	}

	rcut = mlist.full()? mlist.worstRmsd(): rthresh;
	
	for (i = 0; i < qs.getNumSeg(); i++)
	{
		cres[i].setCurMatchCandIdx(-1);
		cres[i].initTsResFlag(ts.getNumRes());
		cres[i].initTsResRmsd(ts.getNumRes());
		ts.readDihedDistr(cres[i].getMatchCand(), cres[i].getPhi(), cres[i].getPsi(), phieps, psieps);
		b = cres[i].getResIdx() - cres[i].getBegRes();
		e = (ts.getNumRes() - 1) - (cres[i].getEndRes() - cres[i].getResIdx());
		for (j = 0; j < cres[i].getMatchCand().size(); j++)
		{
			if ((cres[i].getMatchCand(j) < b) || (e < cres[i].getMatchCand(j)))
			{
				continue;
			}
			cres[i].setTsResFlag(cres[i].getMatchCand(j), MATCH_L1);
		}
		if (cres[i].getRank() >= 0.0)
		{
			cres[i].setRank(float(cres[i].getMatchCand().size()) / float(cres[i].getSegLen()));
		}
	}

	cres.sortByRank();

	for (i = 0; i < qs.getNumSeg(); i++)
	{
		cres[i].setSofarLen(cres[i].getSegLen() + ((i > 0) ? cres[i - 1].getSofarLen() : 0));
		// assume that the remaining part has perfect match
		// NOTE: BB and CA share the same RMSD cutoff, and factor 4 (i.e., number of BB atoms) makes no difference
		cres[i].setSofarRmsdCut(rmsdCut(rmode, rcut, qs.getNumRes(), 0.0, 0, cres[i].getSofarLen(), tune));
	}
	ASSERT(cres[i - 1].getSofarLen() == qs.getNumRes(), "number of residues not consistent for computing RMSD cutoff");

	cri = 0;
	while (1)
	{
		if (cri < 0)
		{
			break;
		}

		chose = false;
		for (c = cres[cri].getCurMatchCandIdx() + 1; c < cres[cri].getMatchCand().size(); c++)
		{
			if (NOT_MATCH == cres[cri].getTsResFlag(cres[cri].getMatchCand(c)))
			{
				continue;
			}
			
			cres[cri].setMatchCandBegRes(cres[cri].getMatchCand(c) - (cres[cri].getResIdx() - cres[cri].getBegRes()));
			cres[cri].setMatchCandEndRes(cres[cri].getMatchCand(c) + (cres[cri].getEndRes() - cres[cri].getResIdx()));
			if (MATCH_L1 == cres[cri].getTsResFlag(cres[cri].getMatchCand(c)))
			{
				f = false;
				for (i = cres[cri].getMatchCandBegRes(); i <= cres[cri].getMatchCandEndRes(); i++)
				{
					if (!ts.getFullBB(i))
					{
						f = true;
						break;
					}
				}
				if (f)
				{
					cres[cri].setTsResFlag(cres[cri].getMatchCand(c), NOT_MATCH);
					continue;
				}
				cres[cri].setTsResFlag(cres[cri].getMatchCand(c), MATCH_L2);
			}
			
			if (cri > 0)
			{
				// guarantee current match candidate satisfy current distance deviation cutoff
				qd = qs.calcDistBB(cres[cri - 1].getResIdx(), cres[cri].getResIdx());
				td = ts.calcDistBB(cres[cri - 1].getMatchCand(cres[cri - 1].getCurMatchCandIdx()), cres[cri].getMatchCand(c));
				if (fabs(td - qd) > cres[cri].getCurDistDevCut())
				{
					continue;
				}
				
				// guarantee matches of different segments not overlap each other
				f = false;
				for (i = 0; i < cri; i++)
				{
					if (!((cres[cri].getMatchCandEndRes() < cres[i].getMatchCandBegRes()) || (cres[cri].getMatchCandBegRes() > cres[i].getMatchCandEndRes())))
					{
						f = true;
						break;
					}
				}
				if (f)
				{
					continue;
				}

				// impose gap length constraints
				if (cres[cri].getGapLen().size() > 0)
				{
					gl = cres[cri].getMatchCandBegRes() - cres[cri - 1].getMatchCandEndRes() - 1;
					if (cres[cri].getGapLen().size() == 1)
					{
						if (gl != cres[cri].getGapLen()[0])
						{
							continue;
						}
					}
					else
					{
						if (cres[cri].getGapLen().size() == 2)
						{
							b = min(cres[cri].getGapLen()[0], cres[cri].getGapLen()[1]);
							e = max(cres[cri].getGapLen()[0], cres[cri].getGapLen()[1]);
							if (!((gl >= b) && (gl <= e)))
							{
								continue;
							}
						}
						else
						{
							error("invalid number of gap length constraints");
						}
					}
				}

				// only for heuristic distance deviation
				if ((dmode > 0) && (cri > 1))
				{
					f = false;
					for (i = 0; i < (cri - 1); i++)
					{
						qd = qs.calcDistBB(cres[i].getResIdx(), cres[cri].getResIdx());
						td = ts.calcDistBB(cres[i].getMatchCand(cres[i].getCurMatchCandIdx()), cres[cri].getMatchCand(c));
						if (fabs(td - qd) > distDevCut(dmode, 0.0, 0, 0.0, 0, deps))
						{
							f = true;
							break;
						}
					}
					if (f)
					{
						continue;
					}
				}
			}

			if (MATCH_L2 == cres[cri].getTsResFlag(cres[cri].getMatchCand(c)))
			{
				qbbai.clear();
				tbbai.clear();
				for (i = 0; (cres[cri].getBegRes() + i) <= cres[cri].getEndRes(); i++)
				{
					qri = cres[cri].getBegRes() + i;
					tri = cres[cri].getMatchCandBegRes() + i;
					for (j = 0; j < NUM_BBA; j++)
					{
						qbbai.push_back(qri * NUM_BBA + j);
						tbbai.push_back(tri * NUM_BBA + j);
					}
				}
				rmsd = calcRmsdKabsch(ts.getBBCoords(), tbbai, qs.getBBCoords(), qbbai, (cres.size() == 1) ? 1 : 0, t, u);
				if (rmsd > cres[cri].getSelfRmsdCut())
				{
					cres[cri].setTsResFlag(cres[cri].getMatchCand(c), NOT_MATCH);
					continue;
				}
				cres[cri].setTsResFlag(cres[cri].getMatchCand(c), MATCH_L3);
				cres[cri].setTsResRmsd(cres[cri].getMatchCand(c), rmsd);
			}

			if (cri > 0)
			{
				if (rmode < 2)
				{
					// NOTE: BB and CA share the same RMSD cutoff, and factor 4 (i.e., number of BB atoms) makes no difference
					if (cres[cri].getTsResRmsd(cres[cri].getMatchCand(c)) > rmsdCut(rmode, rcut, qs.getNumRes(), cres[cri - 1].getSofarRmsd(), cres[cri - 1].getSofarLen(), cres[cri].getSegLen(), tune))
					{
						continue;
					}
				}
				
				qbbai.clear();
				tbbai.clear();
				for (i = 0; i <= cri; i++)
				{
					for (j = 0; (cres[i].getBegRes() + j) <= cres[i].getEndRes(); j++)
					{
						qri = cres[i].getBegRes() + j;
						tri = cres[i].getMatchCandBegRes() + j;
						for (k = 0; k < NUM_BBA; k++)
						{
							qbbai.push_back(qri * NUM_BBA + k);
							tbbai.push_back(tri * NUM_BBA + k);
						}
					}
				}
				rmsd = calcRmsdKabsch(ts.getBBCoords(), tbbai, qs.getBBCoords(), qbbai, (cri == cres.size() - 1) ? 1 : 0, t, u);
				if (rmsd > cres[cri].getSofarRmsdCut())
				{
					continue;
				}
			}
			cres[cri].setSofarRmsd(rmsd);

			chose = true;
			cres[cri].setCurMatchCandIdx(c);
			cri++;
			break;
		}
		if (!chose)
		{
			cres[cri].setCurMatchCandIdx(-1);
			cri--;
			continue;
		}
		if (cri >= cres.size())
		{
			m.setRmsd(rmsd);
			m.setTranslation(t);
			m.setRotation(u);
			for (i = 0; i < cres.size(); i++)
			{
				m.setBegRes(cres[i].getSegIdx(), cres[i].getMatchCandBegRes());
			}
			m.setInsertTime(mlist.getNumMatFnd());
			mlist.insertMatch(m);

			// if really necessary, update RMSD bounds
			if (mlist.full() && (mlist.worstRmsd() < rcut))
			{
				rcut = mlist.worstRmsd();
				for (i = 0; i < cres.size(); i++)
				{
					// NOTE: BB and CA share the same RMSD cutoff, and factor 4 (i.e., number of BB atoms) makes no difference
					cres[i].setSelfRmsdCut(rmsdCut(rmode, rcut, qs.getNumRes(), 0.0, 0, cres[i].getSegLen(), tune));
					cres[i].setSofarRmsdCut(rmsdCut(rmode, rcut, qs.getNumRes(), 0.0, 0, cres[i].getSofarLen(), tune));
				}
			}
			
			cri--;
		}
		else
		{
			if (0 == dmode)
			{
				qbbai.clear();
				tbbai.clear();
				for (i = 0; i < cri; i++)
				{
					for (j = 0; (cres[i].getBegRes() + j) <= cres[i].getEndRes(); j++)
					{
						qri = cres[i].getBegRes() + j;
						tri = cres[i].getMatchCandBegRes() + j;
						for (k = 0; k < NUM_BBA; k++)
						{
							if (((cri - 1) == i) && (cres[i].getResIdx() == qri) && (CA_IDX == k))
							{
								continue;
							}
							qbbai.push_back(qri * NUM_BBA + k);
							tbbai.push_back(tri * NUM_BBA + k);
						}
					}
				}
				rmsd = calcRmsdKabsch(ts.getBBCoords(), tbbai, qs.getBBCoords(), qbbai, 0, t, u);
			}
			// NOTE: BB and CA do NOT share the same distance deviation cutoff
			cres[cri].setCurDistDevCut(distDevCut(dmode, rcut, qs.getNumRes() * NUM_BBA, rmsd, qbbai.size(), deps));
			ts.readDistDistr(cres[cri].getMatchCand(), cres[cri - 1].getMatchCand(cres[cri - 1].getCurMatchCandIdx()), qs.calcDistBB(cres[cri - 1].getResIdx(), cres[cri].getResIdx()), cres[cri].getCurDistDevCut());
		}		
	}	

	for (i = 0; i < 3; ++i)
	{
		delete [] u[i];
	}
	delete [] u;
	delete [] t;
}

void auxSearchByDistDistrCA(SearchResults & mlist, CenResGrp & cres, QueryStruct & qs, TargetStruct & ts, const int & tsidx, const double & phieps, const double & psieps, const int & dmode, const double & deps, const int & rmode, const double & rthresh, const double & tune)
{
	int b, c, cri, e, i, j, k, gl;
	double qd, td, rmsd;
	vector<int> qri, tri;
	bool chose, f;
	Match m(&mlist, tsidx);
	double rcut, *t, **u;
	t = new double [3];
	u = new double* [3];
	for (i = 0; i < 3; ++i)
	{
		u[i] = new double [3];
	}

	rcut = mlist.full()? mlist.worstRmsd(): rthresh;
	
	for (i = 0; i < qs.getNumSeg(); i++)
	{
		cres[i].setCurMatchCandIdx(-1);
		cres[i].initTsResFlag(ts.getNumRes());
		cres[i].initTsResRmsd(ts.getNumRes());
		ts.readDihedDistr(cres[i].getMatchCand(), cres[i].getPhi(), cres[i].getPsi(), phieps, psieps);
		b = cres[i].getResIdx() - cres[i].getBegRes();
		e = (ts.getNumRes() - 1) - (cres[i].getEndRes() - cres[i].getResIdx());
		for (j = 0; j < cres[i].getMatchCand().size(); j++)
		{
			if ((cres[i].getMatchCand(j) < b) || (e < cres[i].getMatchCand(j)))
			{
				continue;
			}
			cres[i].setTsResFlag(cres[i].getMatchCand(j), MATCH_L1);
		}
		if (cres[i].getRank() >= 0.0)
		{
			cres[i].setRank(float(cres[i].getMatchCand().size()) / float(cres[i].getSegLen()));
		}
	}

	cres.sortByRank(); // NOTE: can this be done outside of the auxiliary function? Isn't the order of central resiudes the same for all searches?

	for (i = 0; i < qs.getNumSeg(); i++)
	{
		cres[i].setSofarLen(cres[i].getSegLen() + ((i > 0) ? cres[i - 1].getSofarLen() : 0));
		// assume that the remaining part has perfect match
		cres[i].setSofarRmsdCut(rmsdCut(rmode, rcut, qs.getNumRes(), 0.0, 0, cres[i].getSofarLen(), tune));
	}
	ASSERT(cres[i - 1].getSofarLen() == qs.getNumRes(), "number of residues not consistent for computing appropriate RMSD threshold");

	cri = 0;
	while (1)
	{
		if (cri < 0)
		{
			break;
		}

		chose = false;
		for (c = cres[cri].getCurMatchCandIdx() + 1; c < cres[cri].getMatchCand().size(); c++)
		{
			if (NOT_MATCH == cres[cri].getTsResFlag(cres[cri].getMatchCand(c)))
			{
				continue;
			}
			
			cres[cri].setMatchCandBegRes(cres[cri].getMatchCand(c) - (cres[cri].getResIdx() - cres[cri].getBegRes()));
			cres[cri].setMatchCandEndRes(cres[cri].getMatchCand(c) + (cres[cri].getEndRes() - cres[cri].getResIdx()));			
			if (cri > 0)
			{
				// guarantee current match candidate satisfy current distance deviation cutoff
				qd = qs.calcDistCA(cres[cri - 1].getResIdx(), cres[cri].getResIdx());
				td = ts.calcDistCA(cres[cri - 1].getMatchCand(cres[cri - 1].getCurMatchCandIdx()), cres[cri].getMatchCand(c));
				if (fabs(td - qd) > cres[cri].getCurDistDevCut())
				{
					continue;
				}
				
				// guarantee matches of different segments not overlap each other
				f = false;
				for (i = 0; i < cri; i++)
				{
					if (!((cres[cri].getMatchCandEndRes() < cres[i].getMatchCandBegRes()) || (cres[cri].getMatchCandBegRes() > cres[i].getMatchCandEndRes())))
					{
						f = true;
						break;
					}
				}
				if (f)
				{
					continue;
				}

				// impose gap length constraints
				if (cres[cri].getGapLen().size() > 0)
				{
					gl = cres[cri].getMatchCandBegRes() - cres[cri - 1].getMatchCandEndRes() - 1;
					if (cres[cri].getGapLen().size() == 1)
					{
						if (gl != cres[cri].getGapLen()[0])
						{
							continue;
						}
					}
					else
					{
						if (cres[cri].getGapLen().size() == 2)
						{
							b = min(cres[cri].getGapLen()[0], cres[cri].getGapLen()[1]);
							e = max(cres[cri].getGapLen()[0], cres[cri].getGapLen()[1]);
							if (!((gl >= b) && (gl <= e)))
							{
								continue;
							}
						}
						else
						{
							error("invalid number of gap length constraints");
						}
					}
				}

				// only for heuristic distance deviation
				if ((dmode > 0) && (cri > 1))
				{
					f = false;
					for (i = 0; i < (cri - 1); i++)
					{
						qd = qs.calcDistCA(cres[i].getResIdx(), cres[cri].getResIdx());
						td = ts.calcDistCA(cres[i].getMatchCand(cres[i].getCurMatchCandIdx()), cres[cri].getMatchCand(c));
						if (fabs(td - qd) > distDevCut(dmode, 0.0, 0, 0.0, 0, deps))
						{
							f = true;
							break;
						}
					}
					if (f)
					{
						continue;
					}
				}
			}

			if (MATCH_L1 == cres[cri].getTsResFlag(cres[cri].getMatchCand(c)))
			{
				qri.clear();
				tri.clear();
				for (i = 0; (cres[cri].getBegRes() + i) <= cres[cri].getEndRes(); i++)
				{
					qri.push_back(cres[cri].getBegRes() + i);
					tri.push_back(cres[cri].getMatchCandBegRes() + i);
				}
				rmsd = calcRmsdKabsch(ts.getCACoords(), tri, qs.getCACoords(), qri, (cres.size() == 1) ? 1 : 0, t, u);
				if (rmsd > cres[cri].getSelfRmsdCut())
				{
					cres[cri].setTsResFlag(cres[cri].getMatchCand(c), NOT_MATCH);
					continue;
				}
				cres[cri].setTsResFlag(cres[cri].getMatchCand(c), MATCH_L2);
				cres[cri].setTsResRmsd(cres[cri].getMatchCand(c), rmsd);
			}

			if (cri > 0)
			{
				if (rmode < 2)
				{
					if (cres[cri].getTsResRmsd(cres[cri].getMatchCand(c)) > rmsdCut(rmode, rcut, qs.getNumRes(), cres[cri - 1].getSofarRmsd(), cres[cri - 1].getSofarLen(), cres[cri].getSegLen(), tune))
					{
						continue;
					}
				}
				
				qri.clear();
				tri.clear();
				for (i = 0; i <= cri; i++)
				{
					for (j = 0; (cres[i].getBegRes() + j) <= cres[i].getEndRes(); j++)
					{
						qri.push_back(cres[i].getBegRes() + j);
						tri.push_back(cres[i].getMatchCandBegRes() + j);
					}
				}
				rmsd = calcRmsdKabsch(ts.getCACoords(), tri, qs.getCACoords(), qri, (cri == cres.size() - 1) ? 1 : 0, t, u);
				if (rmsd > cres[cri].getSofarRmsdCut())
				{
					continue;
				}
			}
			cres[cri].setSofarRmsd(rmsd);

			chose = true;
			cres[cri].setCurMatchCandIdx(c);
			cri++;
			break;
		}
		if (!chose)
		{
			cres[cri].setCurMatchCandIdx(-1);
			cri--;
			continue;
		}
		if (cri >= cres.size())
		{
			m.setRmsd(rmsd);
			m.setTranslation(t);
			m.setRotation(u);
			for (i = 0; i < cres.size(); i++)
			{
				m.setBegRes(cres[i].getSegIdx(), cres[i].getMatchCandBegRes());
			}
			m.setInsertTime(mlist.getNumMatFnd());
			mlist.insertMatch(m);

			// if really necessary, update RMSD bounds
			if (mlist.full() && (mlist.worstRmsd() < rcut))
			{
				rcut = mlist.worstRmsd();
				for (i = 0; i < cres.size(); i++)
				{
					cres[i].setSelfRmsdCut(rmsdCut(rmode, rcut, qs.getNumRes(), 0.0, 0, cres[i].getSegLen(), tune));
					cres[i].setSofarRmsdCut(rmsdCut(rmode, rcut, qs.getNumRes(), 0.0, 0, cres[i].getSofarLen(), tune));
				}
			}
			
			cri--;
		}
		else
		{
			if (dmode == 0)
			{
				qri.clear();
				tri.clear();
				for (i = 0; i < cri; i++)
				{
					for (j = 0; (cres[i].getBegRes() + j) <= cres[i].getEndRes(); j++)
					{
						if ((i == (cri - 1)) && ((cres[i].getBegRes() + j) == cres[i].getResIdx()))
						{
							continue;
						}
						qri.push_back(cres[i].getBegRes() + j);
						tri.push_back(cres[i].getMatchCandBegRes() + j);
					}
				}
				rmsd = calcRmsdKabsch(ts.getCACoords(), tri, qs.getCACoords(), qri, 0, t, u);
			}
			cres[cri].setCurDistDevCut(distDevCut(dmode, rcut, qs.getNumRes(), rmsd, qri.size(), deps));
			ts.readDistDistr(cres[cri].getMatchCand(), cres[cri - 1].getMatchCand(cres[cri - 1].getCurMatchCandIdx()), qs.calcDistCA(cres[cri - 1].getResIdx(), cres[cri].getResIdx()), cres[cri].getCurDistDevCut());
		}		
	}	

	for (i = 0; i < 3; ++i)
	{
		delete [] u[i];
	}
	delete [] u;
	delete [] t;
}

void beforeOut(SearchResults & mlist, const string & mifname, const bool & bbrmsd, const string & sofname, const string & sodname, const string & otype, const bool & ddzscore, vector<int>* & oldidx, const vector<vector<int> > & gaplen)
{
	if (mlist.empty())
	{
		return;
	}
	mlist.setPointerVector();

	bool f = false;
	if ((!mifname.empty()) || (!sofname.empty()) || (!sodname.empty()) || ddzscore)
	{
		char fname[500];
		int i, j, k, l;
		fstream ofs;
		vector<string> opdb;
		QueryStruct& qs = mlist.getQs();
		vector<int> qbbai, qcrcai, qri, tbbai, tcrcai, tri;
		double rmsd, *t, **u;
		vector<char*>& seq = mlist.getSeq();
		bool tbbaf;
		TargetStruct ts;
		
		t = new double [3];
		u = new double* [3];
		for (i = 0; i < 3; ++i)
		{
			u[i] = new double [3];
		}

		if (!mifname.empty())
		{
			if (bbrmsd)
			{
				qbbai.clear();
				for (i = 0; i < qs.getNumRes(); i++)
				{
					for (j = 0; j < NUM_BBA; j++)
					{
						qbbai.push_back(i * NUM_BBA + j);
					}
				}
			}
			else
			{
				qs.readCACoords();
				qri.clear();
				for (i = 0; i < qs.getNumRes(); i++)
				{
					qri.push_back(i);
				}
			}
		}

		if (!sofname.empty())
		{
			seq.assign(mlist.getNumTs(), NULL);
		}

		if (!sodname.empty())
		{
			oldidx = new vector<int>;
		}

		if (ddzscore && bbrmsd)
		{
			qcrcai.clear();
			for (i = 0; i < qs.getNumSeg(); i++)
			{
				qcrcai.push_back(qs.getCenRes(i) * NUM_BBA + CA_IDX);
			}
		}

		if (mlist.getMatchLimit() > 0)
		{
			mlist.sortByTarget();
		}
		mlist.initTsNumRes();

		for (i = 0; i < mlist.size(); i++)
		{
			if ((i == 0) || (mlist[i].getTsIdx() != mlist[i - 1].getTsIdx()))
			{
				ts.closeProteinStructFile();
				ts.readProteinStruct(mlist[i].getTsFile());

				if ((!mifname.empty()) || ddzscore)
				{
					if (bbrmsd)
					{
						ts.readBBCoords();
					}
					else
					{
						ts.readCACoords();
					}
				}
				
				if (!sofname.empty())
				{
					ts.readSeq(seq[mlist[i].getTsIdx()]);
				}

				if (!sodname.empty())
				{
					ts.initPdbInfo();
					if (otype.compare("full") == 0)
					{
						ts.readPdbInfo();
					}
				}

				mlist[i].setTsNumRes(ts.getNumRes());
			}

			if (!mifname.empty())
			{
				if (bbrmsd)
				{
					tbbaf = true;
					tbbai.clear();
					for (j = 0; j < mlist[i].getNumSeg(); j++)
					{
						for (k = mlist[i].getBegRes(j); k <= mlist[i].getEndRes(j); k++)
						{
							if (!ts.getFullBB(k))
							{
								tbbaf = false;
								break;
							}
							
							for (l = 0; l < NUM_BBA; l++)
							{
								tbbai.push_back(k * NUM_BBA + l);
							}
						}
						if (!tbbaf)
						{
							break;
						}
					}
					if (tbbaf)
					{
						mlist[i].setRmsd(calcRmsdKabsch(ts.getBBCoords(), tbbai, qs.getBBCoords(), qbbai, 1, t, u));
						mlist[i].setTranslation(t);
						mlist[i].setRotation(u);
					}
					else
					{
						mlist[i].setRmsd(DBL_MAX);
					}
				}
				else
				{
					tri.clear();
					for (j = 0; j < mlist[i].getNumSeg(); j++)
					{
						for (k = mlist[i].getBegRes(j); k <= mlist[i].getEndRes(j); k++)
						{
							tri.push_back(k);
						}
					}
					mlist[i].setRmsd(calcRmsdKabsch(ts.getCACoords(), tri, qs.getCACoords(), qri, 1, t, u));
					mlist[i].setTranslation(t);
					mlist[i].setRotation(u);
				}
			}

			if ((!sodname.empty()) && (mlist[i].getRmsd() != DBL_MAX))
			{
				if (otype.compare("full") == 0)
				{
					tri.clear();
					for (j = 0; j < ts.getNumRes(); j++)
					{
						tri.push_back(j);
					}
				}
				else
				{
					if (otype.compare("wgap") == 0)
					{
						tri.clear();
						for (j = 0; j < mlist[i].getNumSeg(); j++)
						{
							for (k = mlist[i].getBegRes(j); k <= mlist[i].getEndRes(j); k++)
							{
								tri.push_back(k);
							}
							if (j < (mlist[i].getNumSeg() - 1))
							{
								if (gaplen[j].size() > 0)
								{
									for (k = mlist[i].getEndRes(j) + 1; k < mlist[i].getBegRes(j + 1); k++)
									{
										tri.push_back(k);
									}
								}
							}
						}
					}
					else
					{
						if (otype.compare("match") == 0)
						{
							if (!((!mifname.empty()) && (!bbrmsd)))
							{
								tri.clear();
								for (j = 0; j < mlist[i].getNumSeg(); j++)
								{
									for (k = mlist[i].getBegRes(j); k <= mlist[i].getEndRes(j); k++)
									{
										tri.push_back(k);
									}
								}
							}
						}
						else
						{
							error("unidentified output type");
						}
					}
					ts.readPdbInfo(tri);
				}

				transformStructByRMSD(opdb, ts.getPdbInfo(), tri, mlist[i].getTranslation(), mlist[i].getRotation());
				
				sprintf(fname, "%s/tmp%d", sodname.c_str(), i + 1);
				openFileCPP(ofs, string(fname), (ios::out | ios::trunc));
				ofs << "REMARK " << mlist[i] << "\n";
				for (j = 0; j < opdb.size(); j++)
				{
					ofs << opdb[j];
				}
				ofs.close();
			}

			if (ddzscore)
			{
				if (bbrmsd)
				{
					tcrcai.clear();
					for (j = 0; j < mlist[i].getNumSeg(); j++)
					{
						k = mlist[i].getBegRes(j) + (qs.getCenRes(j) - (qs.getResBefBrk(j) + 1));
						tcrcai.push_back(k * NUM_BBA + CA_IDX);
					}
					mlist[i].setMaxDistDev(calcMaxDistDev(ts.getBBCoords(), tcrcai, qs.getBBCoords(), qcrcai));
				}
				else
				{
					tri.clear();
					for (j = 0; j < mlist[i].getNumSeg(); j++)
					{
						k = mlist[i].getBegRes(j) + (qs.getCenRes(j) - (qs.getResBefBrk(j) + 1));
						tri.push_back(k);
					}
					mlist[i].setMaxDistDev(calcMaxDistDev(ts.getCACoords(), tri, qs.getCACoords(), qs.getCenRes()));
				}
			}
		}
		
		for (i = 0; i < 3; ++i)
		{
			delete [] u[i];
		}
		delete [] u;
		delete [] t;

		f = true;
	}
	else
	{
		if (mlist.getMatchLimit() <= 0)
		{
			f = true;
		}		
	}

	if (f)
	{
		mlist.sortByRmsd(oldidx);
	}
}

double calcMaxDistDev(const vector<vector<double> > & align, const vector<int> & aind, const vector<vector<double> > & ref, const vector<int> & rind)
{
	ASSERT(aind.size() == rind.size(), "number of center residues not consistent");
	
	int i, j;
	double ad, rd, dd, maxdd = -1.0;
	for (i = 0; i < rind.size(); i++)
	{
		for (j = i + 1; j < rind.size(); j++)
		{
			rd = sqrt((ref[rind[i]][0] - ref[rind[j]][0]) * (ref[rind[i]][0] - ref[rind[j]][0]) + (ref[rind[i]][1] - ref[rind[j]][1]) * (ref[rind[i]][1] - ref[rind[j]][1]) + (ref[rind[i]][2] - ref[rind[j]][2]) * (ref[rind[i]][2] - ref[rind[j]][2]));
			ad = sqrt((align[aind[i]][0] - align[aind[j]][0]) * (align[aind[i]][0] - align[aind[j]][0]) + (align[aind[i]][1] - align[aind[j]][1]) * (align[aind[i]][1] - align[aind[j]][1]) + (align[aind[i]][2] - align[aind[j]][2]) * (align[aind[i]][2] - align[aind[j]][2]));
			dd = fabs(ad - rd);
			if (dd > maxdd)
			{
				maxdd = dd;
			}
		}
	}

	return maxdd;
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
double calcRmsdKabsch(const vector<vector<double> > & _align, const vector<int> & _aind, const vector<vector<double> > & _ref, const vector<int> & _rind, const int & mode, double* t, double** u)
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
	
	int n=_rind.size();
	ASSERT((n >= 0) && (n == _aind.size()), "invalid protein length for calculating RMSD");

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
		xc[0] += _align[_aind[i]][0];
		xc[1] += _align[_aind[i]][1];
		xc[2] += _align[_aind[i]][2];
		
		yc[0] += _ref[_rind[i]][0];
		yc[1] += _ref[_rind[i]][1];
		yc[2] += _ref[_rind[i]][2];
	}
	for(i=0; i<3; i++){
		xc[i] = xc[i]/(double)n;
		yc[i] = yc[i]/(double)n;        
	}
	
	//compute e0 and matrix r
	for (m=0; m<n; m++) {
		e0 += (_align[_aind[m]][0]-xc[0])*(_align[_aind[m]][0]-xc[0]) \
		  +(_ref[_rind[m]][0]-yc[0])*(_ref[_rind[m]][0]-yc[0]);
		e0 += (_align[_aind[m]][1]-xc[1])*(_align[_aind[m]][1]-xc[1]) \
		  +(_ref[_rind[m]][1]-yc[1])*(_ref[_rind[m]][1]-yc[1]);
		e0 += (_align[_aind[m]][2]-xc[2])*(_align[_aind[m]][2]-xc[2]) \
		  +(_ref[_rind[m]][2]-yc[2])*(_ref[_rind[m]][2]-yc[2]);
		r[0][0] += (_ref[_rind[m]][0] - yc[0])*(_align[_aind[m]][0] - xc[0]);
		r[0][1] += (_ref[_rind[m]][0] - yc[0])*(_align[_aind[m]][1] - xc[1]);
		r[0][2] += (_ref[_rind[m]][0] - yc[0])*(_align[_aind[m]][2] - xc[2]);
		r[1][0] += (_ref[_rind[m]][1] - yc[1])*(_align[_aind[m]][0] - xc[0]);
		r[1][1] += (_ref[_rind[m]][1] - yc[1])*(_align[_aind[m]][1] - xc[1]);
		r[1][2] += (_ref[_rind[m]][1] - yc[1])*(_align[_aind[m]][2] - xc[2]);
		r[2][0] += (_ref[_rind[m]][2] - yc[2])*(_align[_aind[m]][0] - xc[0]);
		r[2][1] += (_ref[_rind[m]][2] - yc[2])*(_align[_aind[m]][1] - xc[1]);
		r[2][2] += (_ref[_rind[m]][2] - yc[2])*(_align[_aind[m]][2] - xc[2]);
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
						u[i][j] = b[i][0]*a[j][0] + b[i][1]*a[j][1]	+ b[i][2]*a[j][2];
					}
				}
			}

			//compute t
			for(i=0; i<3; i++){
				t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1]) - u[i][2]*xc[2];
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

double distDevCut(const int & mode, const double & rtot, const int & ltot, const double & rsofar, const int & lsofar, const double & deps)
{
	switch (mode) {
		case 1:
			return deps;
		case 0:
		default:			
			return sqrt(2.0 * ((rtot * rtot * double(ltot)) - (rsofar * rsofar * double(lsofar))));
	}
}

void extractSeq(string & aa, const Match & m, char* seq, const string & otype, const vector<vector<int> > & gaplen)
{
	aa.clear();
	
	int i, j, pos, seqlen;
	if (otype.compare("full") == 0)
	{
		seqlen = m.getTsNumRes() * LEN_AA_CODE;
		for (i = 0; i < seqlen; i += LEN_AA_CODE)
		{
			aa.append(seq + i, LEN_AA_CODE);
			aa.append(" ");
		}
	}
	else
	{
		if ((otype.compare("match") == 0) || (otype.compare("wgap") == 0))
		{
			for (i = 0; i < m.getNumSeg(); i++)
			{
				pos = m.getBegRes(i) * LEN_AA_CODE;
				seqlen = m.getSegLen(i) * LEN_AA_CODE;
				for (j = 0; j < seqlen; j += LEN_AA_CODE)
				{
					aa.append(seq + pos + j, LEN_AA_CODE);
					aa.append(" ");
				}
				// output gap sequence
				if ((otype.compare("wgap") == 0) && (i < (m.getNumSeg() - 1)))
				{
					if (gaplen[i].size() > 0)
					{
						pos = (m.getEndRes(i) + 1) * LEN_AA_CODE;
						seqlen = (m.getBegRes(i + 1) - m.getEndRes(i) - 1) * LEN_AA_CODE;
						aa.append("[");
						for (j = 0; j < seqlen; j += LEN_AA_CODE)
						{
							aa.append(seq + pos + j, LEN_AA_CODE);
							aa.append(" ");
						}
						aa[aa.size() - 1] = ']';
						aa.append(" ");
					}
				}
			}
		}
	}	
	aa[aa.size() - 1] = '\n';
}

int main(int argc, char* argv[])
{
	double elapse;
	const long int microPerSec = 1000000L;
	SearchResults mlist;
	vector<int> *oldidx = NULL;
	SearchOptions sopts;
	struct timeval start, end;
	clock_t t;
	
	parseCmdLine(argc, argv, sopts);

	// initialize match list
	mlist.setMatchLimit(sopts.getTopN());
	mlist.setQs(sopts.getQsFile());	
	if (sopts.getGapLen().size() > 0)
	{
		ASSERT(mlist.getQs().getNumSeg() == (sopts.getGapLen().size() + 1), "numbers of query segments and gap length constraints not compatible");
	}	

	if (!sopts.getMatchInFile().empty())
	{
		ASSERT(gettimeofday(&start, NULL) == 0, "could not obtain the parse start time");
		t = clock();

		if (sopts.getBBRmsd())
		{
			mlist.getQs().readBBCoords();
			ASSERT(mlist.getQs().fullBB(), "backbone atoms missing in file " + mlist.getQs().getFileName());
		}
		mlist.setMatchList(sopts.getMatchInFile(), sopts.getGapLen());
		
		t = clock() - t;
		ASSERT(gettimeofday(&end, NULL) == 0, "could not obtain the parse finish time");
		elapse = double((end.tv_sec - start.tv_sec) * microPerSec + end.tv_usec - start.tv_usec) / double(microPerSec);
		printf("Parse completed in %.2f seconds, and used %.2f seconds of CPU time.\n", elapse, float(t) / CLOCKS_PER_SEC);
	}
	else
	{
		ASSERT(gettimeofday(&start, NULL) == 0, "could not obtain the search start time");
		t = clock();

		mlist.setTsFiles(sopts.getTsFiles());
		// NOTE: no need to pass all entries of sopts separately
		searchByDistDistr(mlist, sopts.getBBRmsd(), sopts.getPhiDevCut(), sopts.getPsiDevCut(), sopts.getDistDevBoundMode(), sopts.getDistDevCut(), sopts.getRmsdBoundMode(), sopts.getRmsdCut(), sopts.getTuningParam(), sopts.getGapLen());

		t = clock() - t;
		ASSERT(gettimeofday(&end, NULL) == 0, "could not obtain the search finish time");
		elapse = double((end.tv_sec - start.tv_sec) * microPerSec + end.tv_usec - start.tv_usec) / double(microPerSec);
		printf("Search completed in %.2f seconds, and used %.2f seconds of CPU time.\n", elapse, float(t) / CLOCKS_PER_SEC);
	}
	ASSERT(gettimeofday(&start, NULL) == 0, "could not obtain the output start time");
	t = clock();

	beforeOut(mlist, sopts.getMatchInFile(), sopts.getBBRmsd(), sopts.getSeqOutFile(), sopts.getStructOutDir(),  sopts.getOutType(), (sopts.getDistDevZscore() && (mlist.getQs().getNumSeg() > 1)), oldidx, sopts.getGapLen());
	outMatch(mlist, sopts.getMatchOutFile());
	outSeq(mlist, sopts.getSeqOutFile(), sopts.getOutType(), sopts.getGapLen());
	renameStruct(mlist, oldidx, sopts.getStructOutDir(), sopts.getOutType());
	outDistDevZscore(mlist, (sopts.getDistDevZscore() && (mlist.getQs().getNumSeg() > 1)), sopts.getDistDevCut());
	
	t = clock() - t;
	ASSERT(gettimeofday(&end, NULL) == 0, "could not obtain the output finish time");
	elapse = double((end.tv_sec - start.tv_sec) * microPerSec + end.tv_usec - start.tv_usec) / double(microPerSec);
	printf("Output completed in %.2f seconds, and used %.2f seconds of CPU time.\n", elapse, float(t) / CLOCKS_PER_SEC);
	
	return 0;
}

void outDistDevZscore(SearchResults & mlist, const bool & ddzscore, const double & deps)
{
	if (!ddzscore)
	{
		return;
	}
	
	if (mlist.empty())
	{
		cout << "Warning: no match, no distance deviation Z-score.\n";
		return;
	}
	
	int i, j;
	double mean = 0.0, std = 0.0, zscore;

	for (i = 0; i < mlist.size(); ++i)
	{
		mean += mlist[i].getMaxDistDev();
	}
	mean = mean / double(mlist.size());

	for (i = 0; i < mlist.size(); i++)
	{
		std += (mlist[i].getMaxDistDev() - mean) * (mlist[i].getMaxDistDev() - mean);
	}
	std = sqrt(std / double(mlist.size()));

	zscore = (deps - mean) / std;
	fprintf(stdout, "Distance deviation Z-score: %f.\n", zscore);
}

void outMatch(SearchResults & mlist, const string & mofname)
{
	if (mlist.empty())
	{
		cout << "Warning: no match found.\n";
	}
	
#if !defined(TIME_OUTPUT)	
	cout << mlist;
#endif
	if (!mofname.empty())
	{
		fstream ofs;
		openFileCPP(ofs, mofname, (ios::out | ios::trunc));
		ofs << mlist;
		ofs.close();
	}
}

void outSeq(SearchResults & mlist, const string & sofname, const string & otype, const vector<vector<int> > & gaplen)
{
	if (sofname.empty())
	{
		return;
	}
	
	if (mlist.empty())
	{
		cout << "Warning: no match, no sequence.\n";
	}
	
	int i;
	fstream ofs;
	string aa;
	openFileCPP(ofs, sofname, (ios::out | ios::trunc));
	for (i = 0; i < mlist.size(); i++)
	{
		extractSeq(aa, mlist[i], mlist[i].getSeq(), otype, gaplen);
		ofs << std::setw(8);
		if (mlist[i].getRmsd() == DBL_MAX)
		{
			ofs << "NaN";
		}
		else
		{
			ofs << std::setprecision(5) << mlist[i].getRmsd();
		}
		ofs << " " << aa;
	}
	ofs.close();
}

void parseCmdLine(int argc, char** argv, SearchOptions & sopts)
{
	map<string, bool> spec;

	while (1)
	{
		int oind = 0;
		static struct option opts[] =
		{
			{"query", 1, 0, 1},
			{"target", 1, 0, 2},
			{"targetList", 1, 0, 3},
      		{"topN", 1, 0, 4},
      		{"phiEps", 1, 0, 5},
      		{"psiEps", 1, 0, 6},
      		{"dEps", 1, 0, 7},
      		{"rmsdCut", 1, 0, 8},
      		{"structOut", 1, 0, 9},
      		{"matchOut", 1, 0, 10},
      		{"rmsdMode", 1, 0, 11},
      		{"ddZscore", 0, 0, 12},
      		{"tune", 1, 0, 13},
      		{"seqOut", 1, 0, 14},
			{"bbRMSD", 0, 0, 15},
			{"outType", 1, 0, 16},
			{"matchIn", 1, 0, 17},
			{"gapLen", 1, 0, 18},
      		{0, 0, 0, 0}
		};

		int c = getopt_long (argc, argv, "", opts, &oind);
		if (c == -1)
		{
			break;
		}

		switch (c) {
			case 1:
				sopts.setQsFile(string(optarg));
        		spec[string(opts[oind].name)] = true;
        		break;

      		case 2:
				sopts.setTsFile(string(optarg));
        		spec[string(opts[oind].name)] = true;
        		break;

      		case 3:
				sopts.setTsFiles(string(optarg));
        		spec[string(opts[oind].name)] = true;
        		break;

      		case 4:
				sopts.setTopN(optarg);
        		spec[string(opts[oind].name)] = true;
        		break;

	  		case 5:
				sopts.setPhiDevCut(optarg);
	  			spec[string(opts[oind].name)] = true;
				break;

	  		case 6:
				sopts.setPsiDevCut(optarg);
	  			spec[string(opts[oind].name)] = true;
				break;

	  		case 7:
				sopts.setDistDevCut(optarg);
	  			spec[string(opts[oind].name)] = true;
				break;

      		case 8:
				sopts.setRmsdCut(optarg);
        		spec[string(opts[oind].name)] = true;
        		break;

	  		case 9:
				sopts.setStructOutDir(string(optarg));
				spec[string(opts[oind].name)] = true;
	  			break;

	  		case 10:
				sopts.setMatchOutFile(string(optarg));
				spec[string(opts[oind].name)] = true;
	  			break;
			
	  		case 11:
				sopts.setRmsdBoundMode(optarg);
        		spec[string(opts[oind].name)] = true;
        		break;

	  		case 12:
				sopts.setDistDevZscore(true);
        		spec[string(opts[oind].name)] = true;
        		break;
			
	  		case 13:
				sopts.setTuningParam(optarg);
				spec[string(opts[oind].name)] = true;
		        break;
			
	  		case 14:
				sopts.setSeqOutFile(string(optarg));
				spec[string(opts[oind].name)] = true;
			  	break;

			case 15:
				sopts.setBBRmsd(true);
				spec[string(opts[oind].name)] = true;
			  	break;

			case 16:
				sopts.setOutType(string(optarg));
				spec[string(opts[oind].name)] = true;
			  	break;

			case 17:
				sopts.setMatchInFile(string(optarg));
				spec[string(opts[oind].name)] = true;
	  			break;

			case 18:
				sopts.setGapLen(string(optarg));
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

	// make sure all spec options have been specified
  	if (!((spec.find(string("query")) != spec.end())
		&& ((((spec.find(string("target")) != spec.end()) || (spec.find(string("targetList")) != spec.end())) && (spec.find(string("rmsdCut")) != spec.end()))
		|| (spec.find(string("matchIn")) != spec.end()))))
  	{
		usage();
		error("not all required options specified");
	}
	if ((sopts.getOutType().compare("wgap") == 0) && (sopts.getGapLen().size() <= 0))
	{
		usage();
		error("'wgap' output type only works when '--gapLen' option has been specified");
	}
	
#if !defined(DEBUG_ZSCORE)
	if (0 == sopts.getDistDevBoundMode())
	{
		sopts.setDistDevZscore(false);
	}
#endif

#if defined(DEBUG_CMDLINE)
	cout << "output backbone RMSDs? " << sopts.getBBRmsd() << "\n";
	cout << "distance deviation cutoff: " << sopts.getDistDevCut() << "\n";
	cout << "distance deviation mode: " << sopts.getDistDevBoundMode() << "\n";
	cout << "output distance deviation Z-score? " << sopts.getDistDevZscore() << "\n";
	cout << "Phi angle deviation cutoff: " << sopts.getPhiDevCut() << "\n";
	cout << "Psi angle deviation cutoff: " << sopts.getPsiDevCut() << "\n";
	cout << "CA RMSD mode: " << sopts.getRmsdBoundMode() << "\n";
	cout << "CA RMSD cutoff: " << sopts.getRmsdCut() << "\n";
	cout << "output type: " << sopts.getOutType() << "\n";
	cout << "topN: " << sopts.getTopN() << "\n";
	cout << "tuning parameter for CA RMSD bounding: " << sopts.getTuningParam() << "\n";
	exit(-1);
#endif
}

void renameStruct(SearchResults & mlist, vector<int>* & oldidx, const string & sodname, const string & otype)
{
	if (sodname.empty())
	{
		return;
	}
	
	if (mlist.empty())
	{
		cout << "Warning: no match, no structure.\n";
		return;
	}
	
	int i, count = 1, base = 10;
	char format[500], fname[500], tmpfn[500];
	
	while (1)
	{
		if (base > mlist.size())
		{
			break;
		}
		base *= 10;
		count += 1;
	}
	
	if (otype.compare("full") == 0)
	{
		sprintf(format, "%%s/full%%0%dd.pdb", count);
	}
	else
	{
		if (otype.compare("match") == 0)
		{
			sprintf(format, "%%s/match%%0%dd.pdb", count);
		}
		else
		{
			if (otype.compare("wgap") == 0)
			{
				sprintf(format, "%%s/wgap%%0%dd.pdb", count);
			}
			else
			{
				error("unidentified output type");
			}
		}
	}
	
	for (i = 0; i < mlist.size(); i++)
	{
		if (mlist[i].getRmsd() == DBL_MAX)
		{
			continue;
		}
		sprintf(tmpfn, "%s/tmp%d", sodname.c_str(), (*oldidx)[i]);
		ASSERT(exist(tmpfn), "could not find file " + string(tmpfn));
		sprintf(fname, format, sodname.c_str(), i + 1);
		ASSERT(0 == rename(tmpfn, fname), "could not rename file " + string(tmpfn));
	}

	delete oldidx;
}

double rmsdCut(const int & mode, const double & rtot, const int & ltot, const double & rsofar, const int & lsofar, const int & lcur, const double & tune)
{
	switch (mode) {
		case 2:
			return rtot;
		case 1:
			return sqrt(((rtot * rtot * double(ltot)) - (rsofar * rsofar * double(lsofar)) - (tune * rtot * rtot * double(ltot - lsofar - lcur))) / double(lcur));
		case 0:
		default:
			return sqrt(((rtot * rtot * double(ltot)) - (rsofar * rsofar * double(lsofar))) / double(lcur));
	}
}

void searchByDistDistr(SearchResults & mlist, const bool & bbrmsd, const double & phieps, const double & psieps, const int & dmode, const double & deps, const int & rmode, const double & rthresh, const double & tune, const vector<vector<int> > & gaplen)
{
	QueryStruct& qs = mlist.getQs();
	vector<int>& qscr = qs.getCenRes();
	vector<vector<double> >& qscrdihed = qs.getCenResDihed();
	vector<int>& qsrbb = qs.getResBefBrk();
	CenResGrp cres(qs.getNumSeg());
	TargetStruct ts;
	int i, tsi;

	for (i = 0; i < qs.getNumSeg(); i++)
	{
		cres[i].setResIdx(qscr[i]);
		cres[i].setSegIdx(i);		
		cres[i].setBegRes(qsrbb[i] + 1);
		cres[i].setEndRes(qsrbb[i + 1]);
		cres[i].setSegLen(qsrbb[i + 1] - qsrbb[i]);
		cres[i].setPhi(qscrdihed[i][0]);
		cres[i].setPsi(qscrdihed[i][1]);
		// assume that all the other segments have perfect matches
		// NOTE: BB and CA share the same RMSD cutoff, and factor 4 (i.e., number of BB atoms) makes no difference
		cres[i].setSelfRmsdCut(rmsdCut(rmode, rthresh, qs.getNumRes(), 0.0, 0, cres[i].getSegLen(), tune));
		cres[i].setRank(0.0);
		if ((i > 0) && (gaplen.size() > 0))
		{
			cres[i].setGapLen(gaplen[i - 1]);
			if (cres[i].getGapLen().size() > 0)
			{
				cres[i - 1].setRank(float(i - 1 - qs.getNumSeg()));
				cres[i].setRank(float(i - qs.getNumSeg()));
			}
		}
	}

	if (bbrmsd)
	{
		cout << "Searching based on backbone RMSD...\n";
		
		qs.readBBCoords();
		ASSERT(qs.fullBB(), "backbone atoms missing in file " + qs.getFileName());
		for (tsi = 0; tsi < mlist.getNumTs(); tsi++)
		{
#if !defined(TIME_SEARCH)
			printf("[Visiting %d/%d %s...]\n", tsi + 1, mlist.getNumTs(), mlist.getTsFile(tsi).c_str());
#endif
			ts.readTargetStruct(mlist.getTsFile(tsi));
			ts.readBBCoords();
			
			auxSearchByDistDistrBB(mlist, cres, qs, ts, tsi, phieps, psieps, dmode, deps, rmode, rthresh, tune);

			ts.closeProteinStructFile();
		}
	}
	else
	{
		cout << "Searching based on CA RMSD...\n";
		
		qs.readCACoords();
		for (tsi = 0; tsi < mlist.getNumTs(); tsi++)
		{
#if !defined(TIME_SEARCH)
			printf("[Visiting %d/%d %s...]\n", tsi + 1, mlist.getNumTs(), mlist.getTsFile(tsi).c_str());
#endif
			ts.readTargetStruct(mlist.getTsFile(tsi));
			ts.readCACoords();
			
			auxSearchByDistDistrCA(mlist, cres, qs, ts, tsi, phieps, psieps, dmode, deps, rmode, rthresh, tune);

			ts.closeProteinStructFile();
		}
	}
}

void transformStructByRMSD(vector<string> & opdb, const vector<string> & pdbinfo, const vector<int> & ri, double* t, double** u)
{
	opdb.clear();
	
	int i, posn, posc;
	double x, y, z, xn, yn, zn;
	char coords[25];
	coords[24] = '\0'; // the coordinate section is 8*3 = 24 characters long
	string atom;

	for (i = 0; i < ri.size(); i++)
	{
		ASSERT(!pdbinfo[ri[i]].empty(), "could not dump PDB from empty PDB information");
		opdb.push_back(pdbinfo[ri[i]]);
		
		posn = 0;
		posc = 0;
		while ((posn = opdb[i].find("\n", posc)) != string::npos)
		{
			atom = opdb[i].substr(posc, posn - posc + 1);
			ASSERT(sscanf(atom.c_str(), "%*30c%8lf%8lf%8lf", &x, &y, &z) == 3, "could not parse coordinates from " + atom); // coordinates start at character 31 of the line (PDB), so skip the first 30
			// apply transformation
			xn = t[0] + u[0][0] * x + u[0][1] * y + u[0][2] * z;
			yn = t[1] + u[1][0] * x + u[1][1] * y + u[1][2] * z;
			zn = t[2] + u[2][0] * x + u[2][1] * y + u[2][2] * z;
			// write in the new coordinates
		    sprintf(coords, "%8.3f%8.3f%8.3f", xn, yn, zn);
		    opdb[i].replace(posc + 30, 24, coords); // coordinates start at character 31 of the line (PDB), so skip the first 30
			// advance to next atom
		    posc = posn + 1;
		}			
	}
}

void usage()
{
	int w = 80, p1 = 1, p2 = p1 + 18;
	cout << optionUsage("--query", "query PDS file (required).", w, p1, p2) << "\n";
	cout << optionUsage("--target", "target PDS file.", w, p1, p2) << "\n";
	cout << optionUsage("--targetList", "a file with a list of target PDS files. Either --target or --targetList must be given.", w, p1, p2) << "\n";
	cout << optionUsage("--rmsdCut", "RMSD cutoff for defining a match (in Angstrom).", w, p1, p2) << "\n";
	cout << optionUsage("--gapLen", "optional: impose length constraints on gaps between adjacent segments. E.g., --gapLen '1-5;;3-10' will restrain the gap between the first and second segments in the query to be between 1 and 5 residues, no constraint will be placed on the gap between the second and third segments, whereas the gap between the third and fourth segments will need to be between 3 and 10 residues. Note: the order of segments is as they appeared within the PDB file of the query. Also, the number of restraints must match the number of gaps (e.g., in the example above, the query must have four segments).", w, p1, p2) << "\n";
	cout << optionUsage("--matchOut", "optional: file name for storing resulting matches (one line per match); contains all information for defining a match (location of match and RMSD).", w, p1, p2) << "\n";
	cout << optionUsage("--seqOut", "optional: file name for storing the sequences of matching regions (one per line). See --outType for defining what gets output.", w, p1, p2) << "\n";
	cout << optionUsage("--structOut", "optional: name of directory for writing match structures in PDB format (one PDB file per match). See --outType for defining what gets output.", w, p1, p2) << "\n";
	cout << optionUsage("--outType", "optional: specifies what kind of sequences and/or structures to output (only works when --seqOut and/or --structOut have been specified). If set to 'full', will output the entire target sequence and/or structure containing a matching region; if set to 'match', will output just the matching region; if set to 'wgap', will output the matching region with the gap(s) constrained by '--gapLen'. By default, output one PDB file per matching region. In all cases, output structures are aligned to superimpose the matching region onto the query.", w, p1, p2) << "\n";
	cout << optionUsage("--bbRMSD", "optional: search by full-backbone RMSD (default is C-alpha RMSD).", w, p1, p2) << "\n";
	cout << optionUsage("--topN", "optional: keep the best this many matches in terms of the search metric (must be integer); default is 0 (no limit).", w, p1, p2) << "\n";
	cout << optionUsage("--rmsdMode", "optional: RMSD bounding mode. 0 -- provable RMSD bounds will be calculated, guaranteeing that all matches within --rmsdCut will be found (default). 1 -- greedy bound that enforces some uniformity of RMSD residuals (see documentation). 2 -- uses --rmsdCut for both the overall RMSD cutoff as well as the cutoff for partial matches (see documentation).", w, p1, p2) << "\n";
	cout << optionUsage("--tune", "optional: tuning parameter for greedy RMSD cutoff (i.e., when --rmsdMode 1 is specified); 0.5 for default (see documentation).", w, p1, p2) << "\n";
	cout << optionUsage("--dEps", "optional: user-defined greedy distance deviation cutoff (in Angstrom). If given, rather than applying a provable bound on inter-segment distances, this cutoff will be applied.", w, p1, p2) << "\n";
	cout << optionUsage("--phiEps", "optional: phi angle deviation cutoff (in degrees); default is 180.0, meaning no cutoff is applied.", w, p1, p2) << "\n";
	cout << optionUsage("--psiEps", "optional: psi angle deviation cutoff (in degrees); default is 180.0, meaning no cutoff is applied.", w, p1, p2) << "\n";
	cout << optionUsage("--ddZscore", "optional: output a Z-score that describes the distribution of inter-segment distance deviations (between query and matches) relative to the greedy cutoff --dEps. High Z-scores (> 3.5, in our experience), indicating a good choice of --dEps and suggesting that all or nearly all matches were found despite the greedy constraint.", w, p1, p2) << "\n";
	cout << optionUsage("--matchIn","a list of matches from a previously run search (i.e., the result of --matchOut of a previous run). If specified, will skip searching and will produce outputs directly.", w, p1, p2) << "\n";
	
//	cout << optionUsage("--ddOutFile", "optional: file name to output distance deviations of matches to (only works when Z score is output). By default, no distance deviations are output. ", w, p1, p2) << "\n";
}

