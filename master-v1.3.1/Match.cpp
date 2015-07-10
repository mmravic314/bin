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

#include "Match.h"

Match::Match(const Match & other)
{
	Match* m = (Match *)&other;
	initMatch(m->getSearchResults(), m->getTsIdx());
	setInsertTime(m->getInsertTime());
	setRmsd(m->getRmsd());
	setBegRes(m->getBegRes(), m->getNumSeg());
	setRotation(m->getRotation());
	setTranslation(m->getTranslation());
	setMaxDistDev(m->getMaxDistDev());
}

Match::Match(SearchResults* p, const int & idx)
{
	initMatch(p, idx);
}

Match::~Match()
{
	delBegRes();
	delRotation();
	delTranslation();
}

void Match::delBegRes()
{
	delete [] _bres;
}

void Match::delRotation()
{
	int i;
	for (i = 0; i < 3; ++i)
	{
		delete [] _u[i];
	}
	delete [] _u;
}

void Match::delTranslation()
{
	delete [] _t;
}

void Match::initMatch(SearchResults* p, const int & idx)
{
	setTsIdx(idx);
	setSearchResults(p);
	if (NULL != p)
	{
		newBegRes(p->getQsNumSeg());
	}
	newRotation();
	newTranslation();
}

void Match::newRotation()
{
	_u = new double* [3];
	int i;
	for (i = 0; i < 3; ++i)
	{
		_u[i] = new double [3];
	}
}

bool Match::operator<(const Match & other) const
{
	if ((_sr == NULL) || (_sr->getMatchLimit() > 0))
	{
		if (_rmsd != other.getRmsd())
		{
			return (_rmsd < other.getRmsd());
		}
	}
	else
	{
		if (_tsidx != other.getTsIdx())
		{
			return (_tsidx < other.getTsIdx());
		}		
	}
	return (this->getInsertTime() < other.getInsertTime());
}

ostream& operator<<(ostream & os, const Match & m)
{
	int i;
	
	os << std::setw(8);
	if (m.getRmsd() == DBL_MAX)
	{
		os << "NaN";
	}
	else
	{
		os << std::setprecision(5) << m.getRmsd();
	}
	os << " " << m.getTsFile() << " [";
	for (i = 0; i < (m.getNumSeg() - 1); i++)
	{
		os << "(" << m.getBegRes(i) << "," << m.getEndRes(i) << "), ";
	}
	os << "(" << m.getBegRes(i) << "," << m.getEndRes(i) << ")]";

	return os;
}

int Match::parseMatch(char* line, const vector<vector<int> > & gaplen)
{
	ASSERT(_sr != NULL, "'SearchResults' pointer should be specified before parsing match from line");
	
	char *cur;
	int i, eres, b, e, gl;
	string tsfn;

	cur = strtok(line, " [(,)]");
	if (0 == string(cur).compare("NaN"))
	{
		_rmsd = DBL_MAX;
	}
	else
	{
		ASSERT(sscanf(cur, "%lf", &_rmsd) == 1, "could not parse match from line, current state: " + string(line));
	}

	cur = strtok(NULL, " [(,)]");
	tsfn = string(cur);	

	i = 0;
	while (*(cur + strlen(cur) + 1) != ']')
	{
		cur = strtok(NULL, " [(,)]");
		if (i < _sr->getQsNumSeg())
		{
			ASSERT(sscanf(cur, "%d", &(_bres[i])) == 1, "could not parse match from line, current state: " + string(line));
			if ((i > 0) && (gaplen.size() > 0))
			{
				if (gaplen[i - 1].size() > 0)
				{
					gl = _bres[i] - eres - 1;
					if (gaplen[i - 1].size() == 1)
					{
						if (gl != gaplen[i - 1][0])
						{
							return (-1);
						}
					}
					else
					{
						if (gaplen[i - 1].size() == 2)
						{
							b = min(gaplen[i - 1][0], gaplen[i - 1][1]);
							e = max(gaplen[i - 1][0], gaplen[i - 1][1]);
							if (!((gl >= b) && (gl <= e)))
							{
								return (-1);
							}
						}
						else
						{
							error("invalid number of gap length constraints");
						}
					}
				}
			}
		}
		
		cur = strtok(NULL, " [(,)]");
		if (i < _sr->getQsNumSeg())
		{
			ASSERT(sscanf(cur, "%d", &eres) == 1, "could not parse match from line, current state: " + string(line));
		}

		++i;
	}
	ASSERT(_sr->getQsNumSeg() == i, "number of segments not consistent between query and match");

	setTsIdx(_sr->addTsFile(tsfn));
	setInsertTime(_sr->getNumMatFnd());

	return 0;
}

void Match::setRotation(double** u)
{
	int i;
	for (i = 0; i < 3; ++i)
	{
		memcpy((void *)(_u[i]), (const void *)(u[i]), 3 * sizeof(double));
	}
}

bool cmpMatchByTarget(Match* a, Match* b)
{
	return (a->getTsIdx() < b->getTsIdx());
}

bool cmpMatchByRmsd(Match* a, Match* b)
{
	if (a->getRmsd() != b->getRmsd())
	{
		return (a->getRmsd() < b->getRmsd());
	}
	else
	{
		return (a->getInsertTime() < b->getInsertTime());
	}
}

bool cmpMatchPairByRmsd(pair<Match*, int> a, pair<Match*, int> b)
{
	if (a.first->getRmsd() != b.first->getRmsd())
	{
		return (a.first->getRmsd() < b.first->getRmsd());
	}
	else
	{
		return (a.first->getInsertTime() < b.first->getInsertTime());
	}
}

