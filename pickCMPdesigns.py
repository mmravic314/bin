

# look through a list text file summarizing rosetta scores after redesign
import sys
from collections import defaultdict

uniq = defaultdict(list)

scores = []
packs = []
uhbz = []


with open( sys.argv[1] ) as file:
	for i in file:

		path, seq , sc, ps, uhb = tuple( i.split() )

		sc = float( sc )
		ps = float( ps )
		uhb= int( uhb )


		uniq[ seq ].append(  ( ps, sc, uhb, path )  )

		scores.append( sc )
		packs.append(ps)
		uhbz.append(uhb)

print 'best', min( scores ), min( uhbz ), max( packs )

print
print  'packing_score', 'Rosetta_Energy', 'buried_unpair_Hbond_atoms', 'path'

for k, v in uniq.items():

	print k, '#_converged:', len(v)
	print sorted( v , reverse=True)[0]
	print