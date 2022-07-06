import sys
import os
import pandas as pd
import numpy as np

import autostreamtree.functions as ast
from autostreamtree.params import parseArgs

def main():

	params = parseArgs()

	print("Running autostream tree", params.run, "workflow\n")
	G = ast.read_network(params.network, params.shapefile)

	# read point coordinate data
	points = pd.read_csv(params.geodb, sep="\t", header=None)
	(point_coords, pop_coords, popmap) = ast.processSamples(params, points, G)

	# read genetic data from VCF
	seqs = ast.read_vcf(params.vcf, concat=params.concat, popmap=popmap)

	# calculate genetic distance matrix
	if params.run != "STREAMDIST" and params.run != "RUNLOCI":
		gen=None
		if not params.genmat:
			print("\nCalculating genetic distances...")
			(gen, pop_gen) = ast.getGenDistMats(params, point_coords, popmap, seqs)
		else:
			print("\nReading genetic distances from provided matrix:", params.genmat)
			inmat = pd.read_csv(params.genmat, header=0, index_col=0, sep="\t")
			(gen, pop_gen) = ast.parseInputGenMat(params, inmat, point_coords, popmap)

		ast.reportPopGenMats(params, gen, pop_gen, point_coords, pop_coords)

		if params.run == "GENDIST":
			sys.exit(0)


	#IMPORTANT NOTE: This function will silently skip any loci for which calculations aren't possible (e.g., population having no data)
	if params.run == "RUNLOCI":
		genlist=list()
		popgenlist=list()
		for loc in ast.getLocData(seqs):
			try:
				(gen, pop_gen) = ast.getGenDistMats(params, point_coords, popmap, loc)
				genlist.append(gen)
				popgenlist.append(pop_gen)
			except Exception:
				pass
		print("\nCalculating average of locus-wise genetic distance matrices...")
		if genlist[0] is not None:
			gen = np.mean(genlist, axis=0)
		else:
			gen = genlist[0]
		if popgenlist[0] is not None:
			popgen = np.mean(popgenlist, axis=0)
		else:
			popgen = popgenlist[0]

		ast.reportPopGenMatsList(params, genlist, popgenlist, point_coords, pop_coords)


#Call main function
if __name__ == '__main__':
	main()
