import sys
import os
import pandas as pd

import autostreamtree.functions as ast
from autostreamtree.params import parseArgs

def main():

    params = parseArgs()

    print("Running autostream tree", params.run, "workflow\n")
    G = ast.read_network(params.network, params.shapefile)

	# read point coordinate data
    points = pd.read_csv(params.geodb, sep="\t", header="infer")
    (point_coords, pop_coords, popmap) = ast.processSamples(params, points, G)

    print(popmap)



#Call main function
if __name__ == '__main__':
	main()
