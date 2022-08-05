import umap
#import umap.plot
#import json
import sys
#import numpy as np
#import pandas as pd
from pandas import DataFrame
#from matplotlib import pyplot as plt
import pandas
#from csv import reader


def run(input_filename, output_filename, n_neighbors, min_dist, metric, n_components, random_state):
    """this code will perform Univfor Manifild Approximationand Projection (UMAP)
    dimensionality reduction."""
    # opening the CSV file
    #location = r'/Users/marioacuna/_temp_analysis/_MATLAB_CaImaging/temp_data.csv'
    input_name = input_filename
    csvFile = pandas.read_csv(input_name, delimiter=',')
    embedding = umap.UMAP(n_neighbors=int(n_neighbors),
                          min_dist=float(min_dist),
                          metric=metric,
                          n_components=int(n_components),
                          random_state=int(random_state)).fit_transform(csvFile.values)
    mapper = umap.UMAP().fit(csvFile)
    #plt.scatter(embedding[:,0], embedding[:,1]);plt.show()
    df = DataFrame(embedding)
    df.to_csv(output_filename, header=False, index=None)


if __name__ == '__main__':
    # Get user inputs
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
        output_filename = sys.argv[2]
        n_neighbors = sys.argv[3]
        min_dist = sys.argv[4]
        metric = sys.argv[5]
        n_components = sys.argv[6]
        random_state = sys.argv[7]
        run(**dict(arg.split('=') for arg in sys.argv[1:]))
    else:
        input_filename = r'D:\_MATLAB_CaImaging\temp_data.csv'
        output_filename = r'D:\_MATLAB_CaImaging\UMAP.csv'
        n_neighbors = 3
        min_dist = 0.1
        metric = 'euclidean'
        n_components = 2
        random_state = 42
        run(input_filename, output_filename, n_neighbors, min_dist, metric, n_components, random_state)
    # run
    #run(input_filename, output_filename)
