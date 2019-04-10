import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform, pdist

def interatomic_distance(df):
    #   Creating the complete Distance Matrix
    elements = df.index.unique(level='Element')
    
    dist = pd.DataFrame(squareform(pdist(df.iloc[:, 0:])), columns=df.index, index=df.index)
    #print(dist)
    
    #print(pdist(df.iloc[:, 0:]))
    #print(dist.loc[['Au'], ['Au']])
    
    
    #   Creating separate distance matrices for each element_element pair
    #   If the element pair are of the same element type, i.e. Au_Au, then the matrix is symmetric
    #       and you can use squareform
    for first_element in range(len(elements)):
    #    print(elements[first_element])
        for second_element in range(first_element, len(elements)):
    #        print(elements[second_element])
            if first_element == second_element:
                dist2 = dist.loc[[elements[first_element]], [elements[second_element]]]
                dist3 = squareform(dist2)
                dist3.sort()
                filename  = "./dist_" + str(elements[first_element]) + "_" + str(elements[second_element]) + ".dat"            
                with open(filename, "wb") as f:
                    np.savetxt(f, dist3, delimiter=",")
    #            print(elements[first_element] + "_" + elements[second_element])
    #            print(y)
            else:
                dist2 = dist.loc[[elements[first_element]], [elements[second_element]]]
                dist3 = dist2.values.ravel()
                dist3.sort()
                filename  = "./dist_" + str(elements[first_element]) + "_" + str(elements[second_element]) + ".dat"
                with open(filename, "wb") as f:
                    np.savetxt(f, dist3, delimiter=",")
    #            print(filename)
    #            print(len(dist3))