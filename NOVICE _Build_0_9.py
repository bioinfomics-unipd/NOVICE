#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import scipy.spatial.distance as dist
import numpy as np
import sys

#sys.argv[1] is the file you want to use
#sys.argv[2] decide if user want rank genes, PCA or everything
#sys.argv[3] decide if mean or meadian is used to calculate centroids
#sys.argv[4] decide if the distance used is euclidean, ..., ...
#sys.argv[5] provide the name to output files
#sys.argv[6] decide to cut off or not the outliers (yes or no)


if len(sys.argv) < 5:
    print('----------------------------------------------')
    print('----------------------------------------------')
    print('CORRECT SYNTAX: NOVICE_build_0_9.py path/file option1 option2 option3 option4')
    print('available options 1: rank, PCA, all')
    print('available options 2: mean, median')
    print('available options 3: euclidean, braycurtis, correlation')
    print('available options 4: provide the prefix of output files')
    print('available options 5: yes, no')
    print('----------------------------------------------')
    print('----------------------------------------------')
    sys.exit()
          
def ReadFiles(file):
    '''Read tsv, csv and ffn file -- requires Pandas'''
    if file.split('.')[1] == 'tsv':
        file = pd.read_csv(file, sep = '\t')
    elif file.split('.')[1] == 'csv':
        file = pd.read_csv(file)
    elif file.split('.')[1] == 'ffn':
        print('{} is ffn')
    else:
        raise Exception('File not supported')
    return file

data = ReadFiles(sys.argv[1]) #g_c_num is gene coverage number of the MAG
          
data_CLEAN = data.dropna()
          
if sys.argv[2] != 'rank' and sys.argv[2] != 'all' and sys.argv[2] != 'PCA': 

    print('Error, option has been not recognized')
    sys.exit()


print('dimension before NaN removing: {}'.format(data.shape))
print('dimension after NaN removing: {}'.format(data_CLEAN.shape))


if sys.argv[2] == 'rank' or sys.argv[2] == 'all' or sys.argv[2] == 'PCA':

    def centroid(cols, opt = 'mean'):
        '''Calculate centroid for each condition and return them in a list -- requires Numpy'''

        def CalcMeanCond(condition):
            '''Calculate the mean for the condition given and it returns it -- requires Numpy'''
            #filtering for nan and obtaining data available to calculation
            lst = list(data_CLEAN[condition][data_CLEAN[condition].notna()][1:])
            val = [float(el) for el in lst]
            #calculate the mean
            mn = np.mean(val)
            return mn
        
        def CalcMedianCond(condition):
            '''Calculate the mean for the condition given and it returns it -- requires Numpy'''
            #filtering for nan and obtaining data available to calculation
            lst = list(data_CLEAN[condition][data_CLEAN[condition].notna()][1:])
            val = [float(el) for el in lst]
            #calculate the mean
            md = np.median(val)
            return md
        if opt == 'mean':
        	centroids = [CalcMeanCond(e) for e in cols]
        if opt == 'median':
        	centroids = [CalcMedianCond(e) for e in cols] 
        return centroids

    centroids = tuple(zip(data.columns[1:], centroid(data.columns[1:], sys.argv[3])))
    print('centroids created')  


    def CalcEuclDist(obj):
    	eucl_lists = []
    	
    	for n in range(1,len(obj.iloc[:,1])):
    		if sys.argv[4] == 'euclidean':
    			eucl_lists.append(dist.euclidean([float(e) for e in data_CLEAN.iloc[n,:][1:]], centroid(data.columns[1:])))
    		if sys.argv[4] == 'braycurtis':
    			eucl_lists.append(dist.braycurtis([float(e) for e in data_CLEAN.iloc[n,:][1:]], centroid(data.columns[1:])))
    		if sys.argv[4] == 'correlation':
    			eucl_lists.append(dist.correlation([float(e) for e in data_CLEAN.iloc[n,:][1:]], centroid(data.columns[1:])))
    	return eucl_lists        	
        
   
    data_FIN = data_CLEAN.drop(0)

    #creating a df with same dimension of CalcEuclDist
    data_FIN = data_CLEAN.drop(0)
    
    data_FIN['DistanceFromCentroids'] = CalcEuclDist(data_CLEAN)
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2,1)
    ax[0].boxplot(data_FIN.iloc[:,13], showfliers = True)
    ax[1].boxplot(data_FIN.iloc[:,13], showfliers = False)
    ax[0].set_title('Distance Distribution')
    ax[1].set_xlabel('Sample')
    ax[0].set_ylabel('distance - With Outliers')
    ax[1].set_ylabel('distance - Not Ouliers')
    plt.savefig(str(sys.argv[5])+'treshold_barplot.pdf')
    plt.show()
    
    treshold = input('Please provide a threshold based on the barplot: ')
    treshold = float(treshold)
    
    data_FIN['Class'] = ['normal' if e < treshold else 'outlier' for e in data_FIN['DistanceFromCentroids']]
    
    
    #Sort the merged dataframe by distances
    data_SORTED= data_FIN.sort_values(by = 'DistanceFromCentroids',ascending = False)
    
    if sys.argv[6] == 'yes':
    	data_SORTED = data_SORTED[data_SORTED['DistanceFromCentroids'] < treshold]
    
    data_SORTED.to_csv(str(sys.argv[5])+ 'sorted_dataframe.csv')
    
    
    print('dataframe sorted is created')
    print('dataframe sorted is saved as csv')

    plt.clf()     
    print('treshold barplot is created')
    print('treshold barplot is saved')

if sys.argv[2] == 'PCA' or sys.argv[2] == 'all':

    
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline

    scalar = StandardScaler()
    pca = PCA()

    pipeline = make_pipeline(scalar, pca)

    pipeline.fit(data_SORTED.iloc[:, 1:13])

    features = range(pca.n_components_)
    plt.bar(features, pca.explained_variance_)
    plt.xlabel('PCA features')
    plt.ylabel('variance')
    plt.xticks(features)
    plt.savefig(str(sys.argv[5])+'PCA_screenplot.pdf')
    plt.clf()
    print('screenplot is created')
    print('screenplot is saved')

#pca components can be 2, 3, 4
    model = PCA(n_components = 2)
    model.fit(data_SORTED.iloc[:, 1:13])
    transformed = model.transform(data_SORTED.iloc[:, 1:13])
    xs = transformed[:,0]
    ys = transformed[:,1]
    if data_SORTED['Class'].groupby(data_SORTED['Class']).count().shape[0] > 1:
    	number = data_SORTED['Class'].groupby(data_SORTED['Class']).count()[1]
    else:
    	number = 0
    plt.scatter(xs,ys, c = ['black' if e < treshold else 'red' for e in data_SORTED['DistanceFromCentroids']], label = 'outlier numb: '+str(number), alpha = 0.5 )
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend()
    plt.title('PCA plot')
    plt.savefig(str(sys.argv[5])+'PCAplot.pdf')
    plt.clf()
    print('PCA plot is created')
    print('PCA plot is saved')    
    print('')
    print("Now I'll happily die :D") 

