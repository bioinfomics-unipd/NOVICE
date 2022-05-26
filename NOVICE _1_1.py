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
#sys.argv[7] provide filename .ffn
#sys.argv[8] provide kmer size


if len(sys.argv) < 6:
    print('----------------------------------------------')
    print('----------------------------------------------')
    print('CORRECT SYNTAX: NOVICE_1_1.py path/file argument1 argument2 argument3 argument4 argument5 argument6[required only if GC is chosen as option in argument1]')
    print('available options for argument1: rank, PCA, GC, Kmers, all')
    print('available options for argument2: mean, median')
    print('available options for argument3: euclidean, braycurtis, correlation')
    print('available options for argument4: provide the prefix of output files')
    print('available options for argument5: yes, no')
    print('available options for argument6: filename .ffn')
    print('available options for argument7: provide any integer to define the kmers size')
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
    	file = [line for line in open(file)]
    else:
        raise Exception('File not supported')
    return file

data = ReadFiles(sys.argv[1]) #g_c_num is gene coverage number of the MAG
          
data_CLEAN = data.dropna()

print('*******************************************************')
print('>>>>>>>>>>>>>> NOVICE version 1.1 beta <<<<<<<<<<<<<<<<<<<')
print('*******************************************************')
print('>>>>> check our GitHub page! <<<<<')
print('>>>>> Link: https://github.com/bioinfomics-unipd/NOVICE')
print('*******************************************************')
          
if sys.argv[2] != 'rank' and sys.argv[2] != 'all' and sys.argv[2] != 'PCA' and sys.argv[2] != 'GC' and sys.argv[2] != 'Kmers': 

    print('Error, invalid argument1')
    sys.exit()
    
if sys.argv[3] != 'mean' and sys.argv[3] != 'median': 

    print('Error, invalid argument2')
    sys.exit()

if sys.argv[4] != 'euclidean' and sys.argv[4] != 'braycurtis' and sys.argv[4] != 'correlation': 

    print('Error, invalid argument3')
    sys.exit()

if sys.argv[6] != 'yes' and sys.argv[6] != 'no': 

    print('Error, invalid argument5')
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
    plt.savefig(str(sys.argv[5])+'_treshold_barplot.pdf')
    plt.show()
    
    treshold = input('Please provide a threshold based on the barplot: ')
    treshold = float(treshold)
    
    data_FIN['Class'] = ['normal' if e < treshold else 'outlier' for e in data_FIN['DistanceFromCentroids']]
    
    
    #Sort the merged dataframe by distances
    data_SORTED= data_FIN.sort_values(by = 'DistanceFromCentroids',ascending = False)
    
    if sys.argv[6] == 'yes':
    	data_SORTED = data_SORTED[data_SORTED['DistanceFromCentroids'] < treshold]
    
    data_SORTED.to_csv(str(sys.argv[5])+ '_sorted_dataframe.csv')
    
    
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
    plt.savefig(str(sys.argv[5])+'_PCA_screenplot.pdf')
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
    plt.savefig(str(sys.argv[5])+'_PCAplot.pdf')
    plt.clf()
    print('PCA plot is created')
    print('PCA plot is saved')  
    
if sys.argv[2] == 'GC' or sys.argv[2] == 'all':
    
    data = ReadFiles(sys.argv[7])
    
    i = -1
    head = []
    position = []
    chonky = []
    infos = []
    
    for line in data:
    	i += 1
    	if line.find('>') != -1:
    		info = line 
    		header = line[line.find('>')+1:line.find(' ')]
    		head.append(header)
    		position.append(i)
    		infos.append(info)
    	chonky.append(line)
    
    maps = dict(zip(head, position))
    oth_maps = dict(zip(head, infos))
    
    a = -1
    index = {}
    for e in head:
    	a +=1
    	index[e] = a
    
    sequencesToGC = {}
    for elem in data_SORTED['index']: 
    	if elem in maps.keys():
    		indec = index[elem]
    		posnextelem = index[elem]+1
    		nextelem = list(index.keys())[posnextelem]
    		sequence_head = chonky[maps[elem]][1:15]
    		sequence_obj = chonky[maps[elem]+1:maps[nextelem]]
    		sequence_obj = ''.join(sequence_obj)
    		sequence_obj = sequence_obj.replace('\n','')
    		sequencesToGC[sequence_head] = sequence_obj
    
    from collections import Counter
    
    GC_content = {k:Counter(v) for k,v in sequencesToGC.items()}
    
    perc = [(e['G']+e['C'])/(e['A']+e['T']+e['C']+e['G']) for e in GC_content.values()]
    
    data_SORTED['GC'] = perc
    
    data_SORTED.to_csv(str(sys.argv[5])+ '_sorted_dataframe.csv')
    
    
    print('dataframe sorted is overwritten with GC contenent')

    
    x = data_SORTED[data_SORTED['Class'] == 'outlier']['index']
    y = data_SORTED[data_SORTED['Class'] == 'outlier']['GC']
    x1 = data_SORTED[data_SORTED['Class'] == 'normal']['index']
    y1 = data_SORTED[data_SORTED['Class'] == 'normal']['GC']
    plt.scatter(x, y, color = 'red', label = 'outlier', alpha = 0.5)
    plt.scatter(x1, y1, color = 'blue', label = 'normal', alpha = 0.5)
    plt.xticks([])
    plt.legend()
    plt.savefig(str(sys.argv[5])+'_GC_scatterplot.pdf')
    plt.clf()
    
    print('GC scatterplot is saved')
    
if sys.argv[2] == 'Kmers' or sys.argv[2] == 'all':
    print("Kmers analisys in progress...")
    
    def Kmere(sequence):
    	kmrs = []
    	for e in sequence:
    		kmrs.append(sequence[sequence.find(e):sequence.find(e)+int(sys.argv[8])])
    	return kmrs
    
    Kmere_list = {k:Kmere(v) for k,v in sequencesToGC.items()} #for each gene the kmeres for its sequence
    
    unique_kmeres = list(np.unique([np.unique(e) for e in Kmere_list.values()]))
    
    kmere_counts = {k:{e:v.count(e) for e in v} for k,v in Kmere_list.items()} 
    
    data_SORTED = data_SORTED.dropna()
    data_SORTED.index = data_SORTED['index']
    
    for e in unique_kmeres:
    	data_SORTED[e] = [0 for i in data_SORTED['index']]
    
    for e in unique_kmeres:
    	for gene, dct in kmere_counts.items():
    		if e in dct.keys():
    			data_SORTED.loc[gene, e] = dct[e]
    
    
    data_SORTED.to_csv(str(sys.argv[5])+ '_sorted_dataframe.csv')
    print('dataframe sorted is overwritten with Kmers informations')
    
    #PCA on kmere
    #finding how many PCs
    plt.clf()
    
    scalar = StandardScaler()
    pca = PCA()
    
    pipeline = make_pipeline(scalar, pca)
    pipeline.fit(data_SORTED.iloc[:, 18:])
    
    features = range(pca.n_components_)
    plt.bar(features, pca.explained_variance_)
    plt.xlabel('PCA features')
    plt.ylabel('variance')
    plt.xticks([])
    plt.savefig(str(sys.argv[5])+'_Kmers_PCA_screenplot.pdf')
    plt.clf()
    print('Kmers screenplot is created')
    print('Kmers screenplot is saved')
    
    model = PCA(n_components = 4)
    model.fit(data_SORTED.iloc[:, 18:])
    transformed = model.transform(data_SORTED.iloc[:, 18:])
    xs = transformed[:,0]
    ys = transformed[:,1]
    if data_SORTED['Class'].groupby(data_SORTED['Class']).count().shape[0] > 1:
    	number = data_SORTED['Class'].groupby(data_SORTED['Class']).count()[1]
    else:
    	number = 0
    plt.scatter(xs,ys, c = ['black' if e < treshold else 'red' for e in data_SORTED['DistanceFromCentroids']], label = 'outlier numb: '+str(number), alpha = 0.5 )
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA plot')
    plt.savefig(str(sys.argv[5])+'_Kmers_PCAplot.pdf')
    plt.clf()
    print('Kmers PCA plot is created')
    print('Kmers PCA plot is saved')  
    
    
        
print('')
print("Terminating......") 
print("Program closed succesfully")

