## Developed by Fantin Mesny 
## Max Planck Institute For Plant Breeding Research (Cologne, Germany)

# Added line 139 to output phylPCA data - Rowena Hill

import sys
import argparse
import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import pairwise_distances
from itertools import combinations_with_replacement 
import subprocess
import networkx as nx

def get_params(argv):
    parser = argparse.ArgumentParser(description='Analyse the composition in genes of different genome categories, taking phylogeny into account')
    parser.add_argument('-t', '--t', help="Phylogenetic tree, with leaf labels matching the genome names in the data table", required=True)
    parser.add_argument('-i', '--i', help="Dataframe of gene counts (gene families as columns, and genomes as rows). See 'example.csv'", required=True)
    parser.add_argument('-o', '--o', help="Output directory", required=True)
    parser.add_argument('-colors', '--colors', help="Color to match each lifestyle category: 'lifestyleA:blue,lifestyleB:#00FF00,...'", default='')
    a = parser.parse_args()
    return a

def getDistMatrix(tree):
    df=pd.DataFrame()
    nodes=[a.name for a in tree.depths(unit_branch_lengths=True) if a.name!=None]
    comb=list(combinations_with_replacement(nodes, 2))
    for c in comb:
        df.loc[c[0],c[1]]=tree.distance(c[0],c[1])
        df.loc[c[1],c[0]]=tree.distance(c[1],c[0])
    return df

def doPCA(df):
    DF=pd.DataFrame()
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(df)
    DF['PC1'] = pca_result[:,0]
    DF['PC2'] = pca_result[:,1]
    DF['genome']=df.columns
    DF=DF.rename(index=str, columns={'PC1':'PC1 ('+str(pca.explained_variance_ratio_[0])[2:4].lstrip('0')+'%)','PC2':'PC2 ('+str(pca.explained_variance_ratio_[1])[2:4].lstrip('0')+'%)'})
    DF=DF.set_index('genome')
    return DF

def terminal(cmd):
	p = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	p.wait()

def plotPhylPCA(metadata, palet, order, output):
    fig,axes=plt.subplots(1,1, figsize=(12,8))
    sns.scatterplot(x=metadata.columns[0], y=metadata.columns[1], s=70, data=metadata, hue='lifestyle', palette=palet, hue_order=order, ax=axes)
    axes.set_title('PCA of pairwise phylogenetic distances')
    plt.tight_layout()
    plt.savefig(output+'phylogenetic_pca.pdf')    
    
    
def runStat(output):
    script="""
    args <- commandArgs(trailingOnly=TRUE)
    dir <- args[1]

    library(vegan)
    data<-read.csv(paste(dir,'data.csv',sep=''), row.names='genome')
    metadata<-read.csv(paste(dir,'metadata.csv',sep=''), row.names='genome')
    dist <- vegdist(data, method='jaccard')
    distMatrix <- as.data.frame(as.matrix(dist))
    perm <- adonis2(dist~PC1+PC2+lifestyle, data=metadata, permutations = 9999)
    capture.output(perm, file=paste(dir,'permanova.txt',sep=''))
    write.csv(distMatrix,paste(dir,'distMatrix.csv',sep=''), row.names=TRUE)


    library(RVAideMemoire)
    permManova<-pairwise.perm.manova(dist,metadata$lifestyle,,nperm=9999)
    permManova <- as.data.frame(permManova[3])
    write.csv(permManova,paste(dir,'pairwiseComparisons.csv',sep=''), row.names=TRUE)
    
    """
    with open(output+'tmp.R','w+') as tmp:
        tmp.write(script)    
    terminal("Rscript "+output+"tmp.R '"+output+"'")
    terminal("rm "+output+"tmp.R")

def plotPCA(metadata, palet, order, output):
    distMatrix=pd.read_csv(output+'distMatrix.csv').set_index('Unnamed: 0')
    fig,axes=plt.subplots(1,1, figsize=(12,8))
    pca=doPCA(distMatrix).merge(metadata[['lifestyle']], left_index=True,right_index=True)
    pca.rename(index=str, columns={[c for c in pca.columns if 'PC1' in c][0]:'PC1',[c for c in pca.columns if 'PC2' in c][0]:'PC2'}).to_csv(a.o+'pca.csv')
    if palet=={}:
        ls=list(set(metadata['lifestyle']))
        palette=sns.color_palette(n_colors=len(ls))
        palet={ls[l]:palette[l] for l in range(len(ls))}
    sns.scatterplot(x=pca.columns[0], y=pca.columns[1], s=70, data=pca, hue='lifestyle', palette=palet, hue_order=order, ax=axes)
    axes.set_title('PCA of Jaccard distances calculated on genome compositions')
    plt.tight_layout()
    plt.savefig(output+'pca.pdf')
    return palet
    
def plotPvalMatrix(output, palet):
    permManova=pd.read_csv(output+'pairwiseComparisons.csv').rename(index=str, columns={'Unnamed: 0':'Lifestyles'})
    permManova=permManova.rename(index=str, columns={c:c.replace('p.value.','').replace('.',' ') for c in permManova.columns}).set_index('Lifestyles')
    fig,ax=plt.subplots(1,1,figsize=(7,7))
    heatmap = sns.heatmap(permManova, mask=permManova <= 0.05, square = True, linewidths = .5, cmap = 'Blues', cbar=False, vmin = -1000, vmax = 10000)
    heatmap = sns.heatmap(permManova, mask=permManova > 0.05, square = True, linewidths = .5, cmap = 'coolwarm_r', cbar_kws = {'shrink': .4, 'ticks' : [0.5, 0.33, 0]},vmin = -0.1, vmax=1,annot = True,annot_kws = {'size': 12},cbar=False)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_tick_params(rotation=45)
    [label.set_color(palet[label.get_text()] ) for label in ax.get_yticklabels()]
    [label.set_color(palet[label.get_text()] ) for label in ax.get_xticklabels()]
    plt.tight_layout()
    plt.savefig(output+'pvalMatrix.pdf')
    return permManova

def plotNetwork(output, permManova, palet):
    G = nx.Graph()
    for node in set(list(permManova.columns)+list(permManova.index)):
        G.add_node(node)
    for c in permManova.columns:
        for i in permManova.index:
            if permManova.loc[i,c]>0.05:
                G.add_edge(i, c) 
    pos=nx.drawing.spring_layout(G)
    plt.figure(figsize=(5,5))
    nx.draw_networkx(G,pos=pos,with_labels=False,linewidths=1, alpha=1,node_color=[palet[n] for n in G.nodes],font_size=15)
    plt.axis('off')
    plt.savefig(output+'network.pdf')


if __name__ == '__main__':
    a = get_params(sys.argv[1:])
    if a.o[-1]!='/':
        a.o=a.o+'/'
    
    data=pd.read_csv(a.i)
    data=data.set_index('genome')
    
    tree=Phylo.read(a.t, 'newick')
    phylDist=getDistMatrix(tree)
    phylDist.to_csv(a.o+'phyldistmatrix.csv')
    phylPCA=doPCA(phylDist)
    metadata=phylPCA.merge(data[['lifestyle']], left_index=True,right_index=True)
    metadata.rename(index=str, columns={[c for c in metadata.columns if 'PC1' in c][0]:'PC1',[c for c in metadata.columns if 'PC2' in c][0]:'PC2'}).to_csv(a.o+'metadata.csv')
    
    order=sorted(list(set(metadata['lifestyle'])))
    if a.colors=='':
        palet={}
    else:
        palet={c.split(':')[0]:c.split(':')[1] for c in a.colors.split(',')}
        
    
    data=data.drop(columns=['lifestyle'])
    data=data.reindex(metadata.index)
    data.to_csv(a.o+'data.csv')

    runStat(a.o)
    palet=plotPCA(metadata, palet, order, a.o)
    permManova=plotPvalMatrix(a.o, palet)
    plotNetwork(a.o, permManova, palet)
    plotPhylPCA(metadata, palet, order, a.o)




