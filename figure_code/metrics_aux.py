#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.spatial.distance import braycurtis, cosine, pdist, squareform, jensenshannon
from skbio.stats.distance import mantel

h = 10**-8

def confusion_matrix(X, X0, samples, thr=0., thr0=0., to_file=''):
    '''
    Calculation confusion matrix and related statistics

    Arguments:
    X : count table
    X0 : reference count table ('true counts')
    samples : names of samples
    thr : threshold for distinguishing positives and negatives in the count table
    thr0 : threshold for distinguishing positives and negatives in the reference count table
    filepath : path to file for saving the results

    Outputs:
    pandas data frame with the results, line per sample
    '''
    TP = np.sum((X > thr) * (X0 > thr0), axis=0) #true positives
    FN = np.sum((X <= thr) * (X0 > thr0), axis=0) #false negatives
    FP = np.sum((X > thr) * (X0 <= thr0), axis=0) #false positives
    TN = np.sum((X <= thr) * (X0 <= thr0), axis=0) # true negatives

    TPR = TP / (TP + FN + h) #True positive rate = sensitivity
    FPR = FP / (FP + TN + h) #False positive rate = 1 - specificity
    TNR = TN / (TN + FP + h) #True negative rate = specificity
    FNR = FN / (FN + TP + h) #False negative rate

    PPV = TP / (TP + FP + h) #positive predictive value = precision
    NPV = TN / (TN + FN + h) #negative predictive value
    FDR = FP / (FP + TP + h) #false discovery rate
    FOR = FN / (FN + TN + h) #false omission rate

    P_T = np.mean(X > thr, axis=0) #probability of positive outcome
    P_F = np.mean(X <= thr, axis=0) #probability of negative outcome
    
    FPRA = np.sum(X * (X0 <= 0.), axis=0) #false positive relative abundance
    FNRA = np.sum((X <= 0.) * X0, axis=0) #relative abundance of false negatives
    F1 = 2. * TP / (2. * TP + FP + FN + h) # F1-score
    ACC = (TP + TN) / (TP + TN + FP + FN + h) #Accuracy
    df = pd.DataFrame({'TP' : TP,\
                       'FN' : FN,\
                       'FP' : FP,\
                       'TN': TN,\
                       'TPR': TPR,\
                       'FPR': FPR,\
                       'TNR': TNR,\
                       'FNR': FNR,\
                       'PPV': PPV,\
                       'NPV': NPV,\
                       'FDR': FDR,\
                       'FOR': FOR,\
                       'P_T': P_T,\
                       'P_F': P_F,\
                       'FPRA': FPRA,\
                       'FNRA': FNRA,\
                       'F1': F1,\
                       'ACC': ACC,\
                       'Dice': (FN + FP)/(2 * TP + FN + FP)}, index=samples)
    if to_file:
        df.to_csv(to_file, sep='\t', index=False)
    return df

def confusion_matrix_df(X, X0, thr=0., thr0=0., to_file=''):
    '''
    Calculation confusion matrix and related statistics (adapted for pandas data frames)

    Arguments:
    X : count table (pandas DataFrame)
    X0 : reference count table ('true counts', pandas DataFrame)
    thr : threshold for distinguishing positives and negatives in the count table
    thr0 : threshold for distinguishing positives and negatives in the reference count table
    filepath : path to file for saving the results

    Outputs:
    pandas data frame with the results, line per sample
    '''
    TP = np.sum((X > thr) & (X0 > thr0), axis=0) #true positives
    FN = np.sum((X <= thr) & (X0 > thr0), axis=0) #false negatives
    FP = np.sum((X > thr) & (X0 <= thr0), axis=0) #false positives
    TN = np.sum((X <= thr) & (X0 <= thr0), axis=0) # true negatives

    TPR = TP / (TP + FN + h) #True positive rate = sensitivity
    FPR = FP / (FP + TN + h) #False positive rate = 1 - specificity
    TNR = TN / (TN + FP + h) #True negative rate = specificity
    FNR = FN / (FN + TP + h) #False negative rate

    PPV = TP / (TP + FP + h) #positive predictive value = precision
    NPV = TN / (TN + FN + h) #negative predictive value
    FDR = FP / (FP + TP + h) #false discovery rate
    FOR = FN / (FN + TN + h) #false omission rate

    P_T = np.mean(X > thr, axis=0) #probability of positive outcome
    P_F = np.mean(X <= thr, axis=0) #probability of negative outcome
    
    FPRA = np.sum(X * (X0 <= 0.), axis=0) #false positive relative abundance
    F1 = 2. * TP / (2. * TP + FP + FN + h) # F1-score
    ACC = (TP + TN) / (TP + TN + FP + FN + h) #Accuracy

    Dice = (FN + FP) / (2 * TP + FN + FP)

    df = pd.concat([TP, FN, FP, TN, TPR, FPR, TNR, FNR, PPV, NPV, FDR, FOR, P_T, P_F, FPRA, F1, ACC, Dice], axis=1)
    df.columns = ["TP", "FN", "FP", "TN", "TPR", "FPR", "TNR", "FNR", "PPV", "NPV", "FDR", "FOR", "P_T", "P_F", "FPRA", "F1", "ACC", "Dice"]
    if to_file:
        df.to_csv(to_file, sep='\t', index=False)
    return df

def pcoa_dist(dat_dist, to_file='', k=2, clusters={}):
    '''
    Perform MDS/PCoA decomposition
    Arguments:
    df : data frame
    to_file : path to file for saving figure
    k : number of MDS components to calculate
    '''
    n = len(dat_dist)
    labs = dat_dist.columns.tolist()

    H = np.eye(n) - np.ones((n,n))/n
    ZZ = -0.5*H.dot(dat_dist**2).dot(H)
    W, V = eigh(ZZ, subset_by_index = (n-k, n-1))
    V = V[:,::-1]
    W = W[::-1]
    if to_file:
        fig, ax = plt.subplots(1, 1, figsize = (10, 10))
        if clusters:
            ax.scatter(V[:,0], V[:,1], marker='.')
            for (name, jj), c in zip(clusters.items(), cols):
                ax.scatter(V[jj,0], V[jj,1], color=c, label=name)
            ax.legend()
        else:
            ax.scatter(V[:,0], V[:,1])
        ax.set_xlabel('component 1')
        ax.set_ylabel('component 2')
        for j, smpl in enumerate(labs):
            ax.annotate(smpl, (V[j,0], V[j,1]), fontsize=6)

        plt.savefig(to_file)
        plt.close(fig)
    return W, V

def pairwise_nmds(sats, smpls, metric, nPC=3):
    DD = []; VV = []
    for sat in sats:
        if metric == "cosine_sqrt":
            D_bc = pd.DataFrame(squareform(pdist(np.sqrt(sat[smpls]).T, "cosine")), index=smpls, columns=smpls)
        else:
            D_bc = pd.DataFrame(squareform(pdist(sat[smpls].T, metric)), index=smpls, columns=smpls)
        _, V_bc = pcoa_dist(D_bc, k=nPC)
        V_bc = pd.DataFrame(V_bc, index=smpls, columns=["PC{}".format(j+1) for j in range(nPC)])
        DD.append(D_bc)
        VV.append(V_bc)
    return DD, VV

def intermatrix_norms(DD1, DD2, names, filepath=""):
    df = pd.DataFrame([[L1_dist(D1, D2), 
                        L2_dist(D1, D2), 
                        JS_dist(D1, D2),
                        mantel(D1, D2, method="pearson")[0], 
                        mantel(D1, D2, method="spearman")[0]] for D1, D2 in zip(DD1, DD2)], 
                        columns=["L1", "L2", "JSD", "Mantel (Pearson)", "Mantel (Spearman)"],
                        index=names)
    if filepath: df.to_csv(filepath, sep="\t")
    return df

def L1_dist(A, B):
    ii = np.triu_indices(A.shape[0], k=1)
    return np.linalg.norm(A.to_numpy()[ii] - B.to_numpy()[ii], ord=1)
    # return np.abs(A.to_numpy()[ii] - B.to_numpy()[ii]).sum().sum()

def L2_dist(A, B):
    ii = np.triu_indices(A.shape[0], k=1)
    return np.linalg.norm(A.to_numpy()[ii] - B.to_numpy()[ii])

def JS_dist(A, B):
    ii = np.triu_indices(A.shape[0], k=1)
    return jensenshannon(A.to_numpy()[ii], B.to_numpy()[ii])    

def cosine_sqrt(A, B):
    return cosine(np.sqrt(A), np.sqrt(B))

