# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Dimensionality reduction analysis (MDS, PCA) on genotype dataset

import argparse
import itertools
import math
import os
import sys

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# parse command line arguments
parser = argparse.ArgumentParser(description='Filters SNP genotype data.')
parser.add_argument('--analysis', type=str, help='dimensionality reduction analysis type')
parser.add_argument('--genotype', type=str, help='Genotype dataset (tsv file)')
parser.add_argument('--title', type=str, help='dataset title')
parser.add_argument('--output', type=str, help='output file prefix (without ending)')
args = parser.parse_args()

if args.analysis is None or args.genotype is None or args.title is None or args.output is None:
    parser.print_help()
    sys.exit()

if args.analysis != 'mds' and args.analysis != 'pca':
    parser.print_help()
    sys.exit()

# read genotyping data
genotype_df = pd.read_csv(args.genotype, header=0, index_col=0, delimiter=' ')
genotype_df = genotype_df[genotype_df.columns[5:]]

genotype_df.fillna(-2, inplace=True)

# get distance/similarity measure
genotype_measure_df = None
eigen_value = None
eigen_vector = None

genotype_measure_df = pd.DataFrame(index=genotype_df.index, columns=genotype_df.index)

number_markers = len(genotype_df.iloc[0])
if args.analysis == 'mds':
    for genotype1_index in range(len(genotype_measure_df.index)):
        genotype_measure_df.iloc[genotype1_index][genotype1_index] = 0

        for genotype2_index in range(genotype1_index):
            allele_difference = (genotype_df.iloc[genotype1_index] - genotype_df.iloc[genotype2_index])
            allele_difference.replace([3, 4], 2, inplace=True)

            genotype_difference = (allele_difference.abs().sum() / (2 * number_markers))

            genotype_measure_df.iloc[genotype1_index][genotype2_index] = genotype_difference
            genotype_measure_df.iloc[genotype2_index][genotype1_index] = genotype_difference

# output genotype distance data
genotype_measure_df.to_csv(args.output + '/' + args.title + '_genotype_measure.csv', sep='\t', encoding='utf-8')

if args.analysis == 'mds':
    A = genotype_measure_df
    A = A ** 2
    n = A.shape[0]
    J_c = 1. / n * (np.eye(n) - 1 + (n - 1) * np.eye(n))
    B = -0.5 * (J_c.dot(A)).dot(J_c)
    eigen_value = np.linalg.eig(np.array(B, dtype = float))[0]
    eigen_vector = np.linalg.eig(np.array(B, dtype = float))[1].T
else: # pca
    scaler = StandardScaler() 
    scaled_data = scaler.fit_transform(genotype_df)

    pca = PCA()
    pca.fit(scaled_data.T)
    eigen_value = pca.explained_variance_
    eigen_vector = pca.components_

explained_variance = []
dimensionality_reduction_components = []
dimensionality_reduction_eigenvalue = []
dimensionality_reduction_eigenvector = []
pc_indices = []
for pc_index in range(len(eigen_value)):
    pc_indices.append('PC' + str(pc_index + 1))
    explained_variance.append(abs(eigen_value[pc_index].real) / sum(abs(eigen_value.real)) * 100)
    coordinate = (abs(eigen_value[pc_index]) ** 0.5) * eigen_vector[pc_index]
    dimensionality_reduction_components.append(coordinate.real)
    dimensionality_reduction_eigenvalue.append(abs(eigen_value[pc_index].real))
    dimensionality_reduction_eigenvector.append(eigen_vector[pc_index].real)

dimensionality_reduction_explained_variance = pd.DataFrame(data=explained_variance, columns=['percent_explained_variance'], index=pc_indices)
dimensionality_reduction_explained_variance = dimensionality_reduction_explained_variance.sort_values(['percent_explained_variance'], ascending=[False])

dimensionality_reduction_components_transformed = pd.DataFrame(data=np.array(dimensionality_reduction_components).T.tolist(), columns=pc_indices, index=genotype_df.index)
dimensionality_reduction_components_transformed.reindex(dimensionality_reduction_explained_variance.index, axis=1)

dimensionality_reduction_eigenvalue_transformed = pd.DataFrame(data=dimensionality_reduction_eigenvalue, columns=['eigenvalue'], index=pc_indices)
dimensionality_reduction_eigenvalue_transformed = dimensionality_reduction_eigenvalue_transformed.sort_values(['eigenvalue'], ascending=[False])

dimensionality_reduction_eigenvector_transformed = pd.DataFrame(data=np.array(dimensionality_reduction_eigenvector).T.tolist(), columns=pc_indices, index=genotype_df.index)
dimensionality_reduction_eigenvector_transformed.reindex(dimensionality_reduction_explained_variance.index, axis=1)

dimensionality_reduction_explained_variance.index = pc_indices
dimensionality_reduction_components_transformed.columns = pc_indices
dimensionality_reduction_eigenvector_transformed.columns = pc_indices

dimensionality_reduction_explained_variance.to_csv(args.output + '/' + args.title + '_explained_variance.csv', encoding = 'utf-8', sep = '\t')
dimensionality_reduction_components_transformed.to_csv(args.output + '/' + args.title + '_coordinates.csv', encoding = 'utf-8', sep = '\t')
dimensionality_reduction_eigenvalue_transformed.to_csv(args.output + '/' + args.title + '_eigenvalues.csv', encoding = 'utf-8', sep = '\t')
dimensionality_reduction_eigenvector_transformed.to_csv(args.output + '/' + args.title + '_eigenvectors.csv', encoding = 'utf-8', sep = '\t')


