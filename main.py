# Task 3 (the following part of code is finished by Letian Zhang).

import pandas as pd
from tcrdist.repertoire import TCRrep

# read the dataset
df = pd.read_csv(r"vdjdb.csv", sep = ",", usecols=(0,1,2,3,4,5,9,10,16))

# remove vdjdb.score = 0
df1 = df[ df['vdjdb.score'] != 0]
# df1 = df1[ df1['vdjdb.score'] != 1]

# rename some columns for df1
df1 = df1.rename(columns = {'complex.id':'clone_id',
                            'species':'subject',
                            'antigen.gene':'epitope'})

# alpha chain
df1_alpha = df1[ df1['gene'] == 'TRA']
df1_alpha = df1_alpha.dropna(axis=0, how='any')

# beta chain
df1_beta = df1[ df1['gene'] == 'TRB']
df1_beta = df1_beta.dropna(axis=0, how='any')

# rename some columns for df1_alpha and df1_beta
df1_alpha = df1_alpha.rename(columns = {'v.segm':'v_a_gene',
                                        'j.segm':'j_a_gene',
                                        'cdr3':'cdr3_a_aa'})
df1_beta = df1_beta.rename(columns = {'v.segm':'v_b_gene',
                                      'j.segm':'j_b_gene',
                                      'cdr3':'cdr3_b_aa'})

# remove the column 'gene'
df1_alpha = df1_alpha.drop(['gene'], axis=1)
df1_beta = df1_beta.drop(['gene'], axis=1)

# for all the alpha chains
# compute the distance matrix for the alpha chains (mouse)
tr_mouse_a = TCRrep(cell_df = df1_alpha,
                    organism = 'mouse',
                    chains = ['alpha'],
                    db_file = 'alphabeta_gammadelta_db.tsv')
tr_mouse_alpha = tr_mouse_a.pw_alpha
tr_mouse_alpha_color = tr_mouse_a.clone_df
print(tr_mouse_alpha)
#print(tr_mouse_a.pw_cdr3_a_aa)

# compute the distance matrix for the alpha chains (human)
tr_human_a = TCRrep(cell_df = df1_alpha,
                    organism = 'human',
                    chains = ['alpha'],
                    db_file = 'alphabeta_gammadelta_db.tsv')
tr_human_alpha = tr_human_a.pw_alpha
tr_human_alpha_color = tr_human_a.clone_df
print(tr_human_alpha)
#print(tr_mouse_a.pw_cdr3_a_aa)

# for all the beta chains
# compute the distance matrix for the beta chains (mouse)
tr_mouse_b = TCRrep(cell_df = df1_beta,
                    organism = 'mouse',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv')
tr_mouse_beta = tr_mouse_b.pw_beta
tr_mouse_beta_color = tr_mouse_b.clone_df
print(tr_mouse_beta)
#print(tr_mouse_b.pw_cdr3_b_aa)

# compute the distance matrix for the beta chains (human)
tr_human_b = TCRrep(cell_df = df1_beta,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv')
tr_human_beta = tr_human_b.pw_beta
tr_human_beta_color = tr_human_b.clone_df
print(tr_human_beta)
#print(tr_human_b.pw_cdr3_b_aa)

# remove clone_id = 0
df1_alpha = df1_alpha[ df1_alpha['clone_id'] != 0]
df1_beta = df1_beta[ df1_beta['clone_id'] != 0]

# inner connect two dataframe
df1 = pd.merge(df1_alpha, df1_beta)

# only for alpha and beta chains that can be paired
# compute the distance matrix for the alpha and the beta chains (mouse)
tr_mouse = TCRrep(cell_df = df1,
                  organism = 'mouse',
                  chains = ['alpha','beta'],
                  db_file = 'alphabeta_gammadelta_db.tsv')
tr_mouse_alpha1 = tr_mouse.pw_alpha
print(tr_mouse_alpha)
tr_mouse_beta1 = tr_mouse.pw_beta
print(tr_mouse_beta)
#print(tr_mouse.pw_cdr3_a_aa)
#print(tr_mouse.pw_cdr3_b_aa)
tr_mouse_alpha_beta_color = tr_mouse.clone_df
tr_mouse_alpha_beta = tr_mouse.pw_alpha + tr_mouse.pw_beta
print(tr_mouse_alpha_beta)

# compute the distance matrix for the alpha and the beta chains (human)
tr_human = TCRrep(cell_df = df1,
                  organism = 'human',
                  chains = ['alpha','beta'],
                  db_file = 'alphabeta_gammadelta_db.tsv')
tr_human_alpha1 = tr_human.pw_alpha
print(tr_human_alpha)
tr_human_beta1 = tr_human.pw_beta
print(tr_human_beta)
#print(tr_human.pw_cdr3_a_aa)
#print(tr_human.pw_cdr3_b_aa)
tr_human_alpha_beta_color = tr_human.clone_df
tr_human_alpha_beta = tr_human.pw_alpha + tr_human.pw_beta
print(tr_human_alpha_beta)

# transfer matrix to csv file
tr_mouse_alpha_beta_df = pd.DataFrame(tr_mouse_alpha_beta)
tr_mouse_alpha_beta_df.to_csv('tr_mouse_alpha_beta.csv', index = False)
print("Matrix saved as 'tr_mouse_alpha_beta.csv'")

tr_human_alpha_beta_df = pd.DataFrame(tr_human_alpha_beta)
tr_human_alpha_beta_df.to_csv('tr_human_alpha_beta.csv', index = False)
print("Matrix saved as 'tr_human_alpha_beta.csv'")


# Task 4 (the following part of code is finished by Uchit Bhadauriya).


#Task4

import matplotlib.pyplot as plt
import seaborn as sns
from umap import UMAP
from sklearn.preprocessing import StandardScaler
import pandas as pd

# Mouse

##alpha chains

# Scale the data
scaler = StandardScaler()
tr_mouse_alpha_scaled = scaler.fit_transform(tr_mouse_alpha)

# Fit UMAP
reducer = UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding_alpha = reducer.fit_transform(tr_mouse_alpha_scaled)

# Create a DataFrame for the embedding
embedding_alpha_df = pd.DataFrame(embedding_alpha, columns=['UMAP-1', 'UMAP-2'])
embedding_alpha_df['antigen.epitope'] = tr_mouse_alpha_color['antigen.epitope']

# Plot the UMAP embedding
plt.figure(figsize=(12, 10))  # Increase figure size to make room for the legend
scatter_plot = sns.scatterplot(data=embedding_alpha_df, x='UMAP-1', y='UMAP-2', hue='antigen.epitope', s=50, alpha=0.6)

# Here we add the legend outside the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('UMAP: Data Visualization by Antigen Epitope - Mouse Alpha Chains')
plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')

# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

##beta chains


# Scale the data
scaler = StandardScaler()
tr_mouse_beta_scaled = scaler.fit_transform(tr_mouse_beta)

# Fit UMAP
reducer = UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding_beta = reducer.fit_transform(tr_mouse_beta_scaled)

# Create a DataFrame for the embedding
embedding_beta_df = pd.DataFrame(embedding_beta, columns=['UMAP-1', 'UMAP-2'])
embedding_beta_df['antigen.epitope'] = tr_mouse_alpha_beta_color['antigen.epitope']

# Plot the UMAP embedding
plt.figure(figsize=(12, 10))  # Increase figure size to make room for the legend
scatter_plot = sns.scatterplot(data=embedding_df, x='UMAP-1', y='UMAP-2', hue='antigen.epitope', s=50, alpha=0.6)

# Here we add the legend outside the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('UMAP: Data Visualization by Antigen Epitope - Mouse Beta Chains')
plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')
# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

##alpha+beta chains

# Scale the data
scaler = StandardScaler()
tr_mouse_scaled = scaler.fit_transform(tr_mouse_alpha_beta)

# Fit UMAP
reducer = UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding = reducer.fit_transform(tr_mouse_scaled)

# Create a DataFrame for the embedding
embedding_df = pd.DataFrame(embedding, columns=['UMAP-1', 'UMAP-2'])
embedding_df['antigen.epitope'] = tr_mouse_alpha_beta_color['antigen.epitope']

# Plot the UMAP embedding
plt.figure(figsize=(12, 10))  # Increase figure size to make room for the legend
scatter_plot = sns.scatterplot(data=embedding_df, x='UMAP-1', y='UMAP-2', hue='antigen.epitope', s=50, alpha=0.6)

# Here we add the legend outside the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('UMAP: Data Visualization by Antigen Epitope - Mouse Alpha+Beta Chains')
plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')
# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

#Human

##alpha chains

# Scale the data
scaler = StandardScaler()
tr_human_alpha_scaled = scaler.fit_transform(tr_human_alpha)

# Fit UMAP
reducer = UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding_human_alpha = reducer.fit_transform(tr_human_alpha_scaled)

# Create a DataFrame for the embedding
embedding_human_alpha_df = pd.DataFrame(embedding_human_alpha, columns=['UMAP-1', 'UMAP-2'])
embedding_human_alpha_df['antigen.epitope'] = tr_human_alpha_color['antigen.epitope']

# Plot the UMAP embedding
plt.figure(figsize=(12, 10))  # Increase figure size to make room for the legend
scatter_plot = sns.scatterplot(data=embedding_human_alpha_df, x='UMAP-1', y='UMAP-2', hue='antigen.epitope', s=50, alpha=0.6)

# Here we add the legend outside the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('UMAP: Data Visualization by Antigen Epitope - Human Alpha Chains')
plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')
# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

##beta chains

# Scale the data
scaler = StandardScaler()
tr_mouse_human_beta_scaled = scaler.fit_transform(tr_human_alpha_beta)

# Fit UMAP
reducer = UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding_human_beta = reducer.fit_transform(tr_mouse_human_beta_scaled)

# Create a DataFrame for the embedding
embedding_human_beta_df = pd.DataFrame(embedding_human_beta, columns=['UMAP-1', 'UMAP-2'])
embedding_human_beta_df['antigen.epitope'] = tr_human_beta_color['antigen.epitope']

# Plot the UMAP embedding
plt.figure(figsize=(12, 10))  # Increase figure size to make room for the legend
scatter_plot = sns.scatterplot(data=embedding_human_beta_df, x='UMAP-1', y='UMAP-2', hue='antigen.epitope', s=50, alpha=0.6)

# Here we add the legend outside the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('UMAP: Data Visualization by Antigen Epitope - Human Beta Chains')
plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')
# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

##alpha+beta chains

# Scale the data
scaler = StandardScaler()
tr_human_alpha_beta_scaled = scaler.fit_transform(tr_human_alpha_beta)

# Fit UMAP
reducer = UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding_human_alpha_beta = reducer.fit_transform(tr_human_alpha_beta_scaled)

# Create a DataFrame for the embedding
embedding_human_alpha_beta_df = pd.DataFrame(embedding_human_alpha_beta, columns=['UMAP-1', 'UMAP-2'])
embedding_human_alpha_beta_df['antigen.epitope'] = tr_human_alpha_beta_color['antigen.epitope']

# Plot the UMAP embedding
plt.figure(figsize=(12, 10))  # Increase figure size to make room for the legend
scatter_plot = sns.scatterplot(data=embedding_df, x='UMAP-1', y='UMAP-2', hue='antigen.epitope', s=50, alpha=0.6)

# Here we add the legend outside the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('UMAP: Data Visualization by Antigen Epitope - Human Alpha+Beta Chains')
plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')
# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()


# Task 5 (the following part of code is finished by Yu Gu).


import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from sklearn.manifold import MDS

# Mouse

# Read distance matrix.
mouse = pd.read_csv('tr_mouse_alpha_beta.csv')
distance_matrix = mouse.values

# For the distance matrix, we use precomputed as the measure of distance.
dbscan = DBSCAN(eps = 120, min_samples = 35, metric = "precomputed")
clusters = dbscan.fit_predict(distance_matrix)
# print("Mouse Cluster labels:", clusters)

# For visualization, we use the MDS multidimensional scaling technique to reduce the distance matrix to a two-dimensional space.
mds = MDS(n_components = 2, dissimilarity = "precomputed", random_state = 42)
mds_transformed = mds.fit_transform(distance_matrix)

# Calculate silhouette coefficients (only calculate non-noisy clustered data).
valid_clusters = clusters[clusters != -1]
valid_mds_transformed = mds_transformed[clusters != -1]

if len(set(valid_clusters)) > 1:  # There must be at least 2 clusters to calculate the silhouette coefficient.
    silhouette = silhouette_score(valid_mds_transformed, valid_clusters)
    print("The silhouette coefficient of mouse is ", silhouette)
else:
    silhouette = -1  

# Visualize clustering results.
plt.figure(figsize=(10, 8))
plt.scatter(valid_mds_transformed[:, 0], valid_mds_transformed[:, 1], c = valid_clusters, cmap='viridis', marker='o', s=50)
plt.title('DBSCAN Clustering (Silhouette Coefficient: {:.2f})'.format(silhouette))
plt.xlabel('MDS Dimension 1')
plt.ylabel('MDS Dimension 2')
plt.colorbar(label='Cluster Label')
plt.show()

# Human

# Read distance matrix.
human = pd.read_csv('tr_human_alpha_beta.csv')
distance_matrix = human.values

# For the distance matrix, we use precomputed as the measure of distance.
dbscan = DBSCAN(eps = 120, min_samples = 35, metric = "precomputed")
clusters = dbscan.fit_predict(distance_matrix)
# print("Human Cluster labels:", clusters)

# For visualization, we use the MDS multidimensional scaling technique to reduce the distance matrix to a two-dimensional space.
mds = MDS(n_components = 2, dissimilarity = "precomputed", random_state = 42)
mds_transformed = mds.fit_transform(distance_matrix)

# Calculate silhouette coefficients (only calculate non-noisy clustered data).
valid_clusters = clusters[clusters != -1]
valid_mds_transformed = mds_transformed[clusters != -1]

if len(set(valid_clusters)) > 1:  # There must be at least 2 clusters to calculate the silhouette coefficient.
    silhouette = silhouette_score(valid_mds_transformed, valid_clusters)
    print("The silhouette_score of human is ", silhouette)
else:
    silhouette = -1  

# Visualize clustering results.
plt.figure(figsize=(10, 8))
plt.scatter(valid_mds_transformed[:, 0], valid_mds_transformed[:, 1], c = valid_clusters, cmap='viridis', marker='o', s=50)
plt.title('DBSCAN Clustering (Silhouette Coefficient: {:.2f})'.format(silhouette))
plt.xlabel('MDS Dimension 1')
plt.ylabel('MDS Dimension 2')
plt.colorbar(label='Cluster Label')
plt.show()


# Task 6 (the following part of code is finished by Lubin Wan).


import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, recall_score, f1_score, classification_report

# Function definition: perform grid search.
def optimize_knn_parameters(X, y):
    knn = KNeighborsClassifier(metric='precomputed', algorithm='brute')
    param_grid = {
        'n_neighbors': [2, 3, 5, 7, 10, 12, 15], # Different number of neighbors.
        'weights': ['uniform', 'distance'], # Neighbor weight strategy.
        }
    grid_search = GridSearchCV(knn, param_grid, cv=4, scoring='accuracy')
    grid_search.fit(X, y)
    return grid_search.best_estimator_, grid_search.best_params_

# Function definition: Split training and test data and perform model optimization.
def perform_analysis(X, labels):
    indices = np.arange(len(labels))
    train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=1)
    X_train = X[np.ix_(train_indices, train_indices)]
    X_test = X[np.ix_(test_indices, train_indices)]
    y_train = labels[train_indices]
    y_test = labels[test_indices]
    best_knn, best_params = optimize_knn_parameters(X_train, y_train)
    print(f'Best parameters: {best_params}')
    y_pred = best_knn.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred, average='macro')
    f1 = f1_score(y_test, y_pred, average='macro')
    report = classification_report(y_test, y_pred)
    return accuracy, recall, f1, report

# Get label data and perform analysis.
labels = pd.Categorical(tr_human_alpha_beta_color['epitope']).codes
accuracy_human, recall_human, f1_human, report_human = perform_analysis(tr_human_alpha_beta, labels)
print(f'Accuracy-human: {accuracy_human:.2f}')
print(f'Recall-human: {recall_human:.2f}')
print(f'F1 Score-human: {f1_human:.2f}')
print("Classification Report-human:\n", report_human)

labels1 = pd.Categorical(tr_mouse_alpha_beta_color['epitope']).codes
accuracy_mouse, recall_mouse, f1_mouse, report_mouse = perform_analysis(tr_mouse_alpha_beta, labels1)
print(f'Accuracy-mouse: {accuracy_mouse:.2f}')
print(f'Recall-mouse: {recall_mouse:.2f}')
print(f'F1 Score-mouse: {f1_mouse:.2f}')
print("Classification Report-mouse:\n", report_mouse)

# Process distances and plot with the numeric labels
process_tcr_distances(tr_human_alpha, tr_human_alpha_color['epitope_codes'], "Human Alpha Chain")
process_tcr_distances(tr_human_beta, tr_human_beta_color['epitope_codes'], "Human Beta Chain")
process_tcr_distances(tr_human_alpha_beta, tr_human_alpha_beta_color['epitope_codes'], "Human Combined Alpha-Beta Chain")
