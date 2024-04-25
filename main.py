import pandas as pd
from tcrdist.repertoire import TCRrep

# Task 3 (the following part of code is finished by Letian Zhang).

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




# Task 4

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

def plot_dimensionality_reduction(X, labels, method, title):
    """
    Helper function to plot dimensionality reduced data.
    Plots data with colors representing different labels.
    """
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X[:, 0], X[:, 1], c=labels, cmap='viridis', alpha=0.6)
    plt.colorbar(scatter)
    plt.title(f'{method} for {title}')
    plt.grid(True)
    plt.show()

def process_tcr_distances(distance_matrix, labels, title):
    """
    Processes the distance matrix using PCA, t-SNE, and UMAP and plots the results.
    """
    X = distance_matrix

    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    plot_dimensionality_reduction(X_pca, labels, "PCA", title)

    # t-SNE
    tsne = TSNE(n_components=2)
    X_tsne = tsne.fit_transform(X)
    plot_dimensionality_reduction(X_tsne, labels, "t-SNE", title)

    # UMAP
    umap = UMAP(n_components=2)
    X_umap = umap.fit_transform(X)
    plot_dimensionality_reduction(X_umap, labels, "UMAP", title)

# mouse
# Convert 'epitope' labels to categorical numeric codes
labels_numeric = pd.Categorical(df1['epitope']).codes

# Convert 'epitope' labels to categorical numeric codes for coloring
tr_mouse_alpha_color['epitope_codes'] = pd.Categorical(tr_mouse_alpha_color['epitope']).codes
tr_mouse_beta_color['epitope_codes'] = pd.Categorical(tr_mouse_beta_color['epitope']).codes
tr_mouse_alpha_beta_color['epitope_codes'] = pd.Categorical(tr_mouse_alpha_beta_color['epitope']).codes

# Process distances and plot with the numeric labels
process_tcr_distances(tr_mouse_alpha, tr_mouse_alpha_color['epitope_codes'], "Mouse Alpha Chain")
process_tcr_distances(tr_mouse_beta, tr_mouse_beta_color['epitope_codes'], "Mouse Beta Chain")
process_tcr_distances(tr_mouse_alpha_beta, tr_mouse_alpha_beta_color['epitope_codes'], "Mouse Combined Alpha-Beta Chain")

# human
# Convert 'epitope' labels to categorical numeric codes for coloring
tr_human_alpha_color['epitope_codes'] = pd.Categorical(tr_human_alpha_color['epitope']).codes
tr_human_beta_color['epitope_codes'] = pd.Categorical(tr_human_beta_color['epitope']).codes
tr_human_alpha_beta_color['epitope_codes'] = pd.Categorical(tr_human_alpha_beta_color['epitope']).codes

# Process distances and plot with the numeric labels
process_tcr_distances(tr_human_alpha, tr_human_alpha_color['epitope_codes'], "Human Alpha Chain")
process_tcr_distances(tr_human_beta, tr_human_beta_color['epitope_codes'], "Human Beta Chain")
process_tcr_distances(tr_human_alpha_beta, tr_human_alpha_beta_color['epitope_codes'], "Human Combined Alpha-Beta Chain")
