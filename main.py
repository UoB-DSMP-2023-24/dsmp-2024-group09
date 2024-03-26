# We should install pandas and tcrdist3 packages for our python environment first
import pandas as pd
from tcrdist.repertoire import TCRrep

# read the dataset
df = pd.read_csv("vdjdb.csv", sep = ",")

# search TRA data and rename some columns
df_alpha = df[df["gene"] == "TRA"]
df_alpha = df_alpha.rename(columns = {'v.segm':'v_a_gene', 'j.segm':'j_a_gene'})

# search TRB data and rename some columns
df_beta = df[df["gene"] == "TRB"]
df_beta = df_beta.rename(columns = {'v.segm':'v_b_gene', 'j.segm':'j_b_gene'})

# use method in tcrdist3 package to set parameter for alpha chain.
tr_alpha = TCRrep(cell_df = df_alpha,
            organism='human',
            chains = ['alpha'],
            infer_all_genes = True,
            compute_distances = True,
            deduplicate=False,
            db_file = 'alphabeta_gammadelta_db.tsv')

# use method in tcrdist3 package to set parameter for beta chain.
tr_beta = TCRrep(cell_df = df_beta,
            organism='human',
            chains = ['beta'],
            infer_all_genes = True,
            compute_distances = True,
            deduplicate=False,
            db_file = 'alphabeta_gammadelta_db.tsv')

# use default method in tcrdist3 package to calculate the distance
tr_alpha.pw_alpha
tr_beta.pw_beta