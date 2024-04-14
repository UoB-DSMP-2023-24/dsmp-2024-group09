import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt
import scipy.sparse as sp

df = pd.read_csv("vdjdb.csv")
tcr_sequences = df["cdr3"]

label_encoder = LabelEncoder()
encoded_sequences = label_encoder.fit_transform(tcr_sequences)

one_hot_encoder = OneHotEncoder()
one_hot_encoded_sequences = one_hot_encoder.fit_transform(encoded_sequences.reshape(-1, 1))

kmeans = KMeans(n_clusters=3)
kmeans.fit(one_hot_encoded_sequences)
labels = kmeans.labels_

silhouette_avg = silhouette_score(one_hot_encoded_sequences, labels)
print("The silhouette coefficient is ", silhouette_avg)

svd = TruncatedSVD(n_components=2)
reduced_data = svd.fit_transform(one_hot_encoded_sequences)

plt.figure(figsize=(8, 6))
for i in range(3):  
    plt.scatter(reduced_data[labels == i, 0], reduced_data[labels == i, 1], label='Cluster {}'.format(i))

plt.title('TCR Clustering Results')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend()
plt.show()
