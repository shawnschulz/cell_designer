import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import pandas as pd
import umap

import pickle

X_test = # this is inputted by a user requesting it, can be just one cell or many cells
y_train = # labels may or may not need to be loaded separately
# Going to try loading via io a pickle of the umap object for ts dataset,
# To make this faster we might want to have a program that already
# has this loaded
file = open('ts_reducer.pickle', 'rb')
ts_umap = pickle.load(file)

# creates a new embedding based on the ts embeddings
new_embedding = ts_umap.transform(X_test)

# now we need to use knn, SVC or some other unsupervised classifier to check closeness of points to original ts_umap embedding
svc = SVC().fit(new_embedding, ts_umap)
knn = KNeighborsClassifier().fit(new_embedding, ts_umap)

# Regardless of which one we choose we need to map the label numbers to the true labels from ts_umap, then map those to the
# cell names
