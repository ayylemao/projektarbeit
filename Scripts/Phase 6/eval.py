import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
from sklearn.metrics import confusion_matrix
import math
import multiprocessing


#max atoms: 205504
path = "FPVectors/"
lcut_off = 140
n_atoms = 79184
atoms_per_structure = 16
training_cuttoff = math.floor(2*n_atoms/3)


vectors = np.loadtxt("testprints.dat")




gt = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
gt_long = np.tile(gt, math.ceil(n_atoms/atoms_per_structure))
g_truth = np.empty(n_atoms)
g_truth = gt_long[0:n_atoms]

trys = 50
acc = 0
con_matrix = np.zeros([2,2])

print("Training set size: " + str(training_cuttoff))
print("Validation set size: " + str(n_atoms - training_cuttoff))
print(g_truth.shape)
print(vectors.shape)
vectors, g_truth = shuffle(vectors, g_truth)
training = vectors[0:training_cuttoff, :]
training_labels = g_truth[0:training_cuttoff]
test = vectors[training_cuttoff+1:, :]
test_labels = g_truth[training_cuttoff+1:]

clf = RandomForestClassifier(n_estimators=500, max_depth=7, random_state=0, verbose=True)
clf.fit(training, training_labels)

pred = clf.predict(test)
true = np.sum(test_labels == pred)
print("The ratio of correct guesses is: " + str(true/len(pred)))
print(confusion_matrix(test_labels, pred))


# acc += true/len(pred)
# con_matrix += confusion_matrix(test_labels, pred)
# print("The average ratio of correct guesses is: " + str(acc/trys))
# print("Average confusion matrix:")
# print(con_matrix/trys)