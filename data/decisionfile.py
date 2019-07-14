import sys
import os
import re
import pandas as pd
import numpy as np
import pdb
import math
from numpy import linalg as LA
from sklearn.tree import DecisionTreeClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import sklearn.preprocessing
import sklearn.ensemble


# Read training data and partition into features and target label
data = pd.read_csv("data/data/train.csv").as_matrix()
Ytrn = data[:,0]
Xtrn = np.delete(data,0,1)
Xtrn = sklearn.preprocessing.normalize(Xtrn, norm='l2', axis=1, copy=True)

# Read validation data and partition into features and target label
data = pd.read_csv("data/data/validation_1.csv").as_matrix()
Yval1 = data[:,0]
Xval1 = np.delete(data,0,1)
Xval1 = sklearn.preprocessing.normalize(Xval1, norm='l2', axis=1, copy=True)

data = pd.read_csv("data/data/validation_2.csv").as_matrix()
Yval2 = data[:,0]
Xval2 = np.delete(data,0,1)
Xval2 = sklearn.preprocessing.normalize(Xval2, norm='l2', axis=1, copy=True)

data = pd.read_csv("data/data/validation_3.csv").as_matrix()
Yval3 = data[:,0]
Xval3 = np.delete(data,0,1)
Xval3 = sklearn.preprocessing.normalize(Xval3, norm='l2', axis=1, copy=True)

# Read test data and partition into features and target label
data = pd.read_csv("data/data/testfile.csv").as_matrix()
Ytst = data[:,0]
Xtst = np.delete(data,0,1)
Xtst = sklearn.preprocessing.normalize(Xtst, norm='l2', axis=1, copy=True)

Xtrn=np.concatenate((Xtrn,Xval1,Xval2,Xval3))
Ytrn=np.concatenate((Ytrn,Yval1,Yval2,Yval3))


# Build a simple depth-5 decision tree with information gain split criterion
#dtree = DecisionTreeClassifier(criterion="entropy",max_depth=8)
dtree = sklearn.ensemble.RandomForestClassifier(criterion="entropy",max_depth=13,max_features=80)
dtree.fit(Xtrn,Ytrn)

# function score runs prediction on data, and outputs accuracy. 
# If you need predicted labels, use "predict" function
print "\n"

labels = dtree.predict(Xtst)
         
print "training accuracy: ",
print dtree.score(Xtrn,Ytrn)

print "validation accuracy 1 : ",
print dtree.score(Xval1,Yval1)

print "validation accuracy 2 : ",
print dtree.score(Xval2,Yval2)

print "validation accuracy 3 : ",
print dtree.score(Xval3,Yval3)

'''
training accuracy:  0.864357287146
validation accuracy 1 :  0.787526250875
validation accuracy 2 :  0.745924864162
validation accuracy 3 :  0.709856995233
----------------------------


max feature = 40
training accuracy:  0.912058241165
validation accuracy 1 :  0.801660055335
validation accuracy 2 :  0.749991666389
validation accuracy 3 :  0.704223474116

80 is fine 
'''
