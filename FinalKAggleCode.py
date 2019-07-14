
'''
The below implementation is the solution of image comparison problem statement from Kaggle.
The link for the problem statement in the description

'''
import numpy as np
import math
import copy
import sys
import re
import pandas as pd
from sklearn import svm


if __name__ == '__main__':
    
    trainSample = pd.read_csv('data/data/train.csv', header = None).as_matrix()
    Ytrn = trainSample[:, 0]
    
    indices = np.where(trainSample[:,0] == 1)
    Strnmatch = trainSample[indices] 
    
    Xtrn = np.delete(trainSample, 0, 1)
    Strnmatch = np.delete(Strnmatch, 0, 1)

    aMAtrix = Strnmatch[:, range(73)]
    bMatrix = Strnmatch[:, np.array(range(73)) + 73]
    DiffMatrix = abs(aMAtrix - bMatrix) 
    sumDiff = np.sum(DiffMatrix, axis = 0)
    print sumDiff
  
    w = np.argsort(sumDiff)[:61]
    w1 = w + 73

    choosedIndices1 = np.zeros(146)
  
    for i in range(73):
        choosedIndices1[i] = np.corrcoef(Strnmatch[:, i], Strnmatch[:, i + 73])[0][1]
        choosedIndices1[i + 73] = np.corrcoef(Strnmatch[:,i], Strnmatch[:, i + 73])[0][1]
    
    descol = np.where(choosedIndices1 >= 0.4)
  
    choosedIndices = np.concatenate((w, w1))
    indices = np.intersect1d(choosedIndices, descol[0])
    
    
    print choosedIndices
    print len(choosedIndices)
    print descol
    print len(descol[0])

    S = Xtrn[:,indices]
    
    # Test Sample
    testSample = pd.read_csv('data/data/testfile.csv', header = None).as_matrix()
    idx = copy.deepcopy(testSample[:, 0])
    idx = idx.astype(int)
    Xtst = np.delete(testSample, 0, 1)
    Stst = Xtst[:, indices]

    print len(Stst)
    
    # Validation Sample 1
    validSample1 = pd.read_csv('data/data/validation_1.csv', header = None).as_matrix()
    Yval1 = validSample1[:, 0]
    Xval1 = np.delete(validSample1, 0, 1)
    Sval1 = Xval1[:, indices]
    
    # Validation Sample 2
    validSample2 = pd.read_csv('data/data/validation_2.csv', header = None).as_matrix()
    Yval2 = validSample2[:, 0]
    Xval2 = np.delete(validSample2, 0, 1)
    Sval2 = Xval2[:, indices]
    
    # Validation Sample 3
    validSample3 = pd.read_csv('data/data/validation_3.csv', header = None).as_matrix()
    Yval3 = validSample3[:, 0]
    Xval3 = np.delete(validSample3, 0, 1)
    Sval3 = Xval3[:, indices]
    
    Xfinal = np.concatenate((S, Sval1, Sval2, Sval3))
    Yfinal = np.concatenate((Ytrn, Yval1, Yval2, Yval3))

    clf = svm.SVC()
    clf.fit(S, Ytrn)
    
    print "SVM Train Score ",clf.score(S,Ytrn)
    print "SVM validation 1 Score ",clf.score(Sval1,Yval1)
    print "SVM validation 2 Score ",clf.score(Sval2,Yval2)
    print "SVM validation 3 Score ",clf.score(Sval3,Yval3)
    
    labels = clf.predict(Stst)
    print labels
    labels = labels.astype(int)

    print len(idx)  
    print len(labels)
    df = pd.DataFrame(labels, columns = ["TARGET"])
    df.to_csv('svmoutputdiffcov.csv', index_label = "ID")
    
