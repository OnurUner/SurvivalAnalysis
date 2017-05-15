# -*- coding: utf-8 -*-
import numpy as np
from sklearn import datasets
from sklearn.feature_selection import GenericUnivariateSelect, mutual_info_classif, f_classif
from sklearn.ensemble import ExtraTreesClassifier
from load_data import load_dataset, parse_txt
import scipy
from numpy import genfromtxt

scores_path = '..\\data\\scores.mat'
unicox_path = '..\\data\\unicox_scores.csv'

class FeatureScore:
    def __init__(self, label_name, label_index, p_values, p2scores, mut_inf_scores, feature_weights):
        self.label_name = label_name
        self.label_index = label_index
        self.p_values = p_values
        self.p2scores = p2scores
        self.mut_inf_scores = mut_inf_scores
        self.feature_weights = feature_weights
        

def test():
    # The iris dataset
    iris = datasets.load_iris()
    
    # Some noisy data not correlated
    E = np.random.uniform(0, 0.1, size=(len(iris.data), 20))
    
    # Add the noisy data to the informative features
    X = np.hstack((iris.data, E))
    y = iris.target
    
    fselector = GenericUnivariateSelect(f_classif)
    fselector.fit(X, y)
    p2scores = -np.log10(fselector.pvalues_)
    p2scores /= p2scores.max()
    
    mutSelector = GenericUnivariateSelect(mutual_info_classif)
    mutSelector.fit(X, y)
    mutscores = mutSelector.scores_
    
    model = ExtraTreesClassifier()
    model.fit(X, y)
    importance = model.feature_importances_
    
    return p2scores, mutscores, importance

def feature_scores(X, Y):
    fselector = GenericUnivariateSelect(f_classif)
    fselector.fit(X, Y)
    p2scores = -np.log10(fselector.pvalues_)
    p2scores /= p2scores.max()
    
    mutSelector = GenericUnivariateSelect(mutual_info_classif)
    mutSelector.fit(X, Y)
    mutscores = mutSelector.scores_

    return fselector.pvalues_, p2scores, mutscores
    

if __name__ == '__main__':
    dataset = load_dataset()
    X = dataset["data"]
    Y = dataset["labels"]
    feature_score_list = []    
    weights_list = []
    
    for i in range(5):
        print 'Calculating feature scores for label', dataset['y_columns'][i]
        p_values, p2scores, mutscores = feature_scores(X, Y[:,i])
        feature_weights = (p_values < 0.05).astype(int)
        score = FeatureScore(dataset['y_columns'][i], i, p_values, p2scores, mutscores, feature_weights)
        feature_score_list.append(score)
        weights_list.append(feature_weights)
        
    mat_dict = dict()
    unicox = (genfromtxt(unicox_path, delimiter=',')[1] < 0.05).astype(int)
    mat_dict["survival_day"] = unicox
    for i, l in enumerate(weights_list):
        mat_dict[dataset['y_columns'][i]] = l
    scipy.io.savemat(scores_path, mdict=mat_dict)
    
