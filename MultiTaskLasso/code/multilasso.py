import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import MultiTaskLasso
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

def get_fileinfo(filename):
    with open(filename) as f:
        columns = f.readline().split('\t')
    return columns,len(columns)

def parse_subtypes(filename):
    labels = dict()
    with open(filename) as f:
        content = f.readlines()
        for i in xrange(1, len(content)):
            row = content[i].split('\t')
            if len(row[1]) > 1:
                cls = row[0].strip()[0:row[0].strip().rfind('-')]
                labels[cls] = row[1].strip()
    return labels

def get_labels(columns, clinic, survival):
    c_label = []
    s_label = []
    del_ind = []
    for i in xrange(1,len(columns)):
        column = columns[i]
        inSurvival = False
        for patient in survival.keys():
            if column.startswith(patient):
                s_label.append(int(survival[patient]))
                inSurvival = True
                break
            
        if not inSurvival:
            del_ind.append(i-1)
            continue
        
        for patient in clinic.keys():
            if column.startswith(patient):
                c_label.append(clinic[patient])
                break
 
    return sorted(set(c_label)), np.asarray(c_label), np.asarray(s_label), del_ind
        

def plot_coef(coef_multi_task_lasso_):
    fig = plt.figure()
    plt.spy(coef_multi_task_lasso_)
    plt.xlabel('Feature')
    plt.ylabel('Time (or Task)')
    plt.xlim(0, 500)
    fig.suptitle('Coefficient non-zero location')

def get_stats(coef):
    sum_coef = np.sum(coef, axis=0)
    return np.where(sum_coef == 0)[0]
    
    
if __name__ == '__main__':
#    columns, ncols = get_fileinfo('geneExp.txt')
#    X = np.loadtxt(fname='geneExp.txt', delimiter='\t', usecols=range(1,ncols)).T
#    clinic = parse_subtypes('subtype.txt')
#    survival = parse_subtypes('survival.txt')
#    subtypes, sub_labels, survival_labels, del_index = get_labels(columns, clinic, survival)
#    subtype_labels = label_binarize(sub_labels, classes=subtypes)
#    X = np.delete(X, del_index, axis=0)
#    Y = np.concatenate((subtype_labels,np.row_stack(survival_labels)),axis=1)
#    
#    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.4, random_state=0)

    import sys
    sys.path.insert(0, 'C:\\r workspace\\MultiSconES\\py')
    from load_data import load_dataset
    
    dataset = load_dataset()
    X = dataset["data"]
    Y = dataset["labels"]
    
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.33, random_state=42)
    
    clf = MultiTaskLasso(alpha=1)
    print "train start"
    clf.fit(X_train, Y_train)
    print "train end"
    print "coef start"
    coef_multi_task_lasso_ = clf.coef_
    print "coef end"
    plot_coef(coef_multi_task_lasso_)
    zero_coefs = get_stats(coef_multi_task_lasso_)
    print len(zero_coefs)
    
    Y_pred = clf.predict(X_test)
    clf_score = clf.score(X_test, Y_test)
    score = r2_score(Y_test[:,5], Y_pred[:,5])
    
    
    
    
    