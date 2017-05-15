# -*- coding: utf-8 -*-
import numpy as np
import csv
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
from entrez2uniprot import convert_thread
import threading

gene_exp_path = '../data/geneExp.txt'
subtype_path = '../data/subtype.txt'
survival_path = '../data/survival.txt'
ppi_path = '../data/in_bio_map_network_50.mat'

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

def get_genes(filename):
    threads = []
    entrez_ids = []
    uniprot_ids = dict()
    del_index = []
    
    with open(filename, 'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            entrez_ids.append(row[0])
    
    entrez_ids = entrez_ids[1:]
    for i, gene in enumerate(entrez_ids):
        index = gene.find('|')
        if index != -1:
            entrez_ids[i] = gene[index+1:]
        else:
            del_index.append(i)
            
    chunks = [entrez_ids[x:x+100] for x in xrange(0, len(entrez_ids), 100)]
    for i, ids in enumerate(chunks):
        converter = convert_thread(i, ids)
        threads.append(converter)
        converter.start()
        
    for t in threads:
        t.join()
    
    for t in threads:
        uniprot_ids.update(t.uniprot_ids)
    
    return entrez_ids, uniprot_ids, del_index, threads
        
if __name__ == '__main__':
    columns, ncols = get_fileinfo(gene_exp_path)
    entrez_ids, uniprot_ids, del_index, threads = get_genes(gene_exp_path)
#    X = np.loadtxt(fname=gene_exp_path, delimiter='\t', usecols=range(1,ncols)).T
#    clinic = parse_subtypes(subtype_path)
#    survival = parse_subtypes(survival_path)
#    subtypes, sub_labels, survival_labels, del_index = get_labels(columns, clinic, survival)
#    subtype_labels = label_binarize(sub_labels, classes=subtypes)
#    X = np.delete(X, del_index, axis=0)
#    Y = np.concatenate((subtype_labels,np.row_stack(survival_labels)),axis=1)
#    
#    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.4, random_state=0)

    
    
    
    