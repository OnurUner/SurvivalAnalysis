# -*- coding: utf-8 -*-
import numpy as np
import csv
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
from entrez2uniprot import convert_thread
from scipy import io
import pickle
from collections import Counter

gene_exp_path = '../data/geneExp.txt'
subtype_path = '../data/subtype.txt'
survival_path = '../data/survival.txt'
cencoring_path = '../data/cencoring.txt'
ppi_path = '../data/ppi_genes.mat'
ppi_mat_path = '../data/ppi.mat'
dataset_path = '../data/breast_dataset.mat'
modified_ppi_path = '../data/modified_ppi.mat'

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def get_fileinfo(filename):
    with open(filename) as f:
        columns = f.readline().split('\t')
        
    patients = []
    tumors = []
    row_indices = []
    
    for i, column in enumerate(columns[1:]):
        indices = find(column, '-')
        if len(indices) > 0:
            patients.append(column[indices[1]+1: indices[2]].strip())
            tumors.append(int(column[indices[2]+1:indices[2]+3].strip()))


    for i,tumor in enumerate(tumors):
        if tumor == 1:
            row_indices.append(i)
        
    return np.array(patients)[row_indices], row_indices

def parse_txt(filename):
    patients = []
    tumors = []
    values = []
    result = dict()

    with open(filename) as f:
        content = f.readlines()
        for i in xrange(1, len(content)):
            row = content[i].split('\t')
            if len(row[1].strip()) > 0:
                indices = find(row[0], '-')
                patients.append(row[0][indices[1]+1: indices[2]].strip())
                tumors.append(int(row[0][indices[2]+1:].strip()))
                values.append(row[1].strip())
    
    del_ind = []
    for i, tumor in enumerate(tumors):
        if tumor > 1:
            del_ind.append(i)
    
    patients = np.delete(patients, del_ind)
    tumors = np.delete(tumors, del_ind)
    values = np.delete(values, del_ind)
    
    for i, patient in enumerate(patients):
        result[patient] = values[i]
    
    if len(patients) != len(result):
        print 'Warning: The number of patients is not equal to the number of dictionary keys'
        print len(patients)
        print len(result)
    
#    with open(filename) as f:
#        content = f.readlines()
#        for i in xrange(1, len(content)):
#            row = content[i].split('\t')
#            if len(row[1].strip()) > 0:
#                indices = find(row[0], '-')
#                patient = row[0][indices[1]+1: indices[2]]
#                tumor = int(row[0][indices[2]+1:].strip())
#                print tumor
#                if tumor == 1:
#                    if patient not in labels.keys():
#                        labels[patient] = row[1].strip()
#                    else:
#                        same_patient += 1
#                else:
#                    non_cancer += 1

    return result

def get_labels(patients, subtype, survival, cencor):
    subtype_labels = []
    survival_labels = []
    cencor_labels = []
    del_ind = []
    
    for i, patient in enumerate(patients):
        if (patient not in survival.keys()) or (patient not in subtype.keys()) or (patient not in cencor.keys()):
                del_ind.append(i)
    
    rows = np.delete(patients, del_ind)
    
    for i, row in enumerate(rows):
        if (patient not in survival.keys()) or (patient not in subtype.keys()) or (patient not in cencor.keys()):
            print "ERRORRR"
        subtype_labels.append(subtype[row])
        survival_labels.append(int(survival[row]))
        cencor_labels.append(int(cencor[row]))
 
    return sorted(set(subtype_labels)), np.asarray(subtype_labels), np.asarray(survival_labels), np.asarray(cencor_labels), del_ind

def get_genes(filename, only_entrez = True):
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
            
    if not only_entrez:
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
        
def load_entrez2uniprot():
    with open('entrez2uniprot.pkl', 'rb') as f:
        data = dict(pickle.load(f))
    
    result = dict()
    for key,value in data.items():
        if value not in result.values():
            result[key] = value
    return result
    
def prune_dataset(X, entrez_ids, uniprot_ids, ppi_genes):
    delete_genes = list()
    for i, ppi_gene in enumerate(ppi_genes):
        if ppi_gene not in uniprot_ids.values():
            delete_genes.append(i)
    
    print len(ppi_genes)
    ppi_genes = np.delete(ppi_genes, list(delete_genes))
    print len(ppi_genes)
   
    delete_genes = list()

    for i, entrez_id in enumerate(entrez_ids):
        if entrez_id not in uniprot_ids.keys():
           delete_genes.append(i)
        elif uniprot_ids[entrez_id] not in ppi_genes:
           delete_genes.append(i)

    print len(delete_genes)
    new_entrez_ids = [v for i, v in enumerate(entrez_ids) if i not in delete_genes]
    return np.delete(X, list(delete_genes), axis=1), new_entrez_ids, ppi_genes
            
def get_ppi_genes(uniprot_ids):
    ppi_genes = io.loadmat(ppi_path)['ppi_genes']
    delete_genes = list()
    for i, ppi_gene in enumerate(ppi_genes):
        if ppi_gene not in uniprot_ids.values():
            delete_genes.append(i)
    
    return np.delete(ppi_genes, delete_genes), delete_genes

def get_gene_order(entrez_ids, uniprot_ids, ppi_genes):
    order_list = []
    for ppi_gene in ppi_genes:
        for entrez, uniprot in uniprot_ids.iteritems():
            if ppi_gene == uniprot:
                order_list.append(entrez_ids.index(entrez))
                break
    return order_list

def convert_uni2entrez():
    from PyEntrezId import Conversion

    UniProtId = "Q9Y6X1"
    # include your email address
    Id = Conversion('dummyemail@dummybunny.info')
    EntrezID = Id.convert_uniprot_to_entrez(UniProtId)
    # Returns a string
    print EntrezID
    
def convert_entrez2uni(EntrezID):
    from PyEntrezId import Conversion

    # include your email address
    Id = Conversion('dummyemail@dummybunny.info')
    UniProtId = Id.convert_entrez_to_uniprot(EntrezID)
    # Returns a string
    print UniProtId
    return UniProtId

def get_dataset():    
    x_rows, row_indices = get_fileinfo(gene_exp_path)
    entrez_ids = get_genes(gene_exp_path)[0]
    uniprot_ids = load_entrez2uniprot()
    ppi_genes, del_ppi_genes = get_ppi_genes(uniprot_ids)
    ppi_network = io.loadmat(ppi_mat_path)['adjacency_matrix'].toarray()
    ppi_network = np.delete(ppi_network, del_ppi_genes, axis = 1)
    ppi_network = np.delete(ppi_network, del_ppi_genes, axis = 0)
    
    X = np.loadtxt(fname=gene_exp_path, delimiter='\t', usecols=[k+1 for k in row_indices]).T
    gene_orders = get_gene_order(entrez_ids, uniprot_ids, ppi_genes)
    X = X[:,gene_orders]
    x_columns = np.array(entrez_ids)
    x_columns = x_columns[gene_orders]
    
    clinic = parse_txt(subtype_path)
    survival = parse_txt(survival_path)
    cencor = parse_txt(cencoring_path)
    subtypes, sub_labels, survival_labels,cencor_labels, del_index = get_labels(x_rows, clinic, survival, cencor)
    subtype_labels = label_binarize(sub_labels, classes=subtypes)
    X = np.delete(X, del_index, axis=0)
    x_rows = np.delete(x_rows, del_index)
    Y = np.concatenate((subtype_labels,np.row_stack(survival_labels)),axis=1)
    
    zero_std = (np.std(X, axis=0) == 0)
    del_index = np.where(zero_std == True)[0].tolist()
    X = np.delete(X, del_index, axis=1)
    x_columns = np.delete(x_columns, del_index)
    ppi_genes = np.delete(ppi_genes, del_index)
    ppi_network = np.delete(ppi_network, del_index, axis=0)
    ppi_network = np.delete(ppi_network, del_index, axis=1)
#    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.3, random_state=0)
    result = dict()
#    result["train_data"] = X_train
#    result["train_labels"] = Y_train
#    result["test_data"] = X_test
#    result["test_labels"] = Y_test

    result["data"] = X
    result["labels"] = Y
    result["x_rows"] = x_rows
    result["x_columns"] = x_columns
    result["y_columns"] = subtypes + ['survival_day']
    result["cencoring"] = cencor_labels
    result["subtypes"] = sub_labels

          
    ppi = dict()
    ppi["ppi_network"] = ppi_network
    ppi["ppi_genes"] = ppi_genes
    return result, ppi

def save_dataset(dataset):
    io.savemat(dataset_path, mdict=dataset)

def save_ppi(network):
    io.savemat(modified_ppi_path, mdict=network)

def load_dataset():
    return io.loadmat(dataset_path)
    
if __name__ == '__main__':
    dataset, network = get_dataset()
    save_dataset(dataset)
    save_ppi(network)