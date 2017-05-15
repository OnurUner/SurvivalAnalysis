# -*- coding: utf-8 -*-
import lifelines as ll
import numpy as np
import pandas as pd
from numpy import genfromtxt
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 10, 5
from lifelines.estimation import KaplanMeierFitter
from load_data import load_dataset, convert_entrez2uni
from lifelines.utils import k_fold_cross_validation
from lifelines import CoxPHFitter
import random
from sklearn import preprocessing

dataset = load_dataset()
columns = dataset["x_columns"]
features_path = '..\\data\\selected_features.csv'
features = genfromtxt(features_path, delimiter=',')[1:].astype(int)
#features = random.sample(range(0, 15939), 23)
columns = columns[features]
for i, column in enumerate(columns):
    columns[i] = column.strip()

X = dataset["data"]
X = X[:,features]
data = pd.DataFrame(X, columns = columns)
event = dataset["cencoring"]
time = dataset["labels"][:,5]
data["event"] = event.ravel()
data["time"] = time
    
cf = CoxPHFitter()
scores = k_fold_cross_validation(cf, data, 'time', event_col='event', k=3)
print scores
print np.mean(scores)
print np.std(scores)


le = preprocessing.LabelEncoder()
subtypes = le.fit_transform(dataset["subtypes"])
data["subtype"] = subtypes 
T = data["time"]
C = data["event"]

kmf = KaplanMeierFitter()
kmf.fit(T, event_observed=C)
kmf.plot(title = 'Survival Day Profile of Breast Cancer Patients')

# Basal
f1 = data.subtype == 0
T1 = data[f1]['time']
C1 = data[f1]['event']

# Her2
f2 = data.subtype == 1
T2 = data[f2]['time']
C2 = data[f2]['event']

# LumA
f3 = data.subtype == 2
T3 = data[f3]['time']
C3 = data[f3]['event']

# LumB
f4 = data.subtype == 3
T4 = data[f4]['time']
C4 = data[f4]['event']

# Normal
f5 = data.subtype == 4
T5 = data[f5]['time']
C5 = data[f5]['event']

plt.figure()
ax = plt.subplot(111)
kmf.fit(T1, event_observed=C1, label=['Basal'])
kmf.survival_function_.plot(ax=ax)

kmf.fit(T2, event_observed=C2, label=['Her2'])
kmf.survival_function_.plot(ax=ax)

kmf.fit(T3, event_observed=C3, label=['LumA'])
kmf.survival_function_.plot(ax=ax)

kmf.fit(T4, event_observed=C4, label=['LumB'])
kmf.survival_function_.plot(ax=ax)

kmf.fit(T5, event_observed=C5, label=['Normal'])
kmf.survival_function_.plot(ax=ax)

plt.title('Lifespans of different subtypes')

kmf2 = plt.gcf()

uniprots = []

for entrez in columns:
    uniprots.append(convert_entrez2uni(entrez))


#rossi_dataset = load_rossi()
#cf = CoxPHFitter()
#cf.fit(rossi_dataset, 'week', event_col='arrest')
#
#cf.print_summary()  # access the results using cf.summary
#predict = cf.predict_survival_function(rossi_dataset[:10]).plot()