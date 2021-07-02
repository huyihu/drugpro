#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 26.5.2020 10.42
# @Author  : YINYIN
# @Site    :
# @File    : similarity_analysis.py
# @Software: PyCharm

import dataset
import matplotlib.pyplot as plt
import seaborn as sns
from math import isnan
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.preprocessing import MultiLabelBinarizer
import plotly.express as px
import pickle



def split_use(x):
    return x.split(';')


def pre_targets(data):
    data['targets'] = data['all_target'].apply(lambda x: x.split(';') if isinstance(x, str) else None)
    data = data[data['targets'].notna()]
    return data


def pre_tsne():
   simi_pd = dataset.drug_simi_targeted()
   drug_pd = dataset.trial_drug_with_target()
   drug = pre_targets(drug_pd)
   simi = pre_targets(simi_pd)
   name = list(drug['name']) + list(simi['name_cid'])
   target_list = list(drug['targets']) + list(simi['targets'])
   len_drug = drug.shape[0]
   len_simi = simi.shape[0]
   label = ['drug']*len_drug + ['simi']*len_simi
   mlb = MultiLabelBinarizer()
   df_data = pd.DataFrame(mlb.fit_transform(target_list), columns=mlb.classes_, index=name)
   return df_data, label, name


def drug_vector(drug_fingerprint, drug_target, disease_target, disease_id,  inchikey_s):
    # generate target binary matrix
    target_list_list = []
    for inchikey in inchikey_s:
        targets = drug_target.get(inchikey)['all_target']
        targets = [t for t in targets if not isinstance(t, float) ]
        target_list_list.append(targets)
    mlb = MultiLabelBinarizer()
    compound_target_pd = pd.DataFrame(mlb.fit_transform(target_list_list),
                                      columns=mlb.classes_,
                                      index=inchikey_s)

    # generate disease entry overlap matrix
    disease_entrze = disease_target['disease_simple_entrze'][disease_id]
    disease_target_binary_list = []

    for inchikey in inchikey_s:
        target_entrze = drug_target.get(inchikey).get('all_target_entrze')
        if target_entrze == None:
            target_entrze = []
        disease_target_binary = [1 if 'T' + str(dz) in target_entrze else 0 for dz in disease_entrze ]
        disease_target_binary_list.append(disease_target_binary)
    disease_target_binary_pd = pd.DataFrame(disease_target_binary_list, columns=disease_entrze, index=inchikey_s)

    # generate fingerprint  binary
    fingerprint_binary_list = []
    for inchikey in inchikey_s:
        fingerprint_binary = list(drug_fingerprint.get(inchikey).ToBitString())
        fingerprint_binary_list.append(fingerprint_binary)
    fingerprint_binary_pd = pd.DataFrame(fingerprint_binary_list,
                                         columns=['EXTENDED_' + str(i) for i in range(1,1025)],
                                         index=inchikey_s)

    # MERGE ALL THE FINGERPRINT
    feature_vector = pd.concat([disease_target_binary_pd, compound_target_pd, fingerprint_binary_pd],
              axis=1,
              ignore_index=False)

    feature_vector = feature_vector.loc[:,~(feature_vector==0).all(axis=0)]

    return feature_vector


def tsne_plot(data_X,label, com_name):
    tsne = TSNE(n_components=2, random_state=0)
    tsne_obj = tsne.fit_transform(data_X)
    tsne_df = pd.DataFrame({'X': tsne_obj[:, 0],
                            'Y': tsne_obj[:, 1],
                            'source': label, 'com_name':com_name })
    plt.figure(figsize=(15, 15))
    p1 = sns.scatterplot(x="X", y="Y",
                    hue='source',
                    legend='full',
                    data=tsne_df, alpha=0.2, sizes=(500, 300))

    #add annotations one by one with a loop

    for line in range(0, tsne_df.shape[0]):
        p1.text(tsne_df.X[line] + 0.2,
                tsne_df.Y[line],
                tsne_df.com_name[line],
                horizontalalignment='left',
                size='medium', color='black',
                weight='semibold')
    p1.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)',
                      'paper_bgcolor': 'rgba(0, 0, 0, 0)', })
    plt.show()




