import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib_venn
import pickle
import seaborn as sns
import itertools
import collections
from collections import Counter
import plotly.express as px
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
from plotly.offline.offline import plot_mpl
from data import dataset


def prepare_drug_simi_tareget_piars(drugs):
    comp_target_pairs = []
    for i, info in drugs.items():
        if 'simi_mol' in info and 'simi_targeted' in info['simi_mol']:
            for i_2, t in info['simi_mol']['simi_targeted'].items():
                targets = t['target_info']['all_target']
                for t_2 in targets:
                    comp_target_pairs.append(['simi', i_2, t_2])
        if 'target_info' in info and info['target_info'] != None:
            targets_2 = info['target_info']['all_target']
            for t_3 in targets_2:
                comp_target_pairs.append(['drug', i, t_3])
    comp_target_pairs = pd.DataFrame(comp_target_pairs, columns=['source', 'inchikey', 'uniprot'])
    return comp_target_pairs


def plot_distri_tc(score_pd):
    plt.figure(figsize=(8, 8))
    sns.set_style("whitegrid")
    sns.set_context(font_scale=1.4)
    no_simi_tc = list(score_pd['tc'])
    sns.distplot(no_simi_tc, rug=True, rug_kws={"color": "g"},
                 kde_kws={"color": "k", "lw": 3, "label": "KDE"},
                 hist_kws={"histtype": "step", "linewidth": 3,
                           "alpha": 1, "color": "g"})
    plt.xlabel('Tanimoto coefficient between drugs and similar compounds')
    plt.tight_layout()


# 2. plot to similar compounds
def plot_top_inchikey(drugs, path):
    # drugs = pickle.load(open('../result/temp/drug_covid_4', 'rb'))
    simi_inchikeys = [list(info['simi_mol']['simi_targeted'].keys())
                      for i, info in drugs.items() if 'simi_mol' in info and 'simi_targeted' in info['simi_mol']]
    simi_inchikeys = list(itertools.chain.from_iterable(simi_inchikeys))
    simi_count = collections.Counter(simi_inchikeys)

    # The top frequency similar compounds
    plt.figure(figsize=(8, 8))
    sns.set_style('whitegrid')
    sns.set_context(font_scale=1.4)
    simi_top = dict(simi_count.most_common(20))
    sns.barplot(y=list(simi_top.keys()),
                x=list((simi_top.values())))
    plt.xlabel('The number of distruted drugs')
    plt.title('The most 20 frequent similar compounds among clinical drugs')
    plt.tight_layout()


def plot_top_target(drugs):
    comp_target_pairs = prepare_drug_simi_tareget_piars(drugs)

    print('The number of unique inchikeys is {}'.format(len(comp_target_pairs['inchikey'].unique())))
    print('The number of unique uniprot is {}'.format(len(comp_target_pairs['uniprot'].unique())))

    comp_target_pairs_drug = comp_target_pairs[comp_target_pairs['source'] == 'drug']
    drug_target_count = Counter(list(comp_target_pairs_drug['uniprot']))

    # the top frequency targets among drugs
    plt.figure(figsize=(8, 8))
    sns.set_style('whitegrid')
    sns.set_context(font_scale=1.4)
    drug_target_top = dict(drug_target_count.most_common(20))
    sns.barplot(y=list(drug_target_top.keys()),
                x=list((drug_target_top.values())))
    plt.xlabel('The number of distrutibed drugs')
    plt.tight_layout()
    plt.title('The most 20 frequent targets among clinical drugs')
    # plt.show()




