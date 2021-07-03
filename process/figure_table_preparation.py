from data import dataset
import pandas as pd
from data import dataset
import streamlit as st
import collections
from collections import Counter
import pickle


def database_compound_number():
    inchikey_targets_databases = dataset.load_databse_inhcikey()
    changed_name = {'chembl_target': 'ChEMBL',
                    'binding_target': 'BindingDB',
                    'drugbank_target': 'DrugBank',
                    'GtopDB_target': 'GtopDB',
                    'DgiDB_target': 'DGiDB'}
    inchikey_targets_databases = {changed_name[k]: v for k, v in inchikey_targets_databases.items()}
    x = list(inchikey_targets_databases.keys())
    y = [len(v) for k, v in inchikey_targets_databases.items()]
    figure_1_pd = pd.DataFrame({'database':x,
                  '# of Compound':y})

    figure_1_pd.to_csv('data/figure_database_compound_number.csv')


def databases_target_number():
    inchikey_targets_databases = dataset.load_databse_inhcikey()
    changed_name = {'chembl_target': 'ChEMBL',
                    'binding_target': 'BindingDB',
                    'drugbank_target': 'DrugBank',
                    'GtopDB_target': 'GtopDB',
                    'DgiDB_target': 'DGiDB'}
    inchikey_targets_databases = {changed_name[k]: v for k, v in inchikey_targets_databases.items()}
    target_set_list = {}
    target_len_list = {}
    target = []
    for k, v in inchikey_targets_databases.items():
        target_set = []
        for v2 in list(v.values()):
            target_set += v2
            target += v2
        target_set = list(set(target_set))
        target_set_list[k] = target_set
        target_len_list[k] = len(set(target_set))

    # the number of target
    print('the compound_target is {}'.format(len(target)))
    print('The number of target is {}'.format(len(set(target))))



    # generate to x, y

    figure_2_pd = pd.DataFrame({'database': list(target_len_list.keys()),
                                '# of target': list(target_len_list.values())})

    figure_2_pd.to_csv('data/figure_database_target_number.csv')


# information about drug disease
def disease_drug():
    compound_target_dict = dataset.load_inchikey_target_dict()
    disease_drug_dict = dataset.load_disease_drug_dict()
    n = 0
    drug = []

    for k,v in disease_drug_dict.items():
        drug += list(v.keys())
        n +=  len(v)

    drug_unique = set(drug)

    # prepare a dictionary contain drug target
    drug_target_dict_simple = {}
    for d in drug_unique:
        if compound_target_dict.get(d) != None and compound_target_dict[d].get('all_target') !=None:
            drug_target_dict_simple [d] = compound_target_dict[d].get('all_target')
    with open('data/resource_dict/drug_target_simple', 'wb') as handle:
        pickle.dump(drug_target_dict_simple, handle, protocol=pickle.HIGHEST_PROTOCOL)



    print('print disease is {}, drug is {}. pairs is {}'.format(len(disease_drug_dict),
                                                               len(drug_unique),
                                                               n))
    return n

def result_table_standrad():
    drug_annotation = pd.read_csv('data/Drug_annotation.csv')
    drug_name_dict = dict(zip(drug_annotation['standard_inchi_key'],
                              drug_annotation['drug']))
    data = pd.read_csv('data/result_pd_norm_merged_target_min.csv', index_col=0)
    data['drug name'] = data['drug'].apply(lambda x:drug_name_dict[x])

    #
    selected_cols = ['Disease_id', 'Disease_name',
                     'drug', 'drug name',
                     'similar', 'tc',
                     'target_overlaped_rate_with_drug',
                     'disease_distance_norm', 'total_score']

    data_table = data[selected_cols]
    data_table = data_table.rename(columns={
        'Disease_id': 'Disease ID',
        'Disease_name': 'Disease Name',
        'drug': 'Drug Inchikey',
        'drug name': 'Drug Name',
        'similar': 'Compound',
        'tc': 'TC',
        'target_overlaped_rate_with_drug': 'OPTS',
        'disease_distance_norm': 'Disease Distance'})
    data_table.to_csv('data/result_pd.csv')


# prepare drug target set for show the to targets
def prepare_drug_tareget_piars(drugs, drug_target_dict):
    drug_target_pairs = []
    for d in drugs:
        targets = drug_target_dict.get(d)
        if targets != None:
            drug_target_pairs += targets
    common_target_count = Counter(drug_target_pairs)
    return common_target_count


