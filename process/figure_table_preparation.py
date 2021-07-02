from data import dataset
import pandas as pd


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
    disease_drug_dict = dataset.load_disease_drug_dict()
    n = 0
    drug = []
    for k,v in disease_drug_dict.itemvs():
        drug += list(v.keys())
        n +=  len(v)

    drug_unique = set(drug)

    print('print disease is {}, drug is {}. pairs is {}'.format(len(disease_drug_dict),
                                                               len(drug_unique),
                                                               n))
    return n

def add_drug_name():
    drug_annotation = pd.read_csv('data/Drug_annotation.csv')
    drug_name_dict = dict(zip(drug_annotation['standard_inchi_key'],
                              drug_annotation['drug']))
    data = pd.read_csv('data/result_pd_norm_merged_target_min.csv', index_col=0)
    data['drug name'] = data['drug'].apply(lambda x:drug_name_dict[x])
    data.to_csv('data/result_pd_norm_merged_target_min.csv')

