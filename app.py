import streamlit as st
from PIL import Image
from plot import whole_statistic
from matplotlib.font_manager import FontProperties
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
from data.dataset import load_databse_inhcikey
import plotly.express as px
from process.figure_table_preparation import prepare_drug_tareget_piars
import collections
from collections import Counter
from data import dataset


# Title the app
st.title('DrugRepo')

# set the final result table
@st.cache
def get_data():
    path = 'data/result_pd.csv'
    data = pd.read_csv(path, index_col = 0)
    return pd.read_csv(path)

@st.cache
def read_compound_target():
    drug_target = dataset.load_drug_target()
    return drug_target

# compound_target
drug_target_dict = read_compound_target()

#the whole score result pd
data = get_data()

disease_name_id_dict = dict(zip(data['Disease Name'], data['Disease ID']))

# add description for our website
st.markdown("""
 DrugRepo is computational pipeline to repurpose drugs for new indications. 
 The repurposing pipeline has various steps including: Compound-target data analysis, structural analysis, gene-disease 
 relationships and pathway analysis. Pipeline is able to repurpose ~0.8. million compounds across 674 diseases 
 (including various cancers, cardiovascular and kidney diseases)
""")

# show some conceptual figure for work flow
workflow_image = Image.open("images/Workflow.JPG")

st.markdown("### Whole workflow")
st.image(workflow_image)

# show absic data

col1, col2, col3 = st.beta_columns(3)

with col1:
    st.markdown(f"**Compound:** {788078}")
    st.markdown(f"**Disease:** {674}")


with col2:
    st.markdown(f"**Target:** " + \
                f"{9018}")
    st.markdown(f"**Drug:** {1092}")


with col3:
    st.markdown(f"**Com-Tar pairs:** " + \
                f"{1489537}")
    st.markdown(f"**Disease-Drug:** {3775}")


st.text("")


# plot the final statistic result, the number of compound in each databases(left) with taregt (right)
st.write('')
row_space, row, row_space_2, row_2 = st.beta_columns((.1,1, .1, 1))
with row:
    database_compound_number = pd.read_csv('data/figure_database_compound_number.csv', index_col=0)
    st.subheader('Compound library')

    # set figure parameters
    fig_1, ax_1 = plt.subplots(figsize=(6, 4))
    sns.set(font='Arial', style="white", context="paper")
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_name('arial')
    matplotlib.rc('xtick', labelsize=20)
    matplotlib.rc('ytick', labelsize=20)
    font = {
        'weight': 'bold',
        'size': 14}

    matplotlib.rc('font', **font)
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = "sans-serif"
    plt.rcParams['figure.facecolor'] = 'white'
    # plt.style.use('classic')

    splot = sns.barplot(database_compound_number['database'],
                        database_compound_number['# of Compound'], ax=ax_1)

    # add annotation number on figure
    for p in splot.patches:
        splot.annotate(format(p.get_height(), '.0f'),
                       (p.get_x() + p.get_width() / 2,
                        p.get_height()),
                       ha='center', va='center',
                       xytext=(0, 10),
                       textcoords='offset points')
    ax_1.set_xlabel("", fontsize=14, fontname="Arial")
    ax_1.set_ylabel("# of compounds ", fontsize=14, fontname="Arial")
    ax_1.spines['right'].set_visible(False)
    ax_1.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax_1.yaxis.set_ticks_position('left')
    ax_1.xaxis.set_ticks_position('bottom')
    ax_1.tick_params(axis='x', labelsize=14, rotation=45)
    ax_1.tick_params(axis='y', labelsize=14)
    plt.tight_layout()

    st.pyplot(fig_1)


with row_2:
    database_compound_number = pd.read_csv('data/figure_database_target_number.csv', index_col=0)
    st.subheader('Target library')

    # set figure parameters
    fig_2, ax_2 = plt.subplots(figsize=(6, 4))
    sns.set(font='Arial', style="white", context="paper")
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_name('arial')
    matplotlib.rc('xtick', labelsize=20)
    matplotlib.rc('ytick', labelsize=20)
    font = {
        'weight': 'bold',
        'size': 14}

    matplotlib.rc('font', **font)
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = "sans-serif"
    plt.rcParams['figure.facecolor'] = 'white'
    # plt.style.use('classic')

    splot = sns.barplot(database_compound_number['database'],
                        database_compound_number['# of target'], ax=ax_2)

    # add annotation number on figure
    for p in splot.patches:
        splot.annotate(format(p.get_height(), '.0f'),
                       (p.get_x() + p.get_width() / 2,
                        p.get_height()),
                       ha='center', va='center',
                       xytext=(0, 10),
                       textcoords='offset points')
    ax_2.set_xlabel("", fontsize=14, fontname="Arial")
    ax_2.set_ylabel("# of Target ", fontsize=14, fontname="Arial")
    ax_2.spines['right'].set_visible(False)
    ax_2.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax_2.yaxis.set_ticks_position('left')
    ax_2.xaxis.set_ticks_position('bottom')
    ax_2.tick_params(axis='x', labelsize=14, rotation=45)
    ax_2.tick_params(axis='y', labelsize=14)
    plt.tight_layout()

    st.pyplot(fig_2)


# decide the unique selective disease
Sorted_Disease = data['Disease Name'].unique()

# add the box for disease selection
st.markdown("### **Select Disease:**")

select_disease = st.selectbox('', Sorted_Disease)

select_disease_id = disease_name_id_dict.get(select_disease)

# first read feature vector for sleected disease for next TSNE PLOT
feature_vector = pd.read_csv('data/result_stored/{}/feature_vector.csv'.format(select_disease_id))

@st.cache
def get_seperate_disease_data(select_disease):
    # accept the result data of selected disease
    data_table = data[data['Disease Name'].isin([select_disease])]
    return data_table


# show result table of selected disease
data_table = get_seperate_disease_data(select_disease)

# filter by tc, set default is 0.5 1
min_tc, max_tc = st.slider('Select a range of TC', min_value=0.2, max_value=1.0, step=0.1, value=(0.5, 0.9))

# filter by total score, set default is 0.5-1
min_score, max_score = st.slider(f'Select a range of score', min_value=0.0, max_value=1.0, step=0.1, value=(0.5, 1.0))

#based on the tc and total score select data
data_show = data_table[(data_table['TC']>= min_tc) &
                       (data_table['TC'] <= max_tc) &
                       (data_table['total_score']>= min_score)
                       & (data_table['total_score'] <= max_score)]

# rank the show result by order from largest total score to lowest score
data_show = data_show.sort_values(by='total_score', ascending=False)

# use the left drug to get all the drug relatd atrget and their frequency
common_target_count = prepare_drug_tareget_piars(list(data_show['Drug Inchikey'].unique()), drug_target_dict)


# show then selected result table on app
st.dataframe(data_show)

# plot most top target for drugs on app, the top compound on the left, left:right = 1:2
st.write('')
row_space_3, row_3, row_space_4, row_4 = st.beta_columns((.1,1, .1, 1))
with row_3:
    st.subheader('Top targets among drugs')
    # known the most frequent target among the drugs
    largest_count_target = common_target_count.most_common(1)[0][1]
    default_target_n_show = min(largest_count_target, 5)
    # set a bar selection for how many top target to show,
    # with max count is the the most frequent target, default value is also the most frequent target
    top_target_n = st.slider(f'How many top drug target?', min_value=1, max_value=largest_count_target, step=1,
                                     value=default_target_n_show)

    # after receive the top n target that need to be shown
    drug_target_top_dict = dict(common_target_count.most_common(top_target_n))

    # plot the top target selected
    fig_3, ax_3 = plt.subplots(figsize=(6, 4))
    sns.set(font='Arial', style="white", context="paper")
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_name('arial')
    matplotlib.rc('xtick', labelsize=20)
    matplotlib.rc('ytick', labelsize=20)
    font = {
        'weight': 'bold',
        'size': 14}

    matplotlib.rc('font', **font)
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = "sans-serif"
    plt.rcParams['figure.facecolor'] = 'white'
    sns.barplot(y=list(drug_target_top_dict.keys()),
                x=list((drug_target_top_dict.values())), ax = ax_3)

    ax_3.set_title('The most 20 frequent targets among clinical drugs' , fontsize=14, fontname="Arial")
    ax_3.set_xlabel('# of drugs', fontsize=14, fontname="Arial")

    # Only show ticks on the left and bottom spines
    ax_3.spines['right'].set_visible(False)
    ax_3.spines['top'].set_visible(False)
    ax_3.yaxis.set_ticks_position('left')
    ax_3.xaxis.set_ticks_position('bottom')
    ax_3.tick_params(axis='x', labelsize=14)
    ax_3.tick_params(axis='y', labelsize=14)
    plt.tight_layout()

    st.pyplot(fig_3)


# plot top similar compounds by in total score in app
with row_4:

    st.subheader('Top repurposed compounds')

    # knwon the number of al selelcted similar compound
    largest_count_compound = data_show.shape[0]

    # as there might be many similar compound, we set largest is 50
    largest_count_compound = min(largest_count_compound, 50)
    default_compound_n_show = min(largest_count_compound, 5)

    # set a bar selection for how many top compounds to show,
    # with max count is the number of compound, default value is maximum or  20
    top_n_compound = st.slider(f'How many top compounds?', min_value=1, max_value=largest_count_compound, step=1,
                             value=default_compound_n_show)

    # after receive the top n target that need to be shown
    compound_s = data_show['Compound'][0:top_n_compound]
    total_score_s = data_show['total_score'][0:top_n_compound]

    # plot the top compound selected
    fig_4, ax_4 = plt.subplots(figsize=(6, 4))
    sns.set(font='Arial', style="white", context="paper")
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_name('arial')
    matplotlib.rc('xtick', labelsize=20)
    matplotlib.rc('ytick', labelsize=20)
    font = {
        'weight': 'bold',
        'size': 14}

    matplotlib.rc('font', **font)
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = "sans-serif"
    plt.rcParams['figure.facecolor'] = 'white'
    sns.barplot(y=compound_s,
                x=total_score_s, ax=ax_4)

    ax_4.set_title('The most 20 frequent targets among clinical drugs', fontsize=14, fontname="Arial")
    ax_4.set_xlabel('Repurposing score', fontsize=14, fontname="Arial")

    # Only show ticks on the left and bottom spines
    ax_4.spines['right'].set_visible(False)
    ax_4.spines['top'].set_visible(False)
    ax_4.yaxis.set_ticks_position('left')
    ax_4.xaxis.set_ticks_position('bottom')
    ax_4.tick_params(axis='x', labelsize=14)
    ax_4.tick_params(axis='y', labelsize=14)
    plt.tight_layout()

    st.pyplot(fig_4)

# plot tc density , score diatnce

st.write('')
row_space_5, row_5, row_space_6, row_6 = st.beta_columns((.1,1, .1, 1))
with row_5:
    st.subheader('TC density')
    # plot the top target selected
    fig_5, ax_5 = plt.subplots(figsize=(6, 4))
    sns.set(font='Arial', style="white", context="paper")
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_name('arial')
    matplotlib.rc('xtick', labelsize=20)
    matplotlib.rc('ytick', labelsize=20)
    font = {
        'weight': 'bold',
        'size': 14}

    matplotlib.rc('font', **font)
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = "sans-serif"
    plt.rcParams['figure.facecolor'] = 'white'

    no_simi_tc = list(data_show['TC'])
    sns.distplot(no_simi_tc, rug=True, rug_kws={"color": "g"},
                 kde_kws={"color": "k", "lw": 3, "label": "KDE"},
                 hist_kws={"histtype": "step", "linewidth": 3,
                           "alpha": 1, "color": "g"})
    ax_5.set_xlabel('TC', fontsize=14, fontname="Arial")

    ax_5.set_title('Tanimoto coefficient between drugs and similar compounds' , fontsize=14, fontname="Arial")

    # Only show ticks on the left and bottom spines
    ax_5.spines['right'].set_visible(False)
    ax_5.spines['top'].set_visible(False)
    ax_5.yaxis.set_ticks_position('left')
    ax_5.xaxis.set_ticks_position('bottom')
    ax_5.tick_params(axis='x', labelsize=14)
    ax_5.tick_params(axis='y', labelsize=14)
    plt.tight_layout()

    st.pyplot(fig_5)

with row_6:
    st.subheader('Score density')
    # plot the top target selected
    fig_6, ax_6 = plt.subplots(figsize=(6, 4))
    sns.set(font='Arial', style="white", context="paper")
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_name('arial')
    matplotlib.rc('xtick', labelsize=20)
    matplotlib.rc('ytick', labelsize=20)
    font = {
        'weight': 'bold',
        'size': 14}

    matplotlib.rc('font', **font)
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = "sans-serif"
    plt.rcParams['figure.facecolor'] = 'white'

    scores_total = list(data_show['total_score'])
    sns.distplot(scores_total, rug=True, rug_kws={"color": "g"},
                 kde_kws={"color": "k", "lw": 3, "label": "KDE"},
                 hist_kws={"histtype": "step", "linewidth": 3,
                           "alpha": 1, "color": "g"})
    ax_6.set_xlabel('Repurposing score', fontsize=14, fontname="Arial")

    ax_6.set_title('The distribution of total repurposing score', fontsize=14, fontname="Arial")

    # Only show ticks on the left and bottom spines
    ax_6.spines['right'].set_visible(False)
    ax_6.spines['top'].set_visible(False)
    ax_6.yaxis.set_ticks_position('left')
    ax_6.xaxis.set_ticks_position('bottom')
    ax_6.tick_params(axis='x', labelsize=14)
    ax_6.tick_params(axis='y', labelsize=14)
    plt.tight_layout()

    st.pyplot(fig_6)

# plot 3 D score distribution
st.write('')

st.subheader('3D box for drug repurposing score')
fig_7 = px.scatter_3d(data_show,
                      x='TC',
                        y='OPTS',
                        z='Disease Distance',
                        color='Drug Name',
                        hover_data=['Drug Name', 'Compound', 'total_score'],
                        size = 'total_score',
                         opacity=0.8

)

st.plotly_chart(fig_7)


# plot the feature similarity tsne  between drugs and compounds

inchikey_s = list(data_show['Drug Inchikey'].unique()) + list(data_show['Compound'].unique())

feature_vector_selected = feature_vector[feature_vector['inchikey'].isin(inchikey_s)]
feature_vector_selected['size'] = 1.0

st.write('')

st.subheader('The TSNE plot between drug-compound')

plt.figure(figsize=(6, 6))

figure_8 = px.scatter(feature_vector_selected,
                        x='X',
                        y='Y',
                        color='source',
                        size='size',
                        hover_data=['inchikey'],
                        size_max=20,
                         opacity=0.8)

figure_8.update_layout({'plot_bgcolor': 'rgb(255, 255, 255)',
                  'paper_bgcolor': 'rgb(255, 255, 255)'})
st.plotly_chart(figure_8)

# p1 = sns.scatterplot(x="X", y="Y",
#                 hue='source',
#                 legend='full',
#                 data=feature_vector_selected,
#                      alpha=0.2,
#                      sizes=(500, 300))

#add annotations one by one with a loop

# for line in range(0, feature_vector_selected.shape[0]):
#     p1.text(feature_vector_selected.X[line] + 0.2,
#             feature_vector_selected.Y[line],
#             feature_vector_selected.inchikey[line],
#             horizontalalignment='left',
#             size='medium',
#             color='black',
#             weight='semibold')

# p1.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)',
#                   'paper_bgcolor': 'rgba(0, 0, 0, 0)', })
# st.plotly_chart(p1)




