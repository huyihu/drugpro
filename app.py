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

# Title the app
st.title('DrugRepo')

# set the final result table
data = pd.read_csv('data/result_pd_norm_merged_target_min.csv', index_col = 0)

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
    st.subheader(' Compound library')

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
Sorted_Disease = data['Disease_name'].unique()

# add the box for disease selection
st.markdown("### **Select Disease:**")
select_disease = []

select_disease.append(st.selectbox('', Sorted_Disease))

# accept the result data of selected disease
data_df = data[data['Disease_name'].isin(select_disease)]

selected_cols = ['Disease_id', 'Disease_name',
                 'drug', 'drug name',
                 'similar', 'tc',
                 'target_overlaped_rate_with_drug',
                 'disease_distance', 'total_score']

data_table = data_df[selected_cols]
data_table = data_table.rename(columns={'target_overlaped_rate_with_drug':'OPTS'})


# show result table of selected disease
st.dataframe(data_table)

# plot 3 D
st.write('')
row_space, row = st.beta_columns((.1,1))
with row:
    fig_3 = px.scatter_3d(data_table, x='tc',
                            y='OPTS',
                            z='disease_distance',
                            color='drug',
                            hover_data=['drug name', 'similar', 'total_score'],
                            size = 'total_score',
                             opacity=0.6,
                          title = 'The drug repurposing score',
                          width  = 80,
                          height = 80
    )

    st.plotly_chart(fig_3)