from matplotlib.font_manager import FontProperties
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
from data.dataset import load_databse_inhcikey



def plot_database_bar():
    database_compound_number = pd.read_csv('data/figure_database_compound_number.csv', index_col=0)

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
    #plt.style.use('classic')

    splot = sns.barplot(database_compound_number['database'],
                        database_compound_number['# of Compound'], ax=ax_1)

    for p in splot.patches:
        splot.annotate(format(p.get_height(), '.0f'),
                       (p.get_x() + p.get_width() / 2,
                        p.get_height()),
                       ha='center', va='center',
                       xytext=(0, 10),
                       textcoords='offset points')
    ax_1.set_xlabel("",  fontsize=14, fontname="Arial")
    ax_1.set_ylabel("# of compounds ",  fontsize=14, fontname="Arial")
    ax_1.spines['right'].set_visible(False)
    ax_1.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax_1.yaxis.set_ticks_position('left')
    ax_1.xaxis.set_ticks_position('bottom')
    ax_1.tick_params(axis='x', labelsize=14, rotation=45)
    ax_1.tick_params(axis='y', labelsize=14)
    plt.tight_layout()
    return fig_1

def main():
    fig = plot_database_bar()
    plt.savefig()

