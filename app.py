import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



# Title the app
st.title('DrugRepo')

data = pd.read_csv('result_pd_norm_merged_target_min.csv', index_col = 0)
st.markdown("""
 DrugRepo is computational pipeline to repurpose drugs for new indications. 
 The repurposing pipeline has various steps including: Compound-target data analysis, structural analysis, gene-disease 
 relationships and pathway analysis. Pipeline is able to repurpose ~0.8. million compounds across 674 diseases 
 (including various cancers, cardiovascular and kidney diseases)
""")
