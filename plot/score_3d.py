import plotly.express as px

def plot_3_D_score_box(score_pd):

    fig = px.scatter_3d(score_pd, x='target_overlaped_rate_with_drug_norm',
                        y='disease_distance_norm', z='tc_norm',
                  color='drug',
                        hover_data=['drug_name','drug_iupac_name',"simi_name"],
                        size = 'total_score',
                         opacity=0.6)
