import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import numpy as np
import scanpy as sc
import seaborn as sns
from statannot import add_stat_annotation


dict_WT_KO_colors = {'KO1': '#ea693c', 'KO2': '#f28d5c', 'WT1': '#20668d', 'WT2': '#229eb2', 'KO': '#ee7b4c', 'WT': '#2182a0'}


def adata_plot_KOvsWT(adata, list_names, do_return=False, col_cell_type='merged_cell_type'):
    fig, ax = plt.subplots(2, 1, figsize=(9, 12))
    
    df_proportions_KO_WT = pd.DataFrame(columns=['KO1', 'KO2', 'WT1', 'WT2'], index=list_names)
    df_counts_KO_WT = pd.DataFrame(columns=['KO1', 'KO2', 'WT1', 'WT2'], index=list_names)
    df_pval = pd.DataFrame(columns=['p-val'], index=list_names)
    
    counts_KO_all = len(adata[adata.obs['condition'] == 'KO'])
    counts_WT_all = len(adata[adata.obs['condition'] == 'WT'])
    
    for cell_type in list_names:
        adata_sub = adata[adata.obs[col_cell_type] == cell_type]
        counts_KO = len(adata_sub[adata_sub.obs['condition'] == 'KO'])
        counts_WT = len(adata_sub[adata_sub.obs['condition'] == 'WT'])
        
        counts_KO1 = len(adata_sub[adata_sub.obs['batch'] == 'KO1'])
        counts_KO2 = len(adata_sub[adata_sub.obs['batch'] == 'KO2'])
        counts_WT1 = len(adata_sub[adata_sub.obs['batch'] == 'WT1'])
        counts_WT2 = len(adata_sub[adata_sub.obs['batch'] == 'WT2'])
        
        
        df_counts_KO_WT.loc[cell_type] = [counts_KO1, counts_KO2, counts_WT1, counts_WT2]
        df_proportions_KO_WT.loc[cell_type] = [counts_KO1/(counts_KO + counts_WT), counts_KO2/(counts_KO + counts_WT), 
                                               counts_WT1/(counts_KO + counts_WT), counts_WT2/(counts_KO + counts_WT)]
        df_pval.loc[cell_type] = chi2_contingency(np.array([[counts_KO, counts_WT], 
                                                            [(counts_KO + counts_WT) * counts_KO_all/len(adata), (counts_KO + counts_WT) * counts_WT_all/len(adata) ]]))[1]
        


    idx_sort = ((df_proportions_KO_WT['KO1'] + df_proportions_KO_WT['KO2']) - 
                (df_proportions_KO_WT['WT1'] + df_proportions_KO_WT['WT2'])).sort_values(ascending=False).index


    df_proportions_KO_WT.loc[idx_sort].plot(kind='bar', stacked=True, color=list(dict_WT_KO_colors.values()), 
                                            ax=ax[0])
    df_counts_KO_WT.loc[idx_sort].plot(kind='bar', stacked=True, color=list(dict_WT_KO_colors.values()), 
                                       ax=ax[1])
    
    
    for idx, pval in enumerate(df_pval.loc[idx_sort]['p-val'].values):
        if 0.01 < pval < 0.05:
            pval_txt = '*'
        elif 0.001 < pval < 0.01:
            pval_txt = '**'
        elif 0.0001 < pval < 0.001:
            pval_txt = '***'
        elif pval < 0.0001:
            pval_txt = '****'
        else:
            pval_txt = ''
    
        ax[0].text(idx, 1.03, pval_txt, ha='center')
        
        max_val = df_counts_KO_WT.loc[idx_sort].sum(1).max()
        
        ax[1].text(idx, df_counts_KO_WT.loc[idx_sort].sum(1).iloc[idx] + 0.03 * max_val, pval_txt, ha='center')
    
    ax[0].set_ylim([0, 1.1])
    ax[1].set_ylim([0, 1.1 * max_val])
    
    ax[0].grid(False)
    ax[1].grid(False)

    plt.grid(b=None)
    plt.tight_layout()
    
    if do_return:
        return df_proportions_KO_WT, df_counts_KO_WT, df_pval
    
    
    
    
def adata_plot_KOvsWT_2bars(adata, list_names, do_return=False, col_cell_type='merged_cell_type'):
    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
    
    df_counts_KO_WT = pd.DataFrame(columns=['KO1', 'KO2', 'WT1', 'WT2'], index=list_names)
    df_pval = pd.DataFrame(columns=['p-val'], index=list_names)
    
    counts_KO_all = len(adata[adata.obs['condition'] == 'KO'])
    counts_WT_all = len(adata[adata.obs['condition'] == 'WT'])
    
    for cell_type in list_names:
        adata_sub = adata[adata.obs[col_cell_type] == cell_type]
        counts_KO = len(adata_sub[adata_sub.obs['condition'] == 'KO'])
        counts_WT = len(adata_sub[adata_sub.obs['condition'] == 'WT'])
        
        counts_KO1 = len(adata_sub[adata_sub.obs['batch'] == 'KO1'])
        counts_KO2 = len(adata_sub[adata_sub.obs['batch'] == 'KO2'])
        counts_WT1 = len(adata_sub[adata_sub.obs['batch'] == 'WT1'])
        counts_WT2 = len(adata_sub[adata_sub.obs['batch'] == 'WT2'])
        
        
        df_counts_KO_WT.loc[cell_type] = [counts_KO1, counts_KO2, counts_WT1, counts_WT2]
        df_pval.loc[cell_type] = chi2_contingency(np.array([[counts_KO, counts_WT], 
                                                            [(counts_KO + counts_WT) * counts_KO_all/len(adata), (counts_KO + counts_WT) *
                                                             counts_WT_all/len(adata) ]]))[1]
        


    idx_sort = ((df_counts_KO_WT['KO1'] + df_counts_KO_WT['KO2']) + 
                (df_counts_KO_WT['WT1'] + df_counts_KO_WT['WT2'])).sort_values(ascending=False).index

    w, dx = 0.4, 0.2
    for idx, idx_name in enumerate(idx_sort):
        plt.bar(idx + dx, width=w, height=df_counts_KO_WT.loc[idx_name,'KO1'], color=dict_WT_KO_colors['KO1'])
        plt.bar(idx + dx, width=w, height=df_counts_KO_WT.loc[idx_name,'KO2'], bottom=df_counts_KO_WT.loc[idx_name,'KO1'], 
                color=dict_WT_KO_colors['KO2'])
        
        plt.bar(idx - dx, width=w, height=df_counts_KO_WT.loc[idx_name,'WT1'], color=dict_WT_KO_colors['WT1'])
        plt.bar(idx - dx, width=w, height=df_counts_KO_WT.loc[idx_name,'WT2'], bottom=df_counts_KO_WT.loc[idx_name,'WT1'], 
                color=dict_WT_KO_colors['WT2'])
        
    
    list_max_val = []
    
    for idx, pval in enumerate(df_pval.loc[idx_sort]['p-val'].values):
        if 0.01 < pval < 0.05:
            pval_txt = '*'
        elif 0.001 < pval < 0.01:
            pval_txt = '**'
        elif 0.0001 < pval < 0.001:
            pval_txt = '***'
        elif pval < 0.0001:
            pval_txt = '****'
        else:
            pval_txt = ''
            
        max_val = max([df_counts_KO_WT.loc[idx_sort[idx], 'KO1'] + df_counts_KO_WT.loc[idx_sort[idx], 'KO2'], 
                      df_counts_KO_WT.loc[idx_sort[idx], 'WT1'] + df_counts_KO_WT.loc[idx_sort[idx], 'WT2']])
        
        list_max_val.append(max_val)
        
        
        ax.text(idx,  1.03 * max_val, pval_txt, ha='center')
    
    ax.set_ylim([0, 1.1 * max(list_max_val)])
    
    ax.set_xticks(range(len(df_counts_KO_WT)))
    ax.set_xticklabels([i.split('/')[1] for i in idx_sort], rotation=90)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.grid(False)

    
    plt.tight_layout()
    
    if do_return:
        return df_counts_KO_WT, df_pval

    
    
def stat_annot_gene(gene, adata, dict_pops, type_plot='violin', add_stats=True):
    fig = plt.figure()
    df = pd.DataFrame({'x': adata.obs['subtype'].values, 'y': adata[:, gene].X.toarray().ravel(), 'hue': adata.obs['condition'].values})
    
    if type_plot == 'violin':
        g = sns.violinplot(x='x', y='y', hue='hue', data=df, rotation=90, split=True, cut=True, inner="stick", 
                           palette=[dict_WT_KO_colors['KO2'], dict_WT_KO_colors['WT2']])
    elif type_plot == 'box':
        g = sns.boxplot(x='x', y='y', hue='hue', data=df, palette=[dict_WT_KO_colors['KO2'], dict_WT_KO_colors['WT2']])

    g.set_xticklabels(g.get_xticklabels(), rotation=40, ha='right')
    g.set_xlabel(gene)
    g.set_ylabel('')

    
    if add_stats:
        for i in  sorted(dict_pops.keys()):
            try:
                add_stat_annotation(g, data=df, x='x', y='y', hue='hue',
                                box_pairs=[((i, "KO"), (i, "WT"))],
                                test='t-test_ind', loc='inside', comparisons_correction=None, verbose=False)
            except:
                pass

