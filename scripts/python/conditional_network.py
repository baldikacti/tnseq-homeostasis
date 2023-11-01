# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
# from networkx.drawing.nx_agraph import to_agraph
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout, to_agraph, write_dot, from_agraph

df_data = pd.read_csv('data/peter/peter_UC_processed.csv')
data = df_data.iloc[:,7:]
# df_data = df_data.drop(df_data[df_data.time == 0].index)

# droplist = list(data.columns[data.apply(lambda col: col.mean() < 5)])
frank_essential = list(data.columns[data.apply(lambda col: max(col) < 5)])
## remove frank essential genes ###
# df_data = df_data.drop(droplist, axis=1)


conditions = list(df_data.groupby(['strain', 'condition']).groups.keys())
# conditions = [('Carbon Source', 'M2_noCarbon'), ('Nitrogen Source', 'M2_noCarbon'), ('Nitrogen Source', 'M2_noNitrogen'),
#                 ('PH', 'PYE'), ('PYE', 'PYE'), ('Stress', 'PYE')]

# conditions = ['control', 'bile']
# con_size = np.asarray(df_data.groupby(['Group', 'Media']).size())


df_infer = pd.read_csv('inference_oak_1.csv')
df_H = pd.read_csv('H_matrix_oak_1.csv')
locus_tag = df_H['locus_tag']
df_sc = pd.DataFrame(locus_tag)
df_sc.columns=['Name']


# df_H = df_H.set_index('locus_tag')
# df_H = df_H.drop(list(droplist))
df_H = df_H.iloc[:,3:]
H_nrow , H_ncol = df_H.shape
print(df_H.shape)

# H_matrix = np.zeros((5, H_nrow), dtype='int')
# H_matrix_str = np.asarray(df_infer['H'])

H_matrix = np.asarray(df_H.T)

#%%
all_conditon_pred = []
all_keep = []
for i in range(len(conditions)):
    prediction_ls = []
    keep = []
    for j in range(H_nrow):
        pred = H_matrix[:, j]
        prob_ls = list(zip(pred, df_infer['gi_%i' %(i+1)][1:]))
        res = {1:0, 0:0}
        for ess, prob in prob_ls:
            res[ess] += float(prob)

        if res[0] > 0.95:
            keep.append(df_sc['Name'].iloc[j])
            prediction_ls.append(res[0])
    all_conditon_pred.append(prediction_ls)
    all_keep.append(list(set(keep)-set(frank_essential)))

#%%
interection_genes = list(set.intersection(*map(set, all_keep)))
conditional_genes = [list(set(all_keep[i])-set(interection_genes)) for i in range(len(all_keep))]


G = nx.Graph()
G.add_node('Interections')
G.add_nodes_from(conditions)
conditional_11 = list(set(conditional_genes[11])-set(conditional_genes[15])) 
conditional_15 = list(set(conditional_genes[15])-set(conditional_genes[11]))

G.add_nodes_from(conditional_11)
G.add_nodes_from(conditional_15)

edgelist_conditions = [(conditions[i], 'Interections') for i in range(len(conditions))]
edgelist_11 = [(conditions[11], item) for item in conditional_11]
edgelist_15 = [(conditions[15], item) for item in conditional_15]

G.add_edges_from(edgelist_conditions)
G.add_edges_from(edgelist_11)
G.add_edges_from(edgelist_15)
# nx.draw(G)

# write_dot(G, 'cen.dot')

# sizelist = [1000] + [50] * 16

# size_1 = [G.degree(node) * 80 for node in G.nodes() if node in ['Interections']]

#%%
size = [len(all_keep[i]) * 10 for i in range(len(all_keep))]
# size = [5000] * 16
plt.figure(figsize=(20, 20))

# layout_2 = nx.kamada_kawai_layout(G)
# layout_1 = nx.spring_layout(G, iterations=10, scale=1, seed=1, fixed=['Interections'] + conditions, pos = fixed_positions)
layout_1 = nx.circular_layout(conditions, scale=0.1)
layout_1['Interections'] = (0.01, 0.011)

fixed_positions_2 = {conditions[11]:tuple(layout_1[conditions[11]])}
layout_2 = nx.spring_layout(conditional_11 + [conditions[11]], k=1/214**0.5+0.1, iterations=5, scale=10, seed=1, 
fixed=[conditions[11]], pos=fixed_positions_2)
# layout_2=nx.circular_layout(conditional_11, scale=0.3)
# layout_2[conditions[11]] = layout_1[conditions[11]]

for key in layout_2.keys():
    if key != conditions[11]:
        layout_2[key] = np.array([layout_2[key][0], layout_2[key][1]-0.1])   

fixed_positions_3 = {conditions[15]:tuple(layout_1[conditions[15]])}
layout_3 = nx.spring_layout(conditional_15 + [conditions[15]], k=1/214**0.5+0.1, iterations=5, scale=10, seed=1, 
fixed=[conditions[15]], pos=fixed_positions_3)

for key in layout_3.keys():
    if key != conditions[15]:
        layout_3[key] = np.array([layout_3[key][0]+0.17, layout_3[key][1]*0.7])   

# layout_2 = nx.spring_layout(G, k=1/214**0.5-0.1, iterations=21, scale=1)
# layout_2 = nx.spectral_layout(G)
# layout_1 = nx.multipartite_layout(G)

nx.draw_networkx_edges(G, layout_1, edgelist=edgelist_conditions, edge_color='#AAAAAA')
nx.draw_networkx_nodes(G, layout_1, nodelist=['Interections'], node_size=10000, node_color='Red', edgecolors='black')
nx.draw_networkx_nodes(G, layout_1, nodelist=conditions, node_size=size, node_color='lightblue', edgecolors='black')

nx.draw_networkx_edges(G, layout_2, edgelist=edgelist_11, edge_color='#AAAAAA')
nx.draw_networkx_nodes(G, layout_2, nodelist=conditional_11, node_size=100, node_color='#fc8d62')

nx.draw_networkx_edges(G, layout_3, edgelist=edgelist_15, edge_color='#AAAAAA')
nx.draw_networkx_nodes(G, layout_3, nodelist=conditional_15, node_size=100, node_color='#fc8d62')


flat_names = [item[0]+'\n'+item[1] for item in conditions]

condition_dict = dict(zip(conditions, flat_names)) 
gene_dict_11 = dict(zip(conditional_11, conditional_11))
gene_dict_15 = dict(zip(conditional_15, conditional_15))
nx.draw_networkx_labels(G, layout_1, labels={'Interections' : 'Common' +'\n' + 'Essential Gene'}, font_size=20)
nx.draw_networkx_labels(G, layout_1, labels=condition_dict, font_size=15)
nx.draw_networkx_labels(G, layout_2, labels=gene_dict_11, font_size=5)
nx.draw_networkx_labels(G, layout_3, labels=gene_dict_15, font_size=5)

plt.axis('off')

plt.title("Conditionally Essential Genes Network")
plt.savefig('conditional-network.png')
plt.show()


# %%
