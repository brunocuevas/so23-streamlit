import streamlit as st
import networkx as nx
import pandas as pd
import powerlaw as pwl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
import sys
import streamlit.components.v1 as components
from yaml import CLoader as Loader, CDumper as Dumper
from yaml import load, dump
import boto3
import os
import io

s3 = boto3.resource('s3',
  endpoint_url = os.environ['R2_ENDPOINT_URL'],
  aws_access_key_id = os.environ['R2_KEY_ID'],
  aws_secret_access_key = os.environ['R2_KEY']
)

# sns.set_style('whitegrid')


def download_reference_file(file_name):
    if not os.path.isfile(file_name):
        with open(file_name, 'wb') as g:
            s3.Bucket('so23').download_fileobj(file_name, g)



def get_compounds_degrees(G: nx.DiGraph, direction="all"):
    degree = []
    
    for node, _ in filter(lambda x: x[1]['type'] == 'molecule', G.nodes(data=True)):

        if direction == "all":
            degree.append(G.degree[node])
        elif direction == "in":
            degree.append(G.in_degree[node])
        elif direction == "out":
            degree.append(G.out_degree[node])

    return degree

def get_molecule_degree_mass(G, network_label, network_type):
    out = []
    for n, data in filter(lambda x: x[1]['type'] == 'molecule', G.nodes(data=True)):
        try:
            out.append(dict(node=n, log_degree=np.log10(G.degree[n]), log_mw=np.log10(data['mw']), network=network_label, type=network_type))
        except KeyError:
            continue
    return out


def download_network_s3(id):
    f = io.BytesIO()
    # print(f'{id}.pdb')
    # s3.Bucket('nsdb').download_fileobj('nsdb-000001.pdb', f)
    s3.Bucket('so23').download_fileobj(f'{id}.gml', f)
    # path_to_pdb = pathlib.Path(os.getcwd()).parent / 'structures/pdb/{0}.pdb'.format(id)
    f.seek(0)
    return f


download_reference_file('so23.index.yaml')
download_reference_file('so23.micromotifs.json')
download_reference_file('so23.smiles.json')


data = load(open('so23.index.yaml', 'r'), Loader=Loader)
networks = pd.DataFrame.from_records(data['networks'])


st.title("Self Organization across chemical reaction networks")

st.header("Introduction")

"""
While there are only a few ways to be alive, there are countless ways not to be. 
In the context of chemical reaction networks (CRNs), we recognize that certain biologically-related CRNs exhibit specific properties, 
whereas others presumably do not. However, the specifics of these properties remain rather unclear. On this server, 
we offer a catalog of CRNs from various chemical systems—including biotic, prebiotic, and abiotic—aimed at enabling quick 
data visualization and providing raw data for future analysis.

If you want to use this data or server for you research, please cite: 

    TBA
"""


"""
We have curated a set of 16 chemical reaction networks, and we will be happy to include
many more in the future (and to keep updating ours). In the next widget, you can 
choose a chemical reaction network to visualize its properties.

"""

choice = st.selectbox(
    'Choose a network',
    networks['alias'].to_list(), index=None
)

if choice is None:
    sys.exit()

u = download_network_s3(choice)
n = nx.read_gml(u)

"""
This network was obtained from:
"""

st.write(f"*{networks.set_index('alias').loc[choice]['citation']}*")

st.subheader("Download")

"""
The resulting file is a .GML, which can be open in Gephi or with libraries
specialized in graphs (e.g. NetworkX)
"""

bf = io.BytesIO()
nx.write_gml(n, bf)
st.download_button(
    'Download!', data=bf, file_name=choice + '.gml'
)



st.subheader('Network Macrostatistic')

"""
As with every system in physics, we can examine the behavior of each individual part or focus on the system-level behavior. 
Network macrostatistics allow us to understand various aspects of a network, such as the average number of connections each component 
has (mean degree), the number of sub-networks that are disconnected within the network, or the distance between the farthest nodes.
"""


lc = n.subgraph(max(nx.strongly_connected_components(n), key=len))
network_features = [{
    "# nodes": nx.number_of_nodes(n),
    "# edges": nx.number_of_edges(n),
    "# strongly connected components": nx.number_strongly_connected_components(n), 
    "# weakly connected components": nx.number_weakly_connected_components(n),
    "Largest strong component size": nx.number_of_nodes(lc)
}]
network_features = pd.DataFrame.from_records(network_features).T
network_features.columns = ['value']
st.write(network_features)

st.subheader('Degree distribution')

"""
Critical systems are often associated with heavy-tailed distributed attributes. Since the early 2000s, the degree distribution—the probability of each node having a 
specific number of connections—of many networks has been linked to these heavy-tailed distributions, although there remains some controversy over which distribution is most appropriate.

In this section, we compare the fits of a power-law degree distribution to those of an exponential degree distribution. The following figure illustrates 
the complementary cumulative distribution function (CCDF) of the empirical data (in red), the power-law fit (in blue), and the exponential fit (in gray).

"""

degree_distribution = get_compounds_degrees(n)
pwl_fit = pwl.Fit(degree_distribution)
# bins, ccdf = pwl_fit.cdf(survival=True)
# pl_bins = np.unique(pwl.trim_to_range(pwl_fit.data, xmin=pwl_fit.power_law.xmin, xmax=pwl_fit.power_law.xmax))
# pl_ccdf = pwl_fit.power_law.cdf(pl_bins, survival=True)
# exp_bins = np.unique(pwl.trim_to_range(pwl_fit.data, xmin=pwl_fit.exponential.xmin, xmax=pwl_fit.exponential.xmax))
# exp_ccdf = pwl_fit.exponential.cdf(exp_bins, survival=True)

fig, ax = plt.subplots(1)
fig.set_size_inches(6.0, 3.0)
# ax.plot(bins, ccdf, linewidth=1, marker='.', color='#298ACC', label='data')
# ax.plot(pl_bins, pl_ccdf, linewidth=2, linestyle='-', color='#CC2A29', label='power-law')
# ax.plot(exp_bins, exp_ccdf, linewidth=2, linestyle='-', color='gray', label='exponential')
pwl_fit.plot_ccdf(ax=ax, linewidth=1, marker='.', color='#298ACC', label='data')
pwl_fit.power_law.plot_ccdf(ax=ax, linewidth=2, linestyle='-', color='#CC2A29', label='power-law')
pwl_fit.exponential.plot_ccdf(ax=ax, linewidth=2, linestyle='-', color='gray', label='exponential')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('K - degree')
ax.set_ylabel('CCDF(k)')
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()
st.pyplot(fig)

"""
The following table shows the log-ratio and the p values of the comparisson
of the likelihood of both models. A ratio higher than 0 considers the
log-normal model a better fit, but such value is only significative if the p-value
of that comparisson is smaller than a given threshold (e.g. 0.01).
"""

degree_test = [{
    "α":  pwl_fit.alpha,
    "p-value": pwl_fit.loglikelihood_ratio('exponential', 'power_law')[0],
    "log ratio": pwl_fit.loglikelihood_ratio('exponential', 'power_law')[1],
}
]

degree_test = pd.DataFrame.from_records(degree_test).T
degree_test.columns = ['value']
st.write(degree_test)


st.subheader('Mass-Degree')

"""
While we would ideally like to know the distribution of energy and the abundance of each system, this is typically 
intractable without a complete understanding of the partition functions. Alternatively, we can consider that the abundance 
of a compound might correlate with its reactivity—more abundant compounds could be more likely to react—and that the energy of a 
compound might correlate with its mass, as larger molecules generally have more chemical bonds. Therefore, the relationship between 
mass and degree can provide a rough but still insightful view of the nature of the system. For instance, in a system where reactions 
involve minimal changes in enthalpy, elemental species and simpler molecules might be much more interconnected than more complex molecules.


The following figure presents a log-log plot of the degree and mass of each compound in the network.
"""

mass_degree = pd.DataFrame.from_records(get_molecule_degree_mass(n, choice, ''))
try:
    log_mw = mass_degree['log_mw']
    log_degree = mass_degree['log_degree']
    flag = True
except KeyError:
    st.write("**Unfortunately, this network has some issue that prevents us from displaying that data**")
    flag = False

if flag:
    g = sns.lmplot(data=mass_degree, x='log_mw', y='log_degree', height=3.0, aspect=2.0)
    g.set_xlabels('log10(Mass) (logDa)')
    g.set_ylabels('log10(Degree)')
    st.pyplot(g)

    """
    Though this relationship does not have to be linear, we are still 
    interested in finding correlations between these two atributes. The following table
    describes the results of our linear-regression analysis.
    """

    lr = linregress(x=mass_degree['log_mw'], y=mass_degree['log_degree'])

    correlation_test = [{
        'slope':lr.slope,
        'p-value':lr.pvalue,
        'r-value':lr.rvalue,
        
    }]

    correlation_test = pd.DataFrame.from_records(correlation_test).T
    correlation_test.columns = ['value']
    st.write(correlation_test)


st.subheader('Network microstatiscs: Micro-motifs')

col1, col2 = st.columns(2)

with col1:
    """
    Recent research in graph analysis has shown that very different networks can still exhibit similar macro-statistics. By analyzing micro-motifs, we can gain deeper insights into our networks by understanding the frequency of different connection patterns among systems.

    In the case of chemical reaction networks (CRNs), we are particularly interested in understanding how two reactions can relate to each other. We focus on seven specific motifs, which are described in the following figure. For the smallest networks, we generate randomized versions to compare the abundance of each motif against a random network with the same global statistics.
    """

with col2:

    st.image('server-micromotifs-explanation.png', width=300)


motifs = pd.read_json('so23.micromotifs.json')
random_motifs = pd.read_csv('so23.micromotifs.random.csv', index_col=None)
for col in random_motifs.columns[1:-1]:
    random_motifs[col] = random_motifs[col] / random_motifs['total']
random_motifs['mx'] = random_motifs['m5'] + random_motifs['m6']
motifs = motifs.query(f'network == "{choice}"')
random_motifs = random_motifs.query(f'network == "{choice}"')
if len(motifs) > 0:

    motifs_melt = motifs.melt(
        id_vars=['network'], value_vars=['m0', 'm1', 'm2', 'm3', 'm4', 'mx']
    )

    motifs_melt['type'] = 'empirical'
    random_motifs_melt = random_motifs.melt(
        id_vars=['network'], value_vars=['m0', 'm1', 'm2', 'm3', 'm4', 'mx']
    )
    random_motifs_melt['type'] = 'randomized'

    motifs_melt = pd.concat([motifs_melt, random_motifs_melt])
    
    g = sns.catplot(data=motifs_melt, y='variable', x='value', height=3.0, aspect=2.0, hue='type',
                    palette={'empirical': "#51658C", 'randomized': "#D96A42"})
    
    g.set_ylabels('Micro-Motif')
    g.set_xlabels('Frequency')
    st.pyplot(g)

else:
    st.write("**Unfortunately, this network has some issue that prevents us from displaying that data**")

st.subheader("Reactions and molecules")

"""
Chemical reaction networks are collections of reactions and molecules that occur in the real world. Although curating and annotating each of these reactions individually is a daunting task, we are eager to explore which compounds and reactions are most central across the different networks.
"""


some_reactions = pd.read_json('so23.smiles.json')


reaction_subset = some_reactions.query(f'network == "{choice}"')


filter_ = st.text_input(label='filter', value='')
if filter_ != '':
    reaction_subset = reaction_subset.set_index('smiles').filter(like=filter_, axis=0).reset_index()
reaction_subset = reaction_subset.sample(n=min(len(reaction_subset), 20))[['id', 'smiles']]
st.dataframe(reaction_subset)

st.header("About")

"""
This project was conducted at the University of Wisconsin-Madison, in the laboratory of Professor Betül Kaçar. The project was primarily developed by Bruno Cuevas-Zuviría, who is also responsible for maintaining this site. For any inquiries, please contact bruno.czuviria [at] upm.es.
"""
col1, _,  col2, _, col3 = st.columns(5)
col1.image("KacarLab-Logo_Circle-Black.png")
col2.image("black-center-UWlogo-print.png")
col3.image("MUSE-Logo_Patch-Black.png")