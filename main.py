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

sns.set_style('whitegrid')


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
While there are only a few ways to be alive, there are many to not being alive. In
the context of chemical reaction networks (CRNs), we know that certain CRNs related
with biology exhibit certain properties while we assume that other CRNs lack those.
However, the specific of those properties is rather fuzzy. In this server, we aim
to provide a catalogue of CRNs of different chemical systems (biotic, prebiotic and
abiotic) with the goal to enable a quick visualization of the data and to provide
the raw data for future analysis.

If you want to use this data or server for you research, please cite: 

    Paper citation
"""

st.header("Analysis")

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




st.subheader('Network Macrostatistic')

"""
As with every system in physics, we can look at the behaviour of each of the parts,
or we can just look at the system-level behaviour. Network macrostatistics enable us
to understand aspects of the network as how many connections does each component have
(mean degree), how many sub networks are disconnected inside the network, 
or how are the farthest nodes between them. 
"""


lc = n.subgraph(max(nx.strongly_connected_components(n), key=len))
network_features = [dict(
    size_nodes = nx.number_of_nodes(n),
    size_edges = nx.number_of_edges(n),
    n_strong_components = nx.number_strongly_connected_components(n), 
    n_weak_components = nx.number_weakly_connected_components(n),
    largest_component_size = nx.number_of_nodes(lc)
)]
network_features = pd.DataFrame.from_records(network_features).T
network_features.columns = ['value']
st.write(network_features)

st.subheader('Degree distribution')

"""
Critical systems are usually associated with heavy tailed distributed attributes. Since the
early 2000s, the degree distribution (probability of each node to have an X number of 
connections) of many networks have been associated to these heavy tailed degree distributions
— though there is a certain level of controvery about which is the right distribution—. 
"""

"""
In this part, we compare the fits of a power-law degree distribution to the fit of an
exponential degree distribution. The following figure shows the complementary 
cumulative distribution functon (CCDF) of the empirical data (red), the power-law
fit (blue), and the exponential fit (gray).

"""

degree_distribution = get_compounds_degrees(n)
pwl_fit = pwl.Fit(degree_distribution)


fig, ax = plt.subplots(1)
fig.set_size_inches(6.0, 3.0)
pwl_fit.plot_ccdf(ax=ax, linewidth=1, marker='.', color='#298ACC', label='data')
pwl_fit.power_law.plot_ccdf(ax=ax, linewidth=2, linestyle='-', color='#CC2A29', label='power-law')
pwl_fit.exponential.plot_ccdf(ax=ax, linewidth=2, linestyle='-', color='gray', label='exponential')
ax.set_xlabel('K - degree')
ax.set_ylabel('CCDF(k)')
ax.legend()
st.pyplot(fig)

"""
The following table shows the log-ratio and the p values of the comparisson
of the likelihood of both models. A ratio higher than 0 considers the
log-normal model a better fit, but such value is only significative if the p-value
of that comparisson is smaller than a given threshold (e.g. 0.01).
"""

degree_test = [dict(
    alpha = pwl_fit.alpha,
    p_value=pwl_fit.loglikelihood_ratio('exponential', 'power_law')[0],
    log_ratio=pwl_fit.loglikelihood_ratio('exponential', 'power_law')[1],
)]

degree_test = pd.DataFrame.from_records(degree_test).T
degree_test.columns = ['value']
st.write(degree_test)


st.subheader('Mass-Degree')

"""
While we would love to know the distribution of energy and abundance of each system,
that is usually intractable unless we have complete understanding of the partition
functions. Alternatively, we can consider that the abundance of a compound might correlate
with its abundance (more abundanct compounds could be more likely to react), and that
the energy of a compound might correlate with its mass (larger molecules have more chemical
bonds). Therefore, the relationship between mass and degree can provide a very
coarse but still interesting view about the nature of the system. For instance, a system 
where reactions carry few changes in enthalpy could lead to a system where elemental
species and simpler molecules should be much more connected than those more 
complex molecules.
"""

"""
The following figure shows the log-log plot of the degree and the mass of each
of the compounds of the network.
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
    The most recent research in graph analysis is showing that very different networks
    can still provide similar network macro-statistics. The analysis of micro-motifs
    enables us to dive deeper into the understanding of our networks by understanding
    how frequent are the different ways to connect systems among them.
    """

    """
    In the case of CRNs, we are specially intersted in understanding how two reactions
    can relate to each other. We look for 7 motifs that we describe in the following figure.
    In the smallest networks, we attempt to generate randomized versions of the network to
    compare the abundance of each motif against a random network with the same global statistics.
    """

with col2:

    st.image('server-micromotifs-explanation.png', width=300)


motifs = pd.read_json('so23.micromotifs.json')
motifs = motifs.query(f'network == "{choice}"')

if len(motifs) > 0:

    motifs_melt = motifs.melt(
        id_vars=['network'], value_vars=['m0', 'm1', 'm2', 'm3', 'm4', 'mx']
    )

    g = sns.catplot(data=motifs_melt, x='variable', y='value', kind='bar', height=3.0, aspect=2.0, errorbar='sd')
    g.set_xlabels('Micro-Motif')
    g.set_ylabels('Frequency')
    st.pyplot(g)

else:
    st.write("**Unfortunately, this network has some issue that prevents us from displaying that data**")

st.subheader("Reactions and molecules")

"""
This chemical reaction networks are collections of reactions and molecules that
happen to be in the real world. Although curating and annotating individually each of these
reactions is a daunting task, we love to dive deeper into which compounds and reactions are
more central across the different networks.
"""


some_reactions = pd.read_json('so23.smiles.json')


reaction_subset = some_reactions.query(f'network == "{choice}"')


filter_ = st.text_input(label='filter', value='')
if filter_ != '':
    reaction_subset = reaction_subset.set_index('smiles').filter(like=filter_, axis=0).reset_index()
reaction_subset = reaction_subset.sample(n=min(len(reaction_subset), 20))[['id', 'smiles']]
st.dataframe(reaction_subset)


st.subheader("Download")
bf = io.BytesIO()
nx.write_gml(n, bf)
st.download_button(
    'Download!', data=bf, file_name=choice + '.gml'
)