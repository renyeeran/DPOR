"""
@ Ni, C. C., Lin, Y. Y., Gao, J., Gu, X., and Saucan, E. (2015). Ricci curvature of the Internet topology. 2015 IEEE Conference on Computer Communications (INFOCOM). IEEE.
@ Ni, C. C., Lin, Y. Y., Luo, F., and Gao, J. (2019). Community Detection on Networks with Ricci Flow. Scientific reports, 9(1), 9984.

description: This function computes the ORC value for each gene pair with the scRNA-Seq adjacent matrix.

usage: CompORC(adj)

arguments: 
adj: The output value of function Process_Hs or Process_Mm.

value:
curve: The array of Ollivier Ricci Curvature values between pairs of genes.
"""


from GraphRicciCurvature.OllivierRicci import OllivierRicci

import pandas as pd
import networkx as nx
import numpy as np


def CompORC(adj):
    G = nx.Graph()

    n = adj.shape[0]
    adj = np.asarray(adj)
    
    for i in range(n):
        for j in range(n):
            if adj[i,j]==1:
                G.add_edge(i, j) 

    orc = OllivierRicci(G, alpha=0.5, verbose="TRACE")
    orc.compute_ricci_curvature()
    G_orc = orc.G.copy()  # save an intermediate result
    
    curve = np.zeros((n, n))
    for n1, n2 in list(G.edges()):
        curve[n1][n2] = float(G_orc[n1][n2]["ricciCurvature"])
        curve[n2][n1] = float(G_orc[n1][n2]["ricciCurvature"])

    return curve

