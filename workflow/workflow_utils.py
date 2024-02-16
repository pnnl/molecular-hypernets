from IPython.display import Image, HTML
from urllib import parse
import requests

import numpy as np

from hnxwidget import HypernetxWidget
from gnpsdata import taskinfo

""" Summarize and Visualize """

def summarize_component(comp_idx, comps, node_data, s=1):
    H = comps[s][comp_idx]
    data = node_data.loc[H]
    summary = f"{s}-component ({comp_idx}):\n"
    summary += f"\t{H.number_of_nodes()} nodes, {H.number_of_edges()} edges\n"
    summary += f"\tmin m/z: {data['parent mass'].min()}, max m/z: {data['parent mass'].max()}\n"
    summary += f"\tmin RT: {data.RTConsensus.min()}, max RT: {data.RTConsensus.max()}"
    compounds = data.Compound_Name.dropna().unique().tolist()
    if compounds:
        summary += f"\n\tannotated compound(s): {'; '.join(compounds)}"
    adducts = data.Adduct.dropna().unique().tolist()
    if adducts:
        summary += f"\n\tadduct(s): {', '.join(adducts)}"
    print(summary)

def visualize_component(comp_idx, comps, node_labels, node_data, s=1):
    return HypernetxWidget(
        comps[s][comp_idx],
        node_labels=node_labels,
        node_data=node_data,
    )

""" Filter and Search """

def filter_LCMS(node_data, mz=None, mz_rtol=None, rt=None, rt_atol=None):
    filtered_data = node_data
    if mz is not None:
        if isinstance(mz, tuple):
            mz_filter = filtered_data['parent mass'].between(*mz)
        else:
            mz_filter = filtered_data['parent mass'].apply(np.isclose, args=(mz,), rtol=mz_rtol)
        filtered_data = filtered_data[mz_filter]
    
    if rt is not None:
        if isinstance(rt, tuple):
            rt_filter = filtered_data['RTConsensus'].between(*rt)
        else:
            rt_filter = filtered_data['RTConsensus'].apply(np.isclose, args=(rt,), atol=rt_atol)
        filtered_data = filtered_data[rt_filter]
    
    return filtered_data

def filter_annotation(node_data, compound=None, adduct=None):
    filtered_data = node_data
    if compound:
        compound_filter = filtered_data['Compound_Name'] == compound
        filtered_data = filtered_data[compound_filter]
    
    if adduct:
        adduct_filter = filtered_data['Adduct'] == adduct
        filtered_data = filtered_data[adduct_filter]
    return filtered_data

def search_components(comps, filtered_data, s=1, return_index=True):
    filtered_comps = []
    for i, c in enumerate(comps[s]):
        if filtered_data.index.isin(c.nodes).any():
            filtered_comps.append(i if return_index else c)
    return filtered_comps

""" USI and GNPSDashboard """

def get_usi(task, scan_id):
    file = taskinfo.get_task_information(task)['files'][0]
    return f"mzspec:GNPS:TASK-{task}-{file}:scan:{scan_id}"

def get_MS2_dash(usi):
    params = {"usi1": usi}
    return "https://metabolomics-usi.gnps2.org/dashinterface/?{}".format(parse.urlencode(params))

def get_MS2_png(usi):
    params = {"usi1": usi}
    return "https://metabolomics-usi.gnps2.org/png/?{}".format(parse.urlencode(params))

def display_spectra(usi):
    img_url = get_MS2_png(usi)
    img_data = requests.get(img_url).content
    return Image(img_data)

def link_spectra(usi):
    spectra_url = get_MS2_dash(usi)
    scan_id = usi.split(":")[-1]
    return HTML(f"""<a href="{spectra_url}" target="_blank">View MS2 Spectra (scan:{scan_id}) on GNPS</a>""")

def link_network(comp_idx, comps, g, s=1):
    network_url = g.nodes[list(comps[s][comp_idx])[0]]["GNPSLinkout_Network"]
    return HTML(f"""<a href="{network_url}" target="_blank">View MN Component on GNPS</a>""")