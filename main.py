import biogeo_cophylo_simstudy
import treeducken_tools as td
import os
import dendropy
import subprocess
from dendropy import Tree
import pandas as pd
import re

nt = 8
turnover = 0.5
br = 0.07
dr = 0.035
trr = [0.0, 0.035, 0.1, 0.5]
gbr = 0.0
gdr = 0.0
reps = 50
nl = 1
ng = 1
ipp = 1
ne = 1
sd1 = [1859, 1860, 1861, 1862]
sd2 = [2016, 2017, 2018, 2019]

num_settings_regimes = len(trr)

for i in range(0, num_settings_regimes):
    settings_regime_name = str(trr[i]) + "_hsrate"
    settings_dictionary = td.create_settings_dictionary(sbr=br, sdr=dr, lgtr=[i], gbr=gbr, gdr=gdr, num_loci=nl,
                                                        reps=reps, num_genes=1, ipp=1, ne=1, ntax=nt, screen_out=0,
                                                        ofn=settings_regime_name, sd1=sd1[i], sd2=sd2[i])

    settings = td.write_settings_file(settings_dictionary, settings_regime_name)


subprocess.call("src/create_directories_for_simdata.sh", shell=True)


def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file


# TODO: get treeducken function for treeducken_tools running
for s in files("settings/"):
    run_cmd = "treeducken" + " -i " + "settings/" + s
    so = os.popen(run_cmd).read()
    print(so)
    setting_prefix_list = s.split("_")
    setting_prefix = str(setting_prefix_list[0]) + "_" + str(setting_prefix_list[1])
    mv_cmd = "mv " + setting_prefix + "*" + " ./data/" + setting_prefix
    so2 = os.popen(mv_cmd).read()
    print(so2)

# run_sim_regime < - function(sim_dir, prefix_fn, num_reps, use_full_tree=FALSE)

host_tree_suffix_fn = ".sp.tre"
symb_tree_suffix_fn = "_0.loc.tre"

brtimes_suffix_fn = ".times.txt"
ranges_suffix_fn = ".range.nex"

data_dir = ["data/0.0_hsrate/", "data/0.1_hsrate/", "data/0.035_hsrate/", "data/0.5_hsrate/"]
prefix_fn = ["0.0_hsrate_", "0.1_hsrate_", "0.035_hsrate_", "0.5_hsrate_"]

for j in range(0, num_settings_regimes):

    for i in range(0, reps):
        # read in the host and symbiont trees
        host_fn = data_dir[j] + prefix_fn[j] + str(i) + host_tree_suffix_fn
        symb_fn = data_dir[j] + prefix_fn[j] + str(i) + symb_tree_suffix_fn

        host_tree = Tree.get(path=host_fn, schema="nexus", rooting="default-rooted")
        symb_tree = Tree.get(path=symb_fn, schema="nexus", rooting="default-rooted")

        # prune extinct tips on symb tree
        symb_newick_tree = symb_tree.as_string("newick")
        symb_newick_tree = re.sub(r'_\d+', '', symb_newick_tree)
        symb_tree = dendropy.Tree.get_from_string(symb_newick_tree, "newick")
        ext_labels = re.findall("X\d+", symb_newick_tree)
        symb_tree.prune_taxa_with_labels(ext_labels)
        symb_tree.update_bipartitions()
        symb_tree.write_to_path(data_dir[j] + prefix_fn[j] + str(i) + "_pruned" + symb_tree_suffix_fn, schema="nexus")
        symb_newick_tree = symb_tree.as_string('newick')

        # print out epoch times (no uncertainty in age here)
        out_fn = data_dir[j] + prefix_fn[j] + str(i)
        brtimes_out_fn = out_fn + brtimes_suffix_fn
        biogeo_cophylo_simstudy.print_branching_times(host_tree, brtimes_out_fn)

        # print out range files
        nex_range_outfn = out_fn + ranges_suffix_fn
        biogeo_cophylo_simstudy.print_nexus_range_file(nex_range_outfn, host_tree, symb_tree)

        # print out connectivity graphs
        connectivity_prefix_outfn = out_fn + ".connectivity"
        biogeo_cophylo_simstudy.print_connectivity_graph(host_tree=host_tree, of_prefix=connectivity_prefix_outfn)

        # print out "distance" matrix
        dist_mat = host_tree.phylogenetic_distance_matrix()
        distmat_ofn = out_fn + ".distances.txt"
        dm_df = pd.DataFrame(data=dist_mat, index=host_tree.taxon_namespace)
        biogeo_cophylo_simstudy.print_rev_script(data_dir[j], prefix_fn[j], i)

