import dendropy
import numpy as np
import re
import os


def ase_revscript_writer(reps, rev_fn, out_fn_suffix, hs_rates):
    for i in range(0, reps ):
        for j in range(0, len(hs_rates)):
            with open(rev_fn, "r") as f:
                lines = f.readlines()
                lines.append("q()")
                rep_lines = lines
                rep_lines[1] = lines[1].replace("hsrate_[0-9].?[0-9]*_[0-9]*",
                                                "hsrate_" + str(hs_rates[j]) + "_" + str(i))
                out_fn = "data/" + "hsrate_" + str(hs_rates[j]) + "_" + str(i) + str(out_fn_suffix)
                with open(out_fn, "w") as ofn:
                    ofn.writelines(rep_lines)


def print_branching_times(tree, brtimes_out_fn):
    node_ages = np.array(tree.calc_node_ages())
    node_ages = np.ma.masked_equal(node_ages, 0)
    node_ages = node_ages.compressed()
    node_ages = -np.sort(-node_ages)
    node_ages_2 = node_ages - (node_ages * 0.05)
    with open(brtimes_out_fn, "w") as ofn:
        for i in range(0, len(node_ages)):
            print(str(node_ages[i]) + "\t" + str(node_ages_2[i]), file=ofn)



def print_nexus_range_file(ofn, host_tree, symb_tree):
    lines = ["#NEXUS\n\n", "Begin data;\n"]

    total_time = host_tree.max_distance_from_root()
    num_host_lineages = len(host_tree.leaf_edges())
    num_symb_lineages = len(symb_tree.leaf_edges())

    nranges = num_host_lineages
    nsymb = num_symb_lineages

    lines.append("Dimensions ntax=" + str(nsymb) + " nchar=" + str(nranges) + ";\n")
    lines.append("Format datatype=Standard missing=? gap=- labels=\"01\";\n")
    lines.append("Matrix\n")

    host_tree_newick = host_tree.as_string("newick")
    symb_tree_newick = symb_tree.as_string("newick")

    host_tips = re.findall("T\d+", host_tree_newick)
    symb_tips = re.findall("T\d+", symb_tree_newick)
    for i in range(0, len(symb_tips) - 1):
        curr_tip = symb_tips[i].split("_")
        curr_tip = curr_tip[0]
        range_data = ""
        for j in host_tips:
            if curr_tip == j:
                range_data += "1"
            else:
                range_data += "0"
        lines.append("\t" + curr_tip + "_1" + "\t" + range_data + "\n")
    lines.append("\t;")
    lines.append("\nEnd;")
    with open(ofn, "w") as f:
        f.writelines(lines)


def print_connectivity_graph(host_tree, of_prefix):
    total_time = host_tree.max_distance_from_root()
    num_host_lineages = host_tree.num_lineages_at(total_time)

    num_epochs = num_host_lineages
    for i in range(num_epochs):
        connectivity_mat = np.zeros((num_epochs, num_epochs), dtype=int)
        connectivity_mat[0:i, 0:i] = 1
        ofn = of_prefix + "." + str(i) + ".txt"
        with open(ofn, "w") as ff:
            np.savetxt(ff, connectivity_mat, fmt="%i", delimiter=" ", newline="\n")


def print_rev_script(dir_name, prefix_fn, i):

    symb_tree_suffix_fn = "_0.loc.tre"
    rev_script_dir = "src/rev-scripts/"
    brtimes_suffix_fn = ".times.txt"
    ranges_suffix_fn = ".range.nex"

    rev_out_fn = rev_script_dir + prefix_fn[:-1] + "/run_epoch_" + prefix_fn + str(i) + ".Rev"
    os.makedirs(rev_script_dir + prefix_fn[:-1], exist_ok=True)
    with open(rev_out_fn, "w") as rev_file:
        range_fn = dir_name + prefix_fn + str(i) + ranges_suffix_fn
        tree_fn = dir_name + prefix_fn + str(i) + "_pruned" + symb_tree_suffix_fn
        out_fn = dir_name + "output/" + prefix_fn + str(i)
        geo_fn = dir_name + prefix_fn + str(i)
        lines = []
        lines.append("range_fn = \"" + range_fn + "\"\n")
        lines.append("tree_fn = \"" + tree_fn + "\"\n")
        lines.append("out_fn = \"" + out_fn + "\"\n")
        lines.append("geo_fn = \"" + geo_fn + "\"\n")
        lines.append("times_fn = geo_fn + \".times.txt\"\n")
        lines.append("dist_fn  = geo_fn + \".distances.txt\"\n")

        lines.append("moves = VectorMoves()\n")
        lines.append("monitors = VectorMonitors()\n")
        lines.append("n_gen = 5000\n")

        lines.append("dat_range_01 = readDiscreteCharacterData(range_fn)\n")
        lines.append("n_areas <- dat_range_01.nchar()\n")

        lines.append("max_areas <- 2\n")
        lines.append("n_states <- 0\n")
        lines.append("for (k in 0:max_areas) n_states += choose(n_areas, k)\n")

        lines.append("dat_range_n = formatDiscreteCharacterData(dat_range_01, \"DEC\", n_states)\n")
        lines.append("state_desc = dat_range_n.getStateDescriptions()\n")
        lines.append("state_desc_str = \"state,range\"\n")
        lines.append("for (i in 1:state_desc.size())\n")
        lines.append("{\n")
        lines.append("\tstate_desc_str += (i-1) + \",\" + state_desc[i] + \"\\n\"\n")
        lines.append("}\n")

        lines.append("write(state_desc_str, file=out_fn+\".state_labels.txt\")\n")

        lines.append("time_bounds <- readDataDelimitedFile(file=times_fn, delimiter=\" \")\n")
        lines.append("n_epochs <- time_bounds.nrows()\n")

        lines.append("for (i in 1:n_epochs) {\n")
        lines.append("\tepoch_fn = geo_fn + \".connectivity.\" + i + \".txt\"\n")
        lines.append("\tconnectivity[i] <- readDataDelimitedFile(file=epoch_fn, delimiter=\" \")\n")
        lines.append("}\n")

        lines.append("distances <- readDataDelimitedFile(file=dist_fn, delimiter=\" \")\n")
        lines.append("tree <- readTrees(tree_fn)[1]\n")

        lines.append("log10_rate_bg ~ dnUniform(-4,2)\n")
        lines.append("log10_rate_bg.setValue(-2)\n")
        lines.append("rate_bg := 10^log10_rate_bg\n")
        lines.append("moves.append( mvSlide(log10_rate_bg, weight=4) )\n")

        # lines.append("rate_bg <- 1.0\n")

        # lines.append("log_sd <- 0.5\n")
        # lines.append("log_mean <- ln(1) - 0.5*log_sd^2\n")
        # lines.append("dispersal_rate ~ dnExponential(1.0)\n")
        # lines.append("moves.append( mvScale(dispersal_rate, weight=5) )\n")
        lines.append("dispersal_rate <- 1.0\n")
        lines.append("distance_scale ~ dnUnif(0,20)\n")
        lines.append("distance_scale.setValue(0.01)\n")
        lines.append("moves.append( mvScale(distance_scale, weight=3) )\n")

        lines.append("for (i in 1:n_epochs) {\n")
        lines.append("\tfor (j in 1:n_areas) {\n")
        lines.append("\t\tfor (k in 1:n_areas) {\n")
        lines.append("\t\t\tdr[i][j][k] <- 0.0\n")
        lines.append("\t\t\tif (connectivity[i][j][k] > 0) {\n")
        lines.append("\t\t\t\tdr[i][j][k] := dispersal_rate\n")
        lines.append("\t\t\t}\n")
        lines.append("\t\t}\n")
        lines.append("\t}\n")
        lines.append("}\n")

        # lines.append("log_sd <- 0.5\n")
        # lines.append("log_mean <- ln(1) - 0.5*log_sd^2\n")
        # lines.append("extirpation_rate ~ dnLognormal(mean=log_mean, sd=log_sd)\n")
        # lines.append("moves.append( mvScale(extirpation_rate, weight=5) )\n")
        lines.append("extirpation_rate ~ dnExponential(1.0)\n")
        lines.append("moves.append( mvScale(extirpation_rate, weight=5) )\n")

        lines.append("for (i in 1:n_epochs) {\n")
        lines.append("\tfor (j in 1:n_areas) {\n")
        lines.append("\t\tfor (k in 1:n_areas) {\n")
        lines.append("\t\t\ter[i][j][k] <- 0.0\n")
        lines.append("\t\t}\n")
        lines.append("\t\ter[i][j][j] := extirpation_rate\n")
        lines.append("\t}\n")
        lines.append("}\n")

        lines.append("for (i in n_epochs:1) {\n")
        lines.append("Q_DEC[i] := fnDECRateMatrix(dispersalRates=dr[i], extirpationRates=er[i],"
                     " maxRangeSize=max_areas,"
                     "nullRange=\"CondSurv\")\n")
        lines.append("}\n")

        lines.append("for (i in 1:n_epochs) {\n")
        lines.append("\ttime_max[i] <- time_bounds[i][1]\n")
        lines.append("\ttime_min[i] <- time_bounds[i][2]\n")
        lines.append("\tif (i != n_epochs) {\n")
        lines.append("\t\tepoch_times[i] ~ dnUniform(time_min[i], time_max[i])\n")
        lines.append("\t\tmoves.append( mvSlide(epoch_times[i], delta=(time_max[i]-time_min[i])/2) )\n")
        lines.append("\t} else {\n")
        lines.append("\t\tepoch_times[i] <- 0.0\n")
        lines.append("\t}\n")
        lines.append("}\n")

        lines.append("Q_DEC_epoch := fnEpoch(Q=Q_DEC, times=epoch_times, rates=rep(1, n_epochs))\n")

        lines.append("clado_event_types <- [ \"s\", \"a\" ]\n")
        lines.append("p_sympatry ~ dnUniform(0,1)\n")
        lines.append("p_allopatry := abs(1.0 - p_sympatry)\n")
        lines.append("clado_type_probs := simplex(p_sympatry, p_allopatry)\n")
        lines.append("moves.append( mvSlide(p_sympatry, weight=2) )\n")
        lines.append("P_DEC := fnDECCladoProbs(eventProbs=clado_type_probs, eventTypes=clado_event_types,"
                     " numCharacters=n_areas, maxRangeSize=max_areas)\n")

        lines.append("rf_DEC <- rep(0, n_states)\n")
        lines.append("rf_DEC[2] <- 1\n")
        lines.append("rf_DEC_simp <- simplex(rf_DEC)\n")

        lines.append("m_bg ~ dnPhyloCTMCClado(tree=tree, Q=Q_DEC_epoch, cladoProbs=P_DEC, branchRates=rate_bg,"
                     " rootFrequencies=rf_DEC_simp, type=\"NaturalNumbers\", nSites=1)\n")
        lines.append("m_bg.clamp(dat_range_n)\n")
        lines.append("monitors.append( mnScreen(printgen=100, distance_scale, rate_bg, extirpation_rate) )\n")
        lines.append("monitors.append( mnModel(file=out_fn+\".model.log\", printgen=10) )\n")
        lines.append("monitors.append( mnFile(tree, filename=out_fn+\".tre\", printgen=10) )\n")
        lines.append("monitors.append( mnJointConditionalAncestralState(tree=tree, ctmc=m_bg, "
                     "type=\"NaturalNumbers\","
                     " withTips=true, withStartStates=true, filename=out_fn+\".states.log\", printgen=10) )\n")
        lines.append("mymodel = model(m_bg)\n")
        lines.append("mymcmc = mcmc(mymodel, monitors, moves)\n")
        lines.append("mymcmc.run(n_gen)\n")
        rev_file.writelines(lines)
    rev_file.close()