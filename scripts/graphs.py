import argparse
from Bio import SeqIO, Phylo, Seq
from collections import defaultdict, Counter, OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from scipy.stats import kstest
from scipy import stats

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="reconstruct branches from root",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--output', required=True, type=str, help="output png file")
    parser.add_argument('--tsv', required=True, type=str, help="output tsv file")
    parser.add_argument('--length', type=str, help="length of duplication")
    parser.add_argument('--tree', required=True, type=str, help="input newick tree")
    parser.add_argument('--input', type=str, help="input fasta file")
    args = parser.parse_args()

    just_the_duplication = SeqIO.parse(args.input, 'fasta')
    copy1, synonymous_1, nonsynonymous_1, copy2, nonsynonymous_2, synonymous_2, onlyone, synonymous_one, nonsynonymous_one  = (defaultdict(list) for i in range(9))
    tree_ = Phylo.read(args.tree, "newick")
    tree_.root_at_midpoint()
    tree_.find_clades()

    for entry in just_the_duplication:
        if '-' not in entry.seq:
            copy_1 = entry.seq[:int(args.length)][1:-2]
            copy_2 = entry.seq[int(args.length):][1:-2]
            for i in range(0, len(copy_1), 3):
                copy1[entry.id].append(copy_1[i:i+3])
            for i in range(0, len(copy_2), 3):
                copy2[entry.id].append(copy_2[i:i+3])
        else:
            only_one = entry.seq.replace("-", "")
            if len(only_one)== int(args.length):
                for i in range(0, len(entry.seq[1:-2]), 3):
                    onlyone[entry.id].append(only_one[i:i+3])
    for branch in tree_.get_nonterminals(order='postorder'):
        if branch.name in copy1:
            for b in branch:
                if b.name in copy1:
                    index = 0
                    for codon_branch, codon_b in zip(copy1[branch.name], copy1[b.name]):
                        if codon_branch != codon_b:
                            if Seq.translate(codon_branch) == Seq.translate(codon_b):
                                pos = 0
                                for char_branch, char_b in zip(codon_branch, codon_b):
                                    pos +=1
                                    if char_branch != char_b:
                                        synonymous_1[b.name].append(f'{char_b}{char_branch}{pos+(index*3)}')
                            else:

                                pos = 0
                                for char_branch, char_b in zip(codon_branch, codon_b):
                                    pos +=1
                                    if char_branch != char_b:
                                        nonsynonymous_1[b.name].append(f'{char_b}{char_branch}{pos+(index*3)}')
                        index+=1

        if branch.name in copy2:
            for b in branch:
                if b.name in copy2:
                    index = 0
                    for codon_branch, codon_b in zip(copy2[branch.name], copy2[b.name]):
                        if codon_branch != codon_b:
                            if Seq.translate(codon_branch) == Seq.translate(codon_b):
                                pos = 0
                                for char_branch, char_b in zip(codon_branch, codon_b):
                                    pos +=1
                                    if char_branch != char_b:
                                        synonymous_2[b.name].append(f'{char_b}{char_branch}{pos+(index*3)}')
                            else:
                                pos = 0
                                for char_branch, char_b in zip(codon_branch, codon_b):
                                    pos +=1
                                    if char_branch != char_b:
                                        nonsynonymous_2[b.name].append(f'{char_b}{char_branch}{pos+(index*3)}')
                        index+=1

    for branch in tree_.get_nonterminals(order='postorder'):
        if branch.name in synonymous_1:
            for b in branch:
                if b.name in synonymous_1:
                    synonymous_1[b.name] = list(set(synonymous_1[b.name]).difference(set(synonymous_1[branch.name])))
        if branch.name in nonsynonymous_1:
            for b in branch:
                if b.name in synonymous_1:
                    nonsynonymous_1[b.name] = list(set(nonsynonymous_1[b.name]).difference(set(nonsynonymous_1[branch.name])))

    for branch in tree_.get_nonterminals(order='preorder'):
        sort_branch = []
        for e in nonsynonymous_1[branch.name]:
            sort_branch.append(str("".join(sorted(e[:2], key=str.lower))+ e[2:]))

        for b in branch:
            if b.name in nonsynonymous_1:
                sort_b = []
                for entry in nonsynonymous_1[b.name]:
                    sort_b.append(str("".join(sorted(entry[:2], key=str.lower))+ entry[2:]))
                nonsynonymous_1[b.name] = set(set(sort_b).difference(set(sort_branch)))

    for branch in tree_.get_nonterminals(order='postorder'):
        if branch.name in synonymous_2:
            for b in branch:
                if b.name in synonymous_2:
                    synonymous_2[b.name] = list(set(synonymous_2[b.name]).difference(set(synonymous_2[branch.name])))

        if branch.name in nonsynonymous_2:
            for b in branch:
                if b.name in nonsynonymous_2:
                    nonsynonymous_2[b.name] = list(set(nonsynonymous_2[b.name]).difference(set(nonsynonymous_2[branch.name])))

    for branch in tree_.get_nonterminals(order='preorder'):
        sort_branch = []
        for e in nonsynonymous_2[branch.name]:
            sort_branch.append(str("".join(sorted(e[:2], key=str.lower))+ e[2:]))
        for b in branch:
            if b.name in nonsynonymous_2:
                sort_b = []
                for entry in nonsynonymous_2[b.name]:
                    sort_b.append(str("".join(sorted(entry[:2], key=str.lower))+ entry[2:]))
                nonsynonymous_2[b.name] = set(set(sort_b).difference(set(sort_branch)))

    lst_s1, syn_1, lst_n1, nonsyn_1, lst_s2, syn_2, lst_n2, nonsyn_2, lst_sp, syn_one, lst_np, nonsyn_one = ([] for i in range(12))
    scaled_syn_one, scaled_syn_1, scaled_syn_2, scaled_nonsyn_one, scaled_nonsyn_1, scaled_nonsyn_2 = (dict() for i in range(6))

    for i in synonymous_1.values():
        ls = list(i)
        for j in ls: lst_s1.append(j)
    for item_ in lst_s1:
        syn_1.append(int(item_[2:]))

    for i in nonsynonymous_1.values():
        ls = list(i)
        numbers_ = []
        for it in ls:
            numbers_.append(it[2:])
        new_numbers_ = list(set(numbers_))
            
        for j in new_numbers_: lst_n1.append(j)
    for item_ in lst_n1:
        nonsyn_1.append(int(item_))

    for i in synonymous_2.values():
        ls = list(i)
        for j in ls: lst_s2.append(j)
    for item_ in lst_s2:
        syn_2.append(int(item_[2:]))

    for i in nonsynonymous_2.values():
        ls = list(i)
        numbers_ = []
        for it in ls:
            numbers_.append(it[2:])
        new_numbers_ = list(set(numbers_))
            
        for j in new_numbers_: lst_n2.append(j)
    for item_ in lst_n2:
        nonsyn_2.append(int(item_))

    for branch in tree_.get_nonterminals(order='postorder'):
        if branch.name in onlyone:
            for b in branch:
                if b.name in onlyone:
                    index = 0
                    for codon_branch, codon_b in zip(onlyone[branch.name], onlyone[b.name]):
                        if codon_branch != codon_b:
                            if Seq.translate(codon_branch) == Seq.translate(codon_b):
                                pos = 0
                                for char_branch, char_b in zip(codon_branch, codon_b):
                                    pos +=1
                                    if char_branch != char_b:
                                        entry_unsorted = f'{char_b}{char_branch}{pos+(index*3)}'
                                        synonymous_one[b.name].append(str("".join(sorted(entry_unsorted[:2], key=str.lower))+ entry_unsorted[2:]))
                            else:
                                pos = 0
                                for char_branch, char_b in zip(codon_branch, codon_b):
                                    pos +=1
                                    if char_branch != char_b:
                                        entry_unsorted = f'{char_b}{char_branch}{pos+(index*3)}'
                                        nonsynonymous_one[b.name].append(str("".join(sorted(entry_unsorted[:2], key=str.lower))+ entry_unsorted[2:]))
                        index+=1

    for branch in tree_.get_nonterminals(order='postorder'):
        if branch.name in synonymous_one:
            for b in branch:
                if b.name in synonymous_one:
                    synonymous_one[b.name] = list(set(synonymous_one[b.name]).difference(set(synonymous_one[branch.name])))

    for branch in tree_.get_nonterminals(order='preorder'):
        sort_branch = []
        for e in nonsynonymous_one[branch.name]: sort_branch.append(str("".join(sorted(e[:2], key=str.lower))+ e[2:]))

        for b in branch:
            if b.name in nonsynonymous_one:
                sort_b = []
                for entry in nonsynonymous_one[b.name]:
                    sort_b.append(str("".join(sorted(entry[:2], key=str.lower))+ entry[2:]))
                nonsynonymous_one[b.name] = set(set(sort_b).difference(set(sort_branch)))

    for muts in synonymous_one.values():
        ls = list(muts)
        for mut_lst in ls: lst_sp.append(mut_lst)
    for item_ in lst_sp:
        syn_one.append(int(item_[2:]))
    print(syn_one) # locations of each mutation
    for i in nonsynonymous_one.values():
        ls = list(i)
        numbers_ = []
        for it in ls:
            numbers_.append(it[2:])
        new_numbers_ = list(set(numbers_))
        for j in new_numbers_: lst_np.append(j)
    for item_ in lst_np:
        nonsyn_one.append(int(item_))

    total_len = tree_.total_branch_length()
    file_ = SeqIO.parse(args.input,"fasta")
    seq_dict = dict()
    for record in file_:
        seq_dict[record.id] = record.seq
    for branch in tree_.get_nonterminals(order='preorder'):
        if pd.isna(branch.name) == False:
            if '-'*int(args.length) not in seq_dict[branch.name]:
                with_dupl = branch.total_branch_length()
                break
    without_dupl = total_len-with_dupl
    
    predupl_ = [syn_one, nonsyn_one]
    predupl_dicts = [scaled_syn_one, scaled_nonsyn_one]
    post_dupl = [syn_1, syn_2, nonsyn_1, nonsyn_2]
    post_dupl_dicts = [scaled_syn_1, scaled_syn_2, scaled_nonsyn_1, scaled_nonsyn_2]

    for a, b in zip(predupl_, predupl_dicts):
        for key, entry in Counter(a).items():
            b[key] = (entry/without_dupl)/int(args.length)


    for a, b in zip(post_dupl, post_dupl_dicts):
        for key, entry in Counter(a).items():
            b[key] = (entry/with_dupl)/int(args.length) #normalisation by length of gene and tree with duplication


    od = OrderedDict(sorted(scaled_syn_1.items()))
    x = list(od.values())
    res = np.cumsum(x)
    cumulative_syn, cumulative_2_syn, cumulative_one_syn = (dict() for i in range(3))
    for i, j in zip(od.keys(), res): 
        cumulative_syn[i] = j
    od2 = OrderedDict(sorted(scaled_syn_2.items()))
    x2 = list(od2.values())
    res2 = np.cumsum(x2)
    for i, j in zip(od2.keys(), res2):
        cumulative_2_syn[i] = j
    od_ = OrderedDict(sorted(scaled_syn_one.items()))
    x_ = list(od_.values())
    res_ = np.cumsum(x_)
    for i, j in zip(od_.keys(), res_):
        cumulative_one_syn[i] = j

    fig, axs = plt.subplots(1,2)
    axs[0].step(cumulative_syn.keys(), cumulative_syn.values(), label= f'1st copy postduplication')
    axs[0].step(cumulative_2_syn.keys(), cumulative_2_syn.values(), label=f'2nd copy postduplication' )
    axs[0].step(cumulative_one_syn.keys(), cumulative_one_syn.values(), label= f'preduplication')
    
    od = OrderedDict(sorted(scaled_nonsyn_1.items()))
    x = list(od.values())
    res = np.cumsum(x)
    cumulative_, cumulative_one, cumulative_2 = (dict() for i in range(3))
    for i, j in zip(od.keys(), res): cumulative_[i] = j
    od2 = OrderedDict(sorted(scaled_nonsyn_2.items()))
    x2 = list(od2.values())
    res2 = np.cumsum(x2)
    for i, j in zip(od2.keys(), res2): cumulative_2[i] = j
    od_ = OrderedDict(sorted(scaled_nonsyn_one.items()))
    x_ = list(od_.values())
    res_ = np.cumsum(x_)
    for i, j in zip(od_.keys(), res_): cumulative_one[i] = j
    axs[1].step(cumulative_.keys(), cumulative_.values(), label= '1st copy postduplication')
    axs[1].step(cumulative_2.keys(), cumulative_2.values(), label='2nd copy postduplication' )
    axs[1].step(cumulative_one.keys(), cumulative_one.values(), label='preduplication')
    axs[1].legend(loc='lower center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
    
    fig.suptitle('Synonymous and non-synonymous Mutations')
    plt.savefig(args.output)
    #synone, syn_1 and syn_2 are lists (before normalisation by division by tree branch length and length of duplication)
    #they each show the number of times a mutation occurs at a particular location, e.g. 5, 5, 5, 3, 3, 3, 3, 6, 1, 12, ... etc
    # the part of the workflow below calculates the Kolmogorov-Smirnov Statistics for them
    dictionary_ = {'Preduplication statistic': [" ", stats.ks_2samp(syn_one, syn_1),  stats.ks_2samp(syn_one, syn_2)], 'PostDuplication Copy 1 statistic': [stats.ks_2samp(syn_one, syn_1), " ", stats.ks_2samp(syn_1, syn_2)], "PostDuplication Copy 2 statistic": [stats.ks_2samp(syn_one, syn_2), stats.ks_2samp(syn_1, syn_2),  " "]}
    df = pd.DataFrame(dictionary_)
    df.index = ['Preduplication', 'PostDuplication Copy 1', "PostDuplication Copy 2"]
    df.to_csv(args.tsv + "nonsyn.tsv", sep='\t')

    dictionary_ = {'Preduplication statistic': [" ", stats.ks_2samp(nonsyn_one, nonsyn_1),  stats.ks_2samp(nonsyn_one, nonsyn_2)], 'PostDuplication Copy 1 statistic': [stats.ks_2samp(nonsyn_one, nonsyn_1), " ", stats.ks_2samp(nonsyn_1, nonsyn_2)], "PostDuplication Copy 2 statistic": [stats.ks_2samp(nonsyn_one, nonsyn_2), stats.ks_2samp(nonsyn_1, nonsyn_2),  " "]}
    df = pd.DataFrame(dictionary_)
    df.index = ['Preduplication', 'PostDuplication Copy 1', "PostDuplication Copy 2"]
    df.to_csv(args.tsv+"_syn.tsv", sep='\t')


    #Poisson Distribution


number_of_muts_syn = [len(syn_one), len(syn_1), len(syn_2)]
number_of_muts_nonsyn = [len(nonsyn_one), len(nonsyn_1), len(nonsyn_2)]

for nr_of_muts in number_of_muts_syn:
    distr =[]
    for i in np.arange(0.001, 1.001, 0.001):
        poisson = (((i*without_dupl)**(nr_of_muts))*math.exp(-i*without_dupl))/nr_of_muts
        distr.append(poisson)
    print(max(distr))

for nr_of_muts in number_of_muts_nonsyn:
    distr =[]
    for i in np.arange(0.001, 1.001, 0.001):
        poisson = (((i*without_dupl)**(nr_of_muts))*math.exp(-i*without_dupl))/nr_of_muts
        distr.append(poisson)
    print(max(distr))
