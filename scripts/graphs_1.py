import argparse, math
from Bio import SeqIO, Phylo, Seq
from collections import defaultdict, Counter, OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kstest
from scipy import stats

def separate_duplications(duplicationfile, lengthofdupl):
    """
    Divides sequences from duplication file into preduplication, and postduplication copies 1 and 2
    """
    preduplication, postduplication_1, postduplication_2 = (defaultdict(list) for i in range(3))
    for entry in duplicationfile:
        if '-' not in entry.seq:
            copy_1 = entry.seq[:int(lengthofdupl)][1:-2]
            copy_2 = entry.seq[int(lengthofdupl):][1:-2]
            for i in range(0, len(copy_1), 3):
                postduplication_1[entry.id].append(copy_1[i:i+3])
            for i in range(0, len(copy_2), 3):
                postduplication_2[entry.id].append(copy_2[i:i+3])
        else:
            preduplication_ = entry.seq.replace("-", "")
            if len(preduplication_)== int(lengthofdupl):
                for i in range(0, len(entry.seq[1:-2]), 3):
                    preduplication[entry.id].append(preduplication_[i:i+3])
    return(postduplication_1, postduplication_2, preduplication)


def recursive_mutations(treefile, copy):
    """
    Makes basic first reconstruction of common mutations moved to common ancestor - refined by function below
    """
    synonymous_, nonsynonymous_ = (defaultdict(list) for i in range(2))
    for branch in treefile.get_nonterminals(order='postorder'):
        if branch.name in copy:
            for b in branch:
                if b.name in copy:
                    index = 0
                    for codon_branch, codon_b in zip(copy[branch.name], copy[b.name]):
                        if codon_branch != codon_b:
                            if Seq.translate(codon_branch) == Seq.translate(codon_b):
                                pos = 0
                                for char_branch, char_b in zip(codon_branch, codon_b):
                                    pos +=1
                                    if char_branch != char_b:
                                        entry_unsorted = f'{char_b}{char_branch}{pos+(index*3)}'
                                        synonymous_[b.name].append(str("".join(sorted(entry_unsorted[:2], key=str.lower))+ entry_unsorted[2:]))
                            else:
                                pos = 0
                                for char_branch, char_b in zip(codon_branch, codon_b):
                                    pos +=1
                                    if char_branch != char_b:
                                        entry_unsorted = f'{char_b}{char_branch}{pos+(index*3)}'
                                        nonsynonymous_[b.name].append(str("".join(sorted(entry_unsorted[:2], key=str.lower))+ entry_unsorted[2:]))
                        index+=1
    return(nonsynonymous_, synonymous_)

def refine_recursive(treefile, synonymous, nonsynonymous):
    """
    Moves common mutations to common ancestor - up the tree. Returns synonymous and nonsynonymous dictionaries.
    """
    for branch in treefile.get_nonterminals(order='postorder'):
        if branch.name in synonymous:
            for b in branch:
                if b.name in synonymous:
                    synonymous[b.name] = list(set(synonymous[b.name]).difference(set(synonymous[branch.name])))
        if branch.name in nonsynonymous:
            for b in branch:
                if b.name in synonymous:
                    nonsynonymous[b.name] = list(set(nonsynonymous[b.name]).difference(set(nonsynonymous[branch.name])))
    for branch in treefile.get_nonterminals(order='preorder'):
        sort_branch = []
        for e in nonsynonymous[branch.name]:
            sort_branch.append(str("".join(sorted(e[:2], key=str.lower))+ e[2:]))
        for b in branch:
            if b.name in nonsynonymous:
                sort_b = []
                for entry in nonsynonymous[b.name]:
                    sort_b.append(str("".join(sorted(entry[:2], key=str.lower))+ entry[2:]))
                nonsynonymous[b.name] = set(set(sort_b).difference(set(sort_branch)))
    return(nonsynonymous, synonymous)


def mutations(synonymous, nonsynonymous):
    """
    Returns lists of synonymous and nonsynonymous mutation locations in the duplication
    """
    syn_, nonsyn_, lst_s, lst_n = ([] for i in range(4))
    for i in synonymous.values():
        ls = list(i)
        for j in ls: lst_s.append(j)
    for item_ in lst_s: syn_.append(int(item_[2:]))
    for i in nonsynonymous.values():
        ls = list(i)
        numbers_ = []
        for it in ls:numbers_.append(it[2:])
        new_numbers_ = list(set(numbers_))
        for j in new_numbers_: lst_n.append(j)
    for item_ in lst_n: nonsyn_.append(int(item_))
    return(nonsyn_, syn_)

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

    duplication_file = SeqIO.parse(args.input, 'fasta')
    synonymous_one, nonsynonymous_one  = (defaultdict(list) for i in range(2))
    tree_file  = Phylo.read(args.tree, "newick")
    tree_file.root_at_midpoint()
    tree_file.find_clades()

    #length of tree with and without duplication
    total_len = tree_file.total_branch_length()
    seq_dict = dict()
    for record in duplication_file:
        seq_dict[record.id] = record.seq
    for branch in tree_file.get_nonterminals(order='preorder'):
        if pd.isna(branch.name) == False:
            if '-'*int(args.length) not in seq_dict[branch.name]:
                with_dupl = branch.total_branch_length()
                break
    without_dupl = total_len-with_dupl

    copy1, copy2, preduplication = separate_duplications(duplication_file, args.length)
    nonsynonymous_1, synonymous_1 = recursive_mutations(tree_file, copy1)
    nonsynonymous_2, synonymous_2 = recursive_mutations(tree_file, copy2)
    nonsynonymous_pre, synonymous_pre = recursive_mutations(tree_file, preduplication)
    nonsynonymous_1_refined, synonymous_1_refined = refine_recursive(tree_file, nonsynonymous_1, synonymous_1)
    nonsynonymous_2_refined, synonymous_2_refined = refine_recursive(tree_file, nonsynonymous_2, synonymous_2)
    nonsynonymous_pre_refined, synonymous_pre_refined = refine_recursive(tree_file, nonsynonymous_pre, synonymous_pre)
    nonsyn_1, syn_1 = mutations(synonymous_1_refined, nonsynonymous_1_refined)
    nonsyn_2, syn_2 = mutations(synonymous_2_refined, nonsynonymous_2_refined)
    nonsyn_pre, syn_pre = mutations(synonymous_pre_refined, nonsynonymous_pre_refined)

    scaled_syn_one, scaled_syn_1, scaled_syn_2, scaled_nonsyn_one, scaled_nonsyn_1, scaled_nonsyn_2 = (dict() for i in range(6))
    
    predupl_ = [syn_pre, nonsyn_pre]
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
    number_of_muts_syn = [len(syn_one), len(syn_1), len(syn_2)] #number of synonymous mutations in copy of the duplication
    number_of_muts_nonsyn = [len(nonsyn_one), len(nonsyn_1), len(nonsyn_2)] #number of nonsynonymous mutations in each copy of the duplication

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
