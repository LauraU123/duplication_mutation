import argparse, math
from Bio import SeqIO, Phylo, Seq
from collections import defaultdict, Counter, OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

def separate_duplications(lengthofdupl):
    """
    Divides sequences from duplication file into preduplication, and postduplication copies 1 and 2
    By splitting the files into two parts based on duplication length
    
    inputs: length of duplication (72 for RSV-A and 60 for RSV-B)
            duplication file (fasta file containing aligned duplicated regions)
    outputs:
        - preduplication dictionary 
        - postduplication dictionary (first copy)
        -postduplication dictionary (second copy)
    
    """
    duplicationfile = SeqIO.parse(args.input, "fasta")
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
                preduplication_ = preduplication_[1:-2]
                for i in range(0, len(preduplication_), 3):
                    preduplication[entry.id].append(preduplication_[i:i+3])
    return(postduplication_1, postduplication_2, preduplication)


def recursive_mutations(treefile, copy):
    """
    Makes basic first reconstruction of common mutations moved to common ancestor. 
    Mutations are split into synonymous and nonsynonymous based on if they result in amino acid change or not. 

    inputs:
        - open newick tree file (rooted at midpoint)
        - copy: dictionary of sequences in a copy of the duplication
    outputs:
        - dictionary synonymous: rudimentary first reconstruction of the common mutations in the copy in dictionary form (synonymous)
        - dictionary nonsynonymous: same as above but nonsynonymous
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
    Refines the common mutation function described above. 
    Mutations present in multiple branches are moved further up the tree and exceptions due to sequencing errors are taken into account.

    inputs:
        dictionary synonymous (dictionary of branches and synonymous mutations)
        dictionary nonsynonymous (same as above but nonsynonymous)
        read nwk tree file rooted at midpoint, with the relevant branches

    outputs:
        dictionaries for synonymous and nonsynonymous common mutations.
    
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
    Returns locations of synonymous and nonsynonymous mutations in list format
    
    Inputs:
        synonymous mutation dictionary
        nonsynonymous mutation dictionary
    Outputs:
        - list of nonsynonymous mut location
        - list of synonymous mut location
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

def cumulative(dictionary):
    """
    Returns cumulative values for each element of an ordered dictionary
    Input:  
        dictionary with cumulative values for each ordered key 
    """
    od = OrderedDict(sorted(dictionary.items()))
    x = list(od.values())
    res = np.cumsum(x)
    cumul= dict()
    for i, j in zip(od.keys(), res): 
        cumul[i] = j
    return(cumul)

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

    copy1, copy2, preduplication = separate_duplications(args.length)
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
    
    print(syn_pre)
    print(syn_1)
    print(syn_2)


    cumulative_syn_1, cumulative_syn_2, cumulative_one_syn = cumulative(scaled_syn_1), cumulative(scaled_syn_2), cumulative(scaled_syn_one)
    cumulative_nonsyn_1, cumulative_nonsyn_2, cumulative_one_nonsyn = cumulative(scaled_nonsyn_1), cumulative(scaled_nonsyn_2), cumulative(scaled_nonsyn_one)

    print(cumulative_syn_1, cumulative_syn_2)
    #plotting the cumulative distributions
    fig, axs = plt.subplots(1,2)
    axs[0].step(cumulative_syn_1.keys(), cumulative_syn_1.values(), label= f'1st copy postduplication', where='post')
    axs[0].step(cumulative_syn_2.keys(), cumulative_syn_2.values(), label=f'2nd copy postduplication', where='post' )
    axs[0].step(cumulative_one_syn.keys(), cumulative_one_syn.values(), label= f'preduplication', where='post' )
    axs[1].step(cumulative_nonsyn_1.keys(), cumulative_nonsyn_1.values(), label= '1st copy postduplication', where='post' )
    axs[1].step(cumulative_nonsyn_2.keys(), cumulative_nonsyn_2.values(), label='2nd copy postduplication', where='post'  )
    axs[1].step(cumulative_one_nonsyn.keys(), cumulative_one_nonsyn.values(), label='preduplication', where='post' )
    axs[1].legend(loc='lower left', bbox_to_anchor=(1, 1.05))
    fig.suptitle('Synonymous and non-synonymous Mutations')
    plt.savefig(args.output, bbox_inches="tight")




    #plotting synonymous and nonsynonymous
    fig_nonsyn = plt.figure("nonsynonymous duplication")
    plt.title('Nonsynonymous Mutations')
    plt.xlabel("Location within the duplication")
    plt.ylabel("Number of Sequences")
    plt.hist([nonsyn_1, nonsyn_2, nonsyn_pre],stacked=True,  bins=range(int(args.length)+1), label=['first duplication', 'second duplication', 'preduplication'])
    plt.legend(loc='upper left')
    plt.savefig(args.tsv +"nonsynonymous_duplication.png")
    
    fig_syn = plt.figure("synonymous duplication")
    plt.title('Synonymous Mutations')
    plt.xlabel("Location within the duplication")
    plt.ylabel("Number of Sequences")
    plt.hist([syn_1, syn_2, syn_pre],stacked=True,  bins=range(int(args.length)+1), label=['first duplication', 'second duplication', 'preduplication'])
    plt.legend(loc='upper left')
    plt.savefig(args.tsv +"synonymous_duplication.png")

    #KS STATISTIC

    #synone, syn_1 and syn_2 are lists (before normalisation by division by tree branch length and length of duplication)
    #they each show the number of times a mutation occurs at a particular location, e.g. 5, 5, 5, 3, 3, 3, 3, 6, 1, 12, ... etc
    # the part of the workflow below calculates the Kolmogorov-Smirnov Statistics for them

    dictionary_ = {'Preduplication statistic': [" ", stats.ks_2samp(syn_pre, syn_1),  stats.ks_2samp(syn_pre, syn_2)], 'PostDuplication Copy 1 statistic': [stats.ks_2samp(syn_pre, syn_1), " ", stats.ks_2samp(syn_1, syn_2)], "PostDuplication Copy 2 statistic": [stats.ks_2samp(syn_pre, syn_2), stats.ks_2samp(syn_1, syn_2),  " "]}
    df = pd.DataFrame(dictionary_)
    df.index = ['Preduplication', 'PostDuplication Copy 1', "PostDuplication Copy 2"]
    df.to_csv(args.tsv + "nonsyn.tsv", sep='\t')

    dictionary_ = {'Preduplication statistic': [" ", stats.ks_2samp(nonsyn_pre, nonsyn_1),  stats.ks_2samp(nonsyn_pre, nonsyn_2)], 'PostDuplication Copy 1 statistic': [stats.ks_2samp(nonsyn_pre, nonsyn_1), " ", stats.ks_2samp(nonsyn_1, nonsyn_2)], "PostDuplication Copy 2 statistic": [stats.ks_2samp(nonsyn_pre, nonsyn_2), stats.ks_2samp(nonsyn_1, nonsyn_2),  " "]}
    df = pd.DataFrame(dictionary_)
    df.index = ['Preduplication', 'PostDuplication Copy 1', "PostDuplication Copy 2"]
    df.to_csv(args.tsv+"_syn.tsv", sep='\t')

    #POISSON DISTRIBUTION
    
    number_of_muts_syn = [len(syn_pre), len(syn_1), len(syn_2)] #number of synonymous mutations in copy of the duplication
    number_of_muts_nonsyn = [len(nonsyn_pre), len(nonsyn_1), len(nonsyn_2)] #number of nonsynonymous mutations in each copy of the duplication

    """
    print("Mutation rate mu synonymous (pre, post 1, post 2):")
    for nr_of_muts in number_of_muts_syn:
        distr =dict()
        for i in np.arange(0.001, 1.001, 0.001):
            poisson = (((i*without_dupl)**(nr_of_muts))*math.exp(-i*without_dupl))/nr_of_muts #(factorial)
            distr[poisson]=i
        print(max(distr, key=distr.get))
    print("Mutation rate mu nonsynonymous (pre, post 1, post 2):")
    for nr_of_muts in number_of_muts_nonsyn:
        distr =dict()
        for i in np.arange(0.001, 1.001, 0.001):
            poisson = (((i*without_dupl)**(nr_of_muts))*math.exp(-i*without_dupl))/nr_of_muts
            distr[poisson]=i
        print(max(distr, key=distr.get))
        print(nr_of_muts/without_dupl)
    """
    
    mu_pre_syn = len(syn_pre)/without_dupl
    mu_pre_syn_error = math.sqrt(len(syn_pre))/without_dupl

    mu_pre_nonsyn = len(nonsyn_pre)/without_dupl
    mu_pre_nonsyn_error = math.sqrt(len(nonsyn_pre))/without_dupl

    mu_1_syn = len(syn_1)/with_dupl
    mu_1_syn_error = math.sqrt(len(syn_1))/with_dupl

    mu_1_nonsyn = len(nonsyn_1)/with_dupl
    mu_1_nonsyn_error = math.sqrt(len(nonsyn_1))/with_dupl

    mu_2_syn = len(syn_2)/with_dupl
    mu_2_syn_error = math.sqrt(len(syn_2))/with_dupl

    mu_nonsyn_2 = len(nonsyn_2)/with_dupl
    mu_2_nonsyn_error =  math.sqrt(len(nonsyn_2))/with_dupl



    poisson_dataframe = pd.DataFrame({"synonymous preduplication": mu_pre_syn, "synonymous preduplication error": mu_pre_syn_error, "synonymous posduplication 1": mu_1_syn, 
                                      "synonymous postduplication 1 error": mu_1_syn_error,
                                      "synonymous posduplication 2": mu_2_syn, "synonymous postduplication 2 error" : mu_2_syn_error,  "nonsynonymous preduplication": mu_pre_nonsyn,
                                      "nonsynonymous preduplication error" : mu_pre_nonsyn_error, 
                                        "nonsynonymous posduplication 1": mu_1_nonsyn,  "nonsynonymous postduplication 1 error": mu_1_nonsyn_error,
                                        "nonsynonymous posduplication 2": mu_nonsyn_2, "nonsynonymous postduplication 2 error": mu_2_nonsyn_error}, index=[0])
    
    poisson_dataframe.to_csv(args.tsv + "_mutation rate.tsv", sep='\t')
