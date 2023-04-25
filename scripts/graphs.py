import argparse
from Bio import SeqIO, Phylo, Seq
from collections import defaultdict
import matplotlib.pyplot as plt
import re


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="reconstruct branches from root",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--output', required=True, type=str, help="output fasta file")
    parser.add_argument('--tree', required=True, type=str, help="input newick tree")
    parser.add_argument('--input', type=str, help="input fasta file")
    args = parser.parse_args()

    tree_ = Phylo.read(args.tree, "newick")
    tree_.root_at_midpoint()
    tree_.find_clades()

    synonymous_one, nonsynonymous_one, onlyone = (defaultdict(list) for i in range(3))
    lst, syn_one, lst_nonsyn, nonsyn_one  = ([] for i in range(4))

    just_the_duplication = SeqIO.parse(args.input, 'fasta')

    for entry in just_the_duplication:
        if '-' in entry.seq:
            only_one = entry.seq.replace("-", "")
            if len(only_one)== int(args.length):
                for i in range(0, len(entry.seq[1:-2]), 3):
                    onlyone[entry.id].append(only_one[i:i+3])

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
        sort_branch = []
        if branch.name in synonymous_one:
            for b in branch:
                if b.name in synonymous_one:
                    synonymous_one[b.name] = list(set(synonymous_one[b.name]).difference(set(synonymous_one[branch.name])))

    for branch in tree_.get_nonterminals(order='preorder'):
        sort_branch = []
        for e in nonsynonymous_one[branch.name]:
            sort_branch.append(str("".join(sorted(e[:2], key=str.lower))+ e[2:]))

        for b in branch:
            if b.name in nonsynonymous_one:
                sort_b = []
                for entry in nonsynonymous_one[b.name]:
                    sort_b.append(str("".join(sorted(entry[:2], key=str.lower))+ entry[2:]))
                nonsynonymous_one[b.name] = set(set(sort_b).difference(set(sort_branch)))
                
    for i in synonymous_one.values():
        ls = list(i)
        for j in ls: lst.append(j)
    for item_ in lst:
        syn_one.append(int(item_[2:]))

    for i in nonsynonymous_one.values():
        ls = list(i)
        numbers_ = []
        for it in ls:
            numbers_.append(it[2:])
        new_numbers_ = list(set(numbers_))
            
        for j in new_numbers_: lst_nonsyn.append(j)
    for item_ in lst_nonsyn:
        nonsyn_one.append(int(item_))

    plt.title('Mutations in RSV-A no insertion')
    plt.xlabel("Location within the duplication")
    plt.ylabel("Number of Sequences")
    plt.hist([nonsyn_one, syn_one], bins=range(int(args.length)+1), stacked=True, label= ['nonsynonymous mut', 'synonymous mut']) 
    plt.legend(loc='upper left')
    plt.savefig("no_insertions.png")