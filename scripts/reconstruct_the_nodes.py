import argparse, re
from Bio import SeqIO, SeqRecord, Seq, Phylo
from collections import Counter, defaultdict
import pandas as pd

def reconstruct_insertions_(fasta, tree, output):
    """Since insertions are not included in the branches, as they are not included in the input annotated tree,
    It is necessary to reconstruct them from the terminal branches.
    
    Input: fasta file of the reconstructed aligned branches 

    
    """
    tree_ = Phylo.read(tree, "newick")
    f = SeqIO.parse(fasta, "fasta")
    tree_.root_at_midpoint()
    tree_.find_clades()

    mut_dict = defaultdict(list)
    seq_dict = dict()
    all_branch_seq = []

    for record in f:
        seq_dict[record.id] = record.seq
    #writing all the records to their seq in dictionary

    for branch in tree_.get_nonterminals(order='postorder'):
        if pd.isna(branch.name) == False:
            #checking if the branch has a duplication or not
            if '-'*int(int(args.length)/2) in seq_dict[branch.name]:
                substring = '-'*int(int(args.length)/2)
                location = seq_dict[branch.name].find(substring)
                all_branch_seq = [] 
                for b in branch:
                    if '-' not in seq_dict[b.name]:
                        #if there are no gaps in the branch, the relevant sequence is added to the list
                        sequence_in_branch = str(seq_dict[b.name][location:location+int(int(args.length)/2)])
                        all_branch_seq.append(sequence_in_branch)

            if len(all_branch_seq) >= 3:
                common_str = ""
                mutations = []
                for i in range(0, len(all_branch_seq[0])):
                    nuc_at_pos = ""
                    for a in range(0, len(all_branch_seq)): nuc_at_pos+= f'{all_branch_seq[a][i]}'
                    count = Counter(nuc_at_pos)

                    if len(count)>1:
                        # if there is a mutation, most commo character is moved recursively up the tree.
                        most_common_char = Counter(nuc_at_pos).most_common(1)[0][0]
                        common_str+=most_common_char
                        mutations.append(f'{most_common_char}{nuc_at_pos.replace(most_common_char, "")}{location+i}')
                    else:
                        mutations.append(f'{Counter(nuc_at_pos).most_common(1)[0][0]}{Counter(nuc_at_pos).most_common(1)[0][0]}{Counter(nuc_at_pos).most_common(1)[0][0]}{location+i}') 
                        common_str+= nuc_at_pos[0]

                seq_dict[branch.name] =seq_dict[branch.name].replace(substring, common_str)
                mut_dict[branch.name].extend(mutations)
                        
            if len(all_branch_seq) == 2:
                #if branches are of equal length, X is added to the duplication
                common_str = ""
                mutations = []
                for i in range(0, len(all_branch_seq[0])):
                    
                    if all_branch_seq[0][i] == all_branch_seq[1][i]:
                        common_str+=all_branch_seq[0][i]
                        mutations.append(f'{all_branch_seq[0][i]}{all_branch_seq[1][i]}{location+i}')
                    else:
                        common_str+='X'
                        mutations.append(f'{all_branch_seq[0][i]}{all_branch_seq[1][i]}{location+i}')

                seq_dict[branch.name]= seq_dict[branch.name].replace(substring, common_str)
                mut_dict[branch.name].extend(mutations)

    intermediate_file = []
    for i, j in seq_dict.items():
        entry = SeqRecord.SeqRecord(Seq.Seq(j), id=i)
        intermediate_file.append(entry)
    SeqIO.write(intermediate_file, output, "fasta")
    return(mut_dict)

def checkList(lst):
    """this function checks if all elements of a list are equal"""
    if len(lst) < 0:
        val = True
    val = all(ele == lst[0] for ele in lst)
    if(val): return('equal')

def most_frequent(string):
    """finding most frequent character in a string. If most common character is not a nucleotide, but an X,
    the second most common character is returned"""
    c = Counter(string)
    most_frequent_ = c.most_common(1)[0]
    if most_frequent_[0] not in 'ACGT':
        most_frequent_ = c.most_common(2)[0]
    return(most_frequent_)

def find_digit(string):
    num = re.findall(r'\d+', string) 
    return(num)

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="reconstruct branches from root",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-alignment', required=True, type=str, help="input root sequence")
    parser.add_argument('--intermediate', type=str,  help="length of duplication file")
    parser.add_argument('--length', type=str,  help="length of duplication file")
    parser.add_argument('--input-tree', type=str,  help="input newick tree")
    parser.add_argument('--output', type=str, help="output fasta")
    args = parser.parse_args()


    unprocessed_mutations = reconstruct_insertions_(args.input_alignment, args.input_tree, args.intermediate)
    old_file = SeqIO.parse(args.intermediate, "fasta")
    tree_ = Phylo.read(args.input_tree, "newick")
    tree_.root_at_midpoint()
    tree_.find_clades()
    mutations_for_each_node_same = defaultdict(list)

    for i in range(0, int(args.length)):
        dictionary_for_loc = dict()
        for node, entry in unprocessed_mutations.items():
            for e in entry:
                if int(find_digit(e)[0]) == i:
                    dictionary_for_loc[node] = e.replace(find_digit(e)[0], "")

        for branch in tree_.get_nonterminals(order='postorder'):
            if branch.name in dictionary_for_loc:
                for b in branch:
                    if b.name in dictionary_for_loc:
                        old = dictionary_for_loc[branch.name].replace('X', "")
                        n = set(dictionary_for_loc[branch.name]).intersection(dictionary_for_loc[b.name])
                        new = "".join(set(dictionary_for_loc[branch.name]).intersection(dictionary_for_loc[b.name]))+old
                        if n != set() and 'X' in dictionary_for_loc[branch.name]:
                            dictionary_for_loc[branch.name] = new

                        elif n== set() and 'X' in dictionary_for_loc[branch.name]:
                            n_ = "".join(set(dictionary_for_loc[b.name]).difference(dictionary_for_loc[branch.name]))
                            new = n_+ old
                            dictionary_for_loc[branch.name] = new


        for branch in tree_.get_nonterminals(order='postorder'):
                if branch.name in dictionary_for_loc:
                    if len(Counter(dictionary_for_loc[branch.name])) != 1 and branch.is_preterminal() == False :
                        for b in branch:
                                if b.is_terminal():
                                    for b_ in branch: 
                                        if b_.is_terminal() == False and len(Counter(dictionary_for_loc[branch.name])) == 1:
                                            most_common = dictionary_for_loc[b_.name]
                                            dictionary_for_loc[branch.name] = most_common

                    if len(Counter(dictionary_for_loc[branch.name])) != 1 and branch.is_preterminal():
                        parent_= tree_.get_path(branch)[-2].name
                        if parent_ in dictionary_for_loc:
                            parent_seq = (dictionary_for_loc[tree_.get_path(branch)[-2].name])
                            dictionary_for_loc[branch.name] = parent_seq.replace('X', '')
                
        last_one = dict()
        same_ = dict()
        for node, j in dictionary_for_loc.items():
            if len(Counter(dictionary_for_loc[node])) != 1:
                last_one[node] = j
            if len(Counter(dictionary_for_loc[node])) == 1:
                mutations_for_each_node_same[node].append(f'{i}:{j}')

    update_file = []
    for entry in old_file:
        if entry.id in mutations_for_each_node_same:
            for mut_ in mutations_for_each_node_same[entry.id]:
                mut_loc = int(find_digit(mut_)[0])
                if (entry.seq[mut_loc]) == 'X':
                    entry.seq = entry.seq[:mut_loc] + mut_[-1] + entry.seq[mut_loc+1:]
            new_entry = SeqRecord.SeqRecord(Seq.Seq(entry.seq), id=entry.id, description=entry.description)
            update_file.append(new_entry)
        else: 
            update_file.append(entry)
    SeqIO.write(update_file, args.output, "fasta")
