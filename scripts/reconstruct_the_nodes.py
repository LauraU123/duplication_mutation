import argparse, re
from Bio import SeqIO, SeqRecord, Seq, Phylo
from collections import Counter, defaultdict
import pandas as pd


def reconstruct_insertions_(fasta, tree, output):
    tree_ = Phylo.read(tree, "newick")
    f = SeqIO.parse(fasta, "fasta")
    tree_.root_at_midpoint()
    tree_.find_clades()
    mut_dict = defaultdict(list)
    seq_dict = dict()

    for record in f:
        seq_dict[record.id] = record.seq

    for branch in tree_.get_nonterminals(order='postorder'):
        if pd.isna(branch.name) == False:
            print(branch)
            if '-'*72 in seq_dict[branch.name]:
                substring = '-'*72
                location = (seq_dict[branch.name].find(substring))
                all_branch_seq = [] 
                for b in branch:
                    if '-' not in seq_dict[b.name]:
                        all_branch_seq.append(str(seq_dict[b.name][location:location+72]))

            if len(all_branch_seq) >= 3:
                common_str = ""
                mutations = []
                for i in range(0, len(all_branch_seq[0])):
                    nuc_at_pos = ""
                    for a in range(0, len(all_branch_seq)): 
                        nuc_at_pos+= f'{all_branch_seq[a][i]}'
                    count = Counter(nuc_at_pos)

                    if len(count)>1:
                        most_common_char = Counter(nuc_at_pos).most_common(1)[0][0]
                        common_str+=most_common_char
                        mutations.append(f'{most_common_char}{nuc_at_pos.replace(most_common_char, "")}{location+i}')
                    else:
                        mutations.append(f'{Counter(nuc_at_pos).most_common(1)[0][0]}{Counter(nuc_at_pos).most_common(1)[0][0]}{Counter(nuc_at_pos).most_common(1)[0][0]}{location+i}') 
                        common_str+= nuc_at_pos[0]
                seq_dict[branch.name] =seq_dict[branch.name].replace(substring, common_str)
                mut_dict[branch.name].extend(mutations)
                        
            if len(all_branch_seq) == 2:
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

    new_file = []
    for i, j in seq_dict.items():
        entry = SeqRecord.SeqRecord(Seq.Seq(j), id=i)
        new_file.append(entry)
    SeqIO.write(new_file, output, "fasta")
    return(mut_dict)

def checkList(lst):
    if len(lst) < 0:
        res = True
    res = all(ele == lst[0] for ele in lst)
    if(res): return('equal')

def most_frequent(string):
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
    parser.add_argument('--length', type=str,  help="length of duplication file")
    parser.add_argument('--input-tree', type=str,  help="input newick tree")
    parser.add_argument('--output', type=str, help="output fasta")
    args = parser.parse_args()


    unprocessed_mutations = reconstruct_insertions_(args.input_alignment, args.input_tree, "results/{a_or_b}/intermediate_file.fasta")
    tree_ = Phylo.read(args.input_tree, "newick")
    tree_.root_at_midpoint()
    tree_.find_clades()

    mutations_for_each_node_same = defaultdict(list)

    for i in range(0, 144):
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
        for node, j in dictionary_for_loc.items():
            if len(Counter(dictionary_for_loc[node])) == 1:
                mutations_for_each_node_same[node].append(f'{i}:{j}')
            else: last_one[node] = j

    old_file = SeqIO.parse("results/{a_or_b}/intermediate_file.fasta", "fasta")
    update_file_1 = []

    for entry in old_file:
        if entry.id in mutations_for_each_node_same:
            for mut_ in mutations_for_each_node_same[entry.id]:
                mut_loc = int(find_digit(mut_)[0])
                if (entry.seq[mut_loc]) == 'X':
                    entry.seq = entry.seq[:mut_loc] + mut_[-1] + entry.seq[mut_loc+1:]
            new_entry = SeqRecord(Seq(entry.seq), id=entry.id, description=entry.description)
            update_file_1.append(new_entry)
        else: 
            update_file_1.append(entry)

    seq_dict = dict()
    update_file = []
    for record in update_file_1:
        seq_dict[record.id] = record.seq
        if "X" not in record.seq:
            update_file.append(record)
    for branch in tree_.get_nonterminals(order='postorder'):
        all_branch_seq = []
        if pd.isna(branch.name) == False:
            if branch.name in seq_dict:
                if 'X' in seq_dict[branch.name]:
                    indices_of_interest = (find(seq_dict[branch.name], 'X'))
                    branch_names = []
                    for b in branch:
                        if b.is_terminal() == False:
                            branch_names.append(b)
                    if len(branch_names) == 1:
                        for ind in indices_of_interest:
                            if (seq_dict[branch_names[0].name][ind]) != 'X':
                                seq_dict[branch.name] = seq_dict[branch.name][:ind] + seq_dict[branch_names[0].name][ind]+ seq_dict[branch.name][ind+1:]
                            else:
                                which_one = [seq_dict[branch_names[0].get_terminals()[0].name][ind], seq_dict[branch_names[0].get_terminals()[1].name][ind]]
                                for b in branch:
                                    if b.is_terminal() == True:
                                        which_one.append(seq_dict[b.name][ind])
                                new_letter = (max(Counter(which_one).keys()))
                                seq_dict[branch.name] = seq_dict[branch.name][:ind] + new_letter + seq_dict[branch.name][ind+1:]

                        new_entry = SeqRecord(Seq(seq_dict[branch.name]), id=branch.name, description=branch.name)

                    all_term = dict()
                    if len(branch_names) > 1:
                        for i in range(0, len(branch_names)):
                            all_term[branch_names[i].name]= len(branch_names[i].get_terminals())
                        for ind in indices_of_interest:
                            seq_dict[branch.name] = seq_dict[branch.name][:ind] + seq_dict[max(all_term)][ind]+ seq_dict[branch.name][ind+1:]
                        new_entry = SeqRecord(Seq(seq_dict[branch.name]), id=branch.name, description=branch.name)
                    update_file.append(new_entry)

    SeqIO.write(update_file, args.output, 'fasta')