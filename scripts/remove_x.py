import argparse
from Bio import SeqIO, Phylo, SeqRecord, Seq
from collections import Counter
import pandas as pd



def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]




if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="reconstruct branches from root",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--output', required=True, type=str, help="output fasta file")
    parser.add_argument('--tree', required=True, type=str, help="input newick tree")
    parser.add_argument('--input', type=str, help="input fasta file")
    args = parser.parse_args()



    tree_ = Phylo.read(args.tree, "newick")
    f = SeqIO.parse(args.input, "fasta")
    seq_dict = dict()
    update_file = []
    for record in f:
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
                        new_entry = SeqRecord.SeqRecord(Seq.Seq(seq_dict[branch.name]), id=branch.name, description=branch.name)
                    all_term = dict()
                    if len(branch_names) > 1:
                        for i in range(0, len(branch_names)):
                            all_term[branch_names[i].name]= len(branch_names[i].get_terminals())
                        for ind in indices_of_interest:
                            seq_dict[branch.name] = seq_dict[branch.name][:ind] + seq_dict[max(all_term)][ind]+ seq_dict[branch.name][ind+1:]
                        new_entry = SeqRecord.SeqRecord(Seq.Seq(seq_dict[branch.name]), id=branch.name, description=branch.name)
                    update_file.append(new_entry)
    SeqIO.write(update_file, args.output, 'fasta')