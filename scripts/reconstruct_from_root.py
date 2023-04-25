import json, argparse
from Bio import SeqIO, Phylo, SeqRecord, Seq



def new_recursive(node, list_=None, dictionary_=None):

    """ this function does what the previous recursive function did, but all in one go.
     Thus, it returns a dictionary with all the nodes and names with the relevant mutations attached"""

    if list_ is None:
        list_ = []
    if dictionary_ is None:
        dictionary_ = dict()

    if 'mutations' in node['branch_attrs']:
        if 'nuc' in node['branch_attrs']['mutations']:
            list_.extend(node['branch_attrs']['mutations']['nuc'])
            
    if 'name' in node:
            dictionary_[node['name']] = list(list_)

    if 'children' in node:
        for child in node['children']:
            
           new_recursive(child, list(list_), dictionary_)
           
    return(dictionary_)


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="reconstruct branches from root",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-root', required=True, type=str, help="input root sequence")
    parser.add_argument('--input-tree-json', type=str, help="input tree json")
    parser.add_argument('--sequences', type=str,  help="input sequences from tree fasta")
    parser.add_argument('--input-tree', type=str,  help="input newick tree")
    parser.add_argument('--output', type=str, help="output fasta")
    args = parser.parse_args()

    sequences_ = SeqIO.parse(args.sequences, "fasta")
    seq_dict_ = {rec.id : rec.seq for rec in sequences_}
    tree = Phylo.read(args.input_tree, "newick")
    tree.root_at_midpoint()
    tree.find_clades()

    with open (args.input_tree_json) as file_a:
        f_a = json.load(file_a)  
        
    with open(args.input_root) as root_a:
        r_a = json.load(root_a)
        root_sequence = list(r_a['nuc'])

    dictionary_mutations = new_recursive(f_a['tree'])
    all_sequences= dict()
    all_entries = []

    for branch in tree.get_nonterminals(order='postorder'):
        for b in branch:
            root_sequence = list(r_a['nuc'])
            mutations = dictionary_mutations[b.name]
            for mut in mutations:
                root_sequence[int(mut[1:-1])-1] = mut[-1]

            all_sequences[b.name] = "".join(root_sequence) 

    for id_, sequence in all_sequences.items():
        if id_ in seq_dict_:
            entry = SeqRecord.SeqRecord(Seq.Seq(seq_dict_[id_]), id=id_)
        else:
            entry = SeqRecord.SeqRecord(Seq.Seq(all_sequences[id_]), id=id_)
        all_entries.append(entry)

    SeqIO.write(all_entries, args.output, "fasta")