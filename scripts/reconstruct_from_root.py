import json
from Bio import SeqIO, Phylo, SeqRecord
from Bio.Seq import Seq


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





with open ('/home/laura/code/without_G/auspice/rsv_a_genome.json') as file_a:
    f_a = json.load(file_a)  
dictionary_mutations = new_recursive(f_a['tree'])

sequences_ = SeqIO.parse("/home/laura/code/without_G/data/a/sequences.fasta", "fasta")
seq_dict_ = {rec.id : rec.seq for rec in sequences_}
tree = Phylo.read("/home/laura/code/without_G/results/a/genome/tree.nwk", "newick")
tree.root_at_midpoint()
tree.find_clades()

with open ('/home/laura/code/without_G/auspice/rsv_a_genome.json') as file_a:
    f_a = json.load(file_a)  
    
with open('/home/laura/code/without_G/auspice/rsv_a_genome_root-sequence.json') as root_a:
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
        entry = SeqRecord(Seq(seq_dict_[id_]), id=id_)
    else:
        entry = SeqRecord(Seq(all_sequences[id_]), id=id_)
    all_entries.append(entry)

SeqIO.write(all_entries, "allbranches_and_seq_a.fasta", "fasta")