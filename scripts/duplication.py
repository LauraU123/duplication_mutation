from Bio import SeqIO, SeqRecord, Seq
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="reconstruct branches from root",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--output', required=True, type=str, help="output fasta file")
    parser.add_argument('--input', type=str, help="input fasta file")
    parser.add_argument('--start', type=str,  help="input sequences from tree fasta")
    parser.add_argument('--end', type=str,  help="input newick tree")
    args = parser.parse_args()

    alignment_to_cut = SeqIO.parse(args.input, "fasta")
    cut_ids = []
    for entry in alignment_to_cut:
        new_entry = SeqRecord.SeqRecord(Seq.Seq(entry.seq[args.start:args.end]), id=entry.id, description=entry.description)
        cut_ids.append(new_entry)
    SeqIO.write(cut_ids, args.output,"fasta")