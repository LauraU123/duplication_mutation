from Bio import SeqIO, SeqRecord, Seq
import argparse


alignment_to_cut = SeqIO.parse("/home/laura/code/without_G/G_analysis/pairwise_G_a.fasta", "fasta")
cut_ids = []
for entry in alignment_to_cut:
    new_entry = SeqRecord.SeqRecord(Seq.Seq(entry.seq[5465:5609]), id=entry.id, description=entry.description)
    cut_ids.append(new_entry)
SeqIO.write(cut_ids, 'duplication_a.fasta',"fasta")