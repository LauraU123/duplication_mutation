#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene


rule all:
    input:
        graphs = expand("{a_or_b}/cumulative_sum_synonymous.png", "{a_or_b}/cumulative_sum_nonsynonymous.png")

rule branch_from_root:
    input:
        root = ""
        sequences = ""
        tree = ""
    output:
        reconstructed_seq = "reconstructed_sequences.fasta"