#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene
#graphs = expand("{a_or_b}/cumulative_sum_synonymous.png", "{a_or_b}/cumulative_sum_nonsynonymous.png", a_or_b=A_OR_B)

A_OR_B = ["a", "b"]

rule all:
    input:
        "results/{a_or_b}/reconstructed_sequences.fasta"

rule branch_from_root:
    input:
        root = "data/rsv_{a_or_b}_root-sequence.json"
        sequences = "data/{a_or_b}_sequences.fasta"
        tree = "data/{a_or_b}_tree.nwk"
        tree_json = "data/rsv_{a_or_b}_genome.json"
    output:
        reconstructed_seq = "results/{a_or_b}/reconstructed_sequences.fasta"
    shell:
        """
        python3 scripts/reconstruct_from_root.py
        --input_root {input.root} \
        --input_tree_json {input.tree_json}
        --sequences {input.sequences} \
        --tree {input.tree} \
        --output {output.reconstructed_seq} \

        """

rule align_to_ref:
    input:
        reconstructed_seq = rules.branch_from_root.output.reconstructed_seq
        reference = "config/{a_or_b}reference.fasta"
    output:
        aligned = "results/{a_or_b}/pairwise_G.fasta"
    shell:
        """
        nextalign run -j 4 
        --reference /home/laura/code/without_G/config/areference.fasta
        --output-fasta {output.aligned} \
          {input.reconstructed_seq}

        """