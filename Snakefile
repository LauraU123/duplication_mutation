#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene

A_OR_B = ["a", "b"]

rule all:
    input:
        graphs = expand("{a_or_b}/cumulative_sum_synonymous.png", "{a_or_b}/cumulative_sum_nonsynonymous.png", a_or_b=A_OR_B)

rule branch_from_root:
    input:
        root = ""
        sequences = ""
        tree = ""
        tree_json = ""
    output:
        reconstructed_seq = "{a_or_b}_reconstructed_sequences.fasta"
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
        aligned = 
    shell:
        """
        augur align


        
        """