#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene
#graphs = expand("{a_or_b}/cumulative_sum_synonymous.png", "{a_or_b}/cumulative_sum_nonsynonymous.png", a_or_b=A_OR_B)

A_OR_B = ["a", "b"]

rule all:
    input:
        expand("results/{a_or_b}/pairwise_G.fasta", a_or_b=A_OR_B)

rule branch_from_root:
    input:
        root = "data/rsv_{a_or_b}_root_sequence.json",
        sequences = "data/{a_or_b}_sequences.fasta",
        tree = "data/{a_or_b}_tree.nwk",
        treejson = "data/rsv_{a_or_b}_genome.json"
    output:
        reconstructed_seq = "results/{a_or_b}/reconstructed_sequences.fasta"
    shell:
        """
        python3 scripts/reconstruct_from_root.py \
        --input-root {input.root} \
        --input-tree-json {input.treejson} \
        --sequences {input.sequences} \
        --input-tree {input.tree} \
        --output {output.reconstructed_seq} 
        """


rule align_to_ref:
    input:
        reconstructed_seq = rules.branch_from_root.output.reconstructed_seq,
        reference = "config/{a_or_b}reference.fasta"
    output:
        aligned = "results/{a_or_b}/pairwise_G.fasta"
    shell:
        """
        nextalign run -j 4 \
        --reference {input.reference} \
        --output-fasta {output.aligned} \
          {input.reconstructed_seq}
        """

rule just_duplication:
    input:
        reconstructed_seq = rules.align_to_ref.aligned
    params:
        start = config["just_dupl"]["start"][{a_or_b}],
        end = config["just_dupl"]["end"][{a_or_b}]
    output:
        only_dupl = "results/{a_or_b}/only_duplication.fasta"
    shell:
        """
        python3 scripts/duplication.py \
        --input {input.reconstructed_seq} \
        --start {params.start} \
        --end {params.end} \
        --output {output.only_dupl}
        """