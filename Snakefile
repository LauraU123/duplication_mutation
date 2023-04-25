#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene
#graphs = expand("{a_or_b}/cumulative_sum_synonymous.png", "{a_or_b}/cumulative_sum_nonsynonymous.png", a_or_b=A_OR_B)
configfile: "config/configfile.yaml"
A_OR_B = ["a"]


rule all:
    input:
        expand("results/{a_or_b}/reconstructed.fasta", a_or_b=A_OR_B)

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
        reconstructed_seq = rules.align_to_ref.output.aligned
    params:
        start = lambda w: config["just_dupl"]["start"].get(w.a_or_b),
        end = lambda w: config["just_dupl"]["end"].get(w.a_or_b)
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

rule align_G:
    input:
        only_dupl = rules.just_duplication.output.only_dupl,
        reference = rules.align_to_ref.input.reference
    output:
        aligned_dupl = "results/{a_or_b}/aligned_duplication.fasta"
    shell:
        """
        augur align --nthreads 4 \
        --sequences {input.only_dupl} \
        --output {output.aligned_dupl}
        """

rule reconstruct:
    input:
        aligned_seq = rules.align_G.output.aligned_dupl,
        tree = rules.branch_from_root.input.tree,
    output:
        reconstructed_fasta = "results/{a_or_b}/reconstructed.fasta"
    params:
        len_ = lambda w: config["reconstruct"].get(w.a_or_b),
        intermediate = "results/{a_or_b}/intermediate.fasta"
    shell:
        """
        python3 scripts/reconstruct_the_nodes.py \
        --length {params.len_} \
        --intermediate {params.intermediate} \
        --output {output.reconstructed_fasta} \
        --input-tree {input.tree} \
        --input-alignment {input.aligned_seq}
        """

rule find_unknowns:
    input:
        reconstructed_fasta =rules.reconstruct.output.reconstructed_fasta
    output:
        reconstructed_fasta = "results/{a_or_b}/last_reconstruction.fasta"
    shell:
        """
        python3 remove_x.py \
        --input {input.reconstructed_fasta} \
        --output {output.reconstructed_fasta}
        """
