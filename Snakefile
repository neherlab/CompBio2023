reference = "reference/reference_seq.fasta"

import glob, os
all_dataset = [os.path.basename(fname).split('.')[0] for fname in glob.glob('data/*fastq.gz')]

rule all:
    input:
        expand("results/{sample}/consensus.fasta", sample=all_dataset)


rule fetch_primers:
    params:
        primerfile = "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv"
    output:
        "primers.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO, Seq
        raw_primers = pd.read_csv(params.primerfile, sep='\t')
        ref = str(SeqIO.read(reference, 'fasta').seq)
        primers = {}
        for r, row in raw_primers.iterrows():
            start = ref.find(row.seq)
            if start<0:
                start = ref.find(Seq.reverse_complement(row.seq))

            if start>0:
                primers[row.name] = {"segment":"MN908947", "name":row["name"], "seq":row.seq,
                                     "start":start, "end":start+len(row.seq)}
            else:
                print(f"row {row} failed")

        pd.DataFrame(primers).T.to_csv(output[0], sep='\t', index=False)

rule bwa_index:
    input:
        reference
    output:
        reference + '.bwt'
    shell:
        "bwa index {input}"

rule trim_single_read:
    input:
        r1 = "data/{sample}.fastq.gz",
    output:
        r1 = "results/{sample}/{sample}_trimmed.fq.gz"
    params:
        outdir = "results/{sample}",
        min_length = 30
    shell:
        """
        trim_galore --length {params.min_length} --output {params.outdir} {input.r1}
        """

rule map:
    input:
        ref = reference,
        index = rules.bwa_index.output,
        reads = "results/{sample}/{sample}_trimmed.fq.gz",
    output:
        "results/{sample}/reads.bam",
    shell:
        "bwa mem {input.ref} {input.reads} | samtools view -Sb | samtools sort - > {output}"

rule pileup:
    input:
        reads = "results/{sample}/reads.bam",
        primers = "primers.tsv"
    output:
        counts_file = "results/{sample}/counts.tsv",
        insertions_file = "results/{sample}/insertions.json",
    shell:
        """
        python3 scripts/create_allele_counts.py --bam_file {input.reads}
                            --primers {input.primers}\
                            --counts-file {output.counts_file} \
                            --insertions-file {output.insertions_file}
        """


rule consensus_sequence:
    input:
        counts = "results/{sample}/counts.tsv"
    params:
        min_coverage = 10,
        seq_name = "{sample}_consensus"
    output:
        "results/{sample}/consensus.fasta"
    shell:
        """
        python3 scripts/consensus_sequence.py --counts {input.counts} \
                                    --min-coverage {params.min_coverage}\
                                    --seq-name {params.seq_name} --output {output}
        """

