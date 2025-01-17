import json
from sra_data import sra_data
from intervals import interval_list

SRA_IDS = [sra_id for individual in sra_data.values() for sra_id in individual]
individual_sra_pairs = [(individual, sra_id) for individual, sra_ids in sra_data.items() for sra_id in sra_ids]
CHROMOSOMES = [chrm for chrm in interval_list.values()]

rule all:
    input:
        directory(expand("data/genomicsdb/chr_{chrm}", chrm=interval_list.values()))
        # "intervals.list"
        # directory("data/genomicsdb/")
        # expand("data/{individual}/{sra_id}_mergebamalignment.bam", zip, individual=[pair[0] for pair in individual_sra_pairs], sra_id=[pair[1]
        # for pair in individual_sra_pairs])
        # expand("data/{individual}/{individual}_idxstats.txt", individual=sra_data.keys()),
        # expand("data/{individual}/{individual}_flagstat.txt", individual=sra_data.keys()),
        # expand("data/{individual}/{individual}_depth.txt", individual=sra_data.keys()),
        # expand("data/{individual}/{individual}_insert_size_metrics.txt", individual=sra_data.keys()),
        # expand("data/{individual}/{individual}_insert_size_histogram.pdf", individual=sra_data.keys()),
        # expand("data/{individual}/{individual}_alignment_summary_metrics.txt", individual=sra_data.keys()),
        # expand("data/{individual}/{individual}_alignment_summary_histogram.pdf", individual=sra_data.keys())

rule download_sra:
    output:
        "sra/{sra_id}.sra"  
    params:
        sra_id = lambda wildcards: wildcards.sra_id
    shell:
        """
        prefetch {params.sra_id} --output-file sra/{params.sra_id}.sra
        """

rule convert_to_fastq:
    input:
        "sra/{sra_id}.sra"
    output:
        "fastq/{individual}/{sra_id}_1.fastq",  
        "fastq/{individual}/{sra_id}_2.fastq"   
    shell:
        """
        mkdir -p fastq/{wildcards.individual}
        fastq-dump --split-files --outdir fastq/{wildcards.individual}/ {input}
        """

rule download_reference:
    output:
        "reference/GCF_000738735.6_ASM73873v6_genomic.fna.gz"
    shell:
        """
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.6_ASM73873v6/GCF_000738735.6_ASM73873v6_genomic.fna.gz
        """

rule unzip_reference:
    input:
        "reference/GCF_000738735.6_ASM73873v6_genomic.fna.gz"
    output:
        "reference/ref_genome.fasta"
    shell:
        """
        gunzip -c {input} > {output}
        """

rule create_dict:
    input:
        "reference/ref_genome.fasta"
    output:
        "reference/ref_genome.dict"
    shell:
        """
        gatk CreateSequenceDictionary -R {input} -O {output}
        """

rule create_fai:
    input:
        "reference/ref_genome.fasta"
    output:
        "reference/ref_genome.fasta.fai"
    shell:
        """
        samtools faidx {input}
        """
 
rule create_intervals:
    input:
        "reference/ref_genome.fasta.fai"
    output:
        "intervals.list"
    shell:
        """
        awk '{{print $1":"1"-"$2}}' {input} > {output}
        """ 
 
rule fetch_library_names:
    input:
        sra_data="sra_data.json" 
    output:
        library_names="library_names.json" 
    params:
        email="email@gmail.com" 
    run:
        import json
        from Bio import Entrez
        import xml.etree.ElementTree as ET

        def get_library_name(sra_id, email):
            Entrez.email = email 

            try:
                handle = Entrez.efetch(db="sra", id=sra_id, retmode="xml")
                sra_xml = handle.read()
                handle.close()

                root = ET.fromstring(sra_xml)

                library_name = root.find(".//LIBRARY_DESCRIPTOR/LIBRARY_NAME")
                return library_name.text if library_name is not None else "Library name not found"
            except Exception as e:
                print(f"Error fetching data for {sra_id}: {e}")
                return None

        def fetch_library_names(sra_data_file, output_file, email):
            library_names = {}

            with open(sra_data_file, 'r') as f:
                sra_data = json.load(f)

            for srx_id, sra_ids in sra_data.items():
                for sra_id in sra_ids:
                    library_name = get_library_name(sra_id, email)
                    if library_name:
                        library_names[sra_id] = library_name

            with open(output_file, 'w') as f:
                json.dump(library_names, f, indent=4)
                
        fetch_library_names(input.sra_data, output.library_names, params.email)


rule extract_rg_info:
    input:
        "fastq/{individual}/{sra_id}_1.fastq",
        library_names="library_names.json"
    output:
        "data/{individual}/{sra_id}_rg_info.json"
    run:
        import json

        def extract_read_group_info(fastq_file, library_names, sra_id):
            rg_info = None

            with open(library_names, 'r') as lib_file:
                library_data = json.load(lib_file)

            library_name = library_data.get(sra_id, "Unknown")

            with open(fastq_file, 'r') as f:
                for line in f:
                    if line.startswith('@'):
                        parts = line.strip().split(' ')

                        rgid = parts[1]
                        flow_cell_info = rgid.split(':')

                        if len(flow_cell_info) >= 3:
                            sequencer_id = flow_cell_info[0]
                            flow_cell_id = flow_cell_info[1]
                            lane_number = flow_cell_info[2]
                        else:
                            sequencer_id, flow_cell_id, lane_number = "Unknown", "Unknown", "Unknown"

                        rgid = f"{sequencer_id}.{flow_cell_id}.{lane_number}.NA"

                        sample_name = wildcards.individual  
                        platform_unit = f"{flow_cell_id}.{lane_number}"
                        platform = "illumina"

                        rg_info = {
                            "RGID": rgid,
                            "SM": sample_name,
                            "LB": library_name,
                            "PU": platform_unit,
                            "PL": platform
                        }
                        break 

            return rg_info

        rg_info = extract_read_group_info(input[0], input.library_names, wildcards.sra_id)

        with open(output[0], 'w') as out_file:
            import json
            json.dump(rg_info, out_file, indent=4)


rule convert_to_ubam:
    input:
        forward_fastq="fastq/{individual}/{sra_id}_1.fastq",
        reverse_fastq="fastq/{individual}/{sra_id}_2.fastq",
        rg_info="data/{individual}/{sra_id}_rg_info.json"
    output:
        "data/{individual}/{sra_id}_output.bam"
    resources:
        runtime = 120,
        nodes = 20
    run:
        import json

        with open(input.rg_info, 'r') as f:
            rg_info = json.load(f)

        read_group_name = rg_info['RGID']
        sample_name = rg_info['SM']
        library_name = rg_info['LB']
        platform = rg_info['PL']
        
        shell(
            f"""
            java -Xmx64G -jar picard.jar FastqToSam \
                FASTQ={input.forward_fastq} \
                FASTQ2={input.reverse_fastq} \
                OUTPUT={output} \
                READ_GROUP_NAME={read_group_name} \
                SAMPLE_NAME={sample_name} \
                LIBRARY_NAME={library_name} \
                PLATFORM={platform} 
            """
        )

rule mark_illumina_adapters:
    input:
        "data/{individual}/{sra_id}_output.bam"
    output:
        "data/{individual}/{sra_id}_output_marked.bam",
        "data/{individual}/{sra_id}_markilluminaadapters_metrics.txt"
    resources:
        runtime = 120,
        nodes = 5
    params:
        tmp_dir="temp/sra_processing" 
    shell:
        """
        java -Xmx8G -jar picard.jar MarkIlluminaAdapters \
            I={input} \
            O={output[0]} \
            M={output[1]} \
            TMP_DIR={params.tmp_dir}
        """
        
rule sam_to_fastq:
    input:
        "data/{individual}/{sra_id}_output_marked.bam"
    output:
        "data/{individual}/{sra_id}_fastq.fq"
    resources:
        runtime = 60,
        nodes = 5
    params:
        tmp_dir="temp/sra_processing"   
    shell:
        """
        java -Xmx8G -jar picard.jar SamToFastq \
            I={input} \
            FASTQ={output} \
            INTERLEAVE=true \
            CLIPPING_ATTRIBUTE=XT \
            CLIPPING_ACTION=1 \
            NON_PF=true \
            TMP_DIR={params.tmp_dir}
        """
        
rule bwa_mem:
    input:
        reference="reference/ref_genome.fasta",
        reference_dict="reference/ref_genome.dict",
        fastq="data/{individual}/{sra_id}_fastq.fq"
    output:
        "data/{individual}/{sra_id}_bwa_mem.sam"
    resources:
        runtime = 45,
        tasks = 4
    threads: 20
    shell:
         """
        bwa/bwa mem -M -t {threads} -p {input.reference} {input.fastq} > {output}
        """

rule merge_bam_alignment:
    input:
        unmapped_bam="data/{individual}/{sra_id}_output_marked.bam",
        aligned_bam="data/{individual}/{sra_id}_bwa_mem.sam",
        reference="reference/ref_genome.fasta",
        reference_dict="reference/ref_genome.dict"
    output:
        "data/{individual}/{sra_id}_mergebamalignment.bam"
        # "data/{individual}/{sra_id}_mergebamalignment.bai"
    resources:
        runtime = 270
    threads: 20
    params:
        tmp_dir="temp/sra_processing"  
    shell:
        """
        java -Xmx32G -jar picard.jar MergeBamAlignment \
            R={input.reference} \
            UNMAPPED_BAM={input.unmapped_bam} \
            ALIGNED_BAM={input.aligned_bam} \
            O={output[0]} \
            VERBOSITY=DEBUG \
            SORT_ORDER=queryname \
            # CREATE_INDEX=true \ 
            TMP_DIR={params.tmp_dir}
        """

# not used if queryname sorted during merge bam alignment; only used if the .bai files is needed and therefore merge bam alignment will need to be coordinate-sorted.
rule sortbam:
    input: 
        bam="data/{individual}/{sra_id}_mergebamalignment.bam"
    output: 
        bam="data/{individual}/{sra_id}_sortedbam.bam"
    params:
        tmp_dir=lambda wildcards: f"temp/sortbam/{wildcards.individual}" 
    resources:
        runtime = 120
    threads: 20
    shell:
        """
        mkdir -p {params.tmp_dir}
        java -Xmx16G -jar picard.jar SortSam \
          I={input.bam} \
          O={output.bam} \
          SORT_ORDER=queryname \
          TMP_DIR={params.tmp_dir}
        """
        
rule mark_duplicates_spark:
    # input:
      #  bam="data/{individual}/{individual}_merged.bam"
    input:
        lambda wildcards: expand(    
            "data/{individual}/{sra_id}_mergebamalignment.bam", 
            individual=wildcards.individual,
            sra_id=sra_data[wildcards.individual]
        )
    output:
        bam="data/{individual}/{individual}_merged_dedup.bam",
        metrics="data/{individual}/{individual}_merged_dedup_metrics.txt"
    resources:
        runtime = 200,
        nodes = 20
    params:
        tmp_dir=lambda wildcards: f"temp/spark_processing/{wildcards.individual}"  
    run:
        input_bams = " ".join([f"-I {bam}" for bam in input])
        shell(f"""
        mkdir -p {params.tmp_dir}
        rm -rf {output.bam}.parts
        gatk --java-options "-Xmx16G" MarkDuplicatesSpark \
            {input_bams} \
            -O {output.bam} \
            -M {output.metrics} \
            --tmp-dir {params.tmp_dir} 
        """)
        
rule samtools_idxstats:
    input:
        bam="data/{individual}/{individual}_merged_dedup.bam",
    output:
        idxstats="data/{individual}/{individual}_idxstats.txt",
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats}
        """

rule samtools_flagstat:
    input:
        bam="data/{individual}/{individual}_merged_dedup.bam",
    output:
        flagstat="data/{individual}/{individual}_flagstat.txt",
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_depth:
    input:
        bam="data/{individual}/{individual}_merged_dedup.bam",
    output:
        depth="data/{individual}/{individual}_depth.txt",
    shell:
        """
        samtools depth {input.bam} > {output.depth}
        """

rule gatk_insert_size_metrics:
    input:
        bam="data/{individual}/{individual}_merged_dedup.bam",
        genome="reference/ref_genome.fasta",
    output:
        metrics="data/{individual}/{individual}_insert_size_metrics.txt",
        histogram="data/{individual}/{individual}_insert_size_histogram.pdf",
    log:
        "data/{individual}/{individual}_insert_size.log",

    shell:
        """
            gatk --java-options "-Xmx16G" CollectInsertSizeMetrics \
            -R {input.genome} \
            -I {input.bam} \
            -O {output.metrics} \
            -H {output.histogram} \
            2> {log}
        """

rule gatk_alignment_summary_metrics:
    input:
        bam="data/{individual}/{individual}_merged_dedup.bam",
        genome="reference/ref_genome.fasta",
    output:
        metrics="data/{individual}/{individual}_alignment_summary_metrics.txt",
        histogram="data/{individual}/{individual}_alignment_summary_histogram.pdf",
    log:
        "data/{individual}/{individual}_alignment_summary.log",

    shell:
        """
        gatk --java-options "-Xmx16G" CollectAlignmentSummaryMetrics \
            -R {input.genome} \
            -I {input.bam} \
            -O {output.metrics} \
            -H {output.histogram} \
            2> {log}
        """

rule reheader_bam:
    input:
        bam="data/{individual}/{individual}_merged_dedup.bam"
    output:
        bam="data/{individual}/{individual}_final.bam",
        bai="data/{individual}/{individual}_final.bam.bai"
    shell:
        """
        samtools view -H {input.bam} > {output.bam}.header.sam

        sed -i 's/SM:[^ \\t]*/SM:{wildcards.individual}/g' {output.bam}.header.sam

        samtools reheader {output.bam}.header.sam {input.bam} > {output.bam}

        samtools index {output.bam}

        rm {output.bam}.header.sam
        """

rule haplotype_caller:
    input:
        bam="data/{individual}/{individual}_merged_dedup.bam",  
        reference="reference/ref_genome.fasta",
        reference_dict="reference/ref_genome.dict"
    output:
        gvcf="data/{individual}/{individual}_raw_variants.g.vcf",  
        gvcf_index="data/{individual}/{individual}_raw_variants.g.vcf.idx"
    resources:
        runtime = 840
    params:
        tmp_dir="temp/sra_processing",
        threads=80
    shell:
        """
        gatk --java-options "-Xmx16G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bam} \
            -O {output.gvcf} \
            --tmp-dir {params.tmp_dir} \
            -ERC GVCF
        """

rule split_intervals:
    input:
        intervals="intervals.list",
        reference="reference/ref_genome.fasta"
    output:
        "data/intervals/{chunk}-scattered.interval_list"
    params:
        scatter_count=100
    shell:
        """
        gatk SplitIntervals \
            -R {input.reference} \
            -L {input.intervals} \
            -O data/intervals \
            --scatter-count {params.scatter_count} 
        """

rule genomics_db_import:
    input:
        gvcfs=expand("data/{individual}/{individual}_raw_variants.g.vcf", individual=sra_data.keys()),
        gvcf_indexes=expand("data/{individual}/{individual}_raw_variants.g.vcf.idx", individual=sra_data.keys()),
    output:
        db = directory("data/genomicsdb/chr_{chrm}")
    params:
        interval=lambda wildcards: "{wildcards.chrm}",
        tmp_dir=lambda wildcards: f"temp/genomicsdb_import/{wildcards.chrm}",
        threads=80
    resources:
        runtime = 360
    threads: 80
    run:
        input_gvcfs = " ".join([f"-V {gvcf}" for gvcf in input.gvcfs])
        shell(f"""
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xmx64G" GenomicsDBImport \
            --genomicsdb-workspace-path {output.db} \
            -L {params.interval} \
            --tmp-dir {params.tmp_dir} \
            --reader-threads {params.threads} \
            {input_gvcfs}
        """)
    
rule genotype_gvcfs:
    input:
        db="data/genomicsdb/",
        reference="reference/ref_genome.fasta",
        interval_list="intervals.list"
    output:
        vcf="results/cohort_genotyped.vcf.gz",
        vcf_index="results/cohort_genotyped.vcf.gz.tbi"
    params:
        tmp_dir="temp/genotype_gvcfs"
    resources:
        runtime = 500
    shell:
        """
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xmx16G" GenotypeGVCFs \
            -R {input.reference} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            -L {input.interval_list} \
            --tmp-dir {params.tmp_dir}
        """

