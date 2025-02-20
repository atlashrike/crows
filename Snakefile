import json
from sra_data import sra_data
from outgroup_data import outgroup_data
from intervals import interval_list
from intervals2 import interval_list2
from species_samples import species_samples

SRA_IDS = [sra_id for individual in sra_data.values() for sra_id in individual]
OUTGROUP_IDS = [outgroup_id for individual in outgroup_data.values() for outgroup_id in individual]
individual_sra_pairs = [(individual, sra_id) for individual, sra_ids in sra_data.items() for sra_id in sra_ids]
individual_outgroup_pairs = [(individual, outgroup_id) for individual, outgroup_ids in outgroup_data.items() for outgroup_id in outgroup_ids]
CHROMOSOMES = [chrm for chrm in interval_list.values()]
CHROMOSOMES2 = [chrm for chrm in interval_list2.values()] 
# I had two scaffolds that had no data after filtering, so this is the chromosomes + scaffolds list without that 
SPECIES_IDS = [species_id for species in species_samples.values() for species_id in species]
species_pairs = [(species, species_id) for species, species_ids in species_samples.items() for species_id in species_ids]
PARTS = ["1","2","3","4","5","6","7","8"]

rule all:
    input:
        expand("data/haps/chr_{chrm}_haps_orig.haps",
               chrm=interval_list2.values()),
        expand("data/estsfs/output/chr_{chrm}_output-file-sfs_part{part}.txt", chrm=interval_list2.values(), part=PARTS)

rule downloads:
    input:
        expand("sra/{sra_id}.sra", sra_id=SRA_IDS),
        expand("outgroup/sra/{outgroup_id}.sra", outgroup_id=OUTGROUP_IDS)

rule download_sra:
    output:
        "sra/{sra_id}.sra"  
    params:
        sra_id = lambda wildcards: wildcards.sra_id
    shell:
        """
        prefetch {params.sra_id} --output-file sra/{params.sra_id}.sra
        """
    
rule download_outgroup:
    output:
        "outgroup/sra/{outgroup_id}.sra"  
    params:
        outgroup_id = lambda wildcards: wildcards.outgroup_id
    shell:
        """
        prefetch {params.outgroup_id} --output-file outgroup/sra/{params.outgroup_id}.sra
        """

rule converts:
    input:
        expand("fastq/{individual}/{sra_id}_1.fastq", sra_id=SRA_IDS, individual=sra_data.keys()),
        expand("fastq/{individual}/{sra_id}_2.fastq", sra_id=SRA_IDS, individual=sra_data.keys()),
        expand("outgroup/fastq/{individual}/{outgroup_id}_1.fastq", outgroup_id=OUTGROUP_IDS, individual=outgroup_data.keys()),
        expand("outgroup/fastq/{individual}/{outgroup_id}_2.fastq", outgroup_id=OUTGROUP_IDS, individual=outgroup_data.keys())

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

rule convert_to_fastq_out:
    input:
        "outgroup/sra/{outgroup_id}.sra"  
    output:
        "outgroup/fastq/{individual}/{outgroup_id}_1.fastq",  
        "outgroup/fastq/{individual}/{outgroup_id}_2.fastq"   
    shell:
        """
        mkdir -p outgroup/fastq/{wildcards.individual}
        fastq-dump --split-files --outdir outgroup/fastq/{wildcards.individual}/ {input}
        """

rule references:
    input:
        "reference/GCF_000738735.6_ASM73873v6_genomic.fna.gz",
        "reference/ref_genome.fasta",
        "reference/ref_genome.dict",
        "reference/ref_genome.fasta.fai",
        "intervals.list"
        
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

# Note: this rule was actually run seperately as a script off-cluster because the cluster permissions did not allow me to access Entrez

rule library:
    input:
        "library_names.json" 
    
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

rule rg_info:
    input:
        expand("data/{individual}/{sra_id}_rg_info.json", sra_id=SRA_IDS, individual=sra_data.keys()),
        expand("outgroup/{individual}/{outgroup_id}_rg_info.json", outgroup_id=OUTGROUP_IDS, individual=outgroup_data.keys())
        
rule extract_rg_info_in:
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
                            "PL": platform}
                        break 

            return rg_info

        rg_info = extract_read_group_info(input[0], input.library_names, wildcards.sra_id)

        with open(output[0], 'w') as out_file:
            import json
            json.dump(rg_info, out_file, indent=4)

rule extract_rg_info_out:
    input:
        "outgroup/fastq/{individual}/{outgroup_id}_1.fastq",
        library_names="outgroup_library_names.json"
    output:
        "outgroup/{individual}/{outgroup_id}_rg_info.json"
    run:
        import json

        def extract_read_group_info(fastq_file, library_names, outgroup_id):
            rg_info = None

            with open(library_names, 'r') as lib_file:
                library_data = json.load(lib_file)

            library_name = library_data.get(outgroup_id, "Unknown")

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
                            "PL": platform}
                        break 

            return rg_info

        rg_info = extract_read_group_info(input[0], input.library_names, wildcards.outgroup_id)

        with open(output[0], 'w') as out_file:
            import json
            json.dump(rg_info, out_file, indent=4)

rule ubams:
    input:
        expand("data/{individual}/{sra_id}_output.bam", sra_id=SRA_IDS, individual=sra_data.keys()),
        expand("outgroup/data/{individual}/{outgroup_id}_output.bam", outgroup_id=OUTGROUP_IDS, individual=outgroup_data.keys())
           
rule convert_to_ubam_in:
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
        
        shell(f"""
            java -Xmx64G -jar picard.jar FastqToSam \
                FASTQ={input.forward_fastq} \
                FASTQ2={input.reverse_fastq} \
                OUTPUT={output} \
                READ_GROUP_NAME={read_group_name} \
                SAMPLE_NAME={sample_name} \
                LIBRARY_NAME={library_name} \
                PLATFORM={platform} 
            """)

rule convert_to_ubam_out:
    input:
        forward_fastq="outgroup/fastq/{individual}/{outgroup_id}_1.fastq",
        reverse_fastq="outgroup/fastq/{individual}/{outgroup_id}_2.fastq",
        rg_info="outgroup/{individual}/{outgroup_id}_rg_info.json"
    output:
        "outgroup/data/{individual}/{outgroup_id}_output.bam"
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
        
        shell(f"""
            java -Xmx64G -jar picard.jar FastqToSam \
                FASTQ={input.forward_fastq} \
                FASTQ2={input.reverse_fastq} \
                OUTPUT={output} \
                READ_GROUP_NAME={read_group_name} \
                SAMPLE_NAME={sample_name} \
                LIBRARY_NAME={library_name} \
                PLATFORM={platform} 
            """)

rule mark_adapters:
    input:
        expand("data/{individual}/{sra_id}_output_marked.bam", sra_id=SRA_IDS, individual=sra_data.keys()),
        expand("outgroup/data/{individual}/{outgroup_id}_output_marked.bam", outgroup_id=OUTGROUP_IDS, individual=outgroup_data.keys())
        
rule mark_illumina_adapters_in:
    input:
        "data/{individual}/{sra_id}_output.bam"
    output:
        bam = "data/{individual}/{sra_id}_output_marked.bam",
        metrics = "data/{individual}/{sra_id}_markilluminaadapters_metrics.txt"
    resources:
        runtime = 120,
        nodes = 5
    params:
        tmp_dir="temp/sra_processing" 
    shell:
        """
        java -Xmx8G -jar picard.jar MarkIlluminaAdapters \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            TMP_DIR={params.tmp_dir}
        """

rule mark_illumina_adapters_out:
    input:
        "outgroup/data/{individual}/{outgroup_id}_output.bam"
    output:
        bam = "outgroup/data/{individual}/{outgroup_id}_output_marked.bam",
        metrics = "outgroup/data/{individual}/{outgroup_id}_markilluminaadapters_metrics.txt"
    resources:
        runtime = 120,
        nodes = 5
    params:
        tmp_dir="temp/sra_processing" 
    shell:
        """
        java -Xmx8G -jar picard.jar MarkIlluminaAdapters \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            TMP_DIR={params.tmp_dir}
        """
# the parameters used below and in mergebamalignment come from the gatk best practices tutorial 

rule sam_fastq:
    input:
        expand("data/{individual}/{sra_id}_fastq.fq", sra_id=SRA_IDS, individual=sra_data.keys()),
        expand("outgroup/data/{individual}/{outgroup_id}_fastq.fq", outgroup_id=OUTGROUP_IDS, individual=outgroup_data.keys())
        
rule sam_to_fastq_in:
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
        
rule sam_to_fastq_out:
    input:
        "outgroup/data/{individual}/{outgroup_id}_output_marked.bam"
    output:
        "outgroup/data/{individual}/{outgroup_id}_fastq.fq"
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
        
rule bwa_mems:
    input:
        expand("data/{individual}/{sra_id}_bwa_mem.sam", sra_id=SRA_IDS, individual=sra_data.keys()),
        expand("outgroup/data/{individual}/{outgroup_id}_bwa_mem.sam", outgroup_id=OUTGROUP_IDS, individual=outgroup_data.keys())
        
rule bwa_mem_in:
    input:
        reference="reference/ref_genome.fasta",
        reference_dict="reference/ref_genome.dict",
        fastq="data/{individual}/{sra_id}_fastq.fq"
    output:
        "data/{individual}/{sra_id}_bwa_mem.sam"
    resources:
        runtime = 120,
        tasks = 4
    threads: 20
    shell:
         """
        bwa/bwa mem -M -t {threads} -p {input.reference} {input.fastq} > {output}
        """

rule bwa_mem_out:
    input:
        reference="reference/ref_genome.fasta",
        reference_dict="reference/ref_genome.dict",
        fastq="outgroup/data/{individual}/{outgroup_id}_fastq.fq"
    output:
        "outgroup/data/{individual}/{outgroup_id}_bwa_mem.sam"
    resources:
        runtime = 120,
        tasks = 4
    threads: 20
    shell:
         """
        bwa/bwa mem -M -t {threads} -p {input.reference} {input.fastq} > {output}
        """

rule mergebams:
    input:
        expand("data/{individual}/{sra_id}_mergebamalignment.bam", sra_id=SRA_IDS, individual=sra_data.keys()),
        expand("outgroup/data/{individual}/{outgroup_id}_mergebamalignment.bam", outgroup_id=OUTGROUP_IDS, individual=outgroup_data.keys())

rule merge_bam_alignment_in:
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
        java -Xmx16G -jar picard.jar MergeBamAlignment \
            R={input.reference} \
            UNMAPPED_BAM={input.unmapped_bam} \
            ALIGNED_BAM={input.aligned_bam} \
            O={output[0]} \
            CREATE_INDEX=true \
            ADD_MATE_CIGAR=true \
            CLIP_ADAPTERS=false \
            CLIP_OVERLAPPING_READS=true \
            INCLUDE_SECONDARY_ALIGNMENTS=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            ATTRIBUTES_TO_RETAIN=XS \
            TMP_DIR={params.tmp_dir}
        """
        
rule merge_bam_alignment_out:
    input:
        unmapped_bam="outgroup/data/{individual}/{outgroup_id}_output_marked.bam",
        aligned_bam="outgroup/data/{individual}/{outgroup_id}_bwa_mem.sam",
        reference="reference/ref_genome.fasta",
        reference_dict="reference/ref_genome.dict"
    output:
        "outgroup/data/{individual}/{outgroup_id}_mergebamalignment.bam"
        # "outgroup/data/{individual}/{outgroup_id}_mergebamalignment.bai"
    resources:
        runtime = 270
    threads: 20
    params:
        tmp_dir="temp/sra_processing"  
    shell:
        """
        java -Xmx16G -jar picard.jar MergeBamAlignment \
            R={input.reference} \
            UNMAPPED_BAM={input.unmapped_bam} \
            ALIGNED_BAM={input.aligned_bam} \
            O={output[0]} \
            CREATE_INDEX=true \
            ADD_MATE_CIGAR=true \
            CLIP_ADAPTERS=false \
            CLIP_OVERLAPPING_READS=true \
            INCLUDE_SECONDARY_ALIGNMENTS=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            ATTRIBUTES_TO_RETAIN=XS \
            TMP_DIR={params.tmp_dir}
        """

# not used if queryname sorted during merge bam alignment; only used if the .bai files is needed and therefore merge bam alignment will need to be coordinate-sorted.
rule sortbam:
    input:
        bam="outgroup/data/{individual}/{outgroup_id}_mergebamalignment.bam"
        # bam="data/{individual}/{sra_id}_mergebamalignment.bam"
    output: 
        bam="outgroup/data/{individual}/{outgroup_id}_sortedbam.bam"
        # bam="data/{individual}/{sra_id}_sortedbam.bam"
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

rule markdup:
    input:
        expand("data/{individual}/{individual}_merged_dedup.bam", individual=sra_data.keys()),
        expand("outgroup/data/{individual}/{individual}_merged_dedup.bam",  individual=outgroup_data.keys())

rule mark_duplicates_spark_in:
    input:
        lambda wildcards: expand("data/{individual}/{sra_id}_mergebamalignment.bam", individual=wildcards.individual, sra_id=sra_data[wildcards.individual])
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

rule mark_duplicates_spark_out:
    input:
        lambda wildcards: expand("outgroup/data/{individual}/{outgroup_id}_sortedbam.bam", individual=wildcards.individual, outgroup_id=outgroup_data[wildcards.individual])
    output:
        bam="outgroup/data/{individual}/{individual}_merged_dedup.bam",
        metrics="outgroup/data/{individual}/{individual}_merged_dedup_metrics.txt"
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
        
# not used anymore; only used when I realized that I was making the header for fastqs with different SRA ids belonging to the same individual different, which caused issues in haplotype caller, which was fixed. 

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

rule haplotype:
    input:
        expand("data/{individual}/{individual}_raw_variants.g.vcf", individual=sra_data.keys()),
        expand("outgroup/data/{individual}/{individual}_raw_variants.g.vcf",  individual=outgroup_data.keys())

rule haplotype_caller_in:
    input:
        bam="data/{individual}/{individual}_merged_dedup.bam",  
        reference="reference/ref_genome.fasta",
        reference_dict="reference/ref_genome.dict"
    output:
        gvcf="data/{individual}/{individual}_raw_variants.g.vcf",  
        gvcf_index="data/{individual}/{individual}_raw_variants.g.vcf.idx"
    resources:
        runtime = 1440,
        tasks=1
    params:
        tmp_dir="temp/sra_processing"
    threads: 80
    shell:
        """
        gatk --java-options "-Xmx110G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bam} \
            -O {output.gvcf} \
            --tmp-dir {params.tmp_dir} \
            -ERC GVCF
        """

rule haplotype_caller_out:
    input:
        bam="outgroup/data/{individual}/{individual}_merged_dedup.bam",
        reference="reference/ref_genome.fasta",
        reference_dict="reference/ref_genome.dict"
    output:
        gvcf="outgroup/data/{individual}/{individual}_raw_variants.g.vcf",  
        gvcf_index="outgroup/data/{individual}/{individual}_raw_variants.g.vcf.idx"
    resources:
        runtime = 1440,
        tasks=1
    params:
        tmp_dir="temp/sra_processing"
    threads: 80
    shell:
        """
        gatk --java-options "-Xmx110G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bam} \
            -O {output.gvcf} \
            --tmp-dir {params.tmp_dir} \
            -ERC GVCF
        """

# Not used 

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

rule genomicsdb:
    input:
        directory(expand("data/genomicsdb/chr_{chrm}", chrm=CHROMOSOMES)),
        directory(expand("outgroup/genomicsdb/{species}/chr_{chrm}",  chrm=CHROMOSOMES, species=species_samples.keys()))

rule genomics_db_in:
    input:
        gvcfs=expand("data/{individual}/{individual}_raw_variants.g.vcf", individual=sra_data.keys()),
        gvcf_indexes=expand("data/{individual}/{individual}_raw_variants.g.vcf.idx", individual=sra_data.keys()),
    output:
        db = directory("data/genomicsdb/chr_{chrm}")
    params:
        interval=lambda wildcards: "{wildcards.chrm}",
        tmp_dir=lambda wildcards: f"temp/genomicsdb_import/{wildcards.chrm}"
    resources:
        runtime = 120,
        tasks=1
    threads: 80
    run:
        input_gvcfs = " ".join([f"-V {gvcf}" for gvcf in input.gvcfs])
        shell(f"""
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xmx64G" GenomicsDBImport \
            --genomicsdb-workspace-path {output.db} \
            -L {params.interval} \
            --tmp-dir {params.tmp_dir} \
            --reader-threads {threads} \
            {input_gvcfs}
        """)

rule genomics_db_outgroup:
    input:
        gvcfs=lambda wc: expand(
            "outgroup/data/{sample}/{sample}_raw_variants.g.vcf",
            sample=species_samples[wc.species]
        ),
        gvcf_indexes=lambda wc: expand(
            "outgroup/data/{sample}/{sample}_raw_variants.g.vcf.idx",
            sample=species_samples[wc.species]
        )
    output:
        db=directory("outgroup/genomicsdb/{species}/chr_{chrm}")
    params:
        interval=lambda wc: wc.chrm,
        tmp_dir=lambda wc: f"temp/genomicsdb_import/{wc.species}/{wc.chrm}"
    threads: 80
    resources:
        runtime=120,
        tasks=1
    run:
        input_gvcfs = " ".join(f"-V {g}" for g in input.gvcfs)
        shell(f"""
            mkdir -p {params.tmp_dir}
            gatk --java-options "-Xmx64G" GenomicsDBImport \
                --genomicsdb-workspace-path {output.db} \
                -L {params.interval} \
                --tmp-dir {params.tmp_dir} \
                --reader-threads {threads} \
                {input_gvcfs}
        """)

rule genotype:
    input:
        expand("data/vcfs/chr_{chrm}_cohort_genotyped.vcf.gz", chrm=CHROMOSOMES),
        expand("outgroup/vcfs/{species}/chr_{chrm}_cohort_genotyped.vcf.gz",  chrm=CHROMOSOMES, species=species_samples.keys())
        
rule genotype_gvcfs_outgroup:
    input:
        db="outgroup/genomicsdb/{species}/chr_{chrm}",
        reference="reference/ref_genome.fasta"
    output:
        vcf="outgroup/vcfs/{species}/chr_{chrm}_cohort_genotyped.vcf.gz",
        vcf_index="outgroup/vcfs/{species}/chr_{chrm}_cohort_genotyped.vcf.gz.tbi"
    params:
        tmp_dir=lambda wildcards: f"temp/genotype_gvcfs/{wildcards.species}/{wildcards.chrm}"
    resources:
        runtime=90,
        tasks=1
    threads: 80
    shell:
        """
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xmx128G" GenotypeGVCFs \
            -R {input.reference} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            --tmp-dir {params.tmp_dir}
        """

rule genotype_gvcfs:
    input:
        db="data/genomicsdb/chr_{chrm}",
        reference="reference/ref_genome.fasta"
    output:
        vcf="data/vcfs/chr_{chrm}_cohort_genotyped.vcf.gz",
        vcf_index="data/vcfs/chr_{chrm}_cohort_genotyped.vcf.gz.tbi"
    params:
        tmp_dir=lambda wildcards: f"temp/genotype_gvcfs/{wildcards.chrm}"
    resources:
        runtime = 800
    shell:
        """
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xmx128G" GenotypeGVCFs \
            -R {input.reference} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            --tmp-dir {params.tmp_dir}
        """

# not used, ran processing per-chromosome 

rule gather_vcfs:
    input:
        vcfs=expand("outgroup/vcfs/{species}/chr_{chrm}_cohort_genotyped.vcf.gz", species=species_samples.keys(), chrm=interval_list.values()),
        vcf_indexes=expand("outgroup/vcfs/{species}/chr_{chrm}_cohort_genotyped.vcf.gz.tbi", species=species_samples.keys(), chrm=interval_list.values())
        # vcfs=expand("data/vcfs/chr_{chrm}_cohort_genotyped.vcf.gz", chrm=interval_list.values()),
        # vcf_indexes=expand("data/vcfs/chr_{chrm}_cohort_genotyped.vcf.gz.tbi", chrm=interval_list.values())
    output:
        merged_vcf="outgroup/vcfs/{species}/cohort_merged.vcf.gz"
        # merged_vcf="data/vcfs/cohort_merged.vcf.gz"
    params:
        tmp_dir="temp/gather_vcfs"
    resources:
        runtime=10
    threads: 80
    run:
        input_vcfs = " ".join([f"-I {vcf}" for vcf in input.vcfs])
        shell(f"""
        mkdir -p {params.tmp_dir}
        gatk --java-options "-Xmx16G" GatherVcfs \
            {input_vcfs} \
            -O {output.merged_vcf} 
        """)
        
rule index_gvcf:
    input:
        vcf="outgroup/vcfs/{species}/cohort_merged.vcf.gz"
        # vcf="data/vcfs/cohort_merged.vcf.gz"
    output:
        filtered_vcf="outgroup/vcfs/{species}/cohort_merged.vcf.gz.tbi"
        # filtered_vcf="data/vcfs/cohort_merged.vcf.gz.tbi"
    resources:
        runtime=10
    threads: 80
    shell:
        """
        gatk --java-options "-Xmx16G" IndexFeatureFile \
            -I {input.vcf} 
        """

rule filter:
    input:
        expand("data/vcfs/chr_{chrm}_filtered_missing.vcf.gz", chrm=CHROMOSOMES2),
        expand("data/vcfs/chr_{chrm}_filtered.vcf.gz",  chrm=CHROMOSOMES),
        expand("data/vcfs/chr_{chrm}_filtered_missing.vcf.gz.tbi",  chrm=CHROMOSOMES2)
                
rule select_variants_in:
    input:
        vcf="data/vcfs/chr_{chrm}_cohort_genotyped.vcf.gz",
        reference="reference/ref_genome.fasta",
        index="data/vcfs/chr_{chrm}_cohort_genotyped.vcf.gz.tbi"
    output:
        filtered_vcf="data/vcfs/chr_{chrm}_filtered.vcf.gz"
    resources:
        runtime=10
    threads: 80
    shell:
        """
        gatk --java-options "-Xmx16G" SelectVariants \
            -R {input.reference} \
            -V {input.vcf} \
            -O {output.filtered_vcf} \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            --exclude-filtered true
        """
               
rule filter_missing:
    input:
        vcf="data/vcfs/chr_{chrm}_filtered.vcf.gz"
    output:
        vcf="data/vcfs/chr_{chrm}_filtered_missing.vcf.gz"
    shell:
        """
        plink2 \
          --vcf {input.vcf} \
          --allow-extra-chr \
          --geno 0.00 \
          --export vcf bgz \
          --out data/vcfs/chr_{wildcards.chrm}_filtered_missing
          # --out data/vcfs/filtered_missing
        """

rule index_gvcf_two:
    input:
        vcf="data/vcfs/chr_{chrm}_filtered_missing.vcf.gz"
    output:
        index="data/vcfs/chr_{chrm}_filtered_missing.vcf.gz.tbi"
    resources:
        runtime=10
    threads: 80
    shell:
        """
        gatk --java-options "-Xmx16G" IndexFeatureFile \
            -I {input.vcf} 
        """      

rule phase:
    input:
        expand("data/vcfs/chr_{chrm}_phased_vcf.vcf.gz", chrm=CHROMOSOMES)

rule phase_vcfs:
    input:
        vcf="data/vcfs/chr_{chrm}_filtered_missing.vcf.gz",
        path="beagle.17Dec24.224.jar"
    output:
        out="data/vcfs/chr_{chrm}_phased_vcf.vcf.gz"
    shell:
        """
        java -Xmx150g -jar {input.path} gt={input.vcf} out=data/vcfs/chr_{wildcards.chrm}_phased_vcf impute=false
        """
        
rule filter_out:
    input:
        expand("outgroup/vcfs/{species}/chr_{chrm}_filtered_filtered.vcf.gz", chrm=CHROMOSOMES2, species=species_samples.keys())

rule filter_outgroup:
    input:
        out_vcf = "outgroup/vcfs/{species}/chr_{chrm}_cohort_genotyped.vcf.gz",
        out_tbi = "outgroup/vcfs/{species}/chr_{chrm}_cohort_genotyped.vcf.gz.tbi",
        in_vcf = "data/vcfs/chr_{chrm}_filtered_missing.vcf.gz",
        in_tbi = "data/vcfs/chr_{chrm}_filtered_missing.vcf.gz.tbi"
    output:
        vcf = "outgroup/vcfs/{species}/chr_{chrm}_filtered_filtered.vcf.gz",
        tbi = "outgroup/vcfs/{species}/chr_{chrm}_filtered_filtered.vcf.gz.tbi"
    shell:
        """
        bcftools isec -n=2 -w1 {input.out_vcf} {input.in_vcf} | bgzip -c > {output.vcf}
        bcftools index --tbi {output.vcf}
        """

rule estsfss:
    input:
        expand("data/estsfs/chr_{chrm}_estsfs_input.txt", chrm=CHROMOSOMES2),
        expand("data/estsfs/chr_{chrm}_estsfs_input_part{part}.txt", chrm=CHROMOSOMES2, part=PARTS),
        expand("data/estsfs/output/chr_{chrm}_output-file-sfs_part{part}.txt", chrm=CHROMOSOMES2, part=PARTS),
        expand("data/estsfs/output/chr_{chrm}_output-file-pvalues_part{part}.txt", chrm=CHROMOSOMES2, part=PARTS)

rule est_sfs_input:
    input:
        in_vcf = "data/vcfs/chr_{chrm}_filtered_missing.vcf.gz",
        out_vcf1 = "outgroup/vcfs/american_crow/chr_{chrm}_filtered_filtered.vcf.gz",
        out_vcf2 = "outgroup/vcfs/rook/chr_{chrm}_filtered_filtered.vcf.gz",
        out_vcf3 = "outgroup/vcfs/jackdaw/chr_{chrm}_filtered_filtered.vcf.gz"
    output:
        out_txt = "data/estsfs/chr_{chrm}_estsfs_input.txt"
    shell:
        """
        python3.11 est-sfs-input.py {input.in_vcf} {output.out_txt} {input.out_vcf1} {input.out_vcf2} {input.out_vcf3}
        """

############################################################################
        # I tried to split to 1/2, 1/4, but the files are truly massive (for est-sfs). If I split into 1/8, the largest files are ~4 MB, which est-sfs can handle in around ~1.5 hours. the 1/4thed files, which were around 9MB max, took longer than 5 hours...
############################################################################

rule split_estsfs:
    input:
        "data/estsfs/chr_{chrm}_estsfs_input.txt"
    output:
        "data/estsfs/chr_{chrm}_estsfs_input_part1.txt",
        "data/estsfs/chr_{chrm}_estsfs_input_part2.txt",
        "data/estsfs/chr_{chrm}_estsfs_input_part3.txt",
        "data/estsfs/chr_{chrm}_estsfs_input_part4.txt",
        "data/estsfs/chr_{chrm}_estsfs_input_part5.txt",
        "data/estsfs/chr_{chrm}_estsfs_input_part6.txt",
        "data/estsfs/chr_{chrm}_estsfs_input_part7.txt",
        "data/estsfs/chr_{chrm}_estsfs_input_part8.txt"
    shell:
        """
        set -euo pipefail
        INPUT={input}
        P1={output[0]}
        P2={output[1]}
        P3={output[2]}
        P4={output[3]}
        P5={output[4]}
        P6={output[5]}
        P7={output[6]}
        P8={output[7]}
        line_count=$(wc -l < "$INPUT")
        chunk_size=$(( (line_count + 7) / 8 ))
        head -n "$chunk_size" "$INPUT" > "$P1"
        tmp2=$(mktemp)
        tail -n +"$((chunk_size + 1))" "$INPUT" > "$tmp2"
        head -n "$chunk_size" "$tmp2" > "$P2"
        rm -f "$tmp2"
        tmp3=$(mktemp)
        tail -n +"$(((2 * chunk_size) + 1))" "$INPUT" > "$tmp3"
        head -n "$chunk_size" "$tmp3" > "$P3"
        rm -f "$tmp3"
        tmp4=$(mktemp)
        tail -n +"$(((3 * chunk_size) + 1))" "$INPUT" > "$tmp4"
        head -n "$chunk_size" "$tmp4" > "$P4"
        rm -f "$tmp4"
        tmp5=$(mktemp)
        tail -n +"$(((4 * chunk_size) + 1))" "$INPUT" > "$tmp5"
        head -n "$chunk_size" "$tmp5" > "$P5"
        rm -f "$tmp5"
        tmp6=$(mktemp)
        tail -n +"$(((5 * chunk_size) + 1))" "$INPUT" > "$tmp6"
        head -n "$chunk_size" "$tmp6" > "$P6"
        rm -f "$tmp6"
        tmp7=$(mktemp)
        tail -n +"$(((6 * chunk_size) + 1))" "$INPUT" > "$tmp7"
        head -n "$chunk_size" "$tmp7" > "$P7"
        rm -f "$tmp7"
        tmp8=$(mktemp)
        tail -n +"$(((7 * chunk_size) + 1))" "$INPUT" > "$tmp8"
        mv "$tmp8" "$P8"
        """

rule est_sfs:
    input:
        config = "data/estsfs/config.txt",
        data = "data/estsfs/chr_{chrm}_estsfs_input_part{part}.txt"
    output:
        file = "data/estsfs/output/chr_{chrm}_output-file-sfs_part{part}.txt",
        pvalues = "data/estsfs/output/chr_{chrm}_output-file-pvalues_part{part}.txt"
    resources:
        runtime = 150,
        tasks = 80
    threads: 1
    shell:
        """
        est-sfs-release-2.04/est-sfs {input.config} {input.data} \
            data/estsfs/seed-file.txt {output.file} {output.pvalues}
        """

rule make_haps:
    input:
       vcf = "data/vcfs/chr_{chrm}_phased_vcf.vcf.gz" 
    output:
       haps = "data/haps/chr_{chrm}_haps_orig.haps"
    shell:
        """
        plink2 --vcf  {input.vcf} --export haps --out data/haps/chr_{wildcards.chrm}_haps_orig --allow-extra-chr
        """          

rule polarize_haps:
    input:
        est_sfs = "data/estsfs/output/chr_{{chrm}}_output-file-sfs_part{part}.txt",
        est_pvals = "data/estsfs/output/chr_{{chrm}}_output-file-pvalues_part{part}.txt",
        hap = "data/haps/chr_{chrm}_haps_orig.haps",
    output:
        polar = "data/haps/chr_{chrm}_haps_polarized.haps"
    threads: 1
    script:
        "polar.py"
