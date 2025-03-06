import sys
import pysam
from collections import Counter
import random

"""
From est-sfs documentation: 
Input (.txt file)
If there are three outgroups, there are 4 space-separated columns. 

The first column is for the focal species, and the next three columns are for the outgroups. Each column is a comma-separated list of the counts of the four bases in the order A, C, G, T. For example, the first line in the example data file is:

20,0,0,0 0,0,0,1 0,0,0,1 0,0,0,1

At this site, all n = 20 copies sampled in the focal species are A. In the three outgroups, a single copy has been sampled, and in each case it is T. All sites must have the same number of copies sampled in the focal species and up to one copy sampled in each outgroup. If there are missing data in any outgroup, the counts for that outgroup are encoded 0,0,0,0 . Data from polymorphic and non-polymorphic sites are analysed together.
"""

def load_outgroup_vcf(vcf_path):
    """
    Load a multi-sample outgroup VCF, gather all alleles at each site,
    and produce a single base using majority-rule consensus.
    Ties => None (missing).
    
    Returns a dict keyed by (chrom, pos) -> single base 'A','C','G','T' or None.
    """
    outdict = {}
    if not vcf_path:
        return outdict  

    v = pysam.VariantFile(vcf_path)
    all_samples = list(v.header.samples)  

    for rec in v:
        # skips non-biallelic
        if len(rec.alts) != 1:
            continue
        if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
            
        allele_list = []  
        ref_base = rec.ref
        alt_base = rec.alts[0]

        for sample_name in all_samples:
            gt = rec.samples[sample_name]["GT"]
            if gt is None:
                continue  

            for allele_idx in gt:
                if allele_idx == 0:
                    allele_list.append(ref_base)
                elif allele_idx == 1:
                    allele_list.append(alt_base)

        if not allele_list:
            outdict[(rec.chrom, rec.pos)] = None
            continue

        counts = Counter(allele_list)  
        base, top_count = counts.most_common(1)[0]
        top_bases = [b for b,cnt in counts.items() if cnt == top_count]
        if len(top_bases) > 1:
            base = random.choice(top_bases)
        outdict[(rec.chrom, rec.pos)] = base

    return outdict


def base_to_4tuple(base): # converts base to the tuple string needed by est-sfs
    if base == "A":
        return (1,0,0,0)
    elif base == "C":
        return (0,1,0,0)
    elif base == "G":
        return (0,0,1,0)
    elif base == "T":
        return (0,0,0,1)
    else:
        return (0,0,0,0)

def sum_tuples(t1, t2): # element-wise summation
    return tuple(a+b for (a,b) in zip(t1, t2))

def main(): 

    if len(sys.argv) < 3:
        sys.exit(
            "Usage:\n  "
            + f"{sys.argv[0]} INGROUP_VCF OUTPUT_TXT [OUTGROUP1_VCF] [OUTGROUP2_VCF] [OUTGROUP3_VCF]\n"
        ) 
        
    ingroup_vcf = sys.argv[1]
    output_txt = sys.argv[2]
    outgroup1_vcf = sys.argv[3] if len(sys.argv) > 3 else None
    outgroup2_vcf = sys.argv[4] if len(sys.argv) > 4 else None
    outgroup3_vcf = sys.argv[5] if len(sys.argv) > 5 else None

    min_allele_count = 1

    out1_dict = load_outgroup_vcf(outgroup1_vcf)
    out2_dict = load_outgroup_vcf(outgroup2_vcf)
    out3_dict = load_outgroup_vcf(outgroup3_vcf)

    invcf = pysam.VariantFile(ingroup_vcf)
    sample_names = list(invcf.header.samples)

    with open(output_txt, "w") as fout:
        for rec in invcf:
            if len(rec.alts) != 1: #skips non-biallelic, redundant, but might as well
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue

            ref = rec.ref
            alt = rec.alts[0]

            total_4 = (0,0,0,0)
            for sample in sample_names:
                gt = rec.samples[sample]["GT"]
                if gt is None:
                    continue
                for allele_idx in gt:
                    if allele_idx == 0:
                        total_4 = sum_tuples(total_4, base_to_4tuple(ref))
                    elif allele_idx == 1:
                        total_4 = sum_tuples(total_4, base_to_4tuple(alt))

            alt_4 = base_to_4tuple(alt)
            if 1 in alt_4:
                alt_index = alt_4.index(1)
                alt_count = total_4[alt_index]
            else:
                alt_count = 0
            if alt_count < min_allele_count:
                continue

            out1_base = out1_dict.get((rec.chrom, rec.pos), None)
            out2_base = out2_dict.get((rec.chrom, rec.pos), None)
            out3_base = out3_dict.get((rec.chrom, rec.pos), None)

            out1_4 = base_to_4tuple(out1_base)
            out2_4 = base_to_4tuple(out2_base)
            out3_4 = base_to_4tuple(out3_base)

            focal_str = ",".join(map(str, total_4))
            out1_str = ",".join(map(str, out1_4))
            out2_str = ",".join(map(str, out2_4))
            out3_str = ",".join(map(str, out3_4))

            line = f"{focal_str} {out1_str} {out2_str} {out3_str}"
            fout.write(line + "\n")


if __name__ == "__main__":
    main()
