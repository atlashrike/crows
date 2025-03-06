"""

For each interval (left_snp, right_snp), reads LDhelmet’s “mean” (ρ per bp), converts it to cM/Mb via the formula cM/Mb = (ρ × 1e–6 / (4×Ne)) × 100, and multiplies that rate by the base‐pair length of the interval to get an increment in cM, then accumulates those increments to track a running total of the genetic position in cM.  
Finally, writes a three‐column .map file:  
    (1) pos (the right_snp base),  
    (2) the per‐interval cM/Mb (COMBINED_rate),  
    (3) the cumulative genetic map (Genetic_Map).
    
From Paul-Jenkins on GitHub:

The units are rho/bp. So if you wanted to convert the output to cM/Mb you would need to multiply by 10^-6 (to get rho/Mb) then divide by 4Ne where Ne is the diploid effective population size (to get recombination event rate per Mb) then multiply by 100 to get cM/Mb.

- ldhat_out.txt is a file with lines like:
    # left_snp right_snp mean p0.025 p0.500 p0.975
    13 50 1.1925e-04 5.9240e-06 1.1949e-04 2.4561e-04
    50 60 1.1925e-04 5.9240e-06 1.1949e-04 2.4561e-04
    ...

The resulting .map file has 3 columns per line:
  position (bp)   recombination_rate (cM/Mb)   genetic_position (cM)

Relate's documentation states:
This is the standard file format for genetic recombination maps. The three columns are:
Position (b) [integer]
Recombination rate (cM/Mb) [float]
Genetic position (cM) [float]
r[i] = (rdist[i+1] - rdist[i])/(p[i+1] - p[i]) * 1e6
 
The little emoticons are for fun, and will not show unless LC_ALL=en_US.UTF-8
(echo $LC_ALL / export LC_ALL=en_US.UTF-8)
 
"""
import sys
import argparse
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(
        description="꒰ᐢ. .ᐢ꒱₊˚⊹ Convert LDhelmet .txt output (ρ/bp) into a .map file for Relate input."
    )
    parser.add_argument("--input",  "-i", required=True, help="꒰ᐢ. .ᐢ꒱₊˚⊹ Path to LDhelmet text output file.")
    parser.add_argument("--output", "-o", required=True, help="꒰ᐢ. .ᐢ꒱₊˚⊹ Path to the desired .map output.")
    parser.add_argument("--ne",     type=float, default=2e5,
                        help="꒰ᐢ. .ᐢ꒱₊˚⊹ Diploid effective population size (default = 2e5).")
    return parser.parse_args()

def convert_rho_to_cM_Mb(rho_per_bp, ne):
    """
    If LDhelmet outputs ρ = 4*N_e*r (per bp),
    then cM/Mb = (ρ * 1e-6 / (4*Ne)) * 100.
    """
    return rho_per_bp * 1e-6 * 100 / (4 * ne)

def main():
    args = parse_args()
    with open(args.input, "r") as fin:
        all_lines = fin.readlines()
        
    # Collect data lines, skipping comments (#) and blank lines
    data_lines = []
    for line in all_lines:
        line_stripped = line.strip()
        if line_stripped and not line_stripped.startswith("#"):
            data_lines.append(line_stripped)

    cumulative_cM = 0.0
    output_rows   = []

    for row in tqdm(data_lines, desc=" (╭ರ_•́) Converting intervals..."):
        fields = row.split()
        if len(fields) < 3:
            continue  # skip malformed lines

        left_snp  = int(fields[0])
        right_snp = int(fields[1])
        rho_mean  = float(fields[2])  # "mean" column, in ρ per bp

        cM_Mb = convert_rho_to_cM_Mb(rho_mean, args.ne)

        interval_bp = right_snp - left_snp
        cM_increment = cM_Mb * (interval_bp / 1e6)

        cumulative_cM += cM_increment

        # We'll use right_snp as the genomic position in the .map file
        position_bp = right_snp

        output_rows.append((position_bp, cM_Mb, cumulative_cM))

    with open(args.output, "w") as fout:
        fout.write("pos COMBINED_rate Genetic_Map\n")
        for (pos, rate, gpos) in output_rows:
            # Print in scientific notation with 8 decimals
            fout.write(f"{pos} {rate:.8e} {gpos:.8e}\n")

    print(f"ヾ( ˃ᴗ˂ )◞ • *✰ Conversion complete ! Wrote {len(output_rows)} lines to {args.output}.", file=sys.stderr)

if __name__ == "__main__":
    main()
