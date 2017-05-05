# hiccomp
for manuscript purpose only (scripts created by Motakis Efthymios)

# Dependencies:
wavethresh (4.6.8)

 
# Usage:
ChangeIndex(x=inputfile, chrom=chromname, B=999, cut_disagreement=1)
Example:
ChangeIndex(x="Sham_10kb.txt", chrom="chr1", B=999, cut_disagreement=1)

# Details:
inputfile: tab-delimited file with 6 columns [produced by Homer PCA analyis] ( <unique_id> <chromosome> <start_position> <end_position> <strand> <PC1_score>)
1. chromname: name of chromosome, depending on the format written in inputfile (e.g. "chr1")
2. B: The number of subsamples to generate to assess data noise
3. cut_disagreement: the cut-off for the disagreement method


# Output:
A data frame with 9 columns:
1. unique_id : provided in input
2. chromosome : provided in input
3. start_position : provided in input
4. end_position : provided in input
5. strand : provided in input
6. PC1_score : provided in input
7. Disagreement P : quantified the stability of PC1 classification under model-derived data perturbations
8. Instability Index : frequency 
9. Compartment : Active / Inactive compartment (A|B)

