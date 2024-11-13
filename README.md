# Constructing phylogenetic trees
The R script in this repository constructs phylogenetic trees for *known* clonally related cells. It will identify shared and unqiue somatic mutations within clonally related cells and construct a tree based on these mutations. Trunk and branch distances are proportional to the number of shared and unique somatic mutations between cells. The identity of each cell is labeled at the terminus of each branch.

*Note: Cells with clonal relationships should already be established before running this script.* 

Refer to the "Preparation of the Input Sheet" section for a detailed description of this requirement.

### Installation
The necessary R package dependencies will not be installed upon running the script. Please make sure readxl, dplyr, tidyr and ggplot2 have been installed prior to running. This script has been verified to run on R 4.2.2 and higher.

### Preparation of the Input Sheet
The example input file provided here was generated using a number of steps. After aligning sequencing data from genomic DNA to the hg19 version of the genome with BWA and deduplicating with Picard, the reads underwent additional curation to realign indels called with Pindel, and recalibrate base quality using GATK. Then:
1. MuTect2 was used to generate a candidate list of point mutations by comparing the aligned bam files of each cell to the bam files representing the respective patient’s normal DNA.
2. To supplement MuTect calls, variants were called against the reference genome using UnifiedGenotyper and FreeBayes. These variant callers were incorporated to identify any point mutations missed by MuTect. Germline SNPs were removed and variants that were not filtered out were considered somatic mutations and added to the list of somatic variants (mutation list) called by MuTect.
3. Upon generation of the mutation list for each sample, all genome_changes were compiled in an excel sheet with their respective sample_IDs, sorted by Genome_Change and using a simple excel equation, overlapping mutations were identified.
4. Mpileup was then run on the list of unique mutations within subsets of clonally related cells to find mutations that were missed in some samples due to set Mutect2 parameters. These mutations were then added to each sample's final list of mutations.
5. A simple Rscript was then used to concatenate the mutation lists from groups of related samples. Only cells sequenced with same baits are combined.
   - Note: only samples within the same phylogeny can be combined in one input sheet. If there is more than one group of clonally related cells in the sheet (i.e. would be more than one phylogenetic tree), they should be separated out into different sheets and the script should be run separately on each group of related cells.
     - For example; Donor A has cells 1, 2, 3, 4, 5, 6 and 7. Within this donor, cells 1, 2, 3 and 4 are clonally related (let's call this Group A), while cells 6 and 7 are clonally related as well (let's call this Group B). In this case, Group A and B are not related to each other although individually they have clonal relationships. Therefore, you should have 2 input sheets to run separately; one for Group A and one for Group B.

### Running the Script
At minimum, this script requires an excel sheet containing 2 columns named: "Genome_Change" and "ID". The script iterates through the Genome_Change column to find reptitions and creates subsets of clonal relationships based on its occurence within each subset. Therefore, the sheet should be a concatenation of somatic mutations (rows) from the set of clonally related cells (as mentioned in point 5 above), and each mutation (row) should have a corresponding sample ID from which it was extracted. This sheet does not need to be in the same format as the Genome_Change in the sheet provided here; as long as the Genome_Change column contains the repitions between cells and will be able to extract the *ovelapping mutations*.


Using the example input file provided here, the command to run the script would be:

	Rscript phylogenetic_tree_generation.R All_D56_Mutations.xlsx
