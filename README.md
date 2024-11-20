# phylogenetictrees
Create phylogenetic trees for clonally related cells

# Constructing phylogenetic trees
The R script in this repository constructs phylogenetic trees for *known* clonally related cells. It will identify shared and unqiue somatic mutations within clonally related cells and construct a tree based on these mutations. Trunk and branch distances are proportional to the number of shared and unique somatic mutations between cells. The identity of each cell is labeled at the terminus of each branch.

*NOTE: Cells with clonal relationships should already be established before running this script.* Refer to the "Running the Script" section for more details.

### Installation
The necessary R package dependencies will not be installed upon running the script. Please make sure dplyr, tidyr, ggplot2 and readxl have been installed prior to running. This script has been verified to run on R 4.2.2 and higher.

### Running the Script
This script requires an Excel sheet with at least two columns: "Genome_Change" and "ID". It iterates through the "Genome_Change" column to identify repeated entries and creates subsets of clonal relationships based on the occurrence of each entry within each subset. Therefore, the input sheet should be a concatenation of somatic mutations (rows) from the set of clonally related cells, and each mutation (row) should have a corresponding sample ID from which it was inferred. 
 - WARNING: only samples within the same phylogeny can be combined in one input sheet or the script will not work. If there is more than one group of clonally related cells in the sheet (i.e. would be more than one phylogenetic tree), they should be separated out into different sheets and the script should be run separately on each group of related cells.
     - For example; Donor A has cells 1, 2, 3, 4, 5, 6 and 7. Within this donor, cells 1, 2, 3 and 4 are clonally related (let's call this Group A), while cells 6 and 7 are clonally related as well (let's call this Group B). In this case, Group A and B are not related to each other although individually they have clonal relationships. Therefore, you should have 2 input sheets to run separately; one for Group A and one for Group B. However, if were to you run Group A and happened to leave cell 5 in the dataset it will simply get filtered out (since it is neither part of A or B and is unrelated).

*Note*: Your input sheet DOES NOT need to be in the same format as the Genome_Change column in our example sheet; as long as the Genome_Change column are represented in a similar manner for all cells within your dataset, the script should be able to find repitions between cells and extract the *ovelapping and unique mutations*.

Refer to the "Preparation of our Input Sheet" section for a detailed description of how our sheet was built if you are interested in replicating our methods.

Using the example input file provided here, the command to run the script would be:

	Rscript phylogenetic_tree_generation.R All_D56_Mutations.xlsx
 
### Preparation of our Input Sheet
The example input file provided here was generated using a number of steps. After aligning sequencing data from genomic DNA to the hg19 version of the genome with BWA and deduplicating with Picard, the reads underwent additional curation to realign indels called with Pindel, and recalibrate base quality using GATK. Then:
1. MuTect2 was used to generate a candidate list of point mutations by comparing the aligned bam files of each sample to the bam files representing the respective patient’s normal DNA.
2. To enhance rigor, variants were called against the reference genome using both UnifiedGenotyper and FreeBayes. These additional variant callers were included to identify any point mutations potentially missed by MuTect. Germline SNPs were filtered out, and the remaining variants were added to the list of somatic variants identified by MuTect.
3. Indels identified using Pindel were also added to the list of variants.
4. This list of variants was scrutinized as described previously by Tang et al.(Nature, 2020) to remove artifacts.
5. After generating the mutation list for each sample, the mutations from all samples were compiled into an Excel sheet with their respective sample_IDs, organized by "Genome_Change" column.
6. Overlapping mutations were identified between samples based on the recurrence of specific mutations across different samples.
7. The remaining non-overlapping mutations were further scrutinized by running Mpileup on each sample to rule out exclusion due to preset Mutect2 parameters. If identified, these mutations were added to the respective sample's list of mutations.
8. Only samples with known overlaps, as determined by this pipeline, should be used to generate a spreadsheet compiling mutations for those samples. For example, if 10 samples are sequenced and mutation overlaps are identified between samples 1 and 2, as well as between samples 3 and 4, two separate spreadsheets should be created: one for samples 1 and 2, and another for samples 3 and 4. These spreadsheets will then serve as separate inputs for the R script.
9. In the example sheet provided in this repository, all irrelevant columns have been removed for simplicity, however this is not a requirement.
