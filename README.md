# Repository of analysis scripts for the Grosshans 2022 paper on GRH-1

Repo assembled by AATS from own work & previous work by MWMM & SHC.

Main sections:
- GRH-1 ChIP-seq processing
- promoter occupancy (RNA-polII ChIP-seq & mRNA-seq timecourse)
- phase enrichment analysis amongst putative GRH-1 targets



<!-- ----------------------------------------------------------------------- -->
## GRH-1_ChIP-seq

This folder contains RStudio notebooks that process the ChIP-sequencing data to
called peaks, motifs and sites.

### 0-PrepareGenes
Processes WS220 exon definitions for coding transcripts & miRNA non-mitochondrial genes into exon GenomicRanges, pan-isoform gene promotoer & body GenomicRanges, annotated with gene symbols and transcript names, for further use.

### 1-AlignReads
Align ChIP-sequencing SE50 reads from both experiments/runs (2287R and 2633R) to the ce10 genome.

### 2-EnrichedTiles
The ce10 genome was split into 500bp tiles. Tiles were qualified by GC content and by number of reads mapping to them. Input and ChIP samples (and enrichments derived therefrom) were examined for evidence of GC bias or other issues.

### 3-PrepareMACS2Peaks
Peaks were called using MACS2 for sample input/ChIP pairs from both experiments. Peaks overlapping "overmapped" tiles were removed from later analysis. Peaks were then checked for issues (GC content, low mappability, discrepancy between MACS2 and manually-calculated enrichments...).

### 4-Motifs
The ce10 genome was masked to remove unmappable regions. The 10% of leaste enriched peaks were filtered out. HOMER was then run using various setups. Obtained motifs were assessed manually across setups to select candidate GRH-1 motifs for further analysis. Scrambled versions of candidates were generated. The mappable ce10 genome was scanned for sites corresponding to candidate and scrambled motifs. Agreement between peaks and sites was evaluated in various ways. A final candidate motif was retained for further analyses.

### 5-AssignMACS2Peaks2Genes
Binding site information was aggregated to the peak level, and peaks were assigned to genes based on promoter overlap, gene body overlap, or by lying directly upstream (with a distance cutoff).



<!-- ----------------------------------------------------------------------- -->
## Promoter Occupancy

This folder contains scripts analysing the timecourse experiment looking into gene promoter occupancy (via pol-II ChIP-seq) and mRNA expression (RNA-seq).


### polII_timecourse_analysis
This script counts hits to 1kb regions centered on putative TSSs (ie that of the combined "gene body"), from samples of Run 899 aligned to ce10.
Counts are reported per gene, library-size normalised (min lib size,
lib size from project alignment stats), log2-transformed w/ a pseudocount of 8,
filtered to "expressed" genes (mean>3), then quantile normalised.
Only counts for the pull-down samples (ie "polII", ie not input) are kept, 22h-33h.
These are measures of promoter occupancy by RNA pol-II.


### mRNA_timecourse_analysis
This script counts hits to gene exons, from samples of RNA-seq run 1314 aligned to ce10.
Counts are aggregated to gene-level by qCount, library-size normalised
(mean lib size, lib size from summed hit counts), log2-transformed w/ a pseudocount of 8.


### PolII_RNA_corr
This script loads previously-prepared promoter occupancy and mRNA expression level data, calculates correlation between the 2, and describes each gene in terms of each dataset, said correlation and oscillatory status.



<!-- ----------------------------------------------------------------------- -->
## Phase Enrichment Analyses

### GenePeakMotif_Analyses
This script takes site-peak-gene assignments generated above, derives various definitions of "putative GRH-1 target" from this data (eg thresholding on site score, peak enrichment, etc), and investigates whether:
- there is evidence of enrichment in oscillatory genes amongst putative GRH-1 targets
- there is evidence of enrichment in particular phases amongst oscillatory putative GRH-1 targets




<!-- EOF -->
