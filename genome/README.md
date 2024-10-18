# genome

This directory will populate with various STAR and RSEM index files once Part 1 - data preprocessing is complete.
Once the index has been built once, it can be stored and re-used.
Key input files are bolded.

<pre>
.
|______genome 
       |____<b>Species.gtf</b> #Annotation file
       |____<b>Species.fa</b> or <b>Species.fna</b> #Reference genome
       |____<b>species_annotated_id_table_entrez_wholegenome.csv</b> #additional annotations - ensemble_ids, gene names, descriptions, zebrafish homolog ensembl ids and gene names (if applicable)
       |____chrLength.txt #STAR index file - built after part 1 - data preprocessing
       |____chrNameLength.txt #STAR index file - built after part 1 - data preprocessing
       |____chrName.txt #STAR index file - built after part 1 - data preprocessing
       |____chrStart.txt #STAR index file - built after part 1 - data preprocessing
       |____exonGeTrInfo.tab #STAR index file - built after part 1 - data preprocessing
       |____exoninfo.tab #STAR index file - built after part 1 - data preprocessing
       |____geneInfo.tabt #STAR index file - built after part 1 - data preprocessing
       |____genomeParameters.txt #STAR index file - built after part 1 - data preprocessing
       |____Log.out #STAR index file - built after part 1 - data preprocessing
       |____chrStart.txt #STAR index file - built after part 1 - data preprocessing
       |____SjdbInfo.txt #STAR index file - built after part 1 - data preprocessing
       |____SjdbList.fromGTF.out.tab #STAR index file - built after part 1 - data preprocessing
       |____SjdbList.out.tab.txt #STAR index file - built after part 1 - data preprocessing
       |____transcriptInfo.tab #STAR index file - built after part 1 - data preprocessing
       |____Genome #STAR index file - built after part 1 - data preprocessing
       |____SA #STAR index file - built after part 1 - data preprocessing
       |____SAindex #STAR index file - built after part 1 - data preprocessing
       | 
       |____RSEM 
               |____species.grp #RSEM index - built after part 1 - data preprocessing
               |____species.chrlist #RSEM index - built after part 1 - data preprocessing
               |____species.idx.fa #RSEM index - built after part 1 - data preprocessing
               |____species.n2g.idx.fa #RSEM index - built after part 1 - data preprocessing
               |____species.seq #RSEM index - built after part 1 - data preprocessing
               |____species.ti #RSEM index - built after part 1 - data preprocessing
               |____species.transcripts.fa #RSEM index - built after part 1 - data preprocessing
</pre>

## Creating the custom annotation file

First, download the appropriate .FASTA and .GTF files for your species of interest.

Next, in the `/genome/Create_custom_annotation.R` script, define the species for which you want to create an annotation for.
Define the name of your downloaded GTF file in this location.
If your species of interest if Fathead Minnow, you can perform your own BLASTN, or use the `/genome/summarized_annotated_results_FHM_with_ZF_homologs.txt` file for which I have already conducted the BLASTN for your convenience. The Fathead Minnow genome assembly I used was [Genome assembly EPA_FHM_2.0 (NCBI RefSeq) GCF_016745375.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016745375.1/).
Update the version of the ensemble gene list to match the version of your GTF annotation file. It nay be necessary to visit [Ensembl archives](https://useast.ensembl.org/info/website/archives/index.html) and explore [archived version of BioMart](https://may2024.archive.ensembl.org/biomart/martview/4f5066d06ff46f72f5a87e31aff8adca) to find an ensembl gene list for your species.

Run the R Script!