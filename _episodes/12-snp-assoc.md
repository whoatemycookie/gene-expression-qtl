---
title: "SNP association mapping"
teaching: 30
exercises: 30
questions:
- "How do I identify SNPs in a QTL?"
objectives:
- Perform a basic QTL analysis.
- Identify QTL with a genome scan.
- Find SNPs within a QTL.
- Convert founder genotypes to a strain distribution pattern (SDP).
- Infer SNP genotypes for Diversity Outbred mice.
- Perform SNP association analysis 
keypoints:
- "SNP association analysis with DO mice requires SNPs in the QTL regions."
- "."
source: Rmd
---



For multi-parent crosses, it can be useful to collapse the genotype or allele probabilities according to the founder genotypes of the various
SNPs in the region of a QTL.

### QTL analysis in Diversity Outbred mice

To illustrate this sort of SNP association analysis, we'll consider some Diversity Outbred mouse data. The Diversity Outcross (DO) mice are an advanced intercross population derived from the same eight founder strains as the Collaborative Cross (CC). See
[Svenson et al. (2012)](https://www.ncbi.nlm.nih.gov/pubmed/22345611)
and [Gatti et al. (2014)](https://www.ncbi.nlm.nih.gov/pubmed/25237114).

We'll consider a subset of the data from
[Recla et al. (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24700285), available as part of the
[qtl2data github repository](https://github.com/rqtl/qtl2data). (The
full data are in
[`DO_Recla`](https://github.com/rqtl/qtl2data/tree/master/DO_Recla); the directory
[`DOex`](https://github.com/rqtl/qtl2data/tree/master/DOex) contains a reduced set, with just three chromosomes, one phenotype (`OF_immobile_pct`, percent immobile in the open field test), and a
reduced set of markers.

You can download the data from a single zip file, as follows:


~~~
DOex <- read_cross2(file = "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip")
~~~
{: .language-r}

Let's quickly whip through a basic analysis.

We first calculate genotype probabilities and convert them to allele probabilities. We'll just use marker locations and not insert any pseudomarkers.


~~~
pr <- calc_genoprob(DOex, error_prob=0.002)
apr <- genoprob_to_alleleprob(pr)
~~~
{: .language-r}

We calculate kinship matrices (using the "LOCO" method, though with the caveat that here we are only considering genotypes on three chromosomes).


~~~
k <- calc_kinship(apr, "loco")
~~~
{: .language-r}

We create a numeric covariate for sex; be sure to include the individual IDs as names.


~~~
sex <- (DOex$covar$Sex == "male")*1
names(sex) <- rownames(DOex$covar)
~~~
{: .language-r}

Note that you can create the vector and assign the names in one step using the basic R function setNames().


~~~
sex <- setNames( (DOex$covar$Sex == "male")*1, rownames(DOex$covar) )
~~~
{: .language-r}

We perform a genome scan with a linear mixed model (adjusting for a residual polygenic effect), with sex as an additive covariate.


~~~
out <- scan1(apr, DOex$pheno, k, sex)
~~~
{: .language-r}

Here's a plot of the results.


~~~
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out, DOex$gmap)
~~~
{: .language-r}

<img src="../fig/rmd-12-plot_DOex_scan-1.png" title="plot of chunk plot_DOex_scan" alt="plot of chunk plot_DOex_scan" width="612" style="display: block; margin: auto;" />

There's a strong peak on chromosome 2. Let's look at the QTL effects. We estimate them with `scan1coef()`. We need to subset the allele probabilities and the list of kinship matrices.


~~~
coef_c2 <- scan1coef(apr[,"2"], DOex$pheno, k[["2"]], sex)
~~~
{: .language-r}

For the DO, with 8 QTL alleles, we can use the function `plot_coefCC` in the [R/qtl2plot](https://github.com/rqtl/qtl2plot) package, which plots the 8 allele effects in the "official" Collaborative Cross (CC)
colors. (Well, actually _slightly_ modified colors, because I think the official colors are kind of ugly.) The strong locus seems to be mostly
due to the NZO allele. Note that `CCcolors` is a vector of colors included in the qtl2plot package; there's also a `CCorigcolors` object
with the _official_ colors.


~~~
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_coefCC(coef_c2, DOex$gmap["2"], bgcolor="gray95", legend="bottomleft")
~~~
{: .language-r}

<img src="../fig/rmd-12-plot_DOex_effects-1.png" title="plot of chunk plot_DOex_effects" alt="plot of chunk plot_DOex_effects" width="612" style="display: block; margin: auto;" />

If you provide `plot_coefCC()` with the genome scan output, it will display the LOD curve below the coefficient estimates.


~~~
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_coefCC(coef_c2, DOex$gmap["2"], scan1_output=out, bgcolor="gray95", legend="bottomleft")
~~~
{: .language-r}

<img src="../fig/rmd-12-plot_DOex_lod_curve-1.png" title="plot of chunk plot_DOex_lod_curve" alt="plot of chunk plot_DOex_lod_curve" width="612" style="display: block; margin: auto;" />

### Connecting to SNP and gene databases
To perform SNP association analysis in the region of a QTL, we’ll need to grab data on all of the SNPs in the region, including their genotypes in the eight founder strains for the CC and DO populations. As a related task, we’ll also want to identify the genes in the region.

We want R/qtl2 to permit different ways of storing this information. For the CC and DO populations, we’ve prepared SQLite database files for the variants in the CC founders and for the MGI mouse gene annotations. But others might wish to use a different kind of database, or may wish to query an online database.

We provide a template for how to use R/qtl2 to connect to SNP and gene databases, with the functions `create_variant_query_func()` and `create_gene_query_func()`. Each returns a function for querying variant and gene databases, respectively. The query functions that are returned take just three arguments (chr, start, end end), and themselves return a data frame of variants on the one hand and genes on the other.

During [setup](https://smcclatchy.github.io/mapping/setup/) you would have downloaded the SQLite databases from Figshare and placed these in your `data` directory:

+ [SQLite database of variants in Collaborative Cross founder mouse strains (v3)](https://figshare.com/articles/SQLite_database_of_variants_in_Collaborative_Cross_founder_mouse_strains/5280229/3): SNP, indel, and structural variants in the Collaborative Cross founders (3.87 GB)
+ [SQLite database with mouse gene annotations from Mouse Genome Informatics (v6)](https://figshare.com/articles/dataset/SQLite_database_with_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5280238): full set of mouse gene annotations from build 38 mm10 (529.89 MB)
+ [SQLite database with MGI mouse gene annotations from Mouse Genome Informatics (v7)](https://figshare.com/articles/dataset/SQLite_database_with_MGI_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5286019): like the previous, but including only non-duplicate gene records sourced from MGI (11.16 MB)

To create a function for querying the CC variants, call `create_variant_query_func()` with the path to the `cc_variants.sqlite` file:


~~~
query_variants <- create_variant_query_func("../data/cc_variants.sqlite")
~~~
{: .language-r}

To grab the variants in the interval 97-98 Mbp on chromosome 2, you’d then do the following:


~~~
variants_2_97.5 <- query_variants(2, 97, 98)
~~~
{: .language-r}



~~~
Error in query_variants(2, 97, 98): File ../data/cc_variants.sqlite doesn't exist
~~~
{: .error}

Similarly, to create a function for querying the MGI mouse gene annotations, you call `create_gene_query_func()` with the path to the `mouse_genes_mgi.sqlite` file:


~~~
query_genes <- create_gene_query_func("../data/mouse_genes_mgi.sqlite")
~~~
{: .language-r}

To grab the genes overlapping the interval 97-98 Mbp on chromosome 2, you’d then do the following:


~~~
genes_2_97.5 <- query_genes(2, 97, 98)
~~~
{: .language-r}



~~~
Error in query_genes(2, 97, 98): File ../data/mouse_genes_mgi.sqlite doesn't exist
~~~
{: .error}

The way we’ve set this up is a bit complicated, but it allows greatest flexibility on the part of the user. And for our own work, we like to have a local SQLite database, for rapid queries of SNPs and genes.

### SNP associations

Okay, now finally we get to the SNP associations. We have a large peak on chromosome 2, and we want to look at individual SNPs in the region of the locus.

Well, actually, we first need to find the location of the inferred QTL.  The peak LOD score on chromosome 2 occurs at 52.4 cM. But to find nearby SNPs, we really want to know the Mbp position. The calculations were only performed at the marker positions, and so we can use `max()`, giving both the `scan1()` output and the physical map, and then pull out the position from the results.


~~~
peak_Mbp <- max(out, DOex$pmap)$pos
~~~
{: .language-r}

The marker is at 97.5 Mbp. We’ll focus on a 2 Mbp interval centered at 97.5 Mbp.

We can pull out the variants in the 2 Mbp interval centered at 97.5 on chr 2 using the query function we defined above:


~~~
variants <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
~~~
{: .language-r}



~~~
Error in query_variants(2, peak_Mbp - 1, peak_Mbp + 1): File ../data/cc_variants.sqlite doesn't exist
~~~
{: .error}















