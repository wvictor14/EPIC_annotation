---
title: "Annotating DNA methylation array"
date: "07/08/2019"
output:
    html_document:
      keep_md: yes
      theme: spacelab
      toc: yes
      toc_depth: 3
      toc_float:
        collapsed: no
editor_options: 
  chunk_output_type: console
---



This repository contains code for generating an annotation for the Illumina EPIC methylation array.

There are two annotation files, one mapped to **hg19** (`hg19_epic_annotation.rds`) and one mapped to **hg38/GR38** (`hg38_epic_annotation.rds`). They have the same annotation information (columns), but the hg38 annotation is missing 237 probes, since some mappings are lost from converting from hg19 to hg38.

Currently these files are not git tracked because they are too large (~250 mb).





# The process

## 1. Starting annotation

I started with the default annotations provided by Illumina. I used two files, the latest b4 annotation (*MethylationEPIC_v-1-0_B4.csv*), and the list of probes that are missing between b3 and b2 (*MethylationEPIC Missing Legacy CpG (v1.0_B3 vs. v1.0_B2) Annotations.csv*). Both can be found on the product files list from [Illumina's website](https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html).

Using the intersection of these two lists of probes, I used the provided genomic location (chromomsome and position) to map annotations to each cpg. Note that Illumina's provided annotations are based on hg19.

> an example of the starting coordinates from Illumina that this annotation is based on

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> cpg </th>
   <th style="text-align:left;"> chr </th>
   <th style="text-align:right;"> start </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000029 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 53468112 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000103 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 73470186 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000109 </td>
   <td style="text-align:left;"> chr3 </td>
   <td style="text-align:right;"> 171916037 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000155 </td>
   <td style="text-align:left;"> chr7 </td>
   <td style="text-align:right;"> 2590565 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000158 </td>
   <td style="text-align:left;"> chr9 </td>
   <td style="text-align:right;"> 95010555 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000165 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 91194674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000221 </td>
   <td style="text-align:left;"> chr17 </td>
   <td style="text-align:right;"> 54534248 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000236 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 42263294 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000289 </td>
   <td style="text-align:left;"> chr14 </td>
   <td style="text-align:right;"> 69341139 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000292 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 28890100 </td>
  </tr>
</tbody>
</table>

I also kept some probe-specific information that I thought some may find useful. The columns for these variables are all prefixed with "ilmn_".

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> cpg </th>
   <th style="text-align:left;"> chr </th>
   <th style="text-align:right;"> start </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000029 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 53468112 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000103 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 73470186 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000109 </td>
   <td style="text-align:left;"> chr3 </td>
   <td style="text-align:right;"> 171916037 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000155 </td>
   <td style="text-align:left;"> chr7 </td>
   <td style="text-align:right;"> 2590565 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000158 </td>
   <td style="text-align:left;"> chr9 </td>
   <td style="text-align:right;"> 95010555 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000165 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 91194674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000221 </td>
   <td style="text-align:left;"> chr17 </td>
   <td style="text-align:right;"> 54534248 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000236 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 42263294 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000289 </td>
   <td style="text-align:left;"> chr14 </td>
   <td style="text-align:right;"> 69341139 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000292 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 28890100 </td>
  </tr>
</tbody>
</table></div>

## 2. Annotate cpgs

### Transcript-related features, enhancers, cpg islands

I used the R package `annotatr` to access UCSC annotations for cpg islands / transcripts, and FANTOM5 for enhancers.

> UCSC transcript and cpg island -related elements:

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> cpg </th>
   <th style="text-align:left;"> chr </th>
   <th style="text-align:right;"> start </th>
   <th style="text-align:left;"> cpg_id </th>
   <th style="text-align:left;"> cpg_width </th>
   <th style="text-align:left;"> genes_id </th>
   <th style="text-align:left;"> genes_symbol </th>
   <th style="text-align:left;"> genes_tx_id </th>
   <th style="text-align:left;"> genes_width </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000029 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 53468112 </td>
   <td style="text-align:left;"> shore </td>
   <td style="text-align:left;"> 2000 </td>
   <td style="text-align:left;"> promoter, 1to5kb </td>
   <td style="text-align:left;"> RBL2, RBL2 </td>
   <td style="text-align:left;"> uc002ehi.4, uc010vgv.1 </td>
   <td style="text-align:left;"> 1000, 4000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000103 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 73470186 </td>
   <td style="text-align:left;"> sea </td>
   <td style="text-align:left;"> 491623 </td>
   <td style="text-align:left;"> intergenic </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 480899 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000109 </td>
   <td style="text-align:left;"> chr3 </td>
   <td style="text-align:right;"> 171916037 </td>
   <td style="text-align:left;"> sea </td>
   <td style="text-align:left;"> 398648 </td>
   <td style="text-align:left;"> intron, intron, intron </td>
   <td style="text-align:left;"> FNDC3B, FNDC3B, FNDC3B </td>
   <td style="text-align:left;"> uc003fhy.3, uc003fhz.4, uc003fia.3 </td>
   <td style="text-align:left;"> 93324, 93324, 93324 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000155 </td>
   <td style="text-align:left;"> chr7 </td>
   <td style="text-align:right;"> 2590565 </td>
   <td style="text-align:left;"> sea </td>
   <td style="text-align:left;"> 3182 </td>
   <td style="text-align:left;"> intron, intron </td>
   <td style="text-align:left;"> BRAT1, BRAT1 </td>
   <td style="text-align:left;"> uc003smi.3, uc003smj.2 </td>
   <td style="text-align:left;"> 6826, 6826 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000158 </td>
   <td style="text-align:left;"> chr9 </td>
   <td style="text-align:right;"> 95010555 </td>
   <td style="text-align:left;"> sea </td>
   <td style="text-align:left;"> 143935 </td>
   <td style="text-align:left;"> intron, intron, intron, intron, intron </td>
   <td style="text-align:left;"> IARS, IARS, IARS, IARS, IARS </td>
   <td style="text-align:left;"> uc004ars.2, uc004art.2, uc004aru.4, uc010mqr.3, uc010mqt.2 </td>
   <td style="text-align:left;"> 2306, 2306, 2306, 2306, 2306 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000165 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 91194674 </td>
   <td style="text-align:left;"> shore </td>
   <td style="text-align:left;"> 2000 </td>
   <td style="text-align:left;"> intergenic </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 107309 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000221 </td>
   <td style="text-align:left;"> chr17 </td>
   <td style="text-align:right;"> 54534248 </td>
   <td style="text-align:left;"> sea </td>
   <td style="text-align:left;"> 656815 </td>
   <td style="text-align:left;"> exon, intronexonboundary </td>
   <td style="text-align:left;"> ANKFN1, ANKFN1 </td>
   <td style="text-align:left;"> uc002iun.1, uc002iun.1 </td>
   <td style="text-align:left;"> 100, 400 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000236 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 42263294 </td>
   <td style="text-align:left;"> sea </td>
   <td style="text-align:left;"> 10587 </td>
   <td style="text-align:left;"> exon, exon, exon, 3UTR, 3UTR </td>
   <td style="text-align:left;"> VDAC3, VDAC3, VDAC3, VDAC3, VDAC3 </td>
   <td style="text-align:left;"> uc003xpc.3, uc031tay.1, uc022aul.1, uc003xpc.3, uc022aul.1 </td>
   <td style="text-align:left;"> 567, 567, 567, 475, 475 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000289 </td>
   <td style="text-align:left;"> chr14 </td>
   <td style="text-align:right;"> 69341139 </td>
   <td style="text-align:left;"> shore </td>
   <td style="text-align:left;"> 2000 </td>
   <td style="text-align:left;"> exon, exon, exon, exon, exon, 3UTR, 3UTR, 3UTR, 3UTR, 3UTR </td>
   <td style="text-align:left;"> ACTN1, ACTN1, ACTN1, ACTN1, ACTN1, ACTN1, ACTN1, ACTN1, ACTN1, ACTN1 </td>
   <td style="text-align:left;"> uc001xkk.3, uc010ttb.2, uc001xkl.3, uc001xkm.3, uc001xkn.3, uc001xkk.3, uc010ttb.2, uc001xkl.3, uc001xkm.3, uc001xkn.3 </td>
   <td style="text-align:left;"> 895, 895, 895, 895, 895, 736, 736, 736, 736, 736 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000292 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 28890100 </td>
   <td style="text-align:left;"> shore </td>
   <td style="text-align:left;"> 2000 </td>
   <td style="text-align:left;"> 1to5kb, exon, exon, intron </td>
   <td style="text-align:left;"> ATP2A1, ATP2A1, ATP2A1, no_associated_gene </td>
   <td style="text-align:left;"> uc002drp.1, uc002drn.1, uc002dro.1, uc010vct.2 </td>
   <td style="text-align:left;"> 4000, 302, 302, 931314 </td>
  </tr>
</tbody>
</table></div>

> Enhancers

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> cpg </th>
   <th style="text-align:left;"> chr </th>
   <th style="text-align:right;"> start </th>
   <th style="text-align:left;"> enhancers_id </th>
   <th style="text-align:left;"> enhancers_width </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000776 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 156388205 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 116 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00003578 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 12600529 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 328 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00004667 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 16292746 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 536 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00004963 </td>
   <td style="text-align:left;"> chr6 </td>
   <td style="text-align:right;"> 147124996 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 324 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00005325 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 201684967 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 354 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00005461 </td>
   <td style="text-align:left;"> chr3 </td>
   <td style="text-align:right;"> 46131480 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 363 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00007021 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 101819246 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 437 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00007969 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 41633437 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 488 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00009088 </td>
   <td style="text-align:left;"> chr11 </td>
   <td style="text-align:right;"> 60930188 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00009585 </td>
   <td style="text-align:left;"> chr15 </td>
   <td style="text-align:right;"> 33111077 </td>
   <td style="text-align:left;"> enhancer </td>
   <td style="text-align:left;"> 345 </td>
  </tr>
</tbody>
</table></div>

### Placental partially methylated domains (PMDs) from Schroeder et al. 2013:

Taken from the primary article.

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> cpg </th>
   <th style="text-align:left;"> chr </th>
   <th style="text-align:right;"> start </th>
   <th style="text-align:right;"> pmd_width </th>
   <th style="text-align:left;"> pmd_id </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000103 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 73470186 </td>
   <td style="text-align:right;"> 332252 </td>
   <td style="text-align:left;"> chr4:73435322-73767574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000165 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 91194674 </td>
   <td style="text-align:right;"> 81136 </td>
   <td style="text-align:left;"> chr1:91192805-91273941 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000363 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 230560793 </td>
   <td style="text-align:right;"> 68156 </td>
   <td style="text-align:left;"> chr1:230492946-230561102 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000596 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 133098502 </td>
   <td style="text-align:right;"> 77607 </td>
   <td style="text-align:left;"> chr8:133063957-133141564 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000776 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 156388205 </td>
   <td style="text-align:right;"> 162183 </td>
   <td style="text-align:left;"> chr4:156298095-156460278 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000884 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 154609857 </td>
   <td style="text-align:right;"> 74720 </td>
   <td style="text-align:left;"> chr4:154606053-154680773 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000974 </td>
   <td style="text-align:left;"> chr20 </td>
   <td style="text-align:right;"> 6750606 </td>
   <td style="text-align:right;"> 1147 </td>
   <td style="text-align:left;"> chr20:6749547-6750694 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00001099 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 87081553 </td>
   <td style="text-align:right;"> 201811 </td>
   <td style="text-align:left;"> chr8:86879841-87081652 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00001249 </td>
   <td style="text-align:left;"> chr14 </td>
   <td style="text-align:right;"> 60389786 </td>
   <td style="text-align:right;"> 171588 </td>
   <td style="text-align:left;"> chr14:60386751-60558339 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00001520 </td>
   <td style="text-align:left;"> chr14 </td>
   <td style="text-align:right;"> 37666489 </td>
   <td style="text-align:right;"> 24805 </td>
   <td style="text-align:left;"> chr14:37641880-37666685 </td>
  </tr>
</tbody>
</table></div>

### Imprinting regions

These placental imprinted regions were collected from several sources. The merging of these regions into a combined resource is documented at [github.com/wvictor14/human_methylation_imprints](github.com/wvictor14/human_methylation_imprints).

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> cpg </th>
   <th style="text-align:left;"> chr </th>
   <th style="text-align:right;"> start </th>
   <th style="text-align:left;"> imprint_tissue_specificity </th>
   <th style="text-align:left;"> imprint_methylated_allele </th>
   <th style="text-align:left;"> imprint_sources </th>
   <th style="text-align:left;"> imprint_region </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000924 </td>
   <td style="text-align:left;"> chr11 </td>
   <td style="text-align:right;"> 2720463 </td>
   <td style="text-align:left;"> other </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Court 2014, Hanna 2016 </td>
   <td style="text-align:left;"> 11:2719948-2722440 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00050654 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 4576493 </td>
   <td style="text-align:left;"> placental-specific </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Sanchez-Delgado 2016 </td>
   <td style="text-align:left;"> 4:4576220-4577911 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00059930 </td>
   <td style="text-align:left;"> chr13 </td>
   <td style="text-align:right;"> 48894382 </td>
   <td style="text-align:left;"> other </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Court 2014 </td>
   <td style="text-align:left;"> 13:48892341-48895763 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00082664 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 154710796 </td>
   <td style="text-align:left;"> placental-specific </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Sanchez-Delgado 2016, Hamada 2016 </td>
   <td style="text-align:left;"> 4:154709200-154715220 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00082664 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 154710796 </td>
   <td style="text-align:left;"> placental-specific </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Sanchez-Delgado 2016, Hamada 2016 </td>
   <td style="text-align:left;"> 4:154709200-154715220 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00083059 </td>
   <td style="text-align:left;"> chr6 </td>
   <td style="text-align:right;"> 39902348 </td>
   <td style="text-align:left;"> placental-specific </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Hanna 2016 </td>
   <td style="text-align:left;"> 6:39901897-39902693 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00096536 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 154711906 </td>
   <td style="text-align:left;"> placental-specific </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Sanchez-Delgado 2016, Hamada 2016 </td>
   <td style="text-align:left;"> 4:154709200-154715220 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00096536 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 154711906 </td>
   <td style="text-align:left;"> placental-specific </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Sanchez-Delgado 2016, Hamada 2016 </td>
   <td style="text-align:left;"> 4:154709200-154715220 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00098799 </td>
   <td style="text-align:left;"> chr15 </td>
   <td style="text-align:right;"> 99409360 </td>
   <td style="text-align:left;"> other </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Court 2014 </td>
   <td style="text-align:left;"> 15:99408496-99409650 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00155882 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 141110747 </td>
   <td style="text-align:left;"> other </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> Court 2014, Hanna 2016 </td>
   <td style="text-align:left;"> 8:141107717-141111081 </td>
  </tr>
</tbody>
</table></div>

## 3. Map to hg38

Lastly I mapped the annotation to the genome assembly hg38 using UCSC liftover's tool implemented in R. This results in a loss of 237 cpgs.

