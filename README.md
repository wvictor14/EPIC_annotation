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

There are two annotation files, one mapped to **hg19** and one mapped to **hg38/GR38**. They have the same
annotation information (columns), but the hg38 annotation is missing 237 probes, since some mappings
are lost from converting from hg19 to hg38.






# The process

## 1. Starting annotation

I started with the default annotations provided by Illumina. I used two files, the latest b4 annotation (*MethylationEPIC_v-1-0_B4.csv*), and the list of probes that are missing between b3 and b2 (*MethylationEPIC Missing Legacy CpG (v1.0_B3 vs. v1.0_B2) Annotations.csv*). Both can be found on the product files list from [Illumina's website](https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html).

Using the intersection of these two lists of probes, I used the provided genomic location (chromomsome and position) to map annotations to each cpg. Note that Illumina's provided annotations are based on hg19.

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
   <th style="text-align:left;"> ilmn_address_a_id </th>
   <th style="text-align:left;"> ilmn_allele_a_probe_seq </th>
   <th style="text-align:left;"> ilmn_address_b_id </th>
   <th style="text-align:left;"> ilmn_allele_b_probe_seq </th>
   <th style="text-align:left;"> ilmn_infinium_design_type </th>
   <th style="text-align:left;"> ilmn_next_base </th>
   <th style="text-align:left;"> ilmn_color_channel </th>
   <th style="text-align:left;"> ilmn_forward_sequence </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000029 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 53468112 </td>
   <td style="text-align:left;"> 0040668283 </td>
   <td style="text-align:left;"> AACTATACTAACRAAAAAATATCCAAAAAACACTAACRTATAAAAATTTC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> TTTTTTAGATAAGGATATCCAGGCGATGAGGAAGTTTTACTTCTGGGAACAGCCTGGATA[CG]AAACCTTCACACGTCAGTGTCTTTTGGACATTTTCTCGTCAGTACAGCCCTGTTGAATGT </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000103 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 73470186 </td>
   <td style="text-align:left;"> 0076703527 </td>
   <td style="text-align:left;"> AACTTAACTTACAACCAATTATCCCAACRATTTCAACTAATAATATTAAC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> AAATAAGTGGGAAAAAATTCACCCATTACTCATACTCACAAATGGGGTATGTATGGCCAA[CG]TCAACATTATTAGCTGAAATCGTTGGGACAATTGGCTGCAAGCCAAGCCAATGATGAAAC </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000109 </td>
   <td style="text-align:left;"> chr3 </td>
   <td style="text-align:right;"> 171916037 </td>
   <td style="text-align:left;"> 0039798525 </td>
   <td style="text-align:left;"> CAATACTAACAAACACATATACCCCCCCACAAATCTTAACTTCTAAATAC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> GCCTTAGTCCTGAATGAGCCATTTCTCTAAGAAGTCCTGGCTTCTTTTTTAATAGAGAAT[CG]TATTTAGAAGCCAAGATCTGTGGGGGGGTACATGTGCCTGTTAGTATTGCAGTTGTGCCT </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000155 </td>
   <td style="text-align:left;"> chr7 </td>
   <td style="text-align:right;"> 2590565 </td>
   <td style="text-align:left;"> 0022679178 </td>
   <td style="text-align:left;"> AATAAAAAACCACTACACCCAACCTAAACATAATAATTAAATATTCAAAC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> GTCAGACGGGTCACGTACCGACATCATACCAGGGTTCCATGGTGCCACATCTCACACATG[CG]CCTGAACACTCAATTATCATGCTTAGGTTGGGTGCAGTGGCTCTCCACTGTAACCCCAGC </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000158 </td>
   <td style="text-align:left;"> chr9 </td>
   <td style="text-align:right;"> 95010555 </td>
   <td style="text-align:left;"> 0033624521 </td>
   <td style="text-align:left;"> AACTAAACATTCCRTATTATTTACTTCAAACTATTCTCATTTTCCCATCC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> CATGCTTCAGTCATGCTACTGAGGATTACTGATGCGGCAGGGATGGTGATGGCACAGTAG[CG]GATGGGAAAATGAGAACAGCTTGAAGCAAATAATACGGAATGTTTAGTTTTGCCATTATA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000165 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 91194674 </td>
   <td style="text-align:left;"> 0008711382 </td>
   <td style="text-align:left;"> CAAAATCTATTAATACAATAACTTTTAATAAAACAACTAAAACACACATC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> CTAAGTGCAGTCAGGATCTGTTAGTACAGTGGCTTTTGATGGAACAGCTGAGGCACACAT[CG]CCCGTGGCATGGACTCCGGGGCCGAACGCTCACGACCAAGACTTTTGCCCTTTTGAAATG </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000221 </td>
   <td style="text-align:left;"> chr17 </td>
   <td style="text-align:right;"> 54534248 </td>
   <td style="text-align:left;"> 0093795143 </td>
   <td style="text-align:left;"> ATATTTTTTAAAATACATTCCAAAAAACACAAAATTACAAACCACAAACC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> GAATTTATACTGTATTTTTTAAAATACATTCCAGAAAGCACAAAATTACAAACCACAGGC[CG]CAAGCAGTCAGTCTCAAGAAGCCTGAAACACCTGTTCCATTCCTCGAACAAGTTTGTGAA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000236 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 42263294 </td>
   <td style="text-align:left;"> 0063613963 </td>
   <td style="text-align:left;"> TATAACRTCATATTAAAAAAAACRATCTAACCCACCAATTTATACATCAC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> CTCAGCGACAGTGTAGCGTCATGTTAGAGGAGACGATCTGACCCACCAGTTTGTACATCA[CG]TCCTGCATGTCCCACACCATTTTTTCATGACCTTGTAATATACTGGTCTCTGTGCTATAG </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000289 </td>
   <td style="text-align:left;"> chr14 </td>
   <td style="text-align:right;"> 69341139 </td>
   <td style="text-align:left;"> 0012747940 </td>
   <td style="text-align:left;"> ATCTACTATATTCATTTCTCCAATCTCATATCCATTTTAATATAAAAATC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> CAAGTGAGCTAGCAAACACACATGCACCAATGTGCCTTTTGACAAGAGTACCCCCTACCC[CG]ACTCCCACACCAAAATGGACATGAGATTGGAGAAATGAATACAGCAGATGGAACAGATAG </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000292 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 28890100 </td>
   <td style="text-align:left;"> 0094694990 </td>
   <td style="text-align:left;"> AAAACATTAATTACCAACCRCTCTTCCAAAAAACACTTACCATTAAAACC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> TGGGGTGAGTGAGACCACGGGCCTCACCCCGGACCAAGTTAAGCGGAATCTGGAGAAATA[CG]GCCTCAATGGTAAGTGTCCCTTGGAAGAGCGGCTGGTAATTAATGCCCTCCTGCACCCCC </td>
  </tr>
</tbody>
</table></div>

## 2. Annotate cpgs

I used the R package annotatr to access UCSC annotations for cpg islands and transcripts, and then FANTOM5 for enhancers.

### UCSC transcript and cpg island -related elements:

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
   <td style="text-align:left;"> ATP2A1, ATP2A1, ATP2A1, NA </td>
   <td style="text-align:left;"> uc002drp.1, uc002drn.1, uc002dro.1, uc010vct.2 </td>
   <td style="text-align:left;"> 4000, 302, 302, 931314 </td>
  </tr>
</tbody>
</table></div>

### Enhancers

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

### PMDs from Schroeder et al. 2013:

Taken from the primary article.

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> cpg </th>
   <th style="text-align:left;"> chr </th>
   <th style="text-align:right;"> start </th>
   <th style="text-align:left;"> pmd_id </th>
   <th style="text-align:right;"> pmd_width </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000103 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 73470186 </td>
   <td style="text-align:left;"> chr4:73116829-73652646 </td>
   <td style="text-align:right;"> 535817 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000622 </td>
   <td style="text-align:left;"> chr15 </td>
   <td style="text-align:right;"> 23034447 </td>
   <td style="text-align:left;"> chr15:22752148-23105710 </td>
   <td style="text-align:right;"> 353562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000658 </td>
   <td style="text-align:left;"> chr9 </td>
   <td style="text-align:right;"> 139997924 </td>
   <td style="text-align:left;"> chr9:139936610-140037104 </td>
   <td style="text-align:right;"> 100494 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000776 </td>
   <td style="text-align:left;"> chr4 </td>
   <td style="text-align:right;"> 156388205 </td>
   <td style="text-align:left;"> chr4:156349660-156517050 </td>
   <td style="text-align:right;"> 167390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000905 </td>
   <td style="text-align:left;"> chr15 </td>
   <td style="text-align:right;"> 59785306 </td>
   <td style="text-align:left;"> chr15:59309009-59907948 </td>
   <td style="text-align:right;"> 598939 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00000974 </td>
   <td style="text-align:left;"> chr20 </td>
   <td style="text-align:right;"> 6750606 </td>
   <td style="text-align:left;"> chr20:6698958-7897260 </td>
   <td style="text-align:right;"> 1198302 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00001099 </td>
   <td style="text-align:left;"> chr8 </td>
   <td style="text-align:right;"> 87081553 </td>
   <td style="text-align:left;"> chr8:86909966-87150768 </td>
   <td style="text-align:right;"> 240802 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00001364 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:right;"> 214170376 </td>
   <td style="text-align:left;"> chr1:214018414-214462834 </td>
   <td style="text-align:right;"> 444420 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00001506 </td>
   <td style="text-align:left;"> chr18 </td>
   <td style="text-align:right;"> 61050587 </td>
   <td style="text-align:left;"> chr18:59828768-61568754 </td>
   <td style="text-align:right;"> 1739986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00001520 </td>
   <td style="text-align:left;"> chr14 </td>
   <td style="text-align:right;"> 37666489 </td>
   <td style="text-align:left;"> chr14:37161928-37746904 </td>
   <td style="text-align:right;"> 584976 </td>
  </tr>
</tbody>
</table></div>

### Imprinting regions

Both general and placenta specific from GeneImprint and OTAGO databases, plus a court et al. and hanna et al. DMR imprinting papers.

These imprinting information was summarized and provided by GDG

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> cpg </th>
   <th style="text-align:left;"> chr </th>
   <th style="text-align:right;"> start </th>
   <th style="text-align:left;"> imprinted_gene_placenta </th>
   <th style="text-align:left;"> imprinted_gene_general </th>
   <th style="text-align:left;"> imprinted_dmr_general </th>
   <th style="text-align:left;"> imprinted_dmr_placenta </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cg00000924 </td>
   <td style="text-align:left;"> chr11 </td>
   <td style="text-align:right;"> 2720463 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> KCNQ1;KCNQ1OT1 </td>
   <td style="text-align:left;"> chr11:2720229-2722440 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00019511 </td>
   <td style="text-align:left;"> chr14 </td>
   <td style="text-align:right;"> 101201288 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> DLK1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00033209 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 3511787 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NAA60 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00034154 </td>
   <td style="text-align:left;"> chr19 </td>
   <td style="text-align:right;"> 10337163 </td>
   <td style="text-align:left;"> DNMT1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00036440 </td>
   <td style="text-align:left;"> chr16 </td>
   <td style="text-align:right;"> 3507875 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NAA60 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00047185 </td>
   <td style="text-align:left;"> chr19 </td>
   <td style="text-align:right;"> 54196849 </td>
   <td style="text-align:left;"> C19MC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00054741 </td>
   <td style="text-align:left;"> chr19 </td>
   <td style="text-align:right;"> 54201247 </td>
   <td style="text-align:left;"> C19MC </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00058515 </td>
   <td style="text-align:left;"> chr2 </td>
   <td style="text-align:right;"> 80527798 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> LRRTM1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00059930 </td>
   <td style="text-align:left;"> chr13 </td>
   <td style="text-align:right;"> 48894382 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> RB1 </td>
   <td style="text-align:left;"> chr13:48892341-48895763 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cg00068519 </td>
   <td style="text-align:left;"> chr7 </td>
   <td style="text-align:right;"> 94259997 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> SGCE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table></div>

## 3. Map to hg38

Lastly I mapped the annotation to the genome assembly hg38 using UCSC liftover's tool implemented in R. This results in a loss of 237 cpgs.

