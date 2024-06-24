# Dynamic Selection Preference of SARS-CoV-2
This repository provides data and scripts used for the research about dynamic selection preference.
The workflow contains scripts written in C#  and Python and one needs to adjust the file path manually, the [Microsoft VisualStudio](https://visualstudio.microsoft.com/zh-hans/) or other compatible IDEs are needed.

The weekly cases, SARS-CoV-2 sequence mutation, and metadata files were saved under _./Data_
The processed DMS data files were saved under _./Data/DMS_
The scripts used for the research were saved under _./Scripts_
The sequence and mutation data were saved under _./Data_

## **Calculation of weekly cases, haplotype distribution, and mutation frequency matrix**
The weekly cases data were obtained from the global weekly COVID-19 cases reported by the World Health Organization on April 18, 2024 (https://data.who.int/dashboards/covid19/cases). 

The [mutations](https://download.cncb.ac.cn/GVM/Coronavirus/vcf/Mutation.gz ) of SARS-CoV-2 sequences and corresponding [metadata](https://ngdc.cncb.ac.cn/ncov/genome/export/meta) were retrieved from the RCoV19 database.

Calculation of weekly cases, haplotype distribution, and mutation frequency matrix was based on script **CalWeeklyHaplotype_Cases**

The result weekly cases files are saved under _./VOCCases_
The result haplotype diversity and weekly cases files are saved under _./Diversity_
The result weekly mutation frequency files are saved under _./WeeklyMatrix_
The result haplotype files are saved under _./WeeklyMax/Haplotypes_

## **Construction of phylogenetic tree**

The maximum parsimony tree of the eight major clades was inferred using the spike amino acid sequence. The consensus sequence of each clade was obtained from the outbreak.info website with a minimum mutation frequency of 75%. Evolutionary analyses were conducted in [MEGA X (v10.1.8)](https://www.megasoftware.net/).

The result tree files are saved under _./VOCTree_

## ****Calculation of the functional indexes****

The immune escape score and ACE2 binding affinity were calculated based on DMS data following the same methods as in our [previous research](). 
The raw data can be accessed from the following links:
[RBD_ACE2binding_Wuhan-Hu-1](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants/blob/main/results/final_variant_scores/final_variant_scores.csv)
[RBD_ACE2binding_Alpha](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants/blob/main/results/final_variant_scores/final_variant_scores.csv)
[RBD_ACE2binding_Delta](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants/blob/main/results/final_variant_scores/final_variant_scores.csv)
[RBD_ACE2binding_BA.1](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/blob/main/results/final_variant_scores/final_variant_scores.csv)
[RBD_ACE2binding_BA.2](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/blob/main/results/final_variant_scores/final_variant_scores.csv)
[RBD_ACE2binding_BQ.1.1](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ/blob/main/results/final_variant_scores/final_variant_scores.csv)
[RBD_ACE2binding_XBB.1.5](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ/blob/main/results/final_variant_scores/final_variant_scores.csv)

[RBD_Antibodies_Wuhan-Hu-1](https://github.com/jianfcpku/convergent_RBD_evolution/blob/main/antibody_info.xlsx)
[RBD_Antibodies_BA.1](https://github.com/jianfcpku/convergent_RBD_evolution/blob/main/antibody_info.xlsx)
[RBD_Antibodies_BA.2](https://github.com/jianfcpku/convergent_RBD_evolution/blob/main/antibody_info.xlsx)
[RBD_Antibodies_BA.5](https://github.com/jianfcpku/SARS-CoV-2-reinfection-DMS/blob/main/antibody_info.csv)

[RBD_EscapeScore_Wuhan-Hu-1](https://github.com/jianfcpku/convergent_RBD_evolution/blob/main/use_res_clean.csv)
[RBD_EscapeScore_BA.1](https://github.com/jianfcpku/convergent_RBD_evolution/blob/main/use_res_clean.csv)
[RBD_EscapeScore_BA.2](https://github.com/jianfcpku/convergent_RBD_evolution/blob/main/use_res_clean.csv)
[RBD_EscapeScore_BA.5](https://github.com/jianfcpku/SARS-CoV-2-reinfection-DMS/blob/main/antibody_dms_merge.csv.gz)

The spike-mediated cell entry index was obtained from Bloom's research.
The raw data can be accessed from the following links:
[Spike_CellEntry_Delta](https://github.com/dms-vep/SARS-CoV-2_Delta_spike_DMS_REGN10933/blob/main/results/muteffects_functional/muteffects_latent.csv)
[Spike_CellEntry_BA.1](https://github.com/dms-vep/SARS-CoV-2_Omicron_BA.1_spike_DMS_mAbs/tree/main/results/muteffects_functional)
[Spike_CellEntry_BA.2](https://github.com/dms-vep/SARS-CoV-2_Omicron_BA.2_spike_ACE2_binding/blob/main/results/summaries/summary.csv)
[Spike_CellEntry_XBB.1.5](https://github.com/dms-vep/SARS-CoV-2_XBB.1.5_spike_DMS/blob/main/results/summaries/summary.csv)

The antigenic distance index was simplified from the formula of the prediction model [EVEscape](https://evescape.org/), we used the product of accessibility and dissimilarity as the antigenic distance index.

## **Tracking the evolution trajectories of viral population.**

The weekly frequency changes of all mutations were determined using the first week as a reference, based on the script **CalMatrixDistance**. The script also calculates the relative functional indexes of each haplotype.

The result weekly distance files are saved under _./WeeklyMatrix/WeeklyDistance.*.tsv_
The result haplotype distance (to other sequences) files are saved under _./Diversity/HaploSeqs/*.Distance_
The result haplotype distance (to sequences from predecessor) files are saved under _./Diversity/HaploSeqs/*.PreVOCT1Dis_

## **Calculation of the selection preference based on haplotypes.**

The relative growth advantages of each haplotype were calculated using the code provided in [CoV-Spectrum](https://www.sciencedirect.com/science/article/pii/S1755436521000335?via%3Dihub), using all haplotype distance files from the previous step.

The result haplotype distance and fitness files are saved under _./Diversity/HaploSeqs/*.WithGrowthRate_

The Pearson correlation was calculated using [stat_cor](https://www.rdocumentation.org/packages/ggpubr/versions/0.6.0/topics/stat_cor) in R.

## **The prediction of haplotype relative growth advantages.**

The model training and prediction were based on the script **AiDMS.R**, using the result haplotype distance and fitness files from the previous step.
