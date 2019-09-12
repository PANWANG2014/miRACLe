# miRACLe: improving the prediction of miRNA-mRNA interactions by a random contact model
<br>

## Table of contents
1. <a href="#1">Introduction</a><br>
2. <a href="#2">Executing miRACLe</a><br>
	2.1 <a href="#3">Files required</a><br>
	2.2 <a href="#4">Script Execution</a><br>
3. <a href="#5">Benchmarking evaluations</a><br>
4. <a href="#6">References</a>



---
### <a name="1">1. Introduction</a>
**miRACLe (<u>miR</u>NA <u>A</u>nalysis by a <u>C</u>ontact mode<u>L</u>)** is a newly developed miRNA target prediction tool. It combines genome-wide expression profiles and the cumulative weighted context++ score from TargetScan in a random contact model, and then infers miRNA-mRNA interactions (MMIs) by the relative probability of effective contacts. Evaluation by a variety of measures shows that miRACLe consistently outperforms state-of-the-art methods in prediction accuracy, regulatory potential and biological relevance while has a distinct feature of inferring individual-specific miRNA targets. Empirical test suggests that on a laptop Intel Core i7-4712HQ personal computer with a 2.30 GHz CPU and 16 GB of RAM, our source code implementation requires less than 30 seconds of CPU time to complete the prediction for one sample. Importantly, we show that our model can also be applied to other sequence-based algorithms to improve their predictive power, such as DIANA-microT-CDS, miRanda-mirSVR and MirTarget4.  <br>


---
### <a name="2">2. Executing miRACLe</a>
####  <a name="3">2.1 Files required</a>
In order to run the current version of miRACLe, the users should provide two data files that describe the expression levels of each miRNA and mRNA for the same sample. And one additional file that defines the correspondence of samples between the miRNA and mRNA data files. All files are tab-delimited ASCII text files and must comply with the following specifications:

1. **Input miRNA expression file** is organized as follows:<br>

	| miRNA | TCGA-05-4384-01A-01T-1754-13 | TCGA-05-4390-01A-02T-1754-13 | TCGA-05-4396-01A-21H-1857-13 | TCGA-50-5066-01A-01T-1627-13|
	| :-------------: |:-------------:| :-----:| :-----:|:-----:|
	| hsa-let-7a-5p | 19.0144 | 16.2421 | 19.2817 | 18.0721 |
	| hsa-let-7a-3p | 7.31298 | 6.2094 | 7.8392 | 6.2667|
	| hsa-let-7a-2-3p | 6.5235 | 5.4594 | 3.7004 | 7.3837 |
	| hsa-let-7b-5p | 16.9613 | 15.5496 | 17.8444 | 16.9950 |
	| hsa-let-7b-3p | 7.9248 | 5.2094 | 7.6653 | 6.8201 |

	The first line contains the labels Name followed by the identifiers for each sample in the dataset. <br>
	>Line format: `Name(tab)(sample 1 name)(tab)(sample 2 name) (tab) ... (sample N name)`<br>
	>Example: `miRNAName	sample_1	sample_2	...	sample_n`<br>

	The remainder of the file contains data for each of the miRNAs. There is one line for each miRNA. Each line contains the miRNA name and a value for each sample in the dataset.<br>

2. **Input mRNA expression file** is organized as follows:<br>

	| Gene | TCGA-05-4384-01 | TCGA-05-4390-01 | TCGA-05-4396-01 | TCGA-50-5066-01|
	| :-------------: |:-------------:| :-----:| :-----:|:-----:|
	| AARS | 10.7094 | 11.6932 | 12.4282 | 11.0464 |
	| AASDHPPT | 9.9081 | 9.6716 | 10.1113 | 9.98328 |
	| AASDH | 7.9471 | 7.2897 | 8.3216 | 7.6274 |
	| AASS | 9.9649 | 7.7752 | 9.1723 | 5.9506 |
	| AATF | 9.9525 | 9.5380 | 9.3670 | 8.4375 |

	The first line contains the labels Name followed by the identifiers for each sample in the dataset. <br>
	>Line format: `Name(tab)(sample 1 name)(tab)(sample 2 name) (tab) ... (sample N name)`<br>
	>Example: `GeneName	sample_a	sample_b	...	sample_m`<br>

	The remainder of the file contains data for each of the mRNAs. There is one line for each mRNA. Each line contains the mRNA name and a value for each sample in the dataset.<br>

    **Note that** the input miRNA/mRNA expression file should be transformed into a non-negative matrix, in order for the main program to execute correctly. Both microarray profiling and RNA sequencing data are accepted as input. To achieve optimal prediction on the sequencing data, we strongly recommend that users provide log2 transformed normalized counts (e.g. RSEM or RPM) as the input for our program.

3. **Sample matching file** generally contains two columns, which shows the corresponding relationship of the sample identifiers in miRNA expression file and mRNA expression file (miRNA must be the first column and mRNA must be the second column). It also serves as a index to denote which samples we choose to analyze. It is organized as follows:<br>

	| miRNA | Gene | 
	| :-------------: |:-------------:| 
	| TCGA-50-5066-01A-01T-1627-13 | TCGA-50-5066-01 |
	| TCGA-05-4384-01A-01T-1754-13 | TCGA-05-4384-01 |
	| TCGA-05-4390-01A-02T-1754-13 | TCGA-05-4390-01 |
	| TCGA-05-4396-01A-21H-1857-13 | TCGA-05-4396-01 |

	The first line must contain the label Names for samples in each expression dataset with the first column for miRNA and second column for mRNA. <br>
	>Line format: `(sample name in miRNA file)(tab)(sample name in mRNA file)`<br>
	>Example: `sample_1	sample_a`<br>

	The remainder of the file contains sample identifiers used in the miRNA and mRNA expression files. There is one line for each sample. Each line contains the identifiers for that sample.<br>

#### <a name="4">2.2 Script Execution</a><br>
miRACLe is written in R and can be downloaded [here](https://github.com/PANWANG2014/miRACLe-Research/tree/master/Run%20miRACLe) along with test datasets. The source code of miRACLe consists of three parts, namely, 'FUNCTIONS', 'DATA INPUT' and 'MAIN CODE'. The main function "miracle" in "MAIN PROGRAM" calculates the miracle score for each miRNA-mRNA pair at individual and population levels, based on which all putative MMIs are ranked. The essential inputs that the miRACLe algorithm requires to run includes two parts:<br>

The first part contains the sequence-based interaction scores (seqScore) for putative miRNA-mRNA pairs. These scores are originally obtained from TargetSan v7.2 (TargetScan7\_CWCS\_cons and TargetScan7\_CWCS), DIANA-microT-CDS (DIANA\_microT\_CDS), MirTarget v4 (MirTarget4), miRanda-mirSVR (miRanda\_mirSVR) and compiled by the developers to fit the model. Default is **TargetScan7\_CWCS\_cons**. The other scores can be downloaded [here](https://figshare.com/s/0b7c68cd5152da27a191).<br>
	
> seqScore = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head = TRUE, sep = "\t"))<br>

User can also provide their own sequence matching scores, as long as the format of input file meets the requirements. Specifically, the first line must contain the label Names for mRNAs, miRNAs and their associated interaction scores. The remainder of the file contains RNA identifiers corresponding to those used in the expression files and the scores for each miRNA-mRNA pair. Note that the first column must contain identifiers for mRNAs, the second column must contain identifiers for miRNAs with the third column containing the associated scores.<br>

The second part contains paired miRNA-mRNA expression profiles and should be provided by the users. 

>sampleMatch = as.matrix(read.table("Test\_sampleMatch.txt", head = TRUE, sep = "\t"))<br>
>mirExpr = as.matrix(read.table("Test\_miRNA\_expression.txt", head = FALSE, sep = "\t"))<br>
>tarExpr = as.matrix(read.table("Test\_mRNA\_expression.txt", head = FALSE, sep = "\t"))<br>	

The "miracle" function also provides three optional parameters for users, which are: exprFilter (filter of expression profile, miRNAs/mRNAs that are not expressed in more than a given percentage of samples will be removed, default is 1), samSelect (sample selection, users can select a subset of all samples to analyze, default is no selection applied) and OutputSelect (logical variable, select “TRUE” to return the top 10 percent-ranked predictions by scores, and “FALSE” to return the whole prediction result. Default is TRUE).

>miracle(seqScore, sampleMatch, mirExpr, tarExpr)	#default<br>
>miracle(seqScore, sampleMatch, mirExpr, tarExpr, exprFilter = 1, samSelect, OutputSelect = TRUE) #optional parameters added


We also provide an [**R package**](https://github.com/PANWANG2014/miRACLe/tree/master/miRACLe/Sequence%20scores) of the algorithm for ease of use.

---
### <a name="5">3. Benchmarking evaluations</a><br>
1. The codes to reproduce the benchmarking evaluations are written in R.<br> 
2. Generally, all these codes are arranged into three parts as 'FUNCTIONS', 'INPUT DATA' and 'MAIN CODE'. The users need to download and fill in the relevant input files before implementing corresponding analyses.<br>
3. Files required for the reproduction of the evaluations can be broadly classified into three categories:<br>

* Sequence-based predictions (including seqScores for integrative methods)<br>

	| Data file | Description | 
	|:-------------:|:-------------| 
	| TargetScan7\_CWCS\_cons.txt | cumulative weighted context++ scores for conserved targets sites of conserved miRNA families obtained from TargetScan v7.2 |
	| TargetScan7\_CWCS.txt | cumulative weighted context++ scores for all miRNA-mRNA pairs obtained from TargetScan v7.2 |
    | TargetScan7\_qMRE\_cons.txt | number of conserved target sites of conserved miRNA families obtained from TargetScan v7.2 |
	| TargetScan7\_qMRE.txt | number of target sites for all miRNA-mRNA pairs obtained from TargetScan v7.2 |
	| DIANA\_microT\_CDS.txt | human interactions with miTG scores greater than 0.7 obtained from DIANA\-microT\-CDS|
	| miRanda\_mirSVR.txt | human conserved miRNA predictions with good mirSVR score obained from miRanda\-mirSVR|
    | miRmap.txt | predictions from miRmap|
    | miRTar2GO.txt | predictions from the “Highly sensitive” prediction set of miRTar2GO|    
    | miRTar2GO\_HeLa.txt | predictions in HeLa cells from the “Highly sensitive” prediction set of miRTar2GO|
    | MirTarget4.txt | human predictions obtained from miRDB v6.0|
    | miRWalk3.txt | human predictions restricted to 3`UTR obtained from miRWalk v3.0|
    | PITA.txt | the top human predictions with 3/15 flank obtained from PITA|
    | Combine_MMIs.txt | combined predictions from DIANA\-microT\-CDS, miRanda\-mirSVR, MirTarget4, PITA and TargetScan7.CWCS|
    | Symbol\_to\_ID.txt | paired gene symbols and gene entrez IDs downloaded from [HGNC](https://www.genenames.org/download/custom/)|

    These predictions are provided in a compressed file [Sequence\_based\_predictions.7z](https://figshare.com/s/82ffa5da58faf080230e).<br> 

* Input expression data files (mirExpr & tarExpr)<br>

	| Data file| Descriptions |
	|:-------------: |:-------------|
	| HeLa expression data | normalized microarray/RNA-Seq expression data for HeLa cell line |
    | NCI60 data | normalized microarray data for 59 NCI-60 cancer cell lines |
    | TCGA data | log2-transformed RPM/RSEM data for 7991 cancer patients from 32 TCGA cancer types |
	| MCC data | normalized microarray data for 68 tumor tissues and 21 normal tissues |

	These expression data files are provided along with relevant source codes except that the TCGA expression data files are provided in a compressed file [TCGA\_data.7z](https://figshare.com/s/045e07fb7c8278b9b3c2).<br> 

* Validation data (Reference data)<br>
	* Experimentally validated MMIs<br>

	| Data file | Description | validated MMI counts |
	|:-------------:|:-------------|:-----:|
	| Vset\_HeLa.txt | MMIs that are validated in HeLa cells from TarBase v8.0 | 34,263 |
    | Vset\_celllines.txt | MMIs that are validated in cell lines from TarBase v8.0 | 349,726 |
	| Vset\_all.txt | validated MMIs obtained from TarBase v8.0 | 376,205 |
    | Vset\_hc.txt | high-confidence set compiled from TarBase v8.0, miRTarbase v7.0, miRecords and oncomirDB | 10,575 |
	
    * Curated miRNA transfection experiments<br>

    | Data file | Description |
	|:-------------:|:-------------|
	| Transet\_HeLa\_Array.txt | Unified dataset of 5 miRNA transfections in HeLa cell line in which gene exrpession changes are measured by microarray | 
    | Transet\_HeLa\_Seq.txt | Unified dataset of 25 miRNA transfections in HeLa cell line in which gene exrpession changes are measured by RNA-Seq |
    | Transet\_multi.txt | Unified dataset of 105 non\-redundant miRNA transfections that are originally collected from 77 human cell lines or tissues |


	* Known cancer genes<br>

	| Data file | Description | Molecule counts |
	|:-------------:|:-------------|:-----:|
    | Cancer\_gene\_set | cancer genes obtained from cancer gene census | 723 |<br>

    These reference data files are provided along with relevant source codes.

---
### <a name="6">4. References</a><br>

miRACLe: improving the prediction of miRNA-mRNA interactions by a random contact model (in preparation)

