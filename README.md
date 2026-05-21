# RLE retrotransposon mining

In this project, we develop four core pipelines for RLE-related retrotransposon mining and integration site detection.
- 1. **mining.search_RLE_locus**: search RLE-related retrotransposon loci from genomes.
- 2. **mining.build_consensus**: build consensus sequence for multicopy RLE loci with MSA files.
- 3. **mapping.detect_insert_site**: detect insertion sites from mapping NGS data.
- 4. **mapping.enrichment**: calculate and visualize the enrichment status of insertion sites across multiple dimensions, including gene-relative position, chromatin state, and safe harbor features.

# Requirements

## Python
- Python 3.9
    ```bash
    conda env create -f environment.yaml
    conda activate RLE_mining
    ```

## Other software
- openjdk (22)
- picard (3.1.1)
- mafft (v7.310)
- blat (v36x2)

# Usage

## 1. mining.search_RLE_locus
Iterate through the genome files specified by `-genome`, searching for ORFs (open reading frames) and their corresponding peptides. Annotate the peptides using the `-hmm` model, and retain only those ORFs that contain both the RVT_1 and RLE domains. Then, perform a preliminary family classification of these ORFs using the MSA and family files specified by `-RVTmsa` and `-RVTfamily`.

```txt
$  python src/mining.search_RLE_locus/search_RLE_retrotransposon_locus.py -h
usage: search_RLE_retrotransposon_locus.py [-h] -genome GENOME_FILE -hmm RVT_HMM -output OUTPUT_DIR -log
                                           LOG_DIR

optional arguments:
  -h, --help            show this help message and exit
  -genome GENOME_FILE, --genome_file GENOME_FILE
                        CSV describing genome sources (tab-separated columns).
  -hmm RVT_HMM, --RVT_hmm RVT_HMM
                        HMM reference
  -output OUTPUT_DIR, --output_dir OUTPUT_DIR
                        directory to save results
  -log LOG_DIR, --log_dir LOG_DIR
                        directory to save log files
```

## 2. mining.build_consensus
Build consensus sequences for multicopy RLE loci using the MSA files specified by `-msa-list`.

```txt
$ python src/mining.build_consensus/build_consensus_for_multicopy_RLE.py -h
usage: build_consensus_for_multicopy_RLE.py [-h] -msa-list MSA_LIST -output OUTPUT [-conserve] [-mismatch-size MISMATCH_SIZE] [-threads THREADS]

Build consensus from per-group locus MSA files.

optional arguments:
  -h, --help            show this help message and exit
  -msa-list MSA_LIST, --msa-list MSA_LIST
                        Two columns per line: group_ID then MSA path (tab or spaces; use tab if path contains spaces). Lines starting with '#' and empty lines are skipped.
  -output OUTPUT, --output OUTPUT
                        Directory to write <group_ID>.consensus.csv / .warning / .error
  -conserve, --conserve
                        Pass conserve=True into align_consensus_parse_RNA (default: False).
  -mismatch-size MISMATCH_SIZE, --mismatch-size MISMATCH_SIZE
                        mismatch_size for align_consensus_parse_RNA (default: 40).
  -threads THREADS, --threads THREADS
                        Worker processes for Pool (default: 60).
```

## 3. mapping.detect_insert_site
By detecting split alignment events in NGS data, reads that are mapped to both the genome and donor sequences are identified; these reads are considered as supporting evidence for transposition events. During the analysis, various filters are applied to ensure result accuracy, including read quality, the relative mapping position of read1/read2, and whether the primer region is covered.

```txt
$ python src/mapping.detect_insert_site/mapping_analysis_main.py -h
usage: mapping_analysis_main.py [-h] --fastq FASTQ --donor DONOR [--donor_label DONOR_LABEL] [--exist_index EXIST_INDEX] [--bam_ok {True,False}] [--primer_homo PRIMER_HOMO]
                                [--homo_primer HOMO_PRIMER] [--TSD TSD] [--TSD_len TSD_LEN] --UMI UMI --config CONFIG [--sequenator {Illumina,MGI}] --outpath OUTPATH --outlabel OUTLABEL
                                [--species {human,mouse}]

apping analysis for transposons and retrotransposons.

optional arguments:
  -h, --help            show this help message and exit
  --fastq FASTQ         NGS PE raw data; .R1.fq and .R2.fq paths are split by comma.
  --donor DONOR         donor plasmid (.dna).
  --donor_label DONOR_LABEL
                        donor plasmid label used in reference.fa. default same as --outlabel
  --exist_index EXIST_INDEX
                        Existing BWA index prefix.
  --bam_ok {True,False}
                        True: BAM creation steps are already complete. Default False.
  --primer_homo PRIMER_HOMO
                        Coordinates for primer-to-homo (right): primer_5' start, TE_3' end, Homo_3' end (0-based, half-open interval).
  --homo_primer HOMO_PRIMER
                        coordinate of primer to homo(Left), like primer_homo but reverse.
  --TSD TSD             expected TSD sequence, eg. TTAA, TA. Used when specific TSD sequence is known.
  --TSD_len TSD_LEN     [Deprecated] expected TSD length, used when <only> TSD length is known.
  --UMI UMI             Length of UMI. If no UMI is used, it will be set to 0.
  --config CONFIG       config file in yaml format.
  --sequenator {Illumina,MGI}
                        Sequencing platform used to generate NGS data (Illumina or MGI). [default MGI]
  --outpath OUTPATH     specify the output directory.
  --outlabel OUTLABEL   specify the output file label.
  --species {human,mouse}
                        specify the species of genome.
```

## 4. mapping.enrichment

**Step 1:** 
Annotate insertion sites with gene-relative position, chromatin state, and safe harbor features.

```txt
$ python src/mapping.enrichment/insert_site_annotation.py -h
usage: insert_site_annotation.py [-h] --mapping MAPPING_FILE --gtf GTF_FILE --cancergene CANCER_GENELIST --conserved CONSERVED_REGION --miRNA MIRNA --chrom_state CHROM_STATE --outpath
                                 OUTPATH --outlabel OUTLABEL [--thread THREAD] [--fai FAI]

optional arguments:
  -h, --help            show this help message and exit
  --mapping MAPPING_FILE
                        mapping result site, use 0-base coordinate.
  --gtf GTF_FILE        GTF file from UCSC trimmed by table format.
  --cancergene CANCER_GENELIST
                        Cancer gene list from OncoKB.
  --conserved CONSERVED_REGION
                        conserved region in genome from UCSC.
  --miRNA MIRNA         miRNA gene region from UCSC.
  --chrom_state CHROM_STATE
                        chromatin states from epigenomics project.
  --outpath OUTPATH     specify the output directory.
  --outlabel OUTLABEL   specify the output file label.
  --thread THREAD       Number of worker processes for multiprocessing pools.
  --fai FAI             related genome.fai file(samtools faidx..)
```

**Step 2:** 
Calculate the enrichment of insertion sites across multiple dimensions.

```txt
$ python src/mapping.enrichment/calculate_enrichment_ratio_and_show.py -h
usage: calculate_enrichment_ratio_and_show.py [-h] --mapping_result MAPPING_FILE --sample_label SAMPLE_LABEL --random RANDOM --outpath OUTPATH --outlabel OUTLABEL
                                              [--sample_order SAMPLE_ORDER]

Read categorized insertion-site tables produced by insert_site_annotation.py.

optional arguments:
  -h, --help            show this help message and exit
  --mapping_result MAPPING_FILE
                        result of mapping data with category
  --sample_label SAMPLE_LABEL
                        specify sample label matched with mapping_result.
  --random RANDOM       result of simulated data with category.
  --outpath OUTPATH     specify the output directory.
  --outlabel OUTLABEL   specify the output file label.
  --sample_order SAMPLE_ORDER
                        specify the sample order file
```

# PLLM optimization on mutant activity data
 
In this project, we develop two core pipelines to select the best-performing fine-tunning (FT) parameter set for own data and apply it to to our own data.
- 1. **public_data**: Validation of twelve parameter sets based on three widely used public mutation dataset.
- 2. **R4_data**: Benchmarking twelve parameter sets on R4 single-mutation data through repeated holdout validation and using the optimal parameter set to construct a prediction model.

# Requirements

## Python
- Python 3.12
    ```bash
    conda env create -f pllm_optimization_environment.yaml
    conda activate pllm_optimization
    ```

## Base models
- AI4Protein/Prime_690M ([huggingface](https://huggingface.co/AI4Protein/Prime_690M))
- esm2_finetuned_epoch20_R4 (fine-tunned with R4 protein sequences) or facebook/esm2_t33_650M_UR50D ([huggingface](https://huggingface.co/facebook/esm2_t33_650M_UR50D))
- westlake-repl/SaProt_650M_PDB ([huggingface](https://huggingface.co/westlake-repl/SaProt_650M_PDB))

## Other software
- foldseek (version: c1e90b1ebb8559c1e386b7b33667ae51a50745ab)

# Usage

## 1. public_data
- Fine-tunning with a combination from three base protein large language models (PLLMs), two parameter-efficient fine-tunning (FT) strategies, and two loss functions.
- The train/validation/test split for three public datasets are 20,963/52,35/4,004 (avGFP), 3,000/500/850 (KKA2) and 3,489/499/1,417 (TEM-1) respectively.

```txt
usage: 03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py [-h] [--seed SEED] --csv_path CSV_PATH --model_name MODEL_NAME [--batch_size BATCH_SIZE]
                                                                                               [--num_epochs NUM_EPOCHS] [--rank_count RANK_COUNT] [--top_layers TOP_LAYERS] --method
                                                                                               {ft,lora,top_layer} [--learning_rate LEARNING_RATE] [--warming_up WARMING_UP]
                                                                                               [--use_listmle] [--log_path LOG_PATH]

Fine-tuning a regression model using HuggingFace's ESM2_650M model.

options:
  -h, --help            show this help message and exit
  --seed SEED           Random seed for reproducibility, default is 42.
  --csv_path CSV_PATH   Path to the CSV file containing the dataset.
  --model_name MODEL_NAME
                        Path to the pre-trained model.
  --batch_size BATCH_SIZE
                        Batch size for training and evaluation.
  --num_epochs NUM_EPOCHS
                        Number of epochs for training.
  --rank_count RANK_COUNT
                        Number of LoRA ranks to use (only for LoRA method), default is 4
  --top_layers TOP_LAYERS
                        Number of top layers to fine-tune (only for top_layer method), default is 4
  --method {ft,lora,top_layer}
                        Fine-tuning method: "ft", "lora" or "top_layer"
  --learning_rate LEARNING_RATE
                        Learning rate for optimizer, default is 0.0001
  --warming_up WARMING_UP
                        Warming up steps for optimizer, default is 1000
  --use_listmle         Use ListMLE as the loss function
  --log_path LOG_PATH   Path to the log file, default is ./default.log
```

## 2. R4_data

**Step 1:** 
- Fine-tunning with the same twelve parameter sets as in the public data.
- The train/validation split for R4 single-mutation dataset is 80/27. Five independent runs with different sampling seed were conducted.

```txt
usage: 03.Prime_690M_FineTunning_regression_manual_assembled_R4.M1.predictions_TrainValidation.py [-h] [--seed SEED] --csv_path CSV_PATH --model_name MODEL_NAME [--batch_size BATCH_SIZE]
                                                                                                  [--num_epochs NUM_EPOCHS] [--rank_count RANK_COUNT] [--top_layers TOP_LAYERS] --method
                                                                                                  {ft,lora,top_layer} [--learning_rate LEARNING_RATE] [--warming_up WARMING_UP]
                                                                                                  [--use_listmle] [--log_path LOG_PATH]

Fine-tuning a regression model using HuggingFace's Prime_690M model.

options:
  -h, --help            show this help message and exit
  --seed SEED           Random seed for reproducibility, default is 42.
  --csv_path CSV_PATH   Path to the CSV file containing the dataset.
  --model_name MODEL_NAME
                        Path to the pre-trained model.
  --batch_size BATCH_SIZE
                        Batch size for training and evaluation.
  --num_epochs NUM_EPOCHS
                        Number of epochs for training.
  --rank_count RANK_COUNT
                        Number of LoRA ranks to use (only for LoRA method), default is 4
  --top_layers TOP_LAYERS
                        Number of top layers to fine-tune (only for top_layer method), default is 4
  --method {ft,lora,top_layer}
                        Fine-tuning method: "ft", "lora" or "top_layer"
  --learning_rate LEARNING_RATE
                        Learning rate for optimizer, default is 0.0001
  --warming_up WARMING_UP
                        Warming up steps for optimizer, default is 1000
  --use_listmle         Use ListMLE as the loss function
  --log_path LOG_PATH   Path to the log file, default is ./default.log
```

**Step 2:** 
- Using the best-performing parameter to construct a prediction model based on the full traning dataset, that is, 107 data points.
- Saving the fine-tunned model for the use of next round mutation prediction, e.g., doule-mutants.

```text
usage: 06.Prime_690M_FineTunning_regression_manual_assembled_R4.M1.predictions.py [-h] [--seed SEED] --csv_path CSV_PATH --model_name MODEL_NAME [--batch_size BATCH_SIZE]
                                                                                  [--num_epochs NUM_EPOCHS] [--rank_count RANK_COUNT] [--top_layers TOP_LAYERS] --method
                                                                                  {ft,lora,top_layer} [--learning_rate LEARNING_RATE] [--warming_up WARMING_UP] [--use_listmle]
                                                                                  [--log_path LOG_PATH] --output_path OUTPUT_PATH --model_save_path MODEL_SAVE_PATH

Fine-tuning a regression model using HuggingFace's Prime_690M model.

options:
  -h, --help            show this help message and exit
  --seed SEED           Random seed for reproducibility, default is 42.
  --csv_path CSV_PATH   Path to the CSV file containing the dataset.
  --model_name MODEL_NAME
                        Path to the pre-trained model.
  --batch_size BATCH_SIZE
                        Batch size for training and evaluation.
  --num_epochs NUM_EPOCHS
                        Number of epochs for training.
  --rank_count RANK_COUNT
                        Number of LoRA ranks to use (only for LoRA method), default is 4
  --top_layers TOP_LAYERS
                        Number of top layers to fine-tune (only for top_layer method), default is 4
  --method {ft,lora,top_layer}
                        Fine-tuning method: "ft", "lora" or "top_layer"
  --learning_rate LEARNING_RATE
                        Learning rate for optimizer, default is 0.0001
  --warming_up WARMING_UP
                        Warming up steps for optimizer, default is 1000
  --use_listmle         Use ListMLE as the loss function
  --log_path LOG_PATH   Path to the log file, default is ./default.log
  --output_path OUTPUT_PATH
                        Path to save the test data predictions, CSV format.
  --model_save_path MODEL_SAVE_PATH
                        Path to save the PEFT model, directory.
```

- Usage example.
```text
python 06.Prime_690M_FineTunning_regression_manual_assembled_R4.M1.predictions.py --seed 42 --csv_path 01.R4_A43_SingleMutation_chosen_modeling_20250705.csv --model_name model_path_to/Prime_690M --batch_size 5 --num_epochs 400 --rank_count 4 --method lora --learning_rate 0.0001 --log_path 06.Prime_690M_FineTunning_regression_manual_assembled_R4.M1.predictions.Prime_epoch400_lora_MSEloss_lr0.0001.seed42.M1.log --output_path 06.Prime_690M_FineTunning_regression_manual_assembled_R4.M1.predictions.Prime_epoch400_lora_MSEloss_lr0.0001.seed42.M1.csv --model_save_path model_path_to/Prime_690M_peft400
```

# License

Copyright (c) 2026 **Beijing Astragenomics Technology Co., Ltd.**  
All rights reserved.

This source code is released for **non-commercial academic and research purposes only**.
Use, reproduction, and distribution for commercial purposes are **strictly prohibited**
unless expressly authorized by **[Beijing Astragenomics Technology Co., Ltd.](https://www.astragenomics.com)**
