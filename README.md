# Mutational Signatures in B-cell lymphomas

Software repository for our paper "Integration of mutational signature analysis with 3D chromatin data unveils differential AID-related mutagenesis in B-cell lymphoma" (submitted), for reproducibility purposes.

But if you want, you can use you own data too, everything is automated so it will be easy to run if you want a general landscape of mutational signatures in your samples.

## What is included?

- Creation of mutation list from VCF files (optional)
- Collection of variants extra info (context, AID motifs, SNV in Ig loci)
- SBS signature extraction using [SigProfiler](https://github.com/AlexandrovLab/SigProfilerExtractor)
- Sample fitting against COSMIC Signatures using [deconstructSigs](https://github.com/raerose01/deconstructSigs)
- Signature reconstruction using NNLS approach
- A report including plots to graphically visualize the obtained results

## How to use it?

### Requirements
First, you need to have installed [Nextflow](https://www.nextflow.io/) (>=20.07) and [Singularity](https://sylabs.io/guides/3.0/user-guide/).

### Preparation of inputs
You have two options: starting from the VCFs or starting from a list of variants.

- If you want to start from VCFs:
  
  Prepare a CSV file with 3 columns:
  - **name**: will be used as a sample name for the corresponding file
  - **group**: will be used to separe your samples in the general representations of your samples (for example it could be pathology, sample origin, etc)
  - **file**: VCF path, it is recommended to use absolute paths to avoid issues related with that

  It should should look like this:
  ```
  name,group,file
  CLL_01,CLL/MBL,/home/catg/vcf/CLL_01.snp.filter.som.recode.vcf.hg38_multianno.vcf
  CLL_02,CLL/MBL,/home/catg/vcf/CLL_02.snp.filter.som.recode.vcf.hg38_multianno.vcf
  FL_01_1,FL,/home/catg/vcf/FL_01_1.filter.som.recode.vcf.hg38_multianno.vcf
  FL_01_2,FL,/home/catg/vcf/FL_01_2.filter.som.recode.vcf.hg38_multianno.vcf
  ```

- If you already have a list with your variants:
  
  Basically you need to create a CSV file with this format:
  ```
  sample,group,chrom,pos,ref,alt
  CLL_01,CLL/MBL,4,89250352,T,C
  CLL_01,CLL/MBL,5,49600750,T,C
  CLL_01,CLL/MBL,5,49600906,A,C
  ```

### Run the pipeline

To run run the pipeline, execute:
```
nextflow run CATG-UMAG/bcell-lymphomas-mutational-signatures -r main <params>
```

In `<params>`, you need to provide inputs and other options. These are:

| Parameter                       | Required | Default | Description                                                                                                                                                                                                                                |
| ------------------------------- | -------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--vcf_list`                    | yes*     |         | Input CSV if you want to start with the VCFs (according to previous section). Ignored if `--snv_list` is not empty.                                                                                                                |
| `--snv_list`                    | yes*     |         | Input CSV if you want to start with the list of variants (according to previous section).                                                                                                                                                  |
| `--reference`                   | yes      |         | Reference in 2bit format. Must be the same used in the variant calling. For example: [hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit) or [hg38](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit) |
| `--ig_list`                     | yes      |         | Bed file containing the ranges for the Ig loci. Check `data/iglist_hg38.bed` for a example.                                                                                                                                                |
| `--nsignatures_min`             | no       | 2       | Minimum number of signatures to test with sigprofiler.                                                                                                                                                                                     |
| `--nsignatures_max`             | no       | 5       | Maximum number of signatures to test with sigprofiler.                                                                                                                                                                                     |
| `--nsignatures_force`           | no       |         | Ignore the recomendation from SigProfiler regarding the optimal number of signatures, and use a fixed number of signatures as final output. Must be a number between `nsignatures_min` and `nsignatures_max` values (both inclusive). |
| `--cosmic_version`              | no       | 3.2     | Version of COSMIC signatures to use. Check `data/cosmic_signatures_urls.csv` for possible options.                                                                                                                                         |
| `--cosmic_genome`               | no       | GRCh38  | COSMIC signatures genome. Check `data/cosmic_signatures_urls.csv` for possible options.                                                                                                                                                    |
| `--fitting_selected_signatures` | no       |         | Select only a set of reference signatures for the fitting. The value should be a string containing valid signature names from the COSMIC version selected, separated by commas. Example: "SBS1,SBS3,SBS5,SBS6,SBS9,SBS84"                  |
| `--fitting_extra_signatures`    | no       |         | Provide additional (local) signatures for the fitting. Must be a CSV file, check `data/extra_signatures.csv` for the format.                                                                                                               |
| `--results_dir`                 | no       | results | Output directory to store the results.                                                                                                                                                                                                     |
| `--sigprofiler_cpus`            | no       | 8       | Number of CPUs to use with SigProfiler.                                                                                                                                                                                                    |

So, for example, a full execution command should look like this:
```
nextflow run CATG-UMAG/bcell-lymphomas-mutational-signatures -r main \
  --snv_list data/snv_list.csv --reference data/hg38.2bit --ig_list data/iglist_hg38.bed \
  --nsignatures_min 2 --nsignatures_max 10 --fitting_selected_signatures 'SBS1,SBS3,SBS5,SBS6,SBS9,SBS84'
```

Alternatively, you can provide a yaml file containing all the parameters you want to setup (that way you don't have to write everything on the command line). Just download `params.example.yml` and edit it to your needs (you can delete parameters from the file if you don't want to use them). Then execute the pipeline like this:
```
nextflow run CATG-UMAG/bcell-lymphomas-mutational-signatures -r main --params-file params.yml
```

You can also use any option available in Nextflow.

It's also very easy to run on a computing cluster (as long as Singularity is available). I included a profile for SLURM (`-profile slurm`), if your cluster uses a different scheduler, you should look [here](https://www.nextflow.io/docs/latest/executor.html) to find the corresponding configuration.

## Results

Once the pipeline finished running you will find a set of files. These are:
- `snv_list.csv`: a CSV file with all the variants (if you used variant list as input it will be the same file with extra columns)
- `extraction/`
    - `signatures.csv`: the signatures extracted from your samples
    - `contributions.csv`: a list containing the number of mutations contributed by each signature to every one of your samples
    - `statistics.csv`: metrics collected from the extraction of the different number of signatures
    - `sigprofiler_out`: the raw output from SigProfiler
- `fitting.csv`: the results of the sample fitting process using reference signatures
- `reconstruction/`: reconstruction of each one of the extracted denovo signatures using reference signatures 
- `report/`: a summary of all the obtained information with plots, in `.html` for easy visualization and `.ipynb` (Jupyter Notebook) for editing


## Acknowledgements

- Python libraries: cyvcf2, twobitreader, SigProfilerExtractor and all of its dependencies
- R libraries: cluster, cowplot, deconstructSigs, factoextra, IRkernel, NNLS, R.utils, tidyverse
- Others: Jupyter, Nextflow, Singularity

In `containers/` you can find the Singularity recipes used to build the containers for the pipeline (hosted in [Singularity Hub](https://singularity-hub.org/)). These are the ones configured in `nextflow.config`, alongside others from BioContainers.
