#!/usr/bin/env nextflow
nextflow.enable.dsl=2


vcf_list = file(params.vcf_list)
reference = file(params.reference)
ig_list = file(params.ig_list)
cosmic_config = [params.cosmic_version, params.cosmic_genome]
fitting_selected_signatures = params.fitting.cosmic_selected_signatures
extra_signatures = file(params.fitting.extra_signatures)


workflow {
  vcfToCsv(vcf_list)
  annotateVariants(vcfToCsv.out, reference, ig_list)
  sigProfiler(annotateVariants.out)

  runNotebook(
    annotateVariants.out,
    sigProfiler.out.results,
    cosmic_config,
    fitting_selected_signatures,
    extra_signatures
  )
}


process vcfToCsv {
  input:
  path vcf_list

  output:
  path 'snv_list_raw.csv'

  script:
  """
  vcf_to_csv.py -i $vcf_list -o snv_list_raw.csv --autosomes-only 
  """
}


process annotateVariants {
  publishDir 'results/', mode: 'copy'

  input:
  path snv_list
  path reference
  path ig_list

  output:
  path 'snv_list.csv'

  script:
  """
  annotate_variants.py -i $snv_list -o snv_list.csv -r $reference -g $ig_list
  """
}


process sigProfiler {
  publishDir 'results/extraction', mode: 'copy', pattern: 'out/*.csv', saveAs: { it - ~/^.*\// }
  publishDir 'results/extraction', mode: 'copy', pattern: 'sigprofiler_out'
  cpus 8  
  
  input:
  path snv_list
  
  output:
  path 'sigprofiler_out/'
  tuple path('out/signatures.csv'), path('out/statistics.csv'), path('out/contributions.csv'), emit: results
  
  script:
  force_signature = params.signature_extraction.force_n != null
    ? "--force-nsignatures ${params.signature_extraction.force_n}"
    : ''
  """
  make_mutational_matrix.py -i $snv_list -o mutational_matrix.tsv \
    --tab-delimiter

  sig_profiler.py -i mutational_matrix.tsv \
    --output-dir out \
    --min-signatures ${params.signature_extraction.minimum} \
    --max-signatures ${params.signature_extraction.maximum} \
    --n-cpus ${task.cpus} \
    $force_signature
  """
}


process runNotebook {
  publishDir 'results/report', mode: 'copy', pattern: 'SignatureReport.*'
  publishDir 'results/', mode: 'copy', pattern: 'output/*', saveAs: { it - ~/^output\// }

  input:
  path 'data/snv_list.csv'
  tuple path('data/signatures.csv'), path('data/statistics.csv'), path('data/contributions.csv')
  tuple env(COSMIC_VERSION), env(COSMIC_GENOME)
  env FITTING_REFERENCE_SIGNATURES
  path 'data/extra_signatures.csv'

  output:
  path 'SignatureReport.*'
  path 'output/*'
  
  script:
  pre = workflow.containerEngine == 'singularity' ? 'export HOME=""' : ''
  """
  $pre
  export DATA_DIR=\${PWD}/data
  export OUTPUT_DIR=\${PWD}/output
  
  jupyter nbconvert --ExecutePreprocessor.kernel_name=ir --to notebook --execute \
    --output SignatureReport.ipynb --output-dir . ${projectDir}/notebook/Main.ipynb

  jupyter nbconvert --to html --no-input --no-prompt --template lab \
    --output SignatureReport.html --output-dir . SignatureReport.ipynb
  sed -i -e 's/\\(--jp-cell-padding:\\) [0-9]*px/\\1 0px/g' SignatureReport.html
  """
}