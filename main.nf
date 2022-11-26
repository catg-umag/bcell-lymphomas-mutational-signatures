#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { processParams } from './lib/nextflow/utils'


/*
 * Process params and prepare inputs
 */
params = processParams(params)

if (params.vcf_list != null) {
  channel.fromPath(params.vcf_list, checkIfExists: true)
    .splitCsv(header: true, sep: ',')
    .map { [it.name, it.group, file(it.file)] }
    .set { vcfs }

  vcfs
    .map { it[2] }
    .collect()
    .set { vcf_files }

  vcfs
    .map { [it[0], it[1], it[2].name].join(',') }
    .collect()
    .map { 'name,group,file\n' + it.join('\n') }
    .collectFile(name: 'vcf_list.csv')
    .set { vcf_list }
} else {
  vcf_list = null
}
snv_list = (params.snv_list != null) ? file(params.snv_list) : null
reference = file(params.reference)
ig_list = file(params.ig_list)
cosmic_config = [params.cosmic_version, params.cosmic_genome]
fitting_selected_signatures = params.fitting_selected_signatures
extra_signatures = file(params.fitting_extra_signatures)



workflow {
  if (vcf_list != null) {
    vcfToCsv(vcf_list, vcf_files)
    variants_csv =  vcfToCsv.out
  } else {
    variants_csv = snv_list
  }
  
  annotateVariants(variants_csv, reference, ig_list)
  sigProfiler(annotateVariants.out)

  runNotebook(
    annotateVariants.out,
    sigProfiler.out.results,
    cosmic_config,
    fitting_selected_signatures ?: '',
    extra_signatures
  )
}


/*
 * Create a CSV file from the input VCF files
 */
process vcfToCsv {
  input:
  path(vcf_list)
  path(vcf_files)

  output:
  path('snv_list_raw.csv')

  script:
  """
  vcf_to_csv.py -i ${vcf_list} -o snv_list_raw.csv --autosomes-only 
  """
}


/*
 * Add metadata to variants (context, C/T mutationm, AID motifs, ig_locus)
 */
process annotateVariants {
  publishDir "${params.results_dir}", mode: 'copy'

  input:
  path('snv_list_raw.csv')
  path(reference)
  path(ig_list)

  output:
  path('snv_list.csv')

  script:
  """
  annotate_variants.py -i snv_list_raw.csv -o snv_list.csv -r ${reference} -g ${ig_list}
  """
}


/*
 * Run SigProfiler to extract signatures
 */
process sigProfiler {
  clusterOptions params.sigprofiler_gpu ? '--gres=gpu:1' : ''
  containerOptions = params.sigprofiler_gpu ? '--nv' : ''
  publishDir "${params.results_dir}/extraction", mode: 'copy', pattern: 'out/*.csv', saveAs: { it - ~/^.*\// }
  publishDir "${params.results_dir}/extraction", mode: 'copy', pattern: 'sigprofiler_out'
  cpus "${params.sigprofiler_cpus}"
  
  input:
  path(snv_list)
  
  output:
  path('sigprofiler_out/')
  tuple path('out/signatures.csv'), path('out/statistics.csv'), path('out/contributions.csv'), emit: results
  
  script:
  use_gpu = params.sigprofiler_gpu ? "--use-gpu" : ''
  force_signature = params.nsignatures_force != null
    ? "--force-nsignatures ${params.nsignatures_force}"
    : ''
  """
  make_mutational_matrix.py -i $snv_list -o mutational_matrix.tsv --tab-delimiter

  sig_profiler.py -i mutational_matrix.tsv \
    --output-dir out \
    --min-signatures ${params.nsignatures_min} \
    --max-signatures ${params.nsignatures_max} \
    --n-cpus ${task.cpus} \
    ${use_gpu} \
    ${force_signature}
  """
}


/*
 * Execute Jupyter Notebook with the data obtained from the extraction
 */
process runNotebook {
  publishDir "${params.results_dir}/report", mode: 'copy', pattern: 'SignatureReport.*'
  publishDir "${params.results_dir}", mode: 'copy', pattern: 'output/*', saveAs: { it - ~/^output\// }

  input:
  path('data/snv_list.csv')
  tuple path('data/signatures.csv'), path('data/statistics.csv'), path('data/contributions.csv')
  tuple env(COSMIC_VERSION), env(COSMIC_GENOME)
  env(FITTING_REFERENCE_SIGNATURES)
  path('data/extra_signatures.csv')

  output:
  path('SignatureReport.*')
  path('output/*')
  
  script:
  pre = workflow.containerEngine == 'singularity' ? 'export HOME=""' : ''
  """
  ${pre}
  export DATA_DIR=\${PWD}/data
  export OUTPUT_DIR=\${PWD}/output
  
  jupyter nbconvert --ExecutePreprocessor.kernel_name=ir --to notebook --execute \
    --output SignatureReport.ipynb --output-dir . ${projectDir}/notebook/Main.ipynb

  jupyter nbconvert --to html --no-input --no-prompt --template lab \
    --output SignatureReport.html --output-dir . SignatureReport.ipynb
  sed -i -e 's/\\(--jp-cell-padding:\\) [0-9]*px/\\1 0px/g' SignatureReport.html
  """
}