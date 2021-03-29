#!/usr/bin/env nextflow
nextflow.enable.dsl=2


workflow {
  
}


process sigProfiler {
  cpus 8  
  
  input:
  path mutational_matrix
  
  output:
  path "sigprofiler_out/"
  path "out/signatures.csv", emit: signatures
  path "out/statistics.csv", emit: stats
  path "out/contributions.csv", emit: contributions
  
  script:
  force_signature = params.signatures.force_n != null ? "--force-nsignatures ${params.signatures.force_n}" : ""
  """
  sig_profiler.py -i $mutational_matrix \
    --output-dir out \
    --min-signatures ${params.signatures.minimum} \
    --max-signatures ${params.signatures.maximum} \
    --n-cpus ${task.cpus} \
    $force_signature
  """
}