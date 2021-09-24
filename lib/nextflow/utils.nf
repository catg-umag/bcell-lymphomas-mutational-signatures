def processParams(params) {
  // one of two inputs should be mandatory
  if (isEmpty(params, 'vcf_list') && isEmpty(params, 'snv_list')) {
    throwError('You must provide either vcf_list or snv_list as input')
  } else {
    if (!isEmpty(params, 'snv_list')) {
      params.vcf_list = null
    } else {
      params.snv_list = null
    }
  }

  // check obligatory parameters
  ['reference', 'ig_list'].each {
    if (isEmpty(params, it)) { throwError("Param $it is mandatory") }
  }

  // set defaults for optional parameters
  defaults = [
    nsignatures_min: 2,
    nsignatures_max: 6,
    nsignatures_force: null,
    cosmic_version: 3.2,
    cosmic_genome: 'GRCh38',
    fitting_selected_signatures: null,
    fitting_extra_signatures: 'NO_FILE',
    results_dir: 'results',
    sigprofiler_cpus: 8,
    sigprofiler_gpu: false
  ]
  defaults.each { name, default_value ->
    params[name] = getParamValue(params, name, default_value)
  }

  return params
}


def isEmpty(params, name) {
  return !params.containsKey(name) || params[name] == null
}


def getParamValue(params, name, default_value) {
  return !isEmpty(params, name) ? params[name] : default_value
}


def throwError(msg) {
  exit 1, "ERROR: $msg"
}