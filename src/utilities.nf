def print_splash() {
  log.info('----------------------------------------------------------------------------------')
  log.info("""
  ‌‌

                                d8b           d8b
                                88P           88P
                               d88           d88
    88bd88b     d8888b     d888888       d888888       d8888b      d888b8b
    88P'  `    d8b_,dP    d8P' ?88      d8P' ?88      d8P' ?88    d8P' ?88
   d88         88b        88b  ,88b     88b  ,88b     88b  d88    88b  ,88b
  d88'         `?888P'    `?88P'`88b    `?88P'`88b    `?8888P'    `?88P'`88b
                                                                         )88
                                                                        ,88P
                                                                    `?8888P

  ‌‌
  """.stripIndent())
  log.info('----------------------------------------------------------------------------------')
}


// Execute a command and capture standard streams along with return code
def execute_command(command) {
  stdout = new StringBuilder()
  stderr = new StringBuilder()
  process = command.execute()
  process.waitForProcessOutput(stdout, stderr)
  return [process.exitValue(), stdout, stderr]
}


// For an optional stage param variable, check that it is either a Boolean or String
// If it is a string and either 'true' or 'false', return the boolean equivalent
def check_boolean_option(option, name) {
  if (option.getClass() == java.lang.Boolean) {
    return option
  } else if (option.getClass() == java.lang.String) {
    if (option.toLowerCase() == 'true') {
      return true
    } else if (option.toLowerCase() == 'false') {
      return false
    }
  }
  exit 1, "ERROR: ${name} option must be true or false"
}


def check_arguments(params) {
  // Check required input and outputs
  if (! params.reads) {
    exit 1, "ERROR: option 'reads' must be set in nextflow.config"
  }
  if (! params.reference) {
    exit 1, "ERROR: option 'reference' must be set in nextflow.configrequired"
  }
  if (! params.output_dir) {
    exit 1, "ERROR: option 'output_dir' must be set in nextflow.configrequired"
  }
}


def check_output_dir(params) {
  // Do not run if output exists and contains files other than the run info directory (which is created by this point)
  output_dir_files = []
  output_dir = file(params.output_dir)
  output_dir.eachFile { output_dir_files.add(it.name) }
  run_info_dirname = file(params.run_info_dir).simpleName
  output_dir_files.remove(run_info_dirname)
  if (output_dir_files.size() > 0 && ! params.force) {
    exit 1, "ERROR: output directory '${output_dir}' already exists and contains other files, remove or use --force to overwrite"
  }
}


def check_disallowed_arguments(workflow) {
  // Certain options cause issues when specified on the commandline. Prevent user to setting them this way.
  disallowed_args = ['--reads', '--output_dir', '--run_info_dir', '--existing_run_dir']
  commandline_args = workflow.commandLine.tokenize(' ')
  bad_args = []
  disallowed_args.each { if (commandline_args.contains(it)) { bad_args.add(it) } }
  if (bad_args.size() > 0) {
    exit 1, 'ERROR: disallowed argument(s) ' + bad_args.join(', ') + ' were specified on the command line.'
  }
}


def check_host(workflow) {
  // Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
  massive_hostnames = ['m3-login1', 'm3-login2']
  on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
  profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
  if (on_massive && ! profile_explicit) {
    exit 1, "ERROR: to run on MASSIVE you must explicitly set -profile massive"
  }
}


def write_reference_data_to_run_config() {
  script = 'collect_reference_data.py'
  command = "${workflow.projectDir}/bin/${script} --reference_fp ${params.reference} > ${params.run_info_dir}/run_config.tsv"
  command_fq = ["/bin/bash",  '-c', command]
  (return_code, stdout, stderr) = execute_command(command_fq)
  if (return_code != 0) {
    exit 1, "ERROR: failed to collect reference data, producing this error:\n\n${stderr}"
  }
}


def write_param_data_to_run_config() {
  File run_info_fh = new File("${params.run_info_dir}/run_config.tsv")
  run_info_fh.append("bt2_max_frag_len\t${params.bt2_max_frag_len}\n")
  run_info_fh.append("bt2_mode\t${params.bt2_mode}\n")
  run_info_fh.append("var_depth_min\t${params.var_depth_min}\n")
  run_info_fh.append("mapping_cover_min\t${params.mapping_cover_min}\n")
  run_info_fh.append("mapping_depth_min\t${params.mapping_depth_min}\n")
  run_info_fh.append("mapping_mapped_min\t${params.mapping_mapped_min}\n")
  run_info_fh.append("outgroup_mod\t${params.outgroup_mod}\n")
  run_info_fh.append("allele_matrix_cons\t${params.allele_matrix_cons}\n")
}


def validate_merge_data(merge_ignore_errors) {
  script = 'validate_merge_data.py'
  command = "${workflow.projectDir}/bin/${script} --src_dir ${params.existing_run_dir} --dst_dir ${params.output_dir} --reference_fp ${params.reference} --read_globs ${params.reads}"
  (return_code, stdout, stderr) = execute_command(command)
  if (return_code != 0) {
    if (! merge_ignore_errors) {
      msg = "ERROR: validation of merge data failed and exited with the following message:\n\n${stderr}\n"
      msg = msg + 'To ignore these errors at your own peril, re-run with --merge_ignore_errors'
      exit 1, msg
    } else {
      msg = 'validation of merge data failed and exited with the following message '
      msg = msg + '(proceeding due to --merge_ignore_errors):'
      msg = msg + "\n\n${stderr}"
      log.warn(msg)
    }
  }
}
