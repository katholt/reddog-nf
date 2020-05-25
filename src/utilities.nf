def print_splash() {
  log.info """
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
  """.stripIndent()
}


def print_help() {
  log.info """
  ==========================================================================
  reddog-nf: reddog nextflow implementation
  ==========================================================================
  Homepage (original pipeline): https://github.com/katholt/RedDog
  Homepage (nextflow implementation): https://github.com/scwatts/reddog-nf

  Required:
  --------------------
    --reads               Input reads (paired, gzip compressed)

    --reference           Reference genbank

    --output_dir          Output directory
  --------------------

  Other:
  --------------------
    --help                Displays this help message
  --------------------

  ==========================================================================

  """.stripIndent()
}


// Execute a command and capture standard streams along with return code
def execute_command(String command) {
  stdout = new StringBuilder()
  stderr = new StringBuilder()
  process = command.execute()
  process.waitForProcessOutput(stdout, stderr)
  return [process.exitValue(), stdout, stderr]
}


// For an optional stage param variable, check that it is either a Boolean or String
// If it is a string and either 'true' or 'false', return the boolean equivalent
def check_boolean_option(Object option, String name) {
  if (option.getClass() == java.lang.Boolean) {
    return option
  } else if (option.getClass() == java.lang.String) {
    if (option.toLowerCase() == 'true') {
      return true
    } else if (option.toLowerCase() == 'false') {
      return false
    }
  }
  exit 1, "error: ${name} option must be true or false"
}


def check_arguments(Object params) {
  // Check required input and outputs
  if (params.help) {
    print_help()
    exit 0
  }
  if (! params.reads) {
    print_help()
    println "‌‌"
    exit 1, "error: option --reads is required"
  }
  if (! params.reference) {
    print_help()
    println "‌‌"
    exit 1, "error: option --reference is required"
  }
  if (! params.output_dir) {
    print_help()
    println "‌‌"
    exit 1, "error: option --output_dir is required"
  }
}


def check_input_files(Object workflow, Object params) {
  // Validate reference
  script = 'validate_reference.py'
  validate_command = "${workflow.projectDir}/bin/${script} --reference_fp ${params.reference}"
  (return_code, stdout, stderr) = execute_command(validate_command)
  if (return_code != 0) {
    exit 1, "error: validation of reference failed and exited with the following message:\n\n${stderr}"
  }

  // TODO: add process and use local executor
  /*
  // Validate first 10 readsets
  script = 'validate_reads.py'
  ch_read_set_fps = ch_read_sets.take(10).flatMap { it -> it[1..2] }.collect()
  validate_command = "${workflow.projectDir}/bin/${script} --reads_fps ${ch_read_set_fps.val.join(' ')}"
  (return_code, stdout, stderr) = execute_command(validate_command)
  if (return_code != 0) {
    exit 1, "error: validation of first ten reads failed and exited with the following message:\n\n${stderr}"
  }
  */
}


def check_output_dir(Object params) {
  // Do not run if output exists and contains files other than the run info directory (which is created by this point)
  output_dir_files = []
  output_dir = file("${params.output_dir}")
  output_dir.eachFile { output_dir_files.add(it.name) }
  run_info_dirname = file(params.run_info_dir).simpleName
  output_dir_files.remove(run_info_dirname)
  if (output_dir_files.size() > 0 && ! params.force) {
    exit 1, "error: output directory '${output_dir}' already exists and contains other files, remove or use --force to overwrite"
  }
}


def check_host(Object workflow) {
  // Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
  massive_hostnames = ['m3-login1', 'm3-login2']
  on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
  profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
  if (on_massive && ! profile_explicit) {
    exit 1, "error: to run on MASSIVE you must explicitly set -profile"
  }
}
