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
  if (params.help) {
    print_help()
    exit 0
  }
  if (! params.reads) {
    print_help()
    println "‌‌"
    exit 1, "ERROR: option --reads is required"
  }
  if (! params.reference) {
    print_help()
    println "‌‌"
    exit 1, "ERROR: option --reference is required"
  }
  if (! params.output_dir) {
    print_help()
    println "‌‌"
    exit 1, "ERROR: option --output_dir is required"
  }
}


def check_input_files(workflow, params) {
  // Validate reference
  script = 'validate_reference.py'
  validate_command = "${workflow.projectDir}/bin/${script} --reference_fp ${params.reference}"
  (return_code, stdout, stderr) = execute_command(validate_command)
  if (return_code != 0) {
    exit 1, "ERROR: validation of reference failed and exited with the following message:\n\n${stderr}"
  }

  // Disallow non-standard '{' and '}' in globs
  if (params.reads.contains('{') & params.reads.contains('}')) {
    exit 1, 'ERROR: input glob should not contain \'{\' or \'}\''
  }

  // TODO: add process and use local executor
  /*
  // Validate first 10 readsets
  script = 'validate_reads.py'
  ch_read_set_fps = ch_read_sets.take(10).flatMap { it -> it[1..2] }.collect()
  validate_command = "${workflow.projectDir}/bin/${script} --reads_fps ${ch_read_set_fps.val.join(' ')}"
  (return_code, stdout, stderr) = execute_command(validate_command)
  if (return_code != 0) {
    exit 1, "ERROR: validation of first ten reads failed and exited with the following message:\n\n${stderr}"
  }
  */
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


def check_host(workflow) {
  // Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
  massive_hostnames = ['m3-login1', 'm3-login2']
  on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
  profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
  if (on_massive && ! profile_explicit) {
    exit 1, "ERROR: to run on MASSIVE you must explicitly set -profile"
  }
}


def validate_merge_data(merge_ignore_errors) {
  script = 'validate_merge_data.py'
  validate_command = "${workflow.projectDir}/bin/${script} --src_dir ${params.previous_run_dir} --dst_dir ${params.output_dir} --reference_fp ${params.reference} --read_globs ${params.reads}"
  (return_code, stdout, stderr) = execute_command(validate_command)
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


// NOTE: there is some issue where DSL2 creates duplicates of channels and executes this function
//       more than once of those channels, creating a race condition while checking if the file
//       exists. I can't see how this is avoidable right now, so I just catch the exception and
//       continue with execution.
def symlink_merge_data(ch, target_dir) {
    if (! target_dir.exists()) {
        target_dir.mkdirs()
    }
    ch.map { filepath_src ->
        filepath_dst = target_dir / filepath_src.getName()
        try {
            java.nio.file.Files.createSymbolicLink(filepath_dst, filepath_src)
        } catch (java.nio.file.FileAlreadyExistsException ex) {
            // Ignore fails if file already exists
        }
    }
}
