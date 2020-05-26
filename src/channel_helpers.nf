// Read in isolates the replicons that passed and generate a channel to emit [isolate_id, replicon_ids]
// Where replicon_ids is a string with each replicon_id separated by a single space
// Also get count of passing isolates so we can scale resource allocation if using SLURM executor
def collect_passing_isolate_replicons(Object ch) {
  return ch.flatMap { filepath ->
    filepath.readLines().collect { line ->
      tokens = line.tokenize('\t')
      [tokens[0], tokens[1..-1].join(' ')]
    }
  }
}


// Create a channel that emits allele matrices arranged by replicon_id
//   - input of [isolate_id, list(isolate_allele_matrices)]
//     - nextflow returns a list for multiple files or single object for one file
//     - check for different object types and process accordingly
//   - use isolate_id to robustly get replicon_id from allele matrix filename
//   - flat emit [replicon_id, isolate_allele_matrix] for each file
//   - group each matrix by replicon_id to emit [replicon_id, list(isolate_allele_matrices)]
def sort_allele_matrices(Object ch) {
  return ch.flatMap { isolate_id, filepaths ->
    if (! (filepaths instanceof List)) {
      replicon_id = filepaths.getName().minus("_${isolate_id}_alleles.tsv")
      return [[replicon_id, filepaths]]
    } else {
      return filepaths.collect { filepath ->
          replicon_id = filepath.getName().minus("_${isolate_id}_alleles.tsv")
          [replicon_id, filepath]
        }
    }
  }.groupTuple()
}


// Filter matrices that have no alleles so we don't needlessly execute downstream processes
def filter_empty_allele_matrices(Object ch) {
  return ch.filter { replicon_id, fp ->
      // Read first two lines of allele matrix and determine if we have data
      has_alleles = true
      fh = fp.newReader()
      for (int i = 0; i < 2; i++) { has_alleles = fh.readLine() != null }
      return has_alleles
    }
}


def get_replicon_id(Object ch, Object pattern, Object reference_name) {
  return ch.map { filepath ->
    replicon_id = filepath.getName().minus(pattern).minus(reference_name + '_')
    [replicon_id, filepath]
  }
}
