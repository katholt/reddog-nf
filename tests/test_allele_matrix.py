import pathlib
import unittest
import tempfile


import bin.create_allele_matrix
import bin.utility


from . import tests_directory


class AlleleMatrix(unittest.TestCase):

    def test_position_support(self):
        record_info_sets = {
            'DP4=0,0,7,7': (0, 1),
            'DP4=5,5,0,0': (1, 0),
            'DP4=5,5,7,0': (0.5882, 0.4117),
            'DP4=5,0,0,0': (1, 0),
            'DP4=5,5,7,7': (0.4167, 0.5834),
        }
        for record_info, expected in record_info_sets.items():
            result = bin.create_allele_matrix.get_position_support(record_info)
            for a, b in zip(result, expected):
                self.assertAlmostEqual(a, b, delta=1e-4)


class AlleleMatrixFull(unittest.TestCase):

    def setUp(self):
        # Create temporary directory, copy BAM and index
        temp_dir = tempfile.TemporaryDirectory()
        bin.utility.execute_command(f'cp {tests_directory / "data/other/isolate_1.bam"} {temp_dir.name}/')
        bin.utility.execute_command(f'samtools index {temp_dir.name}/isolate_1.bam')
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/contig_1_isolate_1_alleles.tsv').open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.expected_1 = parse_allele_matrix(line_token_gen)
        with (tests_directory / 'data/expected_outputs/contig_2_isolate_1_alleles.tsv').open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.expected_2 = parse_allele_matrix(line_token_gen)
        # Run full script
        program = 'create_allele_matrix.py'
        command_args = {
            '--bam_fp': tests_directory / f'{temp_dir.name}/isolate_1.bam',
            '--sites_fp': tests_directory / 'data/other/snp_sites.tsv',
            '--reference_fp': tests_directory / 'data/other/reference.fasta',
            '--output_dir': temp_dir.name
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        bin.utility.execute_command(command)
        # Read data from output directory and explicitly remove temporary directory
        results_1_fp = pathlib.Path(temp_dir.name, 'contig_1_isolate_1_alleles.tsv')
        results_2_fp = pathlib.Path(temp_dir.name, 'contig_2_isolate_1_alleles.tsv')
        with results_1_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.results_1 = parse_allele_matrix(line_token_gen)
        with results_2_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.results_2 = parse_allele_matrix(line_token_gen)
        temp_dir.cleanup()

    def test_full(self):
        self.assertEqual(self.results_1, self.expected_1)
        self.assertEqual(self.results_2, self.expected_2)


class AlleleMatrixAggregateFull(unittest.TestCase):

    def setUp(self):
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/contig_1_alleles.tsv').open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.expected = parse_allele_matrix(line_token_gen)
        # Run full script
        program = 'aggregate_allele_matrices.py'
        input_fps = [str(fp) for fp in tests_directory.glob('data/allele_matrices/*tsv')]
        command_args = {
            '--allele_fps': ' '.join(input_fps),
            '--sites_fp': tests_directory / 'data/other/replicon_sites.tsv',
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        line_token_gen = (line.rstrip().split('\t') for line in result.stdout.rstrip().split('\n'))
        self.results = parse_allele_matrix(line_token_gen)

    def test_full(self):
        self.assertEqual(self.results, self.expected)


class FilterAlleleMatrixFull(unittest.TestCase):

    def setUp(self):
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/contig_1_alleles_core.tsv').open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.expected = parse_allele_matrix(line_token_gen)
        # Run full script
        program = 'filter_allele_matrix.py'
        command_args = {
            '--allele_fp': tests_directory / 'data/other/contig_1_alleles.tsv',
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        line_token_gen = (line.rstrip().split('\t') for line in result.stdout.rstrip().split('\n'))
        self.results = parse_allele_matrix(line_token_gen)

    def test_full(self):
        self.assertEqual(self.results, self.expected)


def parse_allele_matrix(line_token_gen):
    header_tokens = next(line_token_gen)
    isolates = header_tokens[1:]
    isolate_alleles = {isolate: list() for isolate in isolates}
    for position, *alleles in line_token_gen:
        for isolate, allele in zip(isolates, alleles):
            isolate_alleles[isolate].append((position, allele))
    return isolate_alleles
