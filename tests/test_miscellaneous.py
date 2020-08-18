import pathlib
import tempfile
import unittest


import Bio.SeqIO.FastaIO


import bin.utility


from . import tests_directory


class AggregateSnpSites(unittest.TestCase):

    def setUp(self):
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/snp_sites.tsv').open('r') as fh:
            self.expected = [line.rstrip().split('\t') for line in fh]
        # Run full script
        program = 'aggregate_snp_sites.py'
        input_fps = [str(fp) for fp in tests_directory.glob('data/snp_sites/*tsv')]
        command_args = {
            '--sites_fps': ' '.join(input_fps),
            '--replicons_passing_fp': tests_directory / 'data/other/passing_isolate_replicons.tsv',
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        self.results = [line.rstrip().split('\t') for line in result.stdout.rstrip().split('\n')]

    def test_full(self):
        self.assertEqual(self.results, self.expected)


class CollectReferenceData(unittest.TestCase):

    def setUp(self):
        program = 'collect_reference_data.py'
        command_args = {
            '--reference_fp': tests_directory / 'data/other/reference.gbk',
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        self.results = result.stdout.rstrip().split('\n')

    def test_full(self):
        expected = [
            'reference_name\treference',
            'reference_replicon_names\tcontig_1,contig_2',
            'reference_replicon_sizes\tcontig_1:15836,contig_2:15017',
            (
                'reference_replicon_md5hashes\tcontig_1:fcb0a625261a7572dce701aea09fa813,'
                'contig_2:62ab86f89a5c20dbbe7f8ec9d0e590e1'
            ),
        ]
        self.assertEqual(len(self.results), len(expected))
        for rline, eline in zip(self.results, expected):
            self.assertEqual(rline, eline)


class CreateCoverageDepthMatrix(unittest.TestCase):

    def setUp(self):
        # Get expected
        self.expected_coverage = self.parse_data(tests_directory / 'data/expected_outputs/gene_coverage.csv')
        self.expected_depth = self.parse_data(tests_directory / 'data/expected_outputs/gene_depth.csv')
        # Create temporary directory and define output filepaths
        temp_dir = tempfile.TemporaryDirectory()
        coverage_fp = pathlib.Path(temp_dir.name, 'gene_coverage.csv')
        depth_fp = pathlib.Path(temp_dir.name, 'gene_depth.csv')
        # Run full script
        program = 'create_coverage_depth_matrices.py'
        input_fps = [str(fp) for fp in tests_directory.glob('data/multiway_pileups/*tsv')]
        command_args = {
            '--mpileup_fps': ' '.join(input_fps),
            '--reference_fp': tests_directory / 'data/other/reference.gbk',
            '--output_dir': temp_dir.name,
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        bin.utility.execute_command(command)
        # Read data from output directory and explicitly remove temporary directory
        self.results_coverage = self.parse_data(coverage_fp)
        self.results_depth = self.parse_data(depth_fp)
        temp_dir.cleanup()

    def test_full(self):
        self.assertEqual(self.results_coverage, self.expected_coverage)
        self.assertEqual(self.results_depth, self.expected_depth)

    @staticmethod
    def parse_data(filepath):
        with filepath.open('r') as fh:
            line_token_gen = (line.rstrip().split(',') for line in fh)
            header_tokens = next(line_token_gen)
            isolates = header_tokens[2:]
            results = {isolate: list() for isolate in isolates}
            for replicon, gene, *values in line_token_gen:
                for isolate, value in zip(isolates, values):
                    results[isolate].append((gene, value))
        return results


class CreateSnpAlignment(unittest.TestCase):

    def setUp(self):
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/contig_1_core.mfasta').open('r') as fh:
            self.expected = {desc: seq for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
        # Run full script
        program = 'create_snp_alignment.py'
        command_args = {
            '--allele_fp': tests_directory / 'data/other/contig_1_alleles_cons0.95.csv',
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        fasta_lines = result.stdout.rstrip().split('\n')
        self.results = {desc: seq for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fasta_lines)}

    def test_full(self):
        self.assertEqual(self.results, self.expected)


class DetermineCodingConsequenes(unittest.TestCase):

    def setUp(self):
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/contig_1_consequences_cons0.95.tsv').open('r') as fh:
            self.expected = self.parse_consequences(fh)
        # Run full script
        program = 'determine_coding_consequences.py'
        command_args = {
            '--reference_fp': tests_directory / 'data/other/reference.gbk',
            '--allele_fp': tests_directory / 'data/other/contig_1_alleles_cons0.95.csv',
            '--replicon': 'contig_1',
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        line_iter = iter(result.stdout.rstrip().split('\n'))
        self.results = self.parse_consequences(line_iter)

    def test_full(self):
        self.assertEqual(self.results, self.expected)

    @staticmethod
    def parse_consequences(line_iter):
        records = dict()
        header_tokens = next(line_iter).rstrip().split('\t')
        for line in line_iter:
            line_tokens = line.rstrip().split('\t')
            record = {k: v for k, v in zip(header_tokens, line_tokens)}
            key = (record['position'], record['alt'])
            assert key not in records
            records[key] = record
        return record


class FilterVariants(unittest.TestCase):

    def setUp(self):
        self.result_1 = self.execute_script(tests_directory / 'data/other/isolate_1_contig_1.vcf')
        self.result_2 = self.execute_script(tests_directory / 'data/other/isolate_3_contig_1.vcf')

    def test_full(self):
        self.assertEqual(self.result_1, (19, 0))
        self.assertEqual(self.result_2, (0, 3))

    def execute_script(self, input_fp):
        # Create temporary directory and define output filepaths
        temp_dir = tempfile.TemporaryDirectory()
        snp_homs_q30_fp = pathlib.Path(temp_dir.name, 'homs_q30.tsv')
        snp_hets_fp = pathlib.Path(temp_dir.name, 'hets.tsv')
        # Run full script
        program = 'filter_variants.awk'
        command_args = {
            '-v snp_fp': snp_homs_q30_fp,
            '-v het_fp': snp_hets_fp,
        }
        command_args_fmt = ' '.join(f'{name}={val}' for name, val in command_args.items())
        command = '%s %s %s' % (program, command_args_fmt, input_fp)
        bin.utility.execute_command(command)
        # Read data from output directory and explicitly remove temporary directory
        with snp_homs_q30_fp.open('r') as fh:
            homs = self.get_entries(fh)
        with snp_hets_fp.open('r') as fh:
            hets = self.get_entries(fh)
        temp_dir.cleanup()
        return homs, hets

    @staticmethod
    def get_entries(line_iter):
        i = 0
        for line in line_iter:
            if line.startswith('#'):
                continue
            i += 1
        return i


class GenbankToFasta(unittest.TestCase):

    def setUp(self):
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/reference.fasta').open('r') as fh:
            self.expected = {desc: seq for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
        # Run full script
        program = 'genbank_to_fasta.py'
        command_args = {
            '--input_fp': tests_directory / 'data/other/reference.gbk',
            '--output_fp': '/dev/stdout',
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        line_iter = iter(result.stdout.rstrip().split('\n'))
        self.results = {desc: seq for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(line_iter)}

    def test_full(self):
        self.assertEqual(self.results, self.expected)


class GetPassingRepliconIsolates(unittest.TestCase):

    def setUp(self):
        # Set expected
        with (tests_directory / 'data/expected_outputs/passing_isolate_replicons.tsv').open('r') as fh:
            self.expected = self.parse_data(fh)
        # Run full script
        program = 'get_passing_isolate_replicons.awk'
        input_fp_1 = tests_directory / 'data/other/contig_1_mapping_stats.tsv'
        input_fp_2 = tests_directory / 'data/other/contig_2_mapping_stats.tsv'
        command = f'cd {input_fp_1.parent} && {program} {input_fp_1.name} {input_fp_2.name}'
        result = bin.utility.execute_command(command)
        self.results = self.parse_data(result.stdout.rstrip().split('\n'))

    def test_full(self):
        self.assertEqual(len(self.results), len(self.expected))
        self.assertEqual(self.results.keys(), self.expected.keys())
        self.assertEqual(self.results, self.expected)

    @staticmethod
    def parse_data(line_iter):
        results = dict()
        for line in line_iter:
            isolate, contig = line.rstrip().split('\t')
            results[isolate] = contig
        return results


class MergeGeneStatMatrix(unittest.TestCase):

    def setUp(self):
        # Set expected
        with (tests_directory / 'data/expected_outputs/gene_coverage_merged.csv').open('r') as fh:
            self.expected_coverage = self.parse_data(fh)
        with (tests_directory / 'data/expected_outputs/gene_depth_merged.csv').open('r') as fh:
            self.expected_depth = self.parse_data(fh)
        self.results_coverage = self.execute_script('coverage')
        self.results_depth = self.execute_script('depth')

    def execute_script(self, name):
        program = 'merge_gene_stat_matrix.py'
        command_args = {
            '--fp_1': tests_directory / f'data/other/gene_{name}_premerge_1.csv',
            '--fp_2': tests_directory / f'data/other/gene_{name}_premerge_2.csv',
        }
        # Run script
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        line_iter = iter(result.stdout.rstrip().split('\n'))
        return self.parse_data(line_iter)

    def test_full(self):
        self.assertEqual(self.results_coverage, self.expected_coverage)
        self.assertEqual(self.results_depth, self.expected_depth)

    @staticmethod
    def parse_data(line_iter):
        header_tokens = next(line_iter).rstrip().split(',')
        isolates = header_tokens[2:]
        results = {isolate: list() for isolate in isolates}
        for line in line_iter:
            replicon, gene, *values = line.rstrip().split(',')
            for isolate, value in zip(isolates, values):
                results[isolate].append((gene, value))
        # Order gene, stat tuples
        for isolate, gene_data in results.items():
            results[isolate] = sorted(gene_data)
        return results


class SplitMappedUnmappedReads(unittest.TestCase):

    def setUp(self):
        # Create temporary directory and define output filepaths
        temp_dir = tempfile.TemporaryDirectory()
        input_fp = tests_directory / 'data/other/isolate_10.sam'
        mapped_fp = pathlib.Path(temp_dir.name, 'mapped.bam')
        unmapped_fp = pathlib.Path(temp_dir.name, 'unmapped.bam')
        # Run full script
        program = 'split_mapped_unmapped_reads.awk'
        command_args = {
            '-v mapped_fp': mapped_fp,
            '-v unmapped_fp': unmapped_fp,
        }
        command_args_fmt = ' '.join(f'{name}={val}' for name, val in command_args.items())
        command = '%s %s %s' % (program, command_args_fmt, input_fp)
        bin.utility.execute_command(command)
        # Read data from output directory and explicitly remove temporary directory
        self.mapped_count = self.get_read_count(mapped_fp)
        self.unmapped_count = self.get_read_count(unmapped_fp)
        temp_dir.cleanup()

    def test_full(self):
        self.assertEqual(self.mapped_count, 1312)
        self.assertEqual(self.unmapped_count, 1074)

    @staticmethod
    def get_read_count(filepath):
        result = bin.utility.execute_command(f'samtools view {filepath} | wc -l')
        return int(result.stdout)
