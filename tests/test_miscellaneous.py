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


# TODO: need to first check this is correctly calculated. have concern
#       regarding 1-off index error for gene bounds
class CreateCoverageDepthMatrix(unittest.TestCase):

    def setUp(self):
        pass

    def test_full(self):
        pass


class CreateSnpAlignment(unittest.TestCase):

    def setUp(self):
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/contig_1_core.mfasta').open('r') as fh:
            self.expected = {desc: seq for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
        # Run full script
        program = 'create_snp_alignment.py'
        command_args = {
            '--allele_fp': tests_directory / 'data/other/contig_1_alleles_core.tsv',
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
        with (tests_directory / 'data/expected_outputs/contig_1_consequences_core.tsv').open('r') as fh:
            self.expected = self.parse_consequences(fh)
        # Run full script
        program = 'determine_coding_consequences.py'
        command_args = {
            '--reference_fp': tests_directory / 'data/other/reference.gbk',
            '--allele_fp': tests_directory / 'data/other/contig_1_alleles_core.tsv',
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
            # Order isolates
            record['isolates'] = sorted(record['isolates'].split(','))
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


# TODO: see above re CreateCoverageDepthMatrix
class MergeGeneStatMatrix(unittest.TestCase):

    def setUp(self):
        pass

    def test_full(self):
        pass


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

