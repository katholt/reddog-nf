import pathlib
import tempfile
import unittest


import bin.calculate_mapping_stats
import bin.utility


from . import tests_directory


class TestRecord:

    def __init__(self, total_reads, replicon_coverage, pass_fail):
        self.total_reads = total_reads
        self.replicon_coverage = replicon_coverage
        self.pass_fail = pass_fail

        self.ratio = float()
        self.phylogeny_group = str()


class AssignIngroupOutgroup(unittest.TestCase):

    def setUp(self):
        base_ingroup_records = [TestRecord(1e6, 100, 'p') for i in range(10)]
        fail_records = [TestRecord(1e6, 10, 'f') for i in range(5)]
        self.record_set_1 = [
            TestRecord(4.49e9, 50, 'p'),
            TestRecord(5e9, 50, 'p'),
            *base_ingroup_records,
            *fail_records,
        ]
        self.record_set_2 = [
            TestRecord(5.001e5, 50, 'p'),
            *base_ingroup_records,
        ]
        self.record_set_3 = [
            TestRecord(5e5, 50, 'p'),
            *base_ingroup_records,
        ]

    def test_assign_1(self):
        records_assigned = bin.utility.assign_ingroup_outgroup(self.record_set_1, 2)
        groups_expected = ['i', 'o'] + ['i'] * 10 + [''] * 5
        self.assertEqual([r.phylogeny_group for r in records_assigned], groups_expected)

    def test_assign_2(self):
        records_assigned = bin.utility.assign_ingroup_outgroup(self.record_set_2, 2)
        groups_expected = ['o'] + ['i'] * 10
        self.assertEqual([r.phylogeny_group for r in records_assigned], groups_expected)

    def test_assign_3(self):
        records_assigned = bin.utility.assign_ingroup_outgroup(self.record_set_3, 2)
        self.assertEqual([r.phylogeny_group for r in records_assigned], ['i'] * 11)


class MappingStats(unittest.TestCase):

    def test_unmapped_read_count_1(self):
        filepath = tests_directory / 'data/other/isolate_1_unmapped.bam'
        unmapped_count = bin.calculate_mapping_stats.get_unmapped_read_count(filepath)
        self.assertEqual(unmapped_count, 0)

    def test_unmapped_read_count_2(self):
        filepath = tests_directory / 'data/other/isolate_17_unmapped.bam'
        unmapped_count = bin.calculate_mapping_stats.get_unmapped_read_count(filepath)
        self.assertEqual(unmapped_count, 2436)

    def test_depth_cov_size(self):
        filepath = tests_directory / 'data/other/isolate_1_coverage_depth.tsv'
        depth_averages, coverages, sizes = bin.calculate_mapping_stats.get_depth_coverage_and_sizes(filepath)
        self.assertEqual(depth_averages['contig_1'], 14.49)
        self.assertEqual(depth_averages['contig_2'], 14.46)
        self.assertEqual(coverages['contig_1'], 99.94)
        self.assertEqual(coverages['contig_2'], 99.93)
        self.assertEqual(sizes['contig_1'], 15836)
        self.assertEqual(sizes['contig_2'], 15017)

    def test_read_counts(self):
        filepath = tests_directory / 'data/other/isolate_1.bam'
        mapped_total, mapped_replicons = bin.calculate_mapping_stats.get_read_counts(filepath)
        self.assertEqual(mapped_total, 1786)
        self.assertEqual(mapped_replicons['contig_1'], 918)
        self.assertEqual(mapped_replicons['contig_2'], 868)
        self.assertEqual(sum(mapped_replicons.values()), 1786)

    def test_snp_indel_counts_1(self):
        filepath = tests_directory / 'data/other/isolate_1_q30.vcf'
        replicon_snps, replicon_indels = bin.calculate_mapping_stats.get_snp_indel_counts(filepath)
        self.assertEqual(replicon_snps['contig_1'], 19)
        self.assertEqual(replicon_snps['contig_2'], 2)
        self.assertEqual(replicon_indels['contig_1'], 0)
        self.assertEqual(replicon_indels['contig_2'], 0)

    def test_snp_indel_counts_2(self):
        filepath = tests_directory / 'data/other/isolate_4_q30.vcf'
        replicon_snps, replicon_indels = bin.calculate_mapping_stats.get_snp_indel_counts(filepath)
        self.assertEqual(replicon_snps['contig_1'], 0)
        self.assertEqual(replicon_indels['contig_1'], 1)

    def test_het_counts_1(self):
        filepath = tests_directory / 'data/other/isolate_1_hets.vcf'
        replicon_hets = bin.calculate_mapping_stats.get_het_counts(filepath)
        self.assertEqual(replicon_hets, dict())

    def test_het_counts_2(self):
        filepath = tests_directory / 'data/other/isolate_3_hets.vcf'
        replicon_hets = bin.calculate_mapping_stats.get_het_counts(filepath)
        self.assertEqual(replicon_hets['contig_1'], 3)


class MappingStatsFull(unittest.TestCase):

    def setUp(self):
        with (tests_directory / 'data/expected_outputs/isolate_1_mapping_stats.tsv').open('r') as fh:
            self.expected_lines = [line.rstrip() for line in fh.readlines()]
        program = 'calculate_mapping_stats.py'
        command_args = {
            '--bam_fp': tests_directory / 'data/other/isolate_1.bam',
            '--vcf_q30_fp': tests_directory / 'data/other/isolate_1_q30.vcf',
            '--vcf_hets_fp': tests_directory / 'data/other/isolate_1_hets.vcf',
            '--coverage_depth_fp': tests_directory / 'data/other/isolate_1_coverage_depth.tsv',
            '--bam_unmapped_fp': tests_directory / 'data/other/isolate_1_unmapped.bam'
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        results = bin.utility.execute_command(command)
        self.result_lines = results.stdout.rstrip().split('\n')

    def test_full(self):
        self.assertEqual(len(self.result_lines), 3)
        self.assertEqual(self.result_lines[0], self.expected_lines[0])
        self.assertEqual(self.result_lines[1], self.expected_lines[1])
        self.assertEqual(self.result_lines[2], self.expected_lines[2])


class AggregateMappingStatsFull(unittest.TestCase):

    def setUp(self):
        # Set expected outputs
        with (tests_directory / 'data/expected_outputs/contig_1_mapping_stats.tsv').open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.expected_1 = parse_mapping_stats(line_token_gen)
        with (tests_directory / 'data/expected_outputs/contig_2_mapping_stats.tsv').open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.expected_2 = parse_mapping_stats(line_token_gen)
        # Run full script
        temp_dir = tempfile.TemporaryDirectory()
        program = 'aggregate_mapping_stats.py'
        input_fps = [str(fp) for fp in tests_directory.glob('data/mapping_stats/*tsv')]
        command_args = {
            '--rep_stats_fps': ' '.join(input_fps),
            '--output_dir': temp_dir.name,
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        bin.utility.execute_command(command)
        # Read data from output directory and explicitly remove temporary directory
        results_1_fp = pathlib.Path(temp_dir.name, 'contig_1_mapping_stats.tsv')
        results_2_fp = pathlib.Path(temp_dir.name, 'contig_2_mapping_stats.tsv')
        with results_1_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.results_1 = parse_mapping_stats(line_token_gen)
        with results_2_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.results_2 = parse_mapping_stats(line_token_gen)
        temp_dir.cleanup()

    def test_full(self):
        self.assertEqual(self.results_1, self.expected_1)
        self.assertEqual(self.results_2, self.expected_2)


class MergeMappingStatsFull(unittest.TestCase):

    def setUp(self):
        # Set expected
        with (tests_directory / 'data/expected_outputs/contig_1_mapping_stats_merged.tsv').open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            self.expected = parse_mapping_stats(line_token_gen)
        program = 'merge_mapping_stats.py'
        command_args = {
            '--fp_1': tests_directory / 'data/other/contig_1_mapping_stats_premerge_1.tsv',
            '--fp_2': tests_directory / 'data/other/contig_1_mapping_stats_premerge_2.tsv',
        }
        # Run script
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        result = bin.utility.execute_command(command)
        line_token_gen = (line.rstrip().split('\t') for line in result.stdout.rstrip().split('\n'))
        self.results = parse_mapping_stats(line_token_gen)

    def test_full(self):
        self.assertEqual(self.results, self.expected)


def parse_mapping_stats(line_token_gen):
    # Skip header and return data as dict
    next(line_token_gen)
    return {isolate: data for isolate, *data in line_token_gen}
