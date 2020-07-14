import pathlib
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
        filepath = tests_directory / 'data/isolate_1_unmapped.bam'
        unmapped_count = bin.calculate_mapping_stats.get_unmapped_read_count(filepath)
        self.assertEqual(unmapped_count, 0)

    def test_unmapped_read_count_2(self):
        filepath = tests_directory / 'data/isolate_17_unmapped.bam'
        unmapped_count = bin.calculate_mapping_stats.get_unmapped_read_count(filepath)
        self.assertEqual(unmapped_count, 2436)

    def test_depth_cov_size(self):
        filepath = tests_directory / 'data/isolate_1_coverage_depth.tsv'
        depth_averages, coverages, sizes = bin.calculate_mapping_stats.get_depth_coverage_and_sizes(filepath)
        self.assertEqual(depth_averages['contig_1'], 14.49)
        self.assertEqual(depth_averages['contig_2'], 14.46)
        self.assertEqual(coverages['contig_1'], 99.94)
        self.assertEqual(coverages['contig_2'], 99.93)
        self.assertEqual(sizes['contig_1'], 15836)
        self.assertEqual(sizes['contig_2'], 15017)

    def test_read_counts(self):
        filepath = tests_directory / 'data/isolate_1.bam'
        mapped_total, mapped_replicons = bin.calculate_mapping_stats.get_read_counts(filepath)
        self.assertEqual(mapped_total, 1786)
        self.assertEqual(mapped_replicons['contig_1'], 918)
        self.assertEqual(mapped_replicons['contig_2'], 868)
        self.assertEqual(sum(mapped_replicons.values()), 1786)

    def test_snp_indel_counts_1(self):
        filepath = tests_directory / 'data/isolate_1_q30.vcf'
        replicon_snps, replicon_indels = bin.calculate_mapping_stats.get_snp_indel_counts(filepath)
        self.assertEqual(replicon_snps['contig_1'], 19)
        self.assertEqual(replicon_snps['contig_2'], 2)
        self.assertEqual(replicon_indels['contig_1'], 0)
        self.assertEqual(replicon_indels['contig_2'], 0)

    def test_snp_indel_counts_2(self):
        filepath = tests_directory / 'data/isolate_4_q30.vcf'
        replicon_snps, replicon_indels = bin.calculate_mapping_stats.get_snp_indel_counts(filepath)
        self.assertEqual(replicon_snps['contig_1'], 0)
        self.assertEqual(replicon_indels['contig_1'], 1)

    def test_het_counts_1(self):
        filepath = tests_directory / 'data/isolate_1_hets.vcf'
        replicon_hets = bin.calculate_mapping_stats.get_het_counts(filepath)
        self.assertEqual(replicon_hets, dict())

    def test_het_counts_2(self):
        filepath = tests_directory / 'data/isolate_3_hets.vcf'
        replicon_hets = bin.calculate_mapping_stats.get_het_counts(filepath)
        self.assertEqual(replicon_hets['contig_1'], 3)


class MappingStatsFull(unittest.TestCase):

    def setUp(self):
        self.expected_lines= (
            (
                'replicon\tisolate\treplicon_coverage\treplicon_average_depth\treplicon_reads_mapped\ttotal_mapped_reads'
                '\ttotal_reads\tsnps\tsnps_heterozygous\tindels\tpass_fail'
            ),
            'contig_1\tisolate_1\t99.94\t14.49\t51.3998\t100.0\t1786\t0\t0\t1\tp',
            'contig_2\tisolate_1\t99.93\t14.46\t48.6002\t100.0\t1786\t0\t0\t0\tp'
        )
        program = 'calculate_mapping_stats.py'
        command_args = {
            '--bam_fp':  tests_directory / 'data/isolate_1.bam',
            '--vcf_q30_fp':  tests_directory / 'data/isolate_4_q30.vcf',
            '--vcf_hets_fp':  tests_directory / 'data/isolate_1_hets.vcf',
            '--coverage_depth_fp':  tests_directory / 'data/isolate_1_coverage_depth.tsv',
            '--bam_unmapped_fp':  tests_directory / 'data/isolate_1_unmapped.bam'
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        results = bin.utility.execute_command(command)
        self.result_lines = results.stdout.rstrip().split('\n')

    def test_full(self):
        self.assertEqual(len(self.result_lines), 3)
        self.assertEqual(self.result_lines[0], self.expected_lines[0])
        self.assertEqual(self.result_lines[1], self.expected_lines[1])
        self.assertEqual(self.result_lines[2], self.expected_lines[2])
