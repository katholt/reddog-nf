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
        bin.utility.execute_command(f'cp {tests_directory / "data/isolate_1.bam"} {temp_dir.name}/')
        bin.utility.execute_command(f'samtools index {temp_dir.name}/isolate_1.bam')
        # Set expected outputs
        self.expected_1 = (
            'Position\tReference\tisolate_1\n1158\tG\tA\n1644\tA\tG\n2081\tG\tT\n2086\tT\tA\n2092\tA'
            '\tT\n4807\tG\tA\n4818\tT\tC\n4820\tT\tT\n4850\tC\tC\n6268\tG\tA\n6473\tA\tC\n6484\tT\tA'
            '\n6485\tG\tC\n6506\tC\tA\n6525\tT\tA\n8001\tG\tG\n10954\tT\tT\n10955\tT\tG\n12053\tG\tC'
            '\n13373\tG\tA\n13377\tT\tG\n13392\tA\tG\n14636\tG\tT\n'
        )
        self.expected_2 = 'Position\tReference\tisolate_1\n157\tC\tG\n560\tG\tG\n561\tC\tC\n779\tC\tG\n'
        # Run full script
        program = 'create_allele_matrix.py'
        command_args = {
            '--bam_fp':  tests_directory / f'{temp_dir.name}/isolate_1.bam',
            '--sites_fp':  tests_directory / 'data/snp_sites.tsv',
            '--reference_fp':  tests_directory / 'data/reference.fasta',
            '--output_dir':  temp_dir.name
        }
        command = '%s %s' % (program, ' '.join(f'{name} {val}' for name, val in command_args.items()))
        bin.utility.execute_command(command)
        # Read data from output directory and explicitly remove temporary directory
        results_1_fp = pathlib.Path(temp_dir.name, 'contig_1_isolate_1_alleles.tsv')
        results_2_fp = pathlib.Path(temp_dir.name, 'contig_2_isolate_1_alleles.tsv')
        with results_1_fp.open('r') as fh:
            self.results_1 = fh.read()
        with results_2_fp.open('r') as fh:
            self.results_2 = fh.read()
        temp_dir.cleanup()

    def test_full(self):
        self.assertEqual(self.results_1, self.expected_1)
        self.assertEqual(self.results_2, self.expected_2)
