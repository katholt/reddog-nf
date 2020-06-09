class ReadSet:

    def __init__(self, name):
        self.name = name
        self.variants = dict()

        self.mean_depth = int()
        self.unmapped = float()


class Variants:

    def __init__(self):
        self.homs = list()
        self.hets = list()
        self.indels = list()
        self.lquals = list()


def read_spec_file(filepath):
    # Readsets at top level contain information about broad characteristics
    # Readset variants are stored through readset->contig->variant type
    info = dict()
    readsets = dict()
    with filepath.open('r') as fh:
        # Get defaults
        default_metrics = {'mean_depth', 'outer_length', 'read_length'}
        assert fh.readline().rstrip() == '# Defaults'
        assert fh.readline().rstrip() == 'metric\tvalue\ttype'
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('# Readsets'):
                break
            line_tokens = line.rstrip().split('\t')
            assert len(line_tokens) == 3
            metric, value, metric_type = line_tokens
            if metric not in default_metrics:
                raise ValueError(f'got bad metric {metric}')
            if metric_type == 'int':
                info[metric] = int(value)
            else:
                raise ValueError(f'got back metric type {metric_type}')
        if set(info) ^ default_metrics:
            metric_str = ', '.join(set(info) ^ default_metrics)
            raise ValueError(f'missing default metrics: {metric_str}')
        # Get readsets to simulate
        assert fh.readline().rstrip() == 'isolate_name\tdata'
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('# Variants'):
                break
            line_tokens = line.rstrip().split('\t')
            assert len(line_tokens) == 2
            isolate_name, data_str = line_tokens
            readsets[isolate_name] = ReadSet(isolate_name)
            # Set default mean depth and read_type to pe
            readsets[isolate_name].mean_depth = info['mean_depth']
            readsets[isolate_name].read_type = 'pe'
            if data_str == '-':
                continue
            # Parse readset data
            data = parse_spec_data(data_str)
            for key, val in data.items():
                if key == 'mean_depth':
                    readsets[isolate_name].mean_depth = int(val)
                elif key == 'unmapped':
                    readsets[isolate_name].unmapped = float(val)
                elif key == 'read_type':
                    if val not in {'pe', 'se'}:
                        raise ValueError(f'got bacd read_type {val}')
                    readsets[isolate_name].read_type = val
                else:
                    raise ValueError(f'got bad readset characteristic {key}')
        # Add variant data to readsets
        info['replicons'] = set()
        assert fh.readline().rstrip() == 'isolate_name\treplicon\ttype\tdata'
        for line in fh:
            if line.startswith('##'):
                continue
            line_tokens = line.rstrip().split('\t')
            assert len(line_tokens) == 4
            isolate_name, replicon, variant_type, data_str = line_tokens
            assert isolate_name in readsets
            if replicon not in readsets[isolate_name].variants:
                readsets[isolate_name].variants[replicon] = Variants()
            # Record replicon
            info['replicons'].add(replicon)
            # Parse variant data, cast values, and assign
            data = parse_spec_data(data_str)
            data['position'] = int(data['position'])
            if variant_type == 'hom':
                readsets[isolate_name].variants[replicon].homs.append(data)
            elif variant_type == 'indel':
                assert data['op'] in {'ins', 'del'}
                readsets[isolate_name].variants[replicon].indels.append(data)
            elif variant_type == 'het':
                data['ratio'] = float(data['ratio'])
                readsets[isolate_name].variants[replicon].hets.append(data)
            elif variant_type == 'low_quality':
                readsets[isolate_name].variants[replicon].lquals.append(data)
            else:
                raise ValueError(f'got bad variant type {variant_type}')
    return info, readsets


def parse_spec_data(data_str):
    data = dict()
    for entry in data_str.split(';'):
        key, val = entry.split(':')
        data[key] = val
    return data
