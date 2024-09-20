import pathlib
import sys
from ast import literal_eval

import pandas as pd


donors_key = {
    'B001': 'B001',
    'B002': 'B002',
    'P10': 'P10',
    'P21': 'P21',
    'Control1': 'Control',
    'Control10': 'Control',
    'Control2': 'Control',
    'Control3': 'Control',
    'Control4': 'Control',
    'Control5': 'Control',
    'Control8': 'Control',
    'Control9': 'Control',
}


def extract_donor_id(s):
    donor_pattern = r"([BP]\d+|Ctr\d+|Control\d+)_\d+"
    return s.str.findall(donor_pattern).explode()


def find_paths(sample_data, directory, pat):
    """Return fastq paths for each sample."""
    pat = literal_eval(pat) # Parse from json

    def get_paths(ids):
        data_dir = pathlib.Path(directory)
        s = ids.replace(regex={'_': '-'})
        s = s.apply(lambda x: list(data_dir.glob(f'*{x}*')))
        s = s.apply(lambda x: [str(path) for path in x])
        return s

    def extract_read_number(paths):
        return paths.str.findall(pat).explode()

    sample_paths_data = sample_data.assign(
        paths=get_paths(sample_data['ID'])
    ).explode('paths')

    read_numbers = extract_read_number(sample_paths_data['paths'])
    read_labels = ['fastq_' + num if num in ['1', '2'] else None for num in read_numbers]

    return sample_paths_data.assign(read=read_labels)


def create_samplesheet(input_path, seq_info):
    """Prepare sample sheet for pipeline."""
    data = pd.read_excel(input_path, engine='openpyxl').filter(['ID', 'species'])
    data = data.assign(donor=extract_donor_id(data['ID']))

    return (
        pd.concat([
            find_paths(data, **seq_info[group]).assign(group=group) 
            for group in seq_info
        ])
        .dropna(subset=['read'])
        #
        # NOTE: The aggfunc + explode here handles ID collisions for
        #       aviti	B001	Escherichia coli	B001_4186
        #       aviti	B002	Bifidobacterium spp	B002_503
        #    illumina	B001	Escherichia coli	B001_4186
        #
        # TODO: Fix by...
        #   1. Modify regex for B002_503
        #   2. Use both fastqs for B001_4186 (one sample is low-quality)
        #      ... by tolerating lists in sample sheet?
        # 
        .pivot_table(
            index=['group', 'donor', 'species', 'ID'],
            columns='read',
            values='paths',
            dropna=False,
            aggfunc=lambda x: x
        )
        .dropna(how='all')
        .explode(['fastq_1', 'fastq_2'])
        .reset_index()
        .rename_axis('sample')
        .reset_index()
    )


__all__ = [
    donors_key,
    create_samplesheet,
]