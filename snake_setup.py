#!/usr/bin/env python3

import sys
import pandas as pd
from collections import defaultdict

class ConfigurationError(Exception):
    pass


def unset_defaults(default, outer_key=''):
    """ Return True if default or any value in default
        nested dictionary is an empty string.
    """
    out = 0
    if default == '':
        sys.stderr.write(
            f'\033[31mNo configuration provided for {outer_key} and '
             'no default available.\033[m\n')
        out += 1
    elif isinstance(default, dict):
        for key in default:
            full_key = f'{outer_key} : {key}' if outer_key else key
            out += unset_defaults(default[key], outer_key=full_key)
    return out


def set_default(config, default, key, outer_key):

    RC = 0
    if unset_defaults(default[key], outer_key=outer_key):
        RC = 1
    else:
        config[key] = default[key]
        sys.stderr.write(
            f'\033[33mNo configuration provided for {outer_key}'
            f' - setting to default: {config[key]}.\033[m\n')

    return RC


def set_config(config, default, outer_key=''):

    RC = 0
    for key in default:
        full_key = f'{outer_key} : {key}' if outer_key else key
        try:
            if isinstance(default[key], dict):
                set_config(config[key], default[key], outer_key=full_key)
            elif config[key] is None:
                RC += set_default(config, default, key, outer_key=full_key)
            else:
                sys.stderr.write(
                    f'\033[32mSetting {full_key} to: {config[key]}\033[m\n')
        except KeyError:
            RC += set_default(config, default, key, outer_key=full_key)
        except TypeError:
            RC += set_default(config, default, key, outer_key=full_key)

    if RC > 0:
        raise ConfigurationError(
            '\033[31mInvalid configuration setting.\033[m\n')

    return config


class Samples:

    def __init__(self, samplesFile: str, paired: bool, strand: dict, sep='\s+'):
        self.paired = paired
        self.experimentStrand = strand
        self.table = self.readSamples(samplesFile, sep=sep)


    def readSamples(self, samplesFile, sep='\s+'):
        usecols = ['experiment', 'group', 'rep', 'R1']
        if self.paired:
            usecols = usecols.append('R1')
        table = pd.read_table(
            samplesFile, sep=sep, dtype={'rep': str}, usecols=usecols)

        # Validate read file input with wildcard definitions
        if not table['group'].str.match(r'[^-\.\/]+').all():
            raise ValueError(
                'Groups must not contain the following characters: - . /')
        table['sample'] = (table[['group', 'rep']].apply(lambda x: '-'.join(x), axis=1))

        # Ensure no duplicate names
        duplicateSamples = table['sample'].duplicated()
        if duplicateSamples.any():
            duplicates = list(table['sample'][duplicateSamples])
            raise ValueError(f'Duplicate sample name definitions in {samplesFile} '
                             f'- {duplicates}.')

        if self.paired:
            allFiles = pd.concat([table['R1'], table['R2']])
        else:
            allFiles = table['R1']
        duplicateFiles = allFiles.duplicated()
        if duplicateFiles.any():
            raise ValueError(f'Duplicate file paths in {samplesFile} '
                             f'- {list(allFiles[duplicateFiles])}')

        return table


    def groups(self):
        """ Return unmodified group-rep dictionary """
        return self.table.groupby('group', sort=False)['rep'].apply(list).to_dict()


    def samples(self):
        """ Return unmodified sample list """
        return self.table['sample']


    def singles(self):
        R1 = [f'{sample}-R1' for sample in self.samples()]
        if self.paired:
            R2 = [f'{sample}-R2' for sample in self.samples()]
            return R1 + R2
        else:
            return R1

    def getStrandedness(self):
        strands = {}
        for sample in self.samples():
            # Get experiment associated with sample
            experiment = (self.table
                .loc[self.table['sample'] == sample, 'experiment'].to_list()[0])
            # Get strand associated with experiment
            strand = self.experimentStrand[experiment]
            # Add sample names to dictionary
            strands[sample] = strand
        return strands

    def groupCompares(self):
        """ Return list of pairwise group comparison """
        pairs = itertools.combinations(list(self.groups()), 2)
        return [f'{i[0]}-vs-{i[1]}' for i in pairs]

    def sampleCompares(self):
        """ Return list of pairwise sample comparison """
        pairs = itertools.combinations(self.samples(), 2)
        return [f'{i[0]}-vs-{i[1]}' for i in pairs]


    def path(self, sample, read):
        """ Return path of sample-read """
        paths = self.table.loc[self.table['sample'] == sample, list(read)]
        return paths.values.tolist()[0]
