#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This module provides the class: 'MSIpre.raw_maf.Raw_Maf'
class for adding 'In_repeats' info to each annotation term
of a .maf(mutation annotation format) file.
'''

from intervaltree import IntervalTree
import pandas as pd


def create_repeats_tree(chromosome, ref_repeats_df):
    '''
    Return interval_tree from ref_repeats dataframe of selected chromosome
    '''
    target_ref_repeats = ref_repeats_df[ref_repeats_df['chrom'] == chromosome]
    # axis = 1 means apply to row,
    # @question, why chromEnd + 1, maybe it is 1 based , but python is 0 based
    # @note, bed/BAM is 0 base, GFF、SAM、VCF are 1-base
    interval_tuples = target_ref_repeats.apply(lambda row: (row['chromStart'], row['chromEnd']+1), axis=1)
    repeats_tree = IntervalTree.from_tuples(interval_tuples)   # 线段树, @note, https://blog.csdn.net/SunnyYoona/article/details/43936769
    return repeats_tree


def tag_maf_row(maf_row, repeats_tree):
    '''
    Return a series (a row of mutation annotation format (maf)) tagged with 'In_repeats' info
    '''
    query_tuple = (maf_row['Start_Position'], maf_row['End_Position']+1)
    if len(repeats_tree[query_tuple[0]:query_tuple[1]]):   # @note, maybe use some specificities of IntervalTree
        maf_row['In_repeats'] = 1
    else:
        maf_row['In_repeats'] = 0
    return maf_row


def tag_maf_table(maf_df, ref_repeats_df):

    '''
    Return a dataframe of a maf file tagged with 'In_repeats' info at the end of each row

    '''
    tagged_group_frame = []
    grouped_maf_df = maf_df.groupby('Chromosome')  # groupby
    for chromesome, group_df in grouped_maf_df:
        ref_repeats_tree = create_repeats_tree(chromosome=chromesome, ref_repeats_df=ref_repeats_df)
        tagged_group_df = group_df.apply(lambda row: tag_maf_row(row, ref_repeats_tree), axis =1)
        tagged_group_frame.append(tagged_group_df)
    return pd.concat(tagged_group_frame, ignore_index=True, axis=0)


class Raw_Maf(object):
    """.maf (mutation annotation format) file class"""
    def __init__(self, maf_path):
        '''
        initiate a Raw_Maf class with a path to your .maf file

        '''
        self.maf_path = maf_path

    def create_tagged_maf(self, ref_repeats_file, **tagged_maf_file):
        '''
        Add a column to a .maf file with 'In_repeats' info. In_repeats = 1 denotes that
        this muation annotation term belongs to a simple repeats region, vice versa.

        param ref_repeats_file : path of reference simpleRepeats text file.
        param tagged_maf_file : if not given, a pandas dataframe of a 'In_repeats' info tagged .maf file
        will be returned. If given as: tagged_maf_file = 'your specified output path of tagged_maf file',
        then a tagged .maf file with the given name will be created.
        '''
        maf_file = self.maf_path
        # File read and DataFrame creation
        candidate_chrome = [
            'chr1',
            'chr2',
            'chr3',
            'chr4',
            'chr5',
            'chr6',
            'chr7',
            'chr8',
            'chr9',
            'chr10',
            'chr11',
            'chr12',
            'chr13',
            'chr14',
            'chr15',
            'chr16',
            'chr17',
            'chr18',
            'chr19',
            'chr20',
            'chr21',
            'chr22',
            'chrX',
            'chrY'
        ]

        repeats_column_names = [   # 17 columns
            "bin",
            "chrom",
            "chromStart",
            "chromEnd",
            "name_tag",
            "period_size",
            "copyNUM",
            "consensusSize",
            "perMatch",
            "perIndel",
            "score",
            "A",
            "C",
            "G",
            "T",
            "entropy",
            "unit_sequence"
        ]

        try:
            ref_repeats = pd.read_csv(ref_repeats_file, names=repeats_column_names, sep='\t')
            ref_repeats = ref_repeats[(ref_repeats['chrom'].isin(candidate_chrome)) & (ref_repeats['period_size'] <= 5)][['chrom', 'chromStart', 'chromEnd']]
            # read in maf file by chunks
            chunksize = 10000
            chunks = []
            maf_file_reader = pd.read_csv(maf_file, low_memory=False, comment='#', chunksize = chunksize, sep='\t')
            for chunk in maf_file_reader:
                chunks.append(chunk)
            maf_df = pd.concat(chunks, axis=0)
        except IOError:
            print('Check README for correct usage')
        else:
            # Create tagged('In_repeats' info annotated) maf file dataframe
            annotated_maf = tag_maf_table(maf_df, ref_repeats)
            if not tagged_maf_file:
                return annotated_maf
            else:
                file_path = tagged_maf_file['tagged_maf_file']  # NOTE: 这种写法学习下
                annotated_maf.to_csv(file_path.strip(), sep ='\t', index=False)
