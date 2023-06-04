import numpy as np
import pandas as pd
import subprocess
import scipy.stats as stats


def random_ins(
    indel_len_list: list
) -> str:
    """
    The function generates insertion mutation with 'A', 'T', 'G', 'C'
    uniformly distributed in it
    
    :param indel_len_list with indel lengths to introduce
    :return insertion string
    """
    len_ins = np.random.choice(indel_len_list)
    ins = ''.join(np.random.choice(['A', 'T', 'G', 'C'], len_ins))
    return ins


def gen_ins(
    regions_df: pd.core.frame.DataFrame,
    indel_len_list: list,
    n: int
) -> pd.core.frame.DataFrame:
    """
    The function generates bamsurgeon configure table with insertion mutations
    
    :param regions_df table with columns 'chrom', 'start', 'end'
    :param indel_len_list with indel lengths to introduce
    :param n quantity of mutations to introduce
    :return table with columns 'chrom', 'start', 'end',
        'label' - column filled with 'INS' value,
        'ins' - column with insertion string
    """
    ins_df = regions_df.sample(n).copy()
    ins_df['ins'] = [random_ins(indel_len_list) for _ in range(ins_df.shape[0])]
    ins_df['ins_start'] = (ins_df.start+ins_df.end)/2
    ins_df['ins_start'] = ins_df['ins_start'].astype(int)
    ins_df['ins_end'] = ins_df['ins_start'] + 1
    ins_df['label'] = 'INS'
    
    ins_df = ins_df[['chrom', 'ins_start', 'ins_end', 'label', 'ins']]
    ins_df.columns = ['chrom', 'start', 'end', 'label', 'ins']
    
    return ins_df


def gen_del(
    regions_df: pd.core.frame.DataFrame,
    indel_len_list: list,
    n: int
) -> pd.core.frame.DataFrame:
    """
    The function generates bamsurgeon configure table with deletion mutations
    
    :param regions_df table with columns 'chrom', 'start', 'end'
    :param indel_len_list with indel lengths to introduce
    :param n quantity of mutations to introduce
    :return table with columns 'chrom', 'start', 'end',
        'label' - column filled with 'DEL' value,
    """
    
    del_df = regions_df.sample(n).copy()
    del_df['del_size'] = np.random.choice(indel_len_list, del_df.shape[0])
    del_df['del_start'] = del_df.start
    del_df['del_end'] = del_df['del_start'] + del_df['del_size']
    del_df['label'] = 'DEL'
    
    del_df = del_df[['chrom', 'del_start', 'del_end', 'label']]
    del_df.columns = ['chrom', 'start', 'end', 'label']
    
    return del_df


def gen_snp(
    regions_df: pd.core.frame.DataFrame,
    n: int
) -> pd.core.frame.DataFrame:
    """
    The function generates bamsurgeon configure table with snp mutations
    
    :param regions_df table with columns 'chrom', 'start', 'end'
    :param n quantity of mutations to introduce
    :return table with columns 'chrom', 'start', 'end'
    """
    
    snp_df = regions_df.sample(n).copy()
    snp_df['end'] = snp_df.start
    snp_df = snp_df[['chrom', 'start', 'end']]
    
    return snp_df


def depth_extract(
    path2bam: str,
    chrom: str,
    pos: int
) -> int:
    """
    The function extracts depth from .bam for specified position
    
    :param path2bam path to .bam
    :param chrom chromosome of position
    :param pos start position to extract depth
    :return depth of the position in .bam
    """
    
    cmd = (f'samtools depth -r {chrom}:{pos}-{pos} {path2bam}').split()
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output, error = process.communicate()
    output_list = str(output)[2:-3].split('\\t')
    depth = output_list[-1]
    depth = int(depth) if len(depth)>0 else 0
    
    return depth


def calc_vaf(
    variants_df: pd.core.frame.DataFrame,
    expected_vaf: float,
    path: str,
    bam_name: str
) -> list:
    """
    The function extracts depth from .bam for specified position
    
    :param path2bam path to .bam
    :param expected_vaf expected vaf of mutation
    :param chrom chromosome of position
    :param pos start position to extract depth
    :return depth of the position in .bam
    """
    
    vaf_list = []
    for i in range(variants_df.shape[0]):
        chrom, start = variants_df.iloc[i, :2]
        depth = depth_extract(f'{path}/{bam_name}', chrom, start)
        t_alt_count = stats.binom.rvs(n=depth, p=expected_vaf)
        if depth>0:
            vaf = t_alt_count/depth
        else:
            vaf = 0
        vaf_list.append(vaf)

    return vaf_list
