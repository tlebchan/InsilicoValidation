import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from utils import filtration_mask, maf_processing


def make_reference_set(
    maf_ini_dict: dict, 
    current_filter: bool = False,
    new_filter: bool = False,
    basic: bool = False,
    VAF_threshold: float = 0.0,
    target_variants: list = ['SNP', 'DNP', 'TNP']
) -> dict:
    
    TP_mut_set_dict = {}
    for sample in maf_ini_dict.keys():
        # maf_TP contains only mutations, passed desired filter in initial sample
        maf_TP = maf_ini_dict[sample][
            filtration_mask(maf_ini_dict[sample], current_filter, new_filter, basic, VAF_threshold)
            &maf_ini_dict[sample].Variant_Type.isin(target_variants)
        ]
        
        # save VAFs for mutations in initial sample
        # this will be needed in advanced mode
        TP_mut_set = maf_TP.Tumor_VAF
        
        # save to TP_mut_set_dict
        TP_mut_set_dict[sample] = TP_mut_set
        
    return TP_mut_set_dict


def vaf_filter(
    current_filter: bool = False,
    new_filter: bool = False,
    VAF_threshold: float = 0.0,
    advanced_mode: bool = False
) -> float:
    """
    This function is supporting for plot_reproducability_from_maf_advanced
    
    :param current_filter shows whether to apply current filter
    :param new_filter shows whether to apply new testing filter
    :param VAF_threshold to be applied
    :param advanced_mode indicates whether to adjust reference set
    by VAF or not
    :return final VAF threshold to apply
    """
    if advanced_mode:
        return 0.0
    if current_filter:
        return 0.05
    if new_filter:
        return 0.01
    else:
        return VAF_threshold


def calculate_T_F_P_N(
    TP_mut_set_dict: dict,
    maf_dict: dict,
    current_filter: bool = False,
    new_filter: bool = False,
    basic: bool = False,
    VAF_threshold: float = 0.0,
    target_variants: list = ['SNP', 'DNP', 'TNP'],
    advanced_mode: bool = False
) -> pd.core.frame.DataFrame:
    
    TP_snp_samples = []
    vaf_filter_thresh = vaf_filter(current_filter, new_filter, VAF_threshold, advanced_mode)
    for sample in tqdm.tqdm(maf_dict.keys()):
        sample_ini = sample.split('_')[0]
        fraction = int(sample.split('_')[1])/100
        
        mask1 = (
            (maf_dict[sample].Variant_Type.isin(target_variants))
            &filtration_mask(maf_dict[sample], current_filter, new_filter, basic, VAF_threshold)
        )
        
        # correct reference with vaf_filter_thresh
        TP_mut_set = TP_mut_set_dict[sample_ini][
            TP_mut_set_dict[sample_ini]*fraction >= vaf_filter_thresh
        ].index
        
        mask2 = maf_dict[sample].index.isin(TP_mut_set)
        
        TP = (mask1 & mask2).sum()
        TN = (~mask1 & ~mask2).sum()
        FP = (mask1 & ~mask2).sum()
        FN = (~mask1 & mask2).sum()

        All_TP = len(TP_mut_set)

        TP_snp_samples.append(
            {'TP': TP, 'TN': TN, 'FP': FP, 'FN': FN, 'All_TP': All_TP}
        )
        
    return pd.DataFrame(TP_snp_samples, index = maf_dict.keys())


def plot_reproducability_from_maf(
    maf_ini_dict: dict,
    maf_dict: dict,
    color: str, 
    label: str,
    current_filter: bool = False,
    new_filter: bool = False,
    basic: bool = False,
    VAF_threshold: float = 0.0,
    target_variants: list = ['SNP', 'DNP', 'TNP'],
    coeff: float = 0.1,
    advanced_mode: bool = False
) -> None:
    """
    The function plots reproducability against purity
    for downsampled samples
    
    :param maf_ini_dict dict with .maf of initial samples
    :param maf_dict dict with .maf of downsampled samples
    :param color for the plotting
    :param label name of the filter
    :param current_filter shows whether to apply current filter
    :param new_filter shows whether to apply new testing filter
    :param basic shows whether to apply basic filter by depth
    :param VAF_threshold to be applied
    :param target_variants
    :param coeff decrease scope of the filling std area
    :param advanced_mode indicates whether to adjust reference set
    by VAF or not
    :return None
    """
    
    # reference set
    TP_mut_set_dict = make_reference_set(
        maf_ini_dict, current_filter, new_filter,
        basic, VAF_threshold, target_variants
    )

    # calc TP, TN, FP, FN for each downsampled from maf_dict
    TP_samples_df = calculate_T_F_P_N(
        TP_mut_set_dict, current_filter, new_filter,
        basic, VAF_threshold, target_variants, advanced_mode
    )
        
    # calc reproducability
    TF_samples_df['purity'] = TF_samples_df.index
    TF_samples_df['purity'] = TF_samples_df['purity'].replace(sample2purity)

    TF_samples_df['patient'] = TF_samples_df.index
    TF_samples_df['patient'] = TF_samples_df['patient'].apply(lambda x: x.split('_')[0])
    TF_samples_df['Reproducability'] = TF_samples_df.TP/TF_samples_df.All_TP
    
    # calculate means and stds for plotting
    stats = []
    purities = sorted(TP_samples_df.purity.unique())
    for purity in purities:
        tmp_series = TP_samples_df.Reproducability[
            TP_samples_df.purity.isin([purity])
        ]
        stats.append([tmp_series.mean(), tmp_series.std()])

    stats = np.array(stats)
    means = stats.T[0]
    stds = stats.T[1]
        
    # plot
    plt.plot(
        purities,
        means,
        color = color
    )
    plt.scatter(
        purities,
        means,
        color = color,
        label = label
    )
    plt.fill_between(
        purities,
        means-coeff*stds, means+coeff*stds,
        alpha=0.2,
        color = color
    )
    
    
def plot_metric_for_generated_pats(
    column, metric2plot, color, generated_patients, 
    path2generatedmafs, maf_name, indel_array,
    ref_array
):
    """
    The function plots reproducability against purity
    for downsampled samples
    
    :param column which filter to consider current_filter or new_filter
    :param metric2plot whick metric to plot; works with 'recall' and 'precision'
    :param color for the plotting
    :param generated_patients list with generated patients
    :param path2generatedmafs path to generated samples
    :param maf_name .maf file name
    :param indel_array reference array with indels
    :param ref_array reference array with all mutations
    :return None
    """
    
    tmp_list = []
    for sample in generated_patients:
        maf_path = os.path.join(path2generatedmafs, sample, maf_name)

        # MODIFY maf_processing function from imports
        maf = maf_processing(maf_path)

        purity = file.split('_')[-1]
        maf['ref_index'] = maf.Chromosome.values + '_' + np.where(
            maf.Variant_Type.isin(['INS']), maf.Start_Position+1, maf.Start_Position
        ).astype(str)

        if purity=='75':
            mask1 = maf['ref_index'].isin(indel_array)

        mask1 = maf['ref_index'].isin(ref_array)

        TP = (mask1 & (maf[column] =='PASS')).sum()
        TN = (~mask1 & ~(maf[column] =='PASS')).sum()
        FP = (~mask1 & (maf[column] =='PASS')).sum()
        FN = (mask1 & ~(maf[column] =='PASS')).sum()

        tmp_list.append(
            {'TP': TP, 'TN': TN, 'FP': FP, 'FN': FN, 'purity':purity}
        )

    tmp_df = pd.DataFrame(tmp_list)
    tmp_df['recall'] = tmp_df.TP/(tmp_df.TP + tmp_df.FN)
    tmp_df['precision'] = tmp_df.TP/(tmp_df.TP + tmp_df.FP)
    tmp_df = tmp_df.sort_values(by='purity')

    plt.plot(tmp_df.purity, tmp_df[metric2plot], color = color)
    plt.scatter(tmp_df.purity, tmp_df[metric2plot], color = color)
    plt.ylim(0, 1)
