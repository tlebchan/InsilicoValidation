import pandas as pd


def purity2fraction(
    purity_new: float,
    purity: float,
    ploidy: int
) -> float:
    """
    The function calculates fraction of tumor read needed to obtain
    desired downsampled purity
    
    :param purity_new desired downsample purity
    :param purity initial ngs tumor content of a sample
    :param ploidy ngs rounded ploidy of a sample
    :return fraction of tumor reads to obtain desired purity
    """
    numerator  = purity_new*(1-purity+ploidy*purity/2)
    denumerator = purity*(1-purity_new+ploidy*purity_new/2)
    return numerator/denumerator


def maf_processing(path):
    """
    COMMERCIAL SECRET
    
    make your function, which read maf, preprocess it and
    create filter columns in maf DataFrame with 'PASS'
    value for mutations passed it and '' otherwise
    
    :param path of the .maf file
    :return processed DataFrame
    """
    return maf


def filtration_mask(
    maf: pd.core.frame.DataFrame,
    current_filter: bool = False,
    new_filter: bool = False,
    basic: bool = False,
    vaf: float = 0.0
) -> pd.core.series.Series:
    """
    This function constracts a mask series for appropriate filtering
    with target filter by kit coveraged regions and one of the
    following filters: basic, current or new filter
    
    :param maf Dataframe with mutations and filter columns
    :param current_filter shows whether to apply current filter
    :param new_filter shows whether to apply new testing filter
    :param basic shows whether to apply basic filter by depth
    :param vaf threshold
    :return processed DataFrame
    """
    mask = (maf.Tumor_VAF >= vaf) & ~maf.target_filter.isna()
    
    if current_filter:
        mask = mask & (maf.current_filter=='PASS')
    if new_filter:
        mask = mask & (maf.new_filter=='PASS')
    if basic:
        mask = mask & ~maf.basic.isna()
    
    return mask
