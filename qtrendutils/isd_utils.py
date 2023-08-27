import pandas as pd
import numpy as np
from math import nan

def remove_bad_rows(df):
    """Based on ISD documentation, remove data with the following sources, report types, call signs, and QC.
    ISD docs: https://www1.ncdc.noaa.gov/pub/data/noaa/isd-format-document.pdf
    # remove SOURCE:
    # 2: failed cross checks
    # A: failed cross checks
    # B: failed cross checks
    # O: Summary observation created by NCEI using hourly observations that
    #    may not share the same data source flag
    # 9: missing
    # remove REPORT_TYPE
    # 99999: missing
    # remove CALL_SIGN
    # 99999: missing
    # remove QUALITY_CONTROL
    # V01: no quality control

    From Karen McKinnon 
    https://github.com/karenamckinnon/helpful_utilities/blob/master/helpful_utilities/download_isd.py
    """
    bad_idx = ((df['SOURCE'].astype(str) == '2') |
               (df['SOURCE'].astype(str) == 'A') |
               (df['SOURCE'].astype(str) == 'B') |
               (df['SOURCE'].astype(str) == 'O') |
               (df['SOURCE'].astype(str) == '9') |
               (df['REPORT_TYPE'].astype(str) == '99999') |
               (df['QUALITY_CONTROL'].astype(str) == 'V010'))

    #df[bad_idx] = nan

    return df[~bad_idx]
    #return df

def remove_bad_rows_isla(df):
    bad_idx = ((df['SOURCE'].astype(str) == '2') |
              (df['SOURCE'].astype(str) == 'A') |
              (df['SOURCE'].astype(str) == 'B') |
              (df['SOURCE'].astype(str) == 'O') |
              (df['SOURCE'].astype(str) == '9') |
              (df['REPORT_TYPE'].astype(str) == '99999') |
              (df['QUALITY_CONTROL'].astype(str) == 'V010'))

    df['TMP'][bad_idx] = nan 
    return df


def remove_bad_vals(df, varname):
    """
    Remove values with questionable flags, based on ISD documentation.
    ISD docs: https://www1.ncdc.noaa.gov/pub/data/noaa/isd-format-document.pdf
    """
    flag = np.array([d.split(',')[-1] for d in df[varname]])
    flag = flag.astype(str)
    vals_tmp = np.array([int(d.split(',')[0]) for d in df[varname]])

    if ((varname == 'DEW') | (varname == 'TMP')):
        bad_idx = ((flag == '2') |
                   (flag == '3') |
                   (flag == '6') |
                   (flag == '7') |
                   (flag == 'A') |
                   (flag == 'C') |
                   (vals_tmp == 9999))
    elif varname == 'SLP':
        bad_idx = ((flag == '2') |
                   (flag == '3') |
                   (flag == '6') |
                   (flag == '7') |
                   (vals_tmp == 99999))

    vals = vals_tmp.astype(float)/10

    vals[bad_idx] = np.nan
#    vals = vals[~bad_idx]
    df = df.assign(**{varname: vals})

    return df

