import datetime
from typing import Tuple, List

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
from pathlib import Path
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

from base_classes import PlotConfig, GrowthNumbers, Dataset
import utils
import glob
import os

# --- Core functions ---

def calculate_growth_rates(
        data: pd.DataFrame,
        alpha: float,
        generation_time: float,
        reproduction_number: float
):
    """
    Computes the growth rates from time series data assuming a logistic curve.

    :param data: columns: t:int, b117:int, original:int
    :param alpha:
    :param generation_time:
    :returns:
    """
    RGadvantages = []
    WindowSize = 9

    for i in range(len(data)-WindowSize):  
        extracted_df = data.iloc[i:i+WindowSize].copy().reset_index(drop=True)
        #print(extracted_df.HaploCount)
        window_variant_sum=extracted_df.HaploCount.sum()
        if window_variant_sum > WindowSize*2:
            fd_date = utils.statsmodel_fit(alpha, extracted_df.DateNum, extracted_df.HaploCount, extracted_df.OtherCount + extracted_df.HaploCount, generation_time, reproduction_number)
        else:
            fd_date = 0
        RGadvantages.append(fd_date)

    for i in range(0, WindowSize):
        RGadvantages.append(0)

    return RGadvantages


def main():
    #path_input = 'G://VariationMutation/AiDMS/Spike1205/SequenceData/VOCHaplo/SingleHaplo/*.SingleHaploDateCount'
    path_input1 = 'G://VariationMutation/AiDMS/Spike0417/Diversity/HaploCasesTop1/*.tsv.PreVOCDis'
    input_files = glob.glob(path_input1)  

    for file in input_files:
        if os.path.exists(file+'.WithGrowthRate'): 
            print(file + " existed")
        else:
            print(file)
            variantdata = pd.read_table(file)
            DateRGAdvantages = calculate_growth_rates(variantdata, 0.95, 4.8, 1)
            variantdata.loc[:, 'DailyRelativeGrowthAdvantages'] = DateRGAdvantages
            variantdata.to_csv(file+'.WithGrowthRate', sep='\t', index=False)

if __name__ == "__main__":
    main()
