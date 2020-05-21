
import xlrd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


file = '/Users/yunfei/Desktop/growth curve/03122019 Growth Curve LDTs.xlsx'

workbook = pd.ExcelFile(file)
print(workbook.sheet_names)


def get_dataframe(sheet):
    strain = ["WT", "ΔelsL", "Δldt", "ΔelsLΔldt"]
    replicate = [1, 2]
    index = pd.MultiIndex.from_product([replicate, strain], names=['replicate', 'strain'])
    df = workbook.parse(sheet)
    df.rename(columns={df.columns[0] : 'Time'}, inplace=True)
    time = [0, 2, 4, 6, 8, 10, 12, 24]
    df.loc[:, 'Time'] = time
    df.set_index('Time', inplace=True)
    df = df.T.set_index(index)
    mean = df.mean(level='strain').T
    errors = df.std(level='strain').T
    return mean, errors


def plot_GC(sheet, ax):
    mean, errors = get_dataframe(sheet)
    plt.style.use('ggplot')
    line_format = ['-o', '-s', '-^', '-D']
    for strain in mean.columns:
        idx = mean.columns.tolist().index(strain)
        if strain != 'WT':
            label_name = '$' + strain + '$'
        else:
            label_name = strain
        plt.errorbar(mean.index, mean[strain], yerr=errors[strain],
                    capsize=4, linewidth=2, fmt=line_format[idx], label=label_name)
        ax = plt.gca()

        # get handles
        # remove error bars from handles
    ax.set_yscale("log")
    #ax.set_title("24h Growth Curve", {'fontsize': 20})
    ax.set_xlabel("Time (hrs)", {'fontsize': 20})
    ax.set_xticks(mean.index)
    ax.set_xlim(left=0, right=12)
    if sheet == 'OD':
        ax.set_ylabel("OD$_{600}$", {'fontsize': 20})
    elif sheet == 'CFU':
        ax.set_ylabel("CFU/mL", {'fontsize': 20})
    ax.tick_params(labelsize='large')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)

    handle, label = ax.get_legend_handles_labels()
    handle = [h[0] for h in handle]
    ax.legend(handle, label, loc='lower right', frameon=False, handlelength=2, fontsize=20)


fig, ax = plt.subplots(figsize=(6, 8))
plot_GC('CFU', ax)
plt.show()
