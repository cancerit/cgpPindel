#!/usr/bin/env python3

# core libs
import os  # for mkdir, path stuff
import sys
import argparse

# requirements
import vcfpy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def parse_vcf(df_list, vaf_counts, rl, vcfin, type):
    reader = vcfpy.Reader.from_path(vcfin)
    sample = reader.header.samples.names[0]
    for record in reader:
        if record.INFO['PC'] != type:
            continue
        call = record.call_for_sample[sample]
        if call.data.get('MTP') is None and call.data.get('MTN'):
            continue

        # mismatch data only on Pos
        # diff data only on MT
        df_list.append({'RL': rl, 'SAMPLE': 'MT', 'STRAND': 'POS', 'Reads': call.data.get('MTP'), 'Mismatch': call.data.get('WTM') or 0, 'Diff': call.data.get('MTP') - call.data.get('PP')})
        df_list.append({'RL': rl, 'SAMPLE': 'MT', 'STRAND': 'NEG', 'Reads': call.data.get('MTN'), 'Diff': call.data.get('MTN') - call.data.get('NP')})
        if(call.data.get('VAF') == 0.000):
            if rl in vaf_counts:
                vaf_counts[float(rl)] += 1
            else:
                vaf_counts[float(rl)] = 1
        df_list.append({'RL': rl, 'SAMPLE': 'WT', 'STRAND': 'POS', 'Reads': call.data.get('WTP'), 'Mismatch': call.data.get('MTM') or 0})
        df_list.append({'RL': rl, 'SAMPLE': 'WT', 'STRAND': 'NEG', 'Reads': call.data.get('WTN')})

def process_vcfs(options):
    ### split the labels and create map with vcfs
    r_lengths = options.labels.split(',')
    if len(r_lengths) != len(options.vcfs):
        sys.exit('Error: "-labels" needs to have the same number of elements as the number of vcfs supplied.')
    vcfs = {}
    for l, v in zip(r_lengths, options.vcfs):
        vcfs[float(l)] = v

    df_list = []
    vaf_counts = {}
    for vcf_set in vcfs.items():
        print(f'Processing read length multiplier {vcf_set[0]} for type {options.type}')
        parse_vcf(df_list, vaf_counts, vcf_set[0], vcf_set[1], options.type)
    df = pd.DataFrame.from_records(df_list)

    facet_boxplot(df, 'SAMPLE', 'STRAND', 'RL', 'Reads', options.type+'_reads.png', title='BLAT read depth', ylim=(0, 50))
    facet_boxplot(df, 'SAMPLE', None, 'RL', 'Mismatch', options.type+'_mm.png', title='Mismatch fraction for BLAT reads', aspect=1.2, ylim=(0, 0.05))
    facet_boxplot(df, 'STRAND', None, 'RL', 'Diff', options.type+'_diff.png', title='BLAT reads - Pindel reads', aspect=1.2, ylim=(-20, 10))
    barchart(vaf_counts, 'RL', '0-VAF', options.type+'_0vaf.png', title='Events with VAF=0')

def facet_boxplot(df, row: str, col: str, x_item: str, y_item: str, out_file: str, title=None, aspect=1, ylim=None):
    sns.set()
    grid = sns.FacetGrid(
        df,
        row=row, col=col, margin_titles=True,
        ylim=ylim, aspect=aspect
    )
    grid.map(sns.boxplot, x_item, y_item);
    if title:
        grid.fig.subplots_adjust(top=.9)
        grid.fig.suptitle(title, size=14)
    grid.set_xticklabels(rotation=80)
    grid.savefig(out_file)
    grid.fig.clf() # this clears the figure

def barchart(data_dict, x_label: str, y_label: str, out_file: str, title=None):
    sns.set()
    x_items = []
    y_items = []
    for k in sorted(data_dict):
        x_items.append(k)
        y_items.append(data_dict[k])
    bp = sns.barplot(x=x_items, y=y_items)
    if title:
        bp.set_title(title)
    bp.set_yscale('log')
    bp.set_xlabel(x_label)
    bp.set_ylabel(y_label)
    bp.set_xticklabels(bp.get_xticklabels(), rotation=80)
    plt.savefig(out_file)


parser = argparse.ArgumentParser(description='Generate wisker plots of mismatch rates for a set of vcfs')
parser.add_argument('-d', '--dir', dest='outDir', metavar='outDir', help='Directory for output', required=True)
parser.add_argument('-l', '-labels', dest='labels', metavar='1.0,...', help='CSV of readlength multipliers, same order as vcfs they apply to', required=True)
parser.add_argument('-t', '-type', dest='type', choices=['D','DI','I'], help='Pindel data type (PC=?)', required=True)
parser.add_argument('-f', '--format', dest='imgFormat', choices=['png','pdf','svg'], help='Format to save venn diagram', required=False, default='png')
parser.add_argument('vcfs', nargs='+')

args = parser.parse_args()

# build the output folder before starting work
if os.path.exists(args.outDir) is False:
    os.mkdir(args.outDir, mode=0o700)  # rwx owner only

process_vcfs(args)
