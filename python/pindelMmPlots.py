#!/usr/bin/env python3

# core libs
import os  # for mkdir, path stuff
import sys
import argparse

# requirements
import vcfpy
import pandas as pd
import seaborn as sns

def parse_vcf(df_list, rl, vcfin):
    reader = vcfpy.Reader.from_path(vcfin)
    sample = reader.header.samples.names[0]
    for record in reader:
        #count +=1
        #if count > 1000:
        #    break
        if record.INFO['PC'] != args.type:
            continue
        call = record.call_for_sample[sample]
        if call.data.get('MTP') is None and call.data.get('MTN'):
            continue
        df_list.append({'CHROM': record.CHROM, 'RL': rl, 'DT': 'WT', 'Mismatch': call.data.get('WTM') or 0, 'POS_READS': call.data.get('WTP'), 'NEG_READS': call.data.get('WTN')})
        df_list.append({'CHROM': record.CHROM, 'RL': rl, 'DT': 'MT', 'Mismatch': call.data.get('MTM') or 0, 'POS_READS': call.data.get('MTP'), 'NEG_READS': call.data.get('MTN')})

        if call.data.get('WTP') > 300:
            print(f'{record.CHROM}:{record.POS}..{record.POS}\t{record.CHROM}:{record.POS}-{record.POS}')

def draw_boxplot(df_list: list, x_item: str, y_item: str, hue: str, out_file: str, ylim=None):
    plot = sns.boxplot(
        x=x_item, y=y_item,
        data=pd.DataFrame.from_records(df_list),
        palette="colorblind",
        hue=hue,
    )
    if ylim:
        plot.set(ylim=ylim)
    fig = plot.get_figure()
    fig.savefig(out_file)
    fig.clf() # this clears the figure


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

### split the labels and create map with vcfs
r_lengths = args.labels.split(',')
if len(r_lengths) != len(args.vcfs):
    sys.exit('Error: "-labels" needs to have the same number of elements as the number of vcfs supplied.')
vcfs = {}
for l, v in zip(r_lengths, args.vcfs):
    vcfs[float(l)] = v

df_list = []
for vcf_set in vcfs.items():
    print(f'Processing read length multiplier {vcf_set[0]} for type {args.type}')
    parse_vcf(df_list, vcf_set[0], vcf_set[1])

draw_boxplot(df_list, 'RL', 'Mismatch', 'DT', args.type+'_mm.png') #, ylim=(0, 0.05))
draw_boxplot(df_list, 'RL', 'POS_READS', 'DT', args.type+'_pos.png') #, ylim=(0, 0.05))
draw_boxplot(df_list, 'RL', 'NEG_READS', 'DT', args.type+'_neg.png') #, ylim=(0, 0.05))


