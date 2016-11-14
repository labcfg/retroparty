import pandas as pd
import re
from collections import defaultdict
from collections import Counter
import sys, os, re
from os import listdir
from os.path import isfile, join
from tqdm import tqdm_notebook, tnrange
import distance
from operator import itemgetter
import numpy as np
from datetime import datetime
from joblib import Parallel, delayed

'''
get best re (retroelement part)
1. find max_amount of re_part
2. if max_amount not one return min hamming with target re_part
return(seq, amount, hamming)
'''
def get_best_re(x, target_re):
    x_count = sorted(dict(Counter(x)).items(), key=itemgetter(1), reverse=True)
    x_count_max = [(re,a) for re,a in x_count if a == x_count[0][1]]
    if len(x_count_max) == 1:
        return((x_count_max[0][0],x_count_max[0][1],distance.hamming(x_count_max[0][0], target_re)))
    else:
        x_ham = [(re,a,distance.hamming(re, target_re)) for re,a in x_count_max]
        return(max(x_ham, key=itemgetter(2)))


def collapser(filename, inputdir, outputdir, target_re):

    readsname = os.path.splitext(filename)[0]

    df = pd.read_table(inputdir + filename, '\t')
    humantable = open(outputdir + readsname + '_humanread.txt', 'w')
    humantable.write('\t'.join(['CLUSTER_ID',
                                'READNAME',
                                'CHR',
                                'POS',
                                'STRAND',
                                'RE',
                                'RE_AMOUNT',
                                'RE_HAMMING',
                                'R1',
                                'R2',
                                'TLEN',
                                'CIGAR_R1',
                                'MDFLAG_R1',
                                'NUM_READS',
                                'NUM_BC']) + '\n')
    pctable = open(outputdir + readsname + '_pcread.txt', 'w')
    pctable.write('\t'.join(['CLUSTER_ID',
                             'ID_LIST',
                             'FILENAME',
                             'READNAME',
                             'CHR',
                             'POS',
                             'STRAND',
                             'RE',
                             'RE_AMOUNT',
                             'RE_HAMMING',
                             'RE_LIST',
                             'R1',
                             'R2',
                             'TLEN',
                             'CIGAR_R1',
                             'MDFLAG_R1',
                             'NUM_READS',
                             'NUM_BARCODE',
                             'BARCODE_LIST',
                             'BARCODE_Q_LIST']) + '\n')

    df['MDR1_value'] = df['MDFLAG_R1'].apply(lambda x: sum([int(i) for i in re.findall(r'(\d+)+M', x)]))
    df_group = df.groupby(['CHR', 'INS_STRAND', 'POS'])

    cluster_id = 0
    for (chrom, strand, pos), group in tqdm_notebook(df_group, desc=readsname):
        cluster_id += 1
        best_row = group.loc[group['MDR1_value'] == max(list(group['MDR1_value']))].iloc[0]
        best_re = get_best_re(list(group['RE']), target_re)
        num_bc = len(set(list(group['BARCODE'])))
        humantable.write('\t'.join([str(cluster_id),
                                    best_row['READNAME'],
                                    chrom,
                                    str(pos),
                                    strand,
                                    best_re[0],
                                    str(best_re[1]),
                                    str(best_re[2]),
                                    best_row['R1'],
                                    best_row['R2'],
                                    str(best_row['TLEN']),
                                    best_row['CIGAR_R1'],
                                    best_row['MDFLAG_R1'],
                                    str(np.shape(group)[0]),
                                    str(num_bc)]) + '\n')
        pctable.write('\t'.join([str(cluster_id),
                                 ','.join([str(x) for x in list(group['ID'])]),
                                 best_row['FILENAME'],
                                 best_row['READNAME'],
                                 chrom,
                                 str(pos),
                                 strand,
                                 best_re[0],
                                 str(best_re[1]),
                                 str(best_re[2]),
                                 ','.join(list(group['RE'])),
                                 best_row['R1'],
                                 best_row['R2'],
                                 str(best_row['TLEN']),
                                 best_row['CIGAR_R1'],
                                 best_row['MDFLAG_R1'],
                                 str(np.shape(group)[0]),
                                 str(num_bc),
                                 ''.join(list(group['BARCODE'])),
                                 ''.join(list(group['BARCODE_Q']))]) + '\n')


def main(inputdir, outputdir, target_re, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]

    if len(onlyfiles) == 1:
        stat_series = collapser(filename,
                              inputdir, outputdir, target_re)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(collapser)(filename,
                                            inputdir, outputdir, target_re)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()




