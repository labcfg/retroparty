from Bio import SeqIO
import sys, os, re
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import namedtuple, defaultdict
from tqdm import tqdm, tnrange
from joblib import Parallel, delayed
import multiprocessing
import gzip
from IPython.display import display
from operator import itemgetter
from datetime import datetime
import pyximport
pyximport.install()
from cimple_func import hamming

info = namedtuple('info', 'is_good read alu_barcode errors')


'''
check: is r1 is real r1?
(find primer in start of seq with shift)
'''
def is_r1 (record, primer, shift, mist):
    for i in range(shift):
        if hamming(primer, record.seq[i : len(primer)+i], mist):
            return (True)
    return (False)


'''
ONLY R1 (read1)
trimming primers

1. find primer with hamming in start of seq with shift
2. find primer and adapters in natural_zone* and remove if in
3. find restrict_site (seq like 'AGCT' or 'CTAG' - restriction sites) and remove if in
4. trimm primer and return seq

errors = [primer, adapters, restrict_site, natural_zone, r2_start**]

* natural_zone - zone without primer or adapters
** r2_start - nucleotids in start of trim seq (for example: 'CT') (only r2)
'''
def trim_primers (record, primer, ad1, ad2, shift, mist, restrict_site):
    len_primer = len(primer)
    for i in range(shift):
        if hamming(primer, str(record.seq[i : len_primer + i]), mist):
            for elem in [primer, ad1, ad2]:
                if record.seq[len_primer+6+i :].find(elem, 0) != -1:
                    return (info(is_good=False, read=None,
                                 alu_barcode=None,
                                 errors=np.array([0, 0, 0, 1, 0])))
            if record.seq[len_primer+6+i :].find(restrict_site, 0) != -1:
                return (info(is_good=False, read=None,
                             alu_barcode=None,
                             errors=np.array([0, 0, 1, 0, 0])))
            alu_bar = '__abq:' + str(record.seq[len_primer+i : len_primer+6+i])
            record.description = ''
            record.name = ''
            return (info(is_good=True, read=record[len_primer+6+i :],
                         alu_barcode=alu_bar,
                         errors=np.array([0, 0, 0, 0, 0])))
    return (info(is_good=False, read=None,
                 alu_barcode=None,
                 errors=np.array([1, 0, 0, 0, 0])))


'''
ONLY R2 (read2)
trimming adapters

1. find adapters with hamming in start of seq with shift
    pos(ad1) - pos(ad2) = barlen (length of barcode - UMI)
2. find primer and adapters in natural_zone* and remove if in
3. find restrict_site (seq like 'AGCT' or 'CTAG' - restriction sites) and remove if in
4. find r2_start** and remove if not (r2_start can be empty: '')
5. trimm adapters (with UMI) and return seq (with UMI and its quality)

errors = [primer, adapters, restrict_site, natural_zone, r2_start]

* natural_zone - zone without primer or adapters
** r2_start - nucleotids in start of trim seq (for example: 'CT') (only r2)
'''
def trim_adapters (record, primer, ad1, ad2, blen, shift, mist, restrict_site, r2_start):
    len_ad1 = len(ad1)
    len_ad2 = len(ad2)
    trim_len = len_ad1 + blen + len_ad2
    for i in range(shift):
        seq1 = record.seq[i : len_ad1+i]
        seq2 = record.seq[len_ad1+blen+i : len_ad1+blen+len_ad2+i]
        if hamming(ad1, str(seq1), mist) and hamming(ad2, str(seq2), mist):
            if record.seq[trim_len+i : trim_len+i+len(r2_start)] != r2_start:
                return(info(is_good=False, read=None,
                            alu_barcode=None,
                            errors=np.array([0, 0, 0, 0, 1])))
            if record.seq[len_ad1+blen+len_ad2+i :].find(restrict_site, 0) != -1:
                return (info(is_good=False, read=None,
                             alu_barcode=None,
                             errors=np.array([0, 0, 1, 0, 0])))
            for elem in [primer, ad1, ad2]:
                if record.seq[trim_len+i :].find(elem, 0) != -1:
                    return (info(is_good=False, read=None,
                                 alu_barcode=None,
                                 errors=np.array([0, 0, 0, 1, 0])))
            record.description = ''
            record.name = ''
            barcode = str(record.seq[len_ad1+i : len_ad1+blen+i])
            barcode_q = [chr(x + 33) 
                         for x in record.letter_annotations['phred_quality']]
            barcode_q = ''.join(barcode_q)
            barcode_q = barcode_q[len_ad1+i : len_ad1+blen+i]
            alu_bar = '__abq:' + str(barcode) + '__abq:' + str(barcode_q)
            return (info(is_good=True, read=record[trim_len+i :],
                         alu_barcode=alu_bar,
                         errors=np.array([0, 0, 0, 0, 0])))
    return(info(is_good=False, read=None,
                alu_barcode=None,
                errors=np.array([0, 1, 0, 0, 0])))


'''
Humanreadble for error statistics
'''
def concate_errors(type_err, amount_err):
    result = ','.join([str(t) + '-' + str(a)
                      for t, a in zip(type_err, amount_err)])
    return (result)


'''
Count lines if file
'''
def count_lines(filepath):
    with gzip.open(filepath) as f:
        return (sum(1 for _ in f))

'''
trim reads
1. SeqIO.parse - generator (read 'read' by 'read')
2. Chaos = T | F - if r1 and r2 mix in files (R1, R2)
3. Trim and if good write in goodfile else badfile
4. Print and return statistics
'''
def trim_reads(filename1, filename2, inputdir, outputdir,
               primer, ad1, ad2, blen, shift, mist,
               restrict_site, r2_start, chaos):
    
    readsname = filename1.split('R1')[0].rsplit('.', 1)[0]
    outputfile1, ext = os.path.splitext(filename1)
    outputfile2, ext = os.path.splitext(filename2)
    goodr1 = open(outputdir + outputfile1 + '_good.fastq', 'w')
    goodr2 = open(outputdir + outputfile2 + '_good.fastq', 'w')
    badr1 = open(outputdir + outputfile1 + '_bad.fastq', 'w')
    badr2 = open(outputdir + outputfile2 + '_bad.fastq', 'w')

    print('start: ' + readsname)
    unzip_r1 = gzip.open(inputdir + filename1, 'rt')
    unzip_r2 = gzip.open(inputdir + filename2, 'rt')
    original_R1_reads = SeqIO.parse(unzip_r1, "fastq")
    original_R2_reads = SeqIO.parse(unzip_r2, "fastq")
    
    count = np.array([0, 0, 0, 0, 0])
    elem = ('primer', 'adapters', 'restrict_site', 'natural_zone', 'r2_start')
    count_stat = {'readname':readsname, 'all':0, 'good':0, 'bad':0,
                   'primer':0, 'adapters':0,
                   'restrict_site':0, 'natural_zone':0, 'r2_start':0}
    count_stat_col = ['readname', 'all', 'good', 'bad',
                      'primer', 'adapters',
                      'restrict_site', 'natural_zone', 'r2_start']
    bar = tnrange(int(count_lines(inputdir+filename1)/4), desc=readsname)
    original_R12 = zip(original_R1_reads, original_R2_reads)
    for i in bar:
        count_stat['all'] += 1
        r1, r2 = next(original_R12)
        if chaos:
            if is_r1(r1, primer, shift, mist1):
                rx1 = r1
                rx2 = r2
            else:
                rx1 = r2
                rx2 = r1
        else:
            rx1 = r1
            rx2 = r2
        fr1 = trim_primers(rx1, primer, ad1, ad2,
                           shift, mist, restrict_site)
        if fr1.is_good:
            fr2 = trim_adapters(rx2, primer, ad1, ad2, blen,
                                shift, mist, restrict_site, r2_start)
            if fr2.is_good:
                count_stat['good'] += 1
                fr1.read.id += fr1.alu_barcode + fr2.alu_barcode
                fr2.read.id += fr1.alu_barcode + fr2.alu_barcode
                goodr1.write(fr1.read.format('fastq'))
                goodr2.write(fr2.read.format('fastq'))                  
            else:
                rx2.description += (' reason:' +
                        concate_errors(elem, (np.char.mod('%d', fr2.errors))))
                badr1.write(rx1.format('fastq'))
                badr2.write(rx2.format('fastq'))
                count = np.sum([count, fr2.errors], axis=0)
        else:
            rx2.description += (' reason:' +
                    concate_errors(elem, (np.char.mod('%d', fr1.errors))))
            badr1.write(rx1.format('fastq'))
            badr2.write(rx2.format('fastq'))
            count = np.sum([count, fr1.errors], axis=0)
    
    unzip_r1.close()
    unzip_r2.close()
    goodr1.close()
    goodr2.close()
    badr1.close()
    badr2.close()
    
    count_stat['primer'] = count[0]
    count_stat['adapters'] = count[1]
    count_stat['restrict_site'] = count[2]
    count_stat['natural_zone'] = count[3]
    count_stat['r2_start'] = count[4]  
    
    count_stat['bad'] = round((count_stat['all'] - count_stat['good']) / count_stat['all'], 2)
    count_stat['good'] = round(1 - count_stat['bad'], 2)
    count_stat_pd = pd.Series(count_stat, index = count_stat_col)
    count_elem = concate_errors(elem, (np.char.mod('%d', count)))

    return(count_stat_pd)


'''
main
1. find r1,r2-pairs in input directory (*.fastq.gz)
2. run parallel trimming
3. display stat
'''
def main(inputdir, outputdir,
         primer, ad1, ad2, blen,
         shift, mist,
         restrict_site, r2_start,
         chaos, n_core):
    
    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    # Read files in folder
    onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
    
    r1_files = {}
    r2_files = {}
    
    for filename in onlyfiles:
        filename = filename.rstrip()
        if re.search('R1', filename):
            key_filename = filename.split('R1')[0]
            r1_files[key_filename] = filename
        elif re.search('R2', filename):
            key_filename = filename.split('R2')[0]
            r2_files[key_filename] = filename
    
    conform_files = []
    nonconform_files = []
    
    for key in r1_files:
        if key in r2_files:
            conform_files.append((r1_files[key], r2_files[key]))
            del r2_files[key]
        else: nonconform_files.append(r1_files[key])

    conform_files = sorted(conform_files, key=itemgetter(0))
    nonconform_files = nonconform_files + list(r2_files.values())

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    stat_col = ['readname', 'reads', 'good', 'bad',
                'primer', 'adapters',
                'restrict_site', 'natural_zone', 'r2_start']

    if len(conform_files) == 1:
        stat_series = trim_reads(conform_files[0][0], conform_files[0][1],
                              inputdir, outputdir,
                              primer, ad1, ad2, blen,
                              shift, mist,
                              restrict_site, r2_start, chaos)
        stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(trim_reads)(filename1, filename2,
                                            inputdir, outputdir,
                                            primer, ad1, ad2, blen,
                                            shift, mist,
                                            restrict_site, r2_start, chaos)
                                                for filename1, filename2 in conform_files)
        stat_df = pd.concat(stat_series, axis=1).transpose()

    stat_df_total = stat_df.sum(axis = 0)
    stat_df_total.readname = 'total'
    total_bad = sum([stat_df_total.primer,
                     stat_df_total.adapters,
                     stat_df_total.restrict_site,
                     stat_df_total.natural_zone,
                     stat_df_total.r2_start])
    stat_df_total.bad = round(total_bad / stat_df_total['all'], 2)
    stat_df_total.good = round(1 - stat_df_total.bad, 2)
    stat_df.loc[np.shape(stat_df)[0]] = stat_df_total

    display(stat_df)

    stat_name = ''.join([str(before.year),
                         str(before.month),
                         str(before.day),
                         str(before.hour),
                         str(before.minute),
                         str(before.second)])
    stat_df.to_csv(outputdir + 'statistics_' + stat_name + '.csv', sep=' ', index=False)

    after = datetime.now()
    delta_time = after - before

    # Write logfile
    logfile = open(outputdir + 'logfile_' + stat_name + '.log', 'w')
    logfile.write('#DESC OF PARS:\n'+
                '#MIST - max hamming (for primer and ads)\n'+
                '#SHIFT - for search primer or ads\n'+
                '#PRIMER - seq of primer\n'+
                '#AD1 - seq of adapter1\n'+
                '#AD2 - seq of adapter2\n'+
                '#BLEN - len of UMI(barcode)\n'
                '#RESTRICT_SITE - \'AGCT\' or \'CTAG\'\n'+
                '#R2_START - seq after adapters\n'+
                '#CHAOS - mixed files\n'+
                '#N_CORE - number of active cores\n')
    logfile.write('time_start = ' + str(before) + '\n')
    logfile.write('time_end = ' + str(after) + '\n')
    logfile.write('duration (in sec) = ' + str(round(delta_time.total_seconds(), 2)) + '\n')
    logfile.write('MIST = ' + str(mist) + '\n' +
                  'SHIFT = ' + str(shift) + '\n' +
                  'PRIMER = ' + primer + '\n' +
                  'AD1 = ' + ad1 + '\n' +
                  'AD2 = ' + ad2 + '\n' +
                  'BLEN = ' + str(blen) + '\n' +
                  'RESTRICT_SITE = ' + restrict_site + '\n' +
                  'R2_START = ' + r2_start + '\n' +
                  'CHAOS = ' + str(chaos) + '\n'
                  'N_CORE = ' + str(n_core) + '\n' +
                  'INPUTDIR = ' + inputdir + '\n' +
                  'OUTPUTDIR = ' + outputdir + '\n')
    logfile.close()
    stat_df.to_csv(outputdir + 'logfile_' + stat_name + '.log', index=False, sep=' ', mode='a')

    if len(nonconform_files) != 0:
        print ('I can\'t read this files' + str(nonconform_files))
