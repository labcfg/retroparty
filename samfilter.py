from simplesam import Reader, Writer
import sys, os, re
from os import listdir
from os.path import isfile, join
from tqdm import tqdm, tnrange
from datetime import datetime
from joblib import Parallel, delayed

'''
Count lines if samfile
'''
def count_samlines(filepath):
    with open(filepath) as f:
        i = 0
        for _ in f:
            if _[0] != '@':
                i += 1
        return(i)


def sam2table(filename, inputdir, outputdir):

    autosomeXY = list(range(1, 23))
    autosomeXY.append('X')
    autosomeXY.append('Y')
    autosomeXY = ['chr' + str(x) for x in autosomeXY]

    readsname = os.path.splitext(filename)[0]

    samfile = open(inputdir + filename, 'r')
    tablefile = open(outputdir + readsname + '.txt', 'w')
    errorfile = open(outputdir + readsname + '_error.sam', 'w')
    sam = Reader(samfile)
    error = Writer(errorfile)
    tablefile.write('\t'.join(['ID','FILENAME','READNAME',
                              'CHR','POS','INS_STRAND','RE','R1','R2', 
                              'TLEN','CIGAR_R1','CIGAR_R2','MDFLAG_R1','MDFLAG_R2',
                              'BARCODE','BARCODE_Q']) + '\n')
    id_count = 0
    bar = tnrange(int(count_samlines(inputdir+filename)/2), desc=readsname)
    for i in bar:
        r1 = next(sam)
        r2 = next(sam)
        if r1.rname == r2.rname and r1.rname in autosomeXY:
            if ((r1.flag in [99, 83])
                and (r2.flag in [147, 163])
                and (r1.qname.split('__abq:')[0] == r2.qname.split('__abq:')[0])):
                id_count += 1
                bc = r2.qname.split('__abq:')[2]
                bc_q = r2.qname.split('__abq:')[3]
                re_seq = r1.qname.split('__abq:')[1]
                if r1.reverse:
                    # it's strand of insertion (not read)
                    strand = '+'
                    pos = r2.pos + abs(r1.tlen) - 1
                else:
                    strand = '-'
                    pos = r1.pos
                mdflag_r1 = r1._tags[1]
                mdflag_r2 = r1._tags[1]
                tablefile.write('\t'.join([str(id_count), readsname, r1.qname.split('__abq:')[0], 
                                               r1.rname, str(pos), strand, re_seq,
                                               r1.seq, r2.seq,
                                               str(abs(r1.tlen)),
                                               r1.cigar, r2.cigar, mdflag_r1, mdflag_r2,
                                               bc, bc_q]) + '\n')
            else:
                error.write(r1)
                error.write(r2)
        else:
            error.write(r1)
            error.write(r2)
    samfile.close()
    tablefile.close()
    errorfile.close()
    return(0)


def main(inputdir, outputdir, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.sam')]

    if len(onlyfiles) == 1:
        stat_series = sam2table(filename,
                              inputdir, outputdir)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core)(delayed(sam2table)(filename,
                                            inputdir, outputdir)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()
