#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import gzip
import itertools
import sys,os
import random
import argparse
from io import BufferedReader, TextIOWrapper

parser = argparse.ArgumentParser(description='simulate ancient DNA reads')
parser.add_argument('--start', '-s', dest='start', help='left of the random cut position')
parser.add_argument('--end', '-e', dest='end', help='right of the random cut position')
parser.add_argument('--input_fastq', '-in', dest='inputfq', help='input fastq ')
parser.add_argument('--output_fastq', '-out', dest='output_fq', help='output fastq file')
args = parser.parse_args()


#read fastq to records
def read_fastq(fastq,length1,length2):
    if fastq.endswith('gz'):
        fastq_file = TextIOWrapper(BufferedReader(gzip.open(fastq, mode='rb')))
    else:
        fastq_file = open(fastq, mode='rt')
    while True:
        read=itertools.islice(fastq_file, 4)
        head=next(read)
        #length=length if length else head[head.rfind('='):]
        lt1=int(length1)
        lt2=int(length2) 
        Rand=random.randint(0,lt2-lt1)
        right=lt1+Rand
        element = ''.join([head,next(read)[0:right].rstrip()+os.linesep,next(read),next(read)[0:right].rstrip()+os.linesep])
        if element != '':
            yield element
        else:
            break
    fastq_file.close()
    return element

#write reads to gz file
def write_fastq(reads,out):
    out = gzip.open(out+'.fastq.gz','wb')
    for read in reads:
        out.write(read.encode())

if __name__ == "__main__":
 reads=read_fastq(args.inputfq,args.start,args.end)
 write_fastq(reads,args.outputfq)

## simulate.py ##
# python3 simulate.py -s 80 -e 100 -in data.R1.fastq -out simulate-aDNA.fq.gz
