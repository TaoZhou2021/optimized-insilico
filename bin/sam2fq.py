#!/usr/bin/env python3
# -*- coding:utf-8 -*-


import pysam
import math
import sys
import argparse
from decimal import Decimal

parser = argparse.ArgumentParser(description='select, filtered mate-paris and their align results')
parser.add_argument('--size', '-s', dest='size', help='insert size of the mate-pair libraries')
parser.add_argument('--deviation', '-dev', dest='deviation', help='deviation of the insert size')
parser.add_argument('--readlength', '-length', dest='readlength', help='length of the mate-pair read')
parser.add_argument('--input_sam', '-ins', dest='alignment', help='input sam file of mate-pairs alignmemt file')
parser.add_argument('--select_sam', '-ss', dest='select_sam', help='output sam file of selected mate-pairs')
parser.add_argument('--filter_sam', '-fs', dest='filter_sam', help='output sam file of filtered mate-pairs')
parser.add_argument('--output_fq1', '-out1', dest='output_fastq1', help='output fastq1 file')
parser.add_argument('--output_fq2', '-out2', dest='output_fastq2', help='output fastq2 file')
parser.add_argument('--filtered_fq1', '-f1', dest='filtered_fastq1', help='filtered fastq1 file')
parser.add_argument('--filtered_fq2', '-f2', dest='filtered_fastq2', help='filtered fastq2 file')
parser.add_argument('--align-select', '-as', dest='align_select', help='align results of each selected mate-pairs on the reference')
parser.add_argument('--align-filter', '-af', dest='align_filter', help='align results of each filtered mate-pairs on the reference')
args = parser.parse_args()

#samfile = pysam.AlignmentFile(args.alignment)
# filter_sam = pysam.AlignmentFile(args.filter_sam, "w", template=samfile) 
insert_size=int(args.size)
deviation=Decimal(args.deviation)

def select_mp(alignment,select_sam,output_fastq1,output_fastq2):
 samfile = pysam.AlignmentFile(alignment)
 selected_sam = pysam.AlignmentFile(select_sam, "w", template=samfile)
 with open(output_fastq1,'w') as f1:
  with open(output_fastq2,'w') as f2:
   for read in samfile.fetch():
    a = read.isize
    if read.is_paired:
     if read.is_read1: 
      if (1-deviation)*insert_size <= a <= (1+deviation)*insert_size:
       name= read.query_name
       sequence=read.query_sequence
       quality=read.qual[:]
       selected_sam.write(read)
       f1.writelines("%s\n%s\n+\n%s\n" % ( name,sequence,quality))
     else:
      if (1-deviation)*insert_size <= -a <= (1+deviation)*insert_size:
       name= read.query_name
       sequence=read.query_sequence
       quality=read.qual[:]
       selected_sam.write(read)
       f2.writelines("%s\n%s\n+\n%s\n" % ( name,sequence,quality))
 samfile.close()

def filtered_mp(alignment,filter_sam,filtered_fastq1,filtered_fastq2):
 samfile = pysam.AlignmentFile(alignment)
 filtered_sam = pysam.AlignmentFile(filter_sam, "w", template=samfile)
 with open(filtered_fastq1,'w') as f1:
  with open(filtered_fastq2,'w') as f2:
   for read in samfile.fetch():
    a = read.isize
    if read.is_paired:
     if read.is_read1:
      if a < (1-deviation)*insert_size or a > (1+deviation)*insert_size:
       name= read.query_name
       sequence=read.query_sequence
       quality=read.qual[:]
       filtered_sam.write(read)
       f1.writelines("%s\n%s\n+\n%s\n" % ( name,sequence,quality))
     else:
      if -a < (1-deviation)*insert_size or -a > (1+deviation)*insert_size:
       name= read.query_name
       sequence=read.query_sequence
       quality=read.qual[:]
       filtered_sam.write(read)
       f2.writelines("%s\n%s\n+\n%s\n" % ( name,sequence,quality))
 samfile.close()

def alignresult_mp(mp_sam,summary_file,readlength):
 align_samfile = pysam.AlignmentFile(mp_sam)
 with open(summary_file,'w') as f:
  f.writelines("read_id\tForward\talign\tresults\tread_id\tReverse\talign\tresults\n")
  for read in align_samfile.fetch():
   a = read.isize
   name=read.query_name
   length=read.query_length
   if length == int(readlength):
    if read.is_unmapped:
     if read.is_read1:
      f.writelines("%s R1 unmapped UM\t" % (name))
     else:
      f.writelines("%s R2 unmapped UM\n" % (name))
      # un_map+=1
    else:
     if a == 0:
      if read.is_read1:
       f.writelines("%s R1 mapped-different_scaffolds DS\t" % (name))
      else:
       f.writelines("%s R2 mapped-different_scaffolds DS\n" % (name))
       # map_df+=1
     else:
      if read.is_read1:
       if a == insert_size:
        f.writelines("%s R1 mapped_right RT\t" % (name))
        # map_per+=1
       elif a < 0:
        f.writelines("%s R1 mapped_reversed rev\t" % (name))
        # map_rev+=1
       else:
        deviation1= a - insert_size
        f.writelines("%s R1 mapped_deviation %s\t" % (name,deviation1))
        # map_dev+=1
      else:
       if a == (-1) * insert_size:
        f.writelines("%s R2 mapped_right RT\n" % (name))
        # map_per+=1
       elif a > 0:
        f.writelines("%s R2 mapped_reversed rev\n" % (name))
        # map_rev+=1       
       else:
        deviation2 = a + insert_size
        f.writelines("%s R2 mapped_deviation %s\n" % (name,deviation2))  
        # map_dev+=1  
 align_samfile.close()

if __name__ == "__main__":
 if '-ss' in str(sys.argv):
  select_mp(args.alignment,args.select_sam,args.output_fastq1,args.output_fastq2)
 if '-fs' in str(sys.argv):
  filtered_mp(args.alignment,args.filter_sam,args.filtered_fastq1,args.filtered_fastq2)
 if '-as' in str(sys.argv):
  alignresult_mp(args.select_sam,args.align_select,args.readlength)
 if '-af' in str(sys.argv):
  alignresult_mp(args.filter_sam,args.align_filter,args.readlength)

