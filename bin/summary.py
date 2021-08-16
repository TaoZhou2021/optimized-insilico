#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import pysam
import math
import sys
import argparse

parser = argparse.ArgumentParser(description='input align result file and output summary file')
parser.add_argument('--size', '-s', dest='size', help='insert size of the mate-pair libraries')
parser.add_argument('--input_txt', '-in', dest='input_file', help='input align result file')
parser.add_argument('--output_txt', '-out', dest='output_file', help='output summary file')
args = parser.parse_args()
insert=int(args.size)

def summary_file(input_file,output_file):
 with open(input_file,'r') as f1:
  with open(output_file,'w') as f2:
   line=f1.readlines()[1:]
   (both_UM,onlyone_map,map_df,right,map_dev,map_rev,asm0_10,asm10_20,asm20_30,asm30_40,asm40_100,asm100_x)=[0 for x in range(12)]
   for line in line:
    split=line.split()
    d1=split[3]
    if len(split)>6:
     d2=split[7]
     if d1 == 'UM' and d2 =='UM':
      both_UM+=1
     elif d1 == 'UM' and d2 == 'DS':
      onlyone_map+=1
     elif d1 == 'DS' and d2 == 'UM':
      onlyone_map+=1
     elif d1 == 'DS' and d2 == 'DS':
      map_df+=1
     elif d1 == 'rev' and d2 == 'rev':
      map_rev+=1
     else:
      if d1 =='RT' and d2=='RT':
       right+=1
      else:
       if int(d1)==-1*int(d2):
        percent=abs(int(d1))/insert
        map_dev+=1
        if 0 < percent <= 0.1:
         asm0_10+=1
        elif 0.1 < percent <= 0.2:
         asm10_20+=1
        elif 0.2 < percent <= 0.3:
         asm20_30+=1
        elif 0.3 < percent <= 0.4:
         asm30_40+=1
        elif percent >1:
         asm100_x+=1
        else:
         asm40_100+=1
   with open(output_file,'w') as f2:    
    f2.writelines("insert size: %s mate pairs\t%s\tboth_unmapped\n" % (insert,both_UM))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_right\n" % (insert,right))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_reverse\n" % (insert,map_rev))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_onlyone\n" % (insert,onlyone_map))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_diferent-scaffolds\n" % (insert,map_df))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_deviation\n" % (insert,map_dev))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_deviation within 0 to 10 percent\n" % (insert,asm0_10))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_deviation within 10 to 20 percent\n" % (insert,asm10_20))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_deviation within 20 to 30 percent\n" % (insert,asm20_30))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_deviation within 20 to 30 percent\n" % (insert,asm30_40))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_deviation within 40 to 100 percent\n" % (insert,asm40_100))
    f2.writelines("insert size: %s mate pairs\t%s\tmapped_deviation above 100 percent\n" % (insert,asm100_x))
  f2.close()
 f1.close()

if __name__ == "__main__":
 summary_file(args.input_file,args.output_file)
