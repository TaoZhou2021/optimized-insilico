#! /home/software/Python3/bin/python3
#-*- coding:utf-8 -*-

import os,sys
import csv,operator
import argparse
import pysam
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio import Align
from interval import Interval
from decimal import Decimal
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = argparse.ArgumentParser(description='align sequences in two fasta files')
parser.add_argument('--input_file1', '-input1', dest='inputfile1', help='input txt file1')
parser.add_argument('--output_file', '-output', dest='outputfile', help='output txt file')
parser.add_argument('--input_file2', '-input2', dest='inputfile2', help='input txt file2')
parser.add_argument('--species_name', '-species', dest='speciesname', help='species name')
parser.add_argument('--genome_size', '-gs', dest='genome_size', help='estimate size of genome')
parser.add_argument('--insert_size', '-is', dest='insertsize', help='insert size of mate-pair')
parser.add_argument('--cfg_file', '-cfg', dest='inputcfg', help='input cfg file')
parser.add_argument('--prefix', '-prefix', dest='prefix', help='prefix of assembly')
parser.add_argument('--suffix', '-suf', dest='suffix', help='suffix of file')
parser.add_argument('--step', '-s', dest='step', help='step of programme')
parser.add_argument('--kmer', '-k', dest='kmer', help='kmer to assemble')
parser.add_argument('--thread', '-t', dest='thread', help='thread to assemble')
parser.add_argument('--ref_file', '-ref', dest='inputref', help='input ref file')
parser.add_argument('--ref_target', '-rt', dest='ref_target', help='input ref target')
parser.add_argument('--busco_dataset', '-bd', dest='busco_dataset', help='dataset of busco')
parser.add_argument('--fasta_file', '-fa', dest='fastafile', help='input fasta file')
parser.add_argument('--input_fq1', '-inputfq1', dest='inputfq1', help='input fastq file1')
parser.add_argument('--input_fq2', '-inputfq2', dest='inputfq2', help='input fastq file2')
parser.add_argument('--database_prefix1', '-db1', dest='db_prefix1', help='prefix of database1')
parser.add_argument('--database_prefix2', '-db2', dest='db_prefix2', help='prefix of database2')
parser.add_argument('--E_value', '-evalue', dest='e_value', help='thread of e value')
parser.add_argument('--ref_file1', '-rf1', dest='ref_file1', help='reference file1')
parser.add_argument('--ref_file2', '-rf2', dest='ref_file2', help='reference file2')
parser.add_argument('--unimap_file', '-uf', dest='unimap_file', help='unimap file')
parser.add_argument('--input_folder', '-input_folder', dest='inputfolder', help='folder of the input')
parser.add_argument('--output_folder', '-out_folder', dest='outputfolder', help='folder of the output')
parser.add_argument('--g1_path', '-g1_folder', dest='g1folder', help='path of g1 folder')
parser.add_argument('--g2_path', '-g2_folder', dest='g2folder', help='path of g2 folder')
args = parser.parse_args()

inputfolder=args.inputfolder
outputfolder=args.outputfolder

folder1='/1-preprocess/'
folder2='/2-soapdenovo/'
folder3='/3-unimap_align/'
folder4='/4-insilico/'
folder5='/5-evaluation/'

Non=[]
chr_list2=[]
chr_list=[]
chr_name=[str(i) for i in range(1,100+1)]
deviation=0.2

def make_dir(file_path):
    if os.path.exists(file_path):
        pass
    else:
        os.makedirs(file_path)

def files_name_listdir_local(file_dir,suffix):
    files_local = []
    for files in os.listdir(file_dir):
        name = os.path.splitext(files)[1][1:]
        size = os.path.getsize(file_dir + files)
        if name == suffix and size != 0:
            files_local.append(files)
    return files_local

def files_name_listdir_local2(file_dir,prefix,number):
    files_local = []
    for files in os.listdir(file_dir):
        start = -1 * int(number)
        name = os.path.splitext(files)[0][start:]
        if name == prefix:
            files_local.append(files)
    return files_local

def preprocess(inputfolder):
    preprocess_dir = outputfolder + folder1
    make_dir(preprocess_dir)
    fq_list=files_name_listdir_local(inputfolder,'fq')
    for fq in fq_list:
        outputfq=preprocess_dir + os.path.splitext(fq)[0]+'-rmdup.fq'
        rmdup_cline = "seqkit rmdup " + inputfolder + fq + " -s -o " + outputfq
        subp = subprocess.Popen(str(rmdup_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    rmdup_list=files_name_listdir_local(preprocess_dir,'-rmdup.fq')
    for rmdup_fq in rmdup_list:
        split=os.path.splitext(rmdup_fq)[0]
        identity = split.split('_')[-1]
        if identity == '1-rmdup':
            rmdup_fq2=rmdup_fq.replace('1-rmdup','2-rmdup')       
            trim_cline="trim_galore --phred33 -o " + preprocess_dir + "--paired " + preprocess_dir + rmdup_fq1 + ' ' + preprocess_dir + rmdup_fq2
            subp = subprocess.Popen(str(trim_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
            subp.wait()
            if subp.poll ==0:
                print("trimming OK")
            else:
                print("trimming failed")
        elif identity =='R1-rmdup':
            rmdup_fq2=rmdup_fq.replace('R1-rmdup','R2-rmdup')     
            trim_cline="trim_galore --phred33 -o " + preprocess_dir + "--paired " + preprocess_dir + rmdup_fq1 + ' ' + preprocess_dir + rmdup_fq2
            subp = subprocess.Popen(str(trim_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")     
            subp.wait()
            if subp.poll ==0:
                print("trimming OK")
            else:
                print("trimming failed")

def sparse_contig(cfg_file,kmer,genome_size,thread,prefix):
    #unimap_file=os.path.splitext(query_file)[0] + "-" + os.path.splitext(ref_file)[0] + "_" + str(i+1) + '.unimap'
    outfolder=outputfolder + folder2
    make_dir(outfolder)
    pregraph_cline = "SOAPdenovo2-63mer sparse_pregraph -g 15 -d 4 -e 4 -R -r 0 -s " + cfg_file + " -K " + kmer + " -z " + genome_size + " -p " + thread +" -o " + outfolder + prefix
    subp = subprocess.Popen(str(pregraph_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")   
    contig_cline = "SOAPdenovo2-63mer contig -g " + outfolder + prefix + " -M 1 -p "+thread
    subp = subprocess.Popen(str(contig_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    

def unimap(ref_file,prefix):
    align_folder=outputfolder+folder3
    make_dir(align_folder)
    query_file = outputfolder + folder2 + prefix + '.scafSeq'
    if ',' in ref_file:
        split=ref_file.split(',')
        num=len(split)
        for i in range(0,num):
            cmd_path='/home/software/unimap/uniamp'
            unimap_file=align_folder + os.path.splitext(query_file)[0] + "-" + os.path.splitext(ref_file)[0] + "_" + str(i+1) + '.unimap'
            unimap_cline = cmd_path + " -cxasm5 -t20 " + query_file + " " + ref_file + " > " + unimap_file
            subp = subprocess.Popen(str(unimap_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    else:
        cmd_path='/home/software/unimap/uniamp'
        unimap_file= align_folder + os.path.splitext(query_file)[0] + "-" + os.path.splitext(ref_file)[0] + '.unimap'
        unimap_cline = cmd_path + " -cxasm5 -t20 " + query_file + " " + ref_file + " > " + unimap_file
        subp = subprocess.Popen(str(unimap_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")

def count_contigs(inputfolder,inputfile):
    contig_list=[]
    num=0
    length=0
    with open(inputfolder+inputfile,'r') as f1:
        lines=f1.readlines()
        for line in lines:
            split=line.split()
            contig_name=split[1]
            contig_len=int(split[9])
            if contig_name not in contig_list:
                num+=1
                length+=contig_len
    return num,length    

def compare_ref(inputfolder,suffix,outputfile):
    file_list=files_name_listdir_local(inputfolder,suffix)
    list1=[]
    num_list=[]
    for unimap_file in file_list:
        ref_species=os.path.splitext(unimap_file)[0].split('-')[1]
        num,length=count_contigs(inputfolder,unimap_file)
        percent= length/num
        list1.append((ref_species,num,length,percent))
    a=sorted(list1, key=lambda x: x[3], reverse=True)
    element_num=len(a)    
    with open(outputfile,'w') as f1:
        f1.writelines("ref_name\taligned_contig_number\taligned_length\taligned_length_per_contig\n")
        for i in range(0,element_num):  
            ref_name=a[i][0]
            contig_num=a[i][1]
            contig_len=a[i][2]
            contig_percent=a[i][3]
            f1.writelines("%s\t%s\t%s\t%s\n" % (ref_name,contig_num,contig_len,contig_percent))
    
def cross2mate_pair(insert_size):
    align_folder=outputfolder + folder3
    cross_folder=outputfolder + folder4
    fq_folder = outputfolder + folder1
    for inputfq1 in fq_list:
        inputfq2 = inputfq1.replace('val_1','val_2')
        ref_file=compare_ref(align_folder,'unimap')
        cross_cline='cross-mates -s -c 10 -l 100 -t 20 -i '+ insert_size+' '+ref_file+' '+inputfq1+' '+inputfq2+ ' -o ' + cross_folder
        subp = subprocess.Popen(str(blastn_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")

def map_scaff(cfg_file,kmer,thread,prefix):
    outfolder = outputfolder + folder2
    map_cline = "SOAPdenovo2-63mer map -s "+ cfg_file + "-K "+kmer +"-p "+ thread +" -g " + outfolder + prefix
    subp = subprocess.Popen(str(map_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
    scaff_cline = "SOAPdenovo-63mer scaff " + "-p " + thread +" -g " + outfolder + prefix
    subp = subprocess.Popen(str(scaff_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")

def quast_evaluation(ref_target,busco_dataset,thread,prefix):
    outfolder = outputfolder + folder5
    inputfa = outputfolder + folder2 +prefix+'.scafSeq'
    make_dir(outfolder)
    python3_path ='/home/software/Python3/bin/python3'
    quast_cline = python3_path+ 'quast.py -t '+ thread + ' -R '+ ref_target + ' -o ' + outfolder +'quast-result ' + inputfa 
    subp = subprocess.Popen(str(quast_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")    
    busco_cline = python3_path + '/home/room/users/zhoutao/busco-4.1.4/bin/busco -i ' + inputfa + ' -l ' + busco_dataset + ' -o busco-result' + ' -m genome --cpu 15'
    subp = subprocess.Popen(str(busco_cline),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")  

def paf2list(paf_file):
    chr_list1=[]
    unsorted_csv= os.path.splitext(paf_file)[0]+'-unsorted.csv'
    sorted_csv= os.path.splitext(paf_file)[0]+'-sorted.csv'
    select_region_csv= os.path.splitext(paf_file)[0]+'-select_region.csv'
    csvFile=open(unsorted_csv,'w+',newline='',encoding='utf-8')
    writer = csv.writer(csvFile)
    csvRow = []
    f = open(paf_file,'r')
    for line in f:
        split=line.split()
        chr=split[0]
        csvRow =split[0:11]
        if chr in chr_name:
            csvRow[0] = chr.zfill(12)
            writer.writerow(csvRow)
        else:
            writer.writerow(csvRow)
    f.close()
    csvFile.close()
    data = pd.read_csv(unsorted_csv,header=None,names=['a','b','c','d','e','f','g','h','i','j','k'],engine='python')
    data.sort_values(by=['f', 'h'], inplace=True, ascending=[True, True],kind='quicksort')
    data.to_csv(sorted_csv,encoding='utf-8',index=False)
    f1=open(select_region_csv,'w',encoding='utf-8')
    f=open(sorted_csv,'r',encoding='utf-8')
    data=f.readlines()[1:]
    for line in data:
        split=line.split(',')
        chr=split[5]
        start=int(split[7])
        end=int(split[8])
        field_new=Interval(start,end)
        if chr not in chr_list1:
            chr_list1.append(chr)
            field=Interval(0,1,closed=False)
            if not field_new.overlaps(field):
                f1.writelines(line)
                field1=field_new
                start1=start
                end1=end
        else: 
            if not field_new.overlaps(field1):
                f1.writelines(line)
                field1=field_new
                start1=start
                end1=end
            else:
                if start in field1 and end not in field1:
                   new_start=end1+1
                   length1=end1-start
                   length2=end-end1
                   direction=split[4]
                   if direction == '+':
                       start0=int(split[3])-length2
                       new_list=split[:2]+[str(start0)]+split[3:7]+[str(new_start),str(end),str(length2),str(length2)]
                       new_line=','.join(new_list)                    
                       f1.writelines("%s\n" % (new_line))
                   else:
                       end0=int(split[2])+ length2
                       new_list=split[:3]+[str(end0)]+split[4:7]+ [str(new_start),str(end),str(length2),str(length2)]
                       new_line=','.join(new_list)
                       f1.writelines("%s\n" % (new_line))
                   field1=Interval(start1,end)
                   start1=new_start
                   end1=end 
    f.close()
    f1.close()
    return select_region_csv

def paf2list1(paf_file):
    chr_list1=[]
    unsorted_csv= os.path.splitext(paf_file)[0]+'unsorted.csv'
    sorted_csv=os.path.splitext(paf_file)[0]+'sorted.csv'
    select_region_csv=os.path.splitext(paf_file)[0]+'select_region.csv'
    csvFile=open(unsorted_csv,'w+',newline='',encoding='utf-8')
    writer = csv.writer(csvFile)
    csvRow = []
    f = open(paf_file,'r')
    for line in f:
        split=line.split()
        chr=split[0]
        csvRow =split[0:11]
        if chr in chr_name:
            csvRow[0] = chr.zfill(12)
            writer.writerow(csvRow)
        else:
            writer.writerow(csvRow)
    f.close()
    csvFile.close()
    data = pd.read_csv(unsorted_csv,header=None,names=['a','b','c','d','e','f','g','h','i','j','k'],engine='python')
    data.sort_values(by=['a', 'c'], inplace=True, ascending=[True, True],kind='quicksort')
    data.to_csv(sorted_csv,encoding='utf-8',index=False)    
    f1=open(select_region_csv,'w',encoding='utf-8')
    f=open(sorted_csv,'r',encoding='utf-8')
    data=f.readlines()[1:]
    for line in data:
        split=line.split(',')
        chr=split[0]
        start=int(split[2])
        end=int(split[3])
        field_new=Interval(start,end)
        if chr not in chr_list1:
            chr_list1.append(chr)
            field=Interval(0,1,closed=False)
            if field_new.overlaps(field) is False:
                f1.writelines(line)
                field1=field_new
        else: 
            if field_new.overlaps(field1) is False:
                f1.writelines(line)
                field1=field_new
    f.close()
    f1.close()  
    return sorted_csv    

def get_Length(paf_file,Length_file):
    chr1_list=[] 
    select_region_csv=paf2list(paf_file)    
    with open(Length_file,'w+',encoding='utf-8') as f1:
        with open(select_region_csv,'r',encoding='utf-8') as f:
            line=f.readlines()
            for line in line:
                split=line.split(',')
                chr1=split[0]
                if chr1 not in chr1_list:
                    chr1_list.append(chr1)
                    new_line=line.rstrip('\n')+','+'0'+','+'0'+','+'A'+','+'N'+'\n'
                    f1.writelines(new_line)
                    chr2=split[5]
                    s1=split[2]
                    e1=split[3]
                    s2=split[7]
                    e2=split[8]
                    d=split[4]
                else:
                    chr2_next=split[5]
                    s1_next=split[2]
                    e1_next=split[3]
                    s2_next=split[7]
                    e2_next=split[8]
                    d_next=split[4]
                    if chr2_next==chr2:                 
                        L1=int(s1_next)-int(e1)-1
                        if d_next==d:
                            if d_next=='+':
                                L2=int(s2_next)-int(e2)-1
                                DL=abs((L2-L1)*5)
                                if (1-deviation)*L1 < L2 < (1+deviation)*L1:
                                    new_line=line.rstrip('\n')+','+str(L1)+','+str(L2)+','+'T'+','+str(DL)+'\n'
                                    f1.writelines(new_line)
                                else:
                                    new_line=line.rstrip('\n')+','+str(L1)+','+str(L2)+','+'F'+','+str(DL)+'\n'
                                    f1.writelines(new_line)
                            else:
                                L2=int(s2)-int(e2_next)-1
                                DL=abs((L2-L1)*5)
                                if (1-deviation)*L1 < L2 < (1+deviation)*L1:
                                    new_line=line.rstrip('\n')+','+str(L1)+','+str(L2)+','+'T'+','+str(DL)+'\n'
                                    f1.writelines(new_line)
                                else:
                                    new_line=line.rstrip('\n')+','+str(L1)+','+str(L2)+','+'F'+','+str(DL)+'\n'
                                    f1.writelines(new_line)
                        else:
                            new_line=line.rstrip('\n')+','+str(L1)+','+'NA'+','+'NA'+','+'NA'+'\n'
                            f1.writelines(new_line)
                            s1=s1_next
                            e1=e1_next
                            s2=s2_next
                            e2=e2_next
                            d=d_next
                            chr2=chr2_next
                    else:
                        L1=int(s1_next)-int(e1)-1
                        new_line=line.rstrip('\n')+','+str(L1)+','+'0'+','+'1'+','+'N'+'\n'
                        f1.writelines(new_line)
                        s1=s1_next
                        e1=e1_next
                        s2=s2_next
                        e2=e2_next
                        d=d_next
                        chr2=chr2_next
            f.close()
            f1.close()

def get_Length2(paf_file,Length_file):
    chr1_list=[]
    select_region_csv=paf2list(paf_file)
    with open(Length_file,'w+',encoding='utf-8') as f1:
        with open(select_region_csv,'r',encoding='utf-8') as f:
            line=f.readlines()
            for line in line:
                split=line.split(',')
                chr1=split[5]
                if chr1 not in chr1_list:
                    chr1_list.append(chr1)
                    new_line=line.rstrip('\n')+','+'N'+'\n'
                    f1.writelines(new_line)
                    chr2=split[5]
                    e2=int(split[8])
                else:
                    chr2_next=split[5]
                    s2_next=int(split[7])
                    e2_next=int(split[8])
                    if chr2_next==chr2:
                        L2=s2_next-e2-1
                        new_line=line.rstrip('\n')+','+str(L2)+'\n'
                        f1.writelines(new_line)
                        e2=e2_next
            f.close()
            f1.close()

def get_fa(fasta_file,length_txt,output_file):
    #record_iterator=SeqIO.parse(fasta_file,'fasta')
    fasta=pysam.FastaFile(fasta_file)
    with open(output_file,'w',encoding='utf-8') as f:
        with open(length_txt,'r',encoding='utf-8') as f1:
            line=f1.readlines()
            i=0
            #j=0
            for line in line:
                split=line.rstrip().split(',')
                chr=(split[0])
                start=int(split[2])
                end=int(split[3])
                numbers=split[11]
                de=split[4]
                #if j==0 : 
                    #record=next(record_iterator)
                    #j=j+1
                if chr.startswith('0'):
                    chr=int(chr)
                #while record.id != str(chr):
                    #record=next(record_iterator)
                if numbers == 'N':
                    #head='>'+record.id+' '+str(start)+'-'+'\n'
                    i+=1
                    if i ==1:
                        head='>scaffold'+str(i)+'\n'
                        i+=1
                    else:
                        head='\n>scaffold'+str(i)+'\n'
                        i+=1
                    f.writelines(head)
                    if de == '+':
                        #seq=record.seq[start:end]
                        seq=fasta.fetch(chr,start,end)
                        f.writelines(seq)
                    else:
                        seq=fasta.fetch(chr,start,end)
                        seq_new=revcomp(seq)
                        f.writelines(seq_new)
                else:
                    if int(numbers) >= 0:
                        if de == '+':
                            seq1='N'*int(numbers)+fasta.fetch(chr,start,end)
                            f.writelines(seq1)
                        else:
                            seq='N'*int(numbers)+fasta.fetch(chr,start,end)
                            seq1=revcomp(seq)
                            f.writelines(seq1)
                    #elif int(numbers) < 0:
                        #start_new=start-int(numbers)
                        #if de == '+':
                           # seq2=fasta.fetch(chr,start_new,end)
                           # f.writelines(seq2) 
                        #else:
                           # s1=1-int(numbers)
                           # seq=fasta.fetch(chr,start,end)
                           # seq2=revcomp(seq)[s1:]
                           # f.writelines(seq2)

def complement(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))

def revcomp(seq): 
    return complement(seq)[::-1]

def get_fa1(fasta_file,length_txt,output_file):
    record_iterator=SeqIO.parse(fasta_file,'fasta')
    with open(output_file,'w',encoding='utf-8') as f:
        with open(length_txt,'r',encoding='utf-8') as f1:
            line=f1.readlines()
            i=0
            j=0
            for line in line:
                split=line.split(',')
                chr=(split[0])
                start=int(split[2])
                end=int(split[3])
                numbers=split[11]
                de=split[13]
                if j==0 : 
                    record=next(record_iterator)
                    j=j+1
                if chr.startswith('0'):
                    chr=int(chr)
                while record.id != str(chr):
                    record=next(record_iterator)
                if de != 'A':
                    seq='n'*int(numbers)+record.seq[start:end]
                    f.writelines(seq)
                else:
                    if i==0:
                        head='>'+record.id+' '+str(start)+'-'+'\n'
                        seq1=record.seq[start:end]
                        f.writelines(head)
                        f.writelines(seq1)
                        i=i+1
                    else:
                        head='\n'+'>'+record.id+' '+str(start)+'-'+'\n'
                        seq1=record.seq[start:end]
                        f.writelines(head)               
                        f.writelines(seq1)                              


def optimized_insilico(step):
    if step == 'compare_ref':
        compare_ref(args.inputfolder,args.suffix,args.outputfile)
    elif step == 'paf2list':
        paf2list(args.unimap_file)
    elif step == 'get_fa':
        get_fa(args.fastafile,args.inputfile1,args.outputfile)
    elif step == 'get_length':
        get_Length2(args.unimap_file,args.outputfile)
##
if __name__=="__main__":
    optimized_insilico(args.step)
    #preprocess(args.inputfolder)
    #sparse_contig(args.inputcfg,args.kmer,args.genome_size,args.thread,args.prefix)
    #unimap(args.inputref,args.prefix)
    #compare_ref(args.inputfolder,args.suffix,args.outputfile)
    #cross2mate_pair(args.insertsize)
    #map_scaff(args.inputcfg,args.kmer,args.thread,args.prefix)
    #quast_evaluation(args.ref_target,args.busco_dataset,args.thread,args.prefix)
    #get_Length(paf_file,length_file)
    #get_fa(fasta_file,length_file,output_file)

##Linux_shell##
## python3 paf2fa.py paf_file fasta_file output_file unsorted.csv sorted.csv select_region.csv length.txt
