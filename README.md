optimized-insilico
===

Dependencies:    
---
1. pysam ( sam2fq.py)
2. itertools (simulate.py)

Installation
---
```shell
git clone https://github.com/TaoZhou2021/optimized-insilico.git or git clone git@github.com:TaoZhou2021/optimized-insilico.git
```    

There are three scripts in the optimized-insilico/bin/
---
1. {sam2fq.py}    -------------- select the mate-pairs and also filtered mate-pairs(optional)

2. {simulate.py}  -------------- simulate the ancient DNA reads using one strand of paired-end reads
 
3. {summary.py}   -------------- summary the results of mate-pairs aligned to the reference(select mate-pairs or filtered mate-pairs)

Getting Started
---
1.{sam2fq.py}

*select mate-pairs without output their align results*                                                                                                                             

 (1) only select mate-pairs:    
``` 
python3 sam2fq.py -s 500 -dev 0.2 -length 100 -in 500-mp.sam -ss select-500-mp.sam -out1 select-500-mp.R1.fq -out2 select-500-mp.R2.fq  
```
 (2) select mate-pairs and also output of filtered mate-pairs without their align results:    
```
python3 sam2fq.py -s 500 -dev 0.2 -length 100 -in 500-mp.sam -ss select-500-mp.sam -out1 select-500-mp.R1.fq -out2 select-500-mp.R2.fq -fs filter-500-mp.sam -f1 filter-500-mp.R1.fq -f2 filter-500-mp.R2.fq
```
*select mate-paris with output their align results*

 (1) only select mate-pairs and its summary file:     
```
python3 sam2fq.py -s 500 -dev 0.2 -length 100 -in 500-mp.sam -ss select-500-mp.sam -out1 select-500-mp.R1.fq -out2 select-500-mp.R2.fq -as align_select-500-mp.txt
```
 (2) select mate-pairs, output of filtered mate-pairs with their align results:    
```
python3 sam2fq.py -s 500 -dev 0.2 -length 100 -in 500-mp.sam -ss select-500-mp.sam -out1 select-500-mp.R1.fq -out2 select-500-mp.R2.fq -as align_select-500-mp.txt -fs filter-500-mp.sam -f1 filter-500-mp.R1.fq -f2 filter-500-mp.R2.fq -af algin_filter-500-mp.txt
```

2.{simulate.py}

*simulate ancient DNA reads*

 (1) simulate ancient DNA reads (80 to 100 bp long)
``` 
python3 simulate.py -s 80 -e 100 -in my.fastq -out simulate-aDNA.fq.gz
```

3.{summary.py}

*summary the align results of mate-pairs*

 (1) summary the align results of selected mate-pairs
```
python3 summary.py -s 500 -in align_select-500-mp.txt -out summary-select-500-mp.txt
```
 (2) summary the align results of filtered mate-pairs
```
python3 summary.py -s 500 -in align_filter-500-mp.txt -out summary-filter-500-mp.txt
```
