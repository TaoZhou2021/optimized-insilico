#! /bin/bash

path=$(pwd)

export PATH=$PATH:/home/software/SOAPdenovo2-master
export PATH=$PATH:/home/software/TrimGalore-0.6.4
export PATH=$PATH:/home/software/Python3/bin
export PATH=$PATH:$path
export PATH=$PATH:/home/software/unimap
export PERL5LIB=$PERL5LIB:/home/software/insilico/lib
#export PATH=$PATH:/home/room/users/zhoutao/insilico/cros/bin
#export BUSCO_CONFIG_FILE=/home/room/users/zhoutao/busco-4.1.4/config/myconfig.ini
#export AUGUSTUS_CONFIG_PATH=/home/room/users/zhoutao/conda/envs/busco_env/config
#export PATH=$PATH:/home/software/quast-master

./opt-insilico -t 50 -i 300,500,1000,1500,2000,5000,10000,20000,50000,100000,200000 -c 30 -P genome -d trimmed -C my_config -m takifugu_flavidus.fa.masked,takifugu_rubripes.fa.masked,tetradon_nigroviridis_rm.fa,Mola_rm.fa -f corrected_srr2_1.fq -k 37 -G 400000000 -Q quast_result -B busco_result -b /home/software/vertebrata_odb10 

#./opt-insilico test_1.fq test_2.fq

exit 0
