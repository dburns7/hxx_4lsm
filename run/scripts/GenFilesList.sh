#!/usr/bin/env bash
/mnt/hadoop/store/user/dburns/data_LPC/pu

# $1 = source dir
# $2 = name

TAR1='/mnt/hadoop/store/user/dburns/DAS_Full'
TAR2='/mnt/hadoop/store/user/dburns/data_LPC'

python filesList.py $TAR1/GluGluToHToZZTo4L GluGluToHToZZTo4L.txt
mv GluGluToHToZZTo4L.txt $TAR2

python filesList.py $TAR2/pu pu.txt
mv pu.txt $TAR2

python filesList.py $TAR2/no_pu no_pu.txt
mv no_pu.txt $TAR2

python filesList.py $TAR2/pu_10k pu_10k.txt
mv pu_10k.txt $TAR2

python filesList.py $TAR2/pu_nocuts pu_nocuts.txt
mv pu_nocuts.txt $TAR2

python filesList.py $TAR2/nopu_nocuts nopu_nocuts.txt
mv nopu_nocuts.txt $TAR2

python filesList.py $TAR2/pu_nocuts_START53_V19 pu_nocuts_START53_V19.txt
mv pu_nocuts_START53_V19.txt $TAR2

python filesList.py $TAR2/pu_nocuts_100k pu_nocuts_100k.txt
mv pu_nocuts_100k.txt $TAR2

python filesList.py $TAR2/hxx_1GeV_nocuts hxx_1GeV_nocuts.txt
mv hxx_1GeV_nocuts.txt $TAR2

python filesList.py $TAR2/hxx_10GeV_nocuts hxx_10GeV_nocuts.txt
mv hxx_10GeV_nocuts.txt $TAR2

python filesList.py $TAR2/hxx_100GeV_nocuts hxx_100GeV_nocuts.txt
mv hxx_100GeV_nocuts.txt $TAR2

python filesList.py $TAR2/hxx_500GeV_nocuts hxx_500GeV_nocuts.txt
mv hxx_500GeV_nocuts.txt $TAR2

python filesList.py $TAR2/hxx_1000GeV_nocuts hxx_1000GeV_nocuts.txt
mv hxx_1000GeV_nocuts.txt $TAR2

python filesList.py $TAR2/4emuta 4emuta.txt
mv 4emuta.txt $TAR2
