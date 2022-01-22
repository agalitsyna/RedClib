#!/bin/bash

cd ./bin
wget https://github.com/agalitsyna/rklib/archive/refs/heads/main.zip; unzip main.zip; cd rklib-main
make; cd ../
cp rklib-main/bin/* ./
rm main.zip; rm -rf rklib-main
cd ../

wget https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/hisat2_extract_splice_sites.py -O ./bin/hisat2_extract_splice_sites.py
chmod 777 ./bin/hisat2_extract_splice_sites.py

echo "Binaries prepared!"
