code_num=$1

mkdir tmp/
cd tmp/
  
git clone https://github.com/statgen/METAL
cd METAL
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make test
make install
cd ../..

git clone https://github.com/seanjosephjurgens/UKBB_200KWES_CVD/
chmod +x ./UKBB_200KWES_CVD/metal_meta.sh
chmod +x ./UKBB_200KWES_CVD/genotype_count_based_meta_analysis_dominant.R

git clone --branch patch-1 https://github.com/seanjosephjurgens/TOPMed_AFib_pipeline/
chmod +x ./TOPMed_AFib_pipeline/DNANexus/PheWAS/GC_based_meta/
