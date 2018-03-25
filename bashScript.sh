#!/bin/bash/

mkdir -p "logs"
mkdir -p fpkm/data-raw
mkdir -p fpkm/data
mkdir -p fpkm/data-paired
echo "Preparing Directory Structure"
for type in `cat types.txt
	mkdir -p "fpkm/data-raw/$type"
	mkdir -p "fpkm/data/$type/mir"
	mkdir -p "fpkm/data/$type/rna"
	cp "supplementary_files/gdc-client.exe" "fpkm/data-raw/$type/"
	cp "supplementary_files/"$type"_manifest.txt" "fpkm/data-raw/$type/"
	cp "supplementary_files/"$type"_metadata.json" "fpkm/data/"
done
printf "\n"
echo "Starting Download"
for type in `cat types.txt`
do
	echo "Downloading $type"
	cd "fpkm/data-raw/$type"
	./gdc-client.exe download -m $type"_manifest.txt"
	echo "$type complete"
	cd ../../..
done
echo "Downloads complete"
echo "Organising files for pairing"
for type in `cat types.txt`
do
	cp fpkm/data-raw/"$type"/*/*quantification.txt fpkm/data/"$type"/mir/
	cp fpkm/data-raw/"$type"/*/*.gz fpkm/data/"$type"/rna/
	printf "$type \t"
done
printf "\n"
echo "Unzipping"
gunzip fpkm/data/ -r 
echo "RNA Files Unzipped\n"
echo "Proceed to run 0_data_Prep.R and 1_Generate_Correlation_Matrices.R through R"