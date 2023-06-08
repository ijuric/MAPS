#!/usr/bin/bash

filename="/home/jurici/MAPS/PLAC-Seq_datasets/test_dataset2/feather_output/test_current/test.all.bedpe"
cat $filename | awk '{
if($9=="+" && $10=="+") printf "%s\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", $7, 0, $1, $2, 0, 0, $4, $5, 1, 60, 60
if($9=="+" && $10=="-") printf "%s\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", $7, 0, $1, $2, 0, 16, $4, $5, 1, 60, 60
if($9=="-" && $10=="+") printf "%s\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", $7, 16, $1, $2, 0, 0, $4, $5, 1, 60, 60
if($9=="-" && $10=="-") printf "%s\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", $7, 16, $1, $2, 0, 16, $4, $5, 1, 60, 60
}' | sort -k3,3 -k7,7 - > sample_name.long.txt
