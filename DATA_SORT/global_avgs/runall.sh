#!/bin/bash

for ifile in `ls *GOGA2*.py` ; do
   python $ifile
done

for ifile in `ls *BEST*.py` ; do
   python $ifile
done


