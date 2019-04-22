#!/bin/bash

#  Without log
rm -f -r scatterplots/without-log/*
./scatterplots.R --sample combo7 >scatterplots/without-log/scatterplots-without-log.txt 2>&1
scatterplots/clean-logfile.pl <scatterplots/without-log/scatterplots-without-log.txt >scatterplots/without-log/correlations.tsv
mv *.eps scatterplots/without-log/

#  With log
rm -f -r scatterplots/with-log/*
./scatterplots.R --sample combo7 --log >scatterplots/with-log/scatterplots-with-log.txt 2>&1
scatterplots/clean-logfile.pl <scatterplots/with-log/scatterplots-with-log.txt >scatterplots/with-log/correlations.tsv
mv *.eps scatterplots/with-log/

