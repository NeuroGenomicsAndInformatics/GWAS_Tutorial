#!/bin/bash
# This script splits same position warnings into 2 columns. This is useful for removing multiallelics.

LOGFILE=$1

perl -p -e 's/\n/ /' $LOGFILE \
| perl -p -e 's/Warning:/\n/g' \
| perl -ne "print if /'/" \
| awk '{print $2, $4}' \
| sed 's/'\''//g' \
| awk '{print $1}' > ${LOGFILE%.log}_Col1.txt

perl -p -e 's/\n/ /' $LOGFILE \
| perl -p -e 's/Warning:/\n/g' \
| perl -ne "print if /'/" \
| awk '{print $2, $4}' \
| sed 's/'\''//g' \
| awk '{print $2}' > ${LOGFILE%.log}_Col2.txt
