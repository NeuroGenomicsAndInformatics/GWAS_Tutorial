#!/bin/bash
LOGFILE=$1
echo -n "" > ${LOGFILE}.MAcheck
perl -p -e 's/\n/ /' $LOGFILE | perl -p -e 's/(Warning:)/\n/g'| grep Variants | while read LINE; do
	VAR1=$(echo $LINE | cut -d \' -f2)
	VAR2=$(echo $LINE | cut -d \' -f4)
	REF1=$(echo $VAR1 | cut -d: -f3 | cut -c1-1)
	REF2=$(echo $VAR2 | cut -d: -f3 | cut -c1-1)
	[[ $REF1 != $REF2 ]] && echo "$VAR1 $VAR2" >> ${LOGFILE}.MAcheck
done	