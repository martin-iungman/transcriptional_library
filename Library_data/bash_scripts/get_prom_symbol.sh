#!/bin/sh
# Saca el codigo FP de los que lo tienen. Y despues te quedas con los "convencionales": solo symbol_numero. 
sed -E "s/^FP[[:digit:]]{6}_//" seq_id.txt | grep -E "^[[:alnum:]]+_[[:digit:]]{1}$" > promoter_id_symbol.txt

