#!/usr/bin/bash
inputfile=$1
gvcPath="/usr/local/ANN"
#gvcPath="/disk/Tools/code/ANN/GVC-ANN"
sort -k 1,1 -k 2,2n $inputfile|uniq > $inputfile.sort.uniq
less -S $inputfile.sort.uniq | awk -F '\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$48,"Tumor-1",$40,$49,$50,$60,$7,$16,$17,$27,$56,$57,$58}' > $inputfile.sort.uniq.apart
cat head.title $inputfile.sort.uniq.apart > $inputfile.sort.uniq.apart.head
$gvcPath/Annotation  -p 1,2,2,3,4 -A -b -G -C -t -l -r -i $inputfile.sort.uniq.apart.head  -o $inputfile.WS.snv.simp


