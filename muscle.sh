#!/usr/bin/env bash






filename=orthoIDsuniq
exec 4<$filename
echo $filename
echo Start
while read -u4 k ; do
	pidlist_sampe=""
	echo $k
	fas='.fas'
	outfas='al.fas'
	out=$k$outfas
	ins=$k$fas
	out1='1'
	out1=$out$out1
	outap='apis'
	outap=$outap$out
	muscle -in $ins -out $out
	awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' $out > $out1
	rm $out
	mv $out1 $out
	sed -i 's/ /_/g' $out
	grep -A 1 7460 $out > $outap
done




