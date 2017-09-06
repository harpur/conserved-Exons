#!/usr/bin/env bash









filename=orthoIDsuniq
exec 4<$filename
echo $filename
echo Start
while read -u4 k ; do
	pidlist_sampe=""
	echo $k
	fas='.fas'
	out=$k$fas
	grep $k  ortho_IDs | cut -f2 > $k
	grep -A 1 -f $k focal.fs > $out
done




