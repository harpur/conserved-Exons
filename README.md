# Conserved Exon Identification 

The Conserved Exons (CE) pipeline identifies conserved- and non-conserved sites within a protein alignment for a focal species against closely and more distantly-related taxa. The piepline makes use of [OrthoDB v 9.1](http://www.orthodb.org/v9/?page=downloads). This database contains orthologous proteins for  172 vertebrates, 133 arthropods, 227 fungi, 25 basal metazoans, 3663 bacteria and 31 plants. While a great resource, the database does change IDs and formats between iterations, so use extreme caution if you are applying this pipeline to any database beyond v9.1.

By downloading the peptide data, extracting relevant taxa, aligning the output, and outputing conservation scores, this pipeline identifies conserved peptides across a given phylogeny within known orthologs. 


This workflow is rough (v0.01) at the moment and has some scripting steps I need to remove and assumptions I need to test.


# Requirements 
* R 3.1 or higher, with the following packages installed:
	* taxize
* [muscle 3.8.3](http://www.drive5.com/muscle/downloads.htm)
* Database files from OrthoDB v9.1:
	* odb9_species.tab
	* odb9_OG2genes.tab 
	* odb9_fasta_metazoa.tgz (very large data set)

# Extract ORTHODB IDs 
I've provided an R script that extracts the relevant taxonomic levels of all species within orthodb v9.1 using `taxize`. It is important to check the output of this file by eye before proceeding. Occasionally, `taxize` can misplace taxonomic levels (e.g. *Loa loa* is apparently a Hymenopteran). 

```R
output_taxa.r
```
Depending on how you've run this script, it will output species lists for you to check and continue with that looks like this (focal_spp)

```
7460 Apis mellifera 12101 83040 C
7461 Apis cerana 8915 64537 C
7462 Apis dorsata 10452 75545 C
7463 Apis florea 11816 80899 C
7493 Trichogramma pretiosum 11988 43779 C
```
The first column contains the OrthoDB species ID and the second is the species name and can be cut out of the file (in this case, by file is 'lower_spp'). 


```
cut -f1 -d' ' lower_spp > lower_spp_IDs

```


# ID Species of interest, extract genes 

Here, I am extracting all *Apis mellifera* genes
```
grep 7460: odb9_OG2genes.tab > ortho_test_genes 
cut -f1 ortho_test_genes > orthoIDs
sort orthoIDs | uniq -u > orthoIDsuniq

```


# Extract gene(s) of interest and align

## Extract sequences
I first make a new .fasta with only species of interest
```
grep -A 1 -f  lower_spp_IDs /home/blencowe/blencowe31/harpurbr/orthodb/metazoa/metas.fs > focal.fs
```

I then extract sequences of interest
```			
grep -f orthoIDsuniq  odb9_OG2genes.tab  > ortho_IDs	
./extract_fasta.sh
```

## Align 
I use `muscle v3.8` to align the orthologs. As of right now, I use default options, but this will be need to be optimized. 

```
./muscle.sh
```

This results in an alignment file in .fasta format.



# Quantify conservation 
I use [Capra and Singh's 2007 method](http://compbio.cs.princeton.edu/conservation/score.html) to quantify conservation across the alignment 



<!---
python score_conservation.py -s js_divergence -w 3 -d
      swissprot.distribution -o alignment.scores alignment.clustal

python2 score_conservation.py -s js_divergence -w 3 -d swissprot.distribution -o alignment.scores /home/blencowe/blencowe31/harpurbr/orthodb/out_test_focal.afa 





























#Trim these files do get species DB------------------------

cut -f1 -d' ' focal_spp > focal_spp_IDs
grep -A 1 -f  focal_spp_IDs out_test.fas > out_test_focal.fas

cut -f1 -d' ' lower_spp > lower_spp_IDs
cut -f1 -d' ' non_ap_spp > non_ap_spp_IDs




				#dl orthodb - probably best to start with honey bee peptides and work out from there 

				#7460	Apis mellifera	12101	83040	C


				grep 7460: odb9_OG2genes.tab > ortho_test_genes #NOTE: in this file are duplicated within genomes come back to this.

				cut -f1 ortho_test_genes > orthoIDs

				sort orthoIDs | uniq -u > orthoIDsuniq

				grep -f orthoIDsuniq  odb9_OG2genes.tab  > ortho_IDs


#test case ---
head -n426 ortho_IDs | cut -f2 > GID
grep '^EOG091B00HZ'  ortho_IDs | cut -f2 > GID #GB53150

#Hymenopterans----
grep -A 1 -f  GID /home/blencowe/blencowe31/harpurbr/orthodb/metazoa/metas.fs > out_test.fas

cut -f1 -d' ' focal_spp > focal_spp_IDs

grep -A 1 -f  focal_spp_IDs out_test.fas > out_test_focal.fas

muscle -in out_test_focal.fas -out out_test_focal.afa



sed -i 's/ /_/g' out_test_focal.afa


#remove the stupid extra lines:
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' out_test_focal.afa > out_test_focal_1.afa


grep -A 1 -f  non_ap_spp_IDs out_test_focal_1.afa > out_test_non_app.fas


#Bees----

#extract alignment of bees only ----
cut -f1 -d' ' lower_spp > lower_spp_IDs
grep -A 1 -f  lower_spp_IDs out_test_focal_1.afa > out_test_local.fas



#align bees only-------

muscle -in out_test_local.fas -out 
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' bee_align.fas > bee_align1.fas
grep -A 1 7460 bee_align1.fas> apis.fas

#I think this is the way to go for now....demonstrate that bee-specific regions are more/less likely to be conserved






#cons score---
Capra JA and Singh M. Predicting functionally important residues from
sequence conservation. Bioinformatics. 23(15):1875-82, 2007.  





python2 score_conservation.py -s js_divergence -w 100 -d swissprot.distribution -g 0.8 -p FALSE -o alignment_non_ap.scores /home/blencowe/blencowe31/harpurbr/orthodb/out_test_non_app.fas


python2 score_conservation.py -s js_divergence -w 100 -d swissprot.distribution -g 0.8 -p FALSE -o alignment_local.scores /home/blencowe/blencowe31/harpurbr/orthodb/out_test_local.fas





system("python2 score_conservation.py -s js_divergence -w 3 -d swissprot.distribution -g 0.4 -p FALSE -o alignment_non_ap.scores /home/blencowe/blencowe31/harpurbr/orthodb/bee_align1.fas.fas")



# ID exons with high values -----
R here






#################
#cons score---
Capra JA and Singh M. Predicting functionally important residues from
sequence conservation. Bioinformatics. 23(15):1875-82, 2007.  

python score_conservation.py -s js_divergence -w 3 -d
      swissprot.distribution -o alignment.scores alignment.clustal

python2 score_conservation.py -s js_divergence -w 3 -d swissprot.distribution -o alignment.scores /home/blencowe/blencowe31/harpurbr/orthodb/out_test_focal.afa 





cons_score -s js_divergence -w 3 -n TRUE -d /home/blencowe/blencowe31/harpurbr/git/cons_scores/conservation_code/distributions/swissprot.distribution  -o test_scores out_test_focal.afa



grep 7370:000927 odb9_genes.tab
...
7370:000927	7370	17007358
...

odb9_genes.tab
1.	Ortho DB unique gene id (not stable between releases)
2.	organism tax id
3.	protein sequence id, as downloaded together with the sequence
4.	Uniprot id, evaluated by mapping
5.	ENSEMBL gene name, evaluated by mapping
6.	NCBI gid, evaluated by mapping
7.	description, evaluated by mapping



grep -A 2 7370:000927 /home/blencowe/blencowe31/harpurbr/orthodb/metazoa/7370.fs

ls /home/blencowe/blencowe31/harpurbr/orthodb/metazoa/7370*
/home/blencowe/blencowe31/harpurbr/orthodb/metazoa/7370.fs


 cat *.fs > metas.fs


grep -A 2 7370:000927 /home/blencowe/blencowe31/harpurbr/orthodb/metazoa/metas.fs


-->
