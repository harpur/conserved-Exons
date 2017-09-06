###
# Output Taxa IDs 
###



#extract hymenopteran sequences:
library("taxize")

odb = read.table('odb9_species.tab',sep="\t",fill=T,header=F)

x = classification(odb$V2, db = 'itis')

fams = unlist(lapply(x,function(x) return(unlist(x)[11])))
fams.super = unlist(lapply(x,function(x) return(unlist(x)[14])))

odb.ap = odb[which(fams.super=="Apoidea"),]
odb.hym = odb[which(fams=="Hymenoptera"),]
odb.hym.non.ap = odb.hym[-which(odb.hym$V2 %in% odb.ap$V2),]
odb.hym.non.ap = odb.hym.non.ap[-1,]


write.table(odb[which(fams=="Hymenoptera"),], file="focal_spp",col.names=F,row.names=F,quote=F)
write.table(odb[which(fams.super=="Apoidea"),], file="lower_spp",col.names=F,row.names=F,quote=F)
write.table(odb.hym.non.ap, file="non_ap_spp",col.names=F,row.names=F,quote=F)

