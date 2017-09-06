###
# Find high and low conservation
###



#loads input of an aligned conservation score and identifies the highest and lowest regions






#Libraries ---------





#Functions ---------

MultiWayOverlapper = function(win.start,win.end,gene.start,gene.end,gene.list) {
  #this is a monster, but basically, looks within each row for genes  overlapping with whatever you want

  blah=outer(as.numeric(unlist(win.start)), as.numeric(unlist(gene.end)), "<=") 
  blah1=outer(as.numeric(unlist(win.end)), as.numeric(unlist(gene.start)), ">=") 
  blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	if(!is.null(nrow(blah))){return(blah)}	
  

  }


#load df -----------
ap.fas = readLines('/home/blencowe/blencowe31/harpurbr/orthodb/apis.fas')
ap.nms = ap.fas[1]
ap.nms = gsub(".*_","",ap.nms)
ap.fas = ap.fas[-1]
ap.fas = unlist(strsplit(as.character(ap.fas),"")) #extract Apis sequence 
ap.fas.n = ap.fas[which(ap.fas != "-")]
cons = read.table('alignment_non_ap.scores',skip=1)
cons$seq = ap.fas
cons = cons[which(cons$seq!="-"),]
cons$ap.pos = c(1:nrow(cons))
cons = cons[which(cons$V2 > 0.001),]


#get CDS positions ------------
cons$cds = (3*cons$ap.pos)-2


#load GFF -----------
system('grep GB53150 /home/blencowe/blencowe31/harpurbr/genomes/amel/ensemb/amelgtf_chr.gtf > temp')
x = read.table('temp',header=F,sep='\t')
x = x[which(x$V3=="CDS"),]
x = x[c(1,4,5,7)]


#for pos in the peptide ----
	#each is in correct position

f_evens = x$V5
f_odds = x$V4

f_odds[1] = 1
f_evens[1] = (x$V5[1]-x$V4[1]) +1

for(i in 2:length(f_odds)){
	diff = (f_evens[i] - f_odds[i]) 
	f_odds[i] = f_evens[i-1] +1
	f_evens[i] = f_odds[i] + diff

}

x$st = f_odds
x$en = f_evens

test = MultiWayOverlapper(x$st, x$en, cons$cds, cons$cds, cons)
cons = cbind(x[test[,1],], cons[test[,2],])
names(cons)[c(1)] = c("chr")

#ID high and low regions --------------
lw1 <- loess(V2 ~ V1,data=cons,span = 0.01)
#plot(V2 ~ V1, data=cons,pch=19,cex=0.1)
j <- order(cons$V1)
#lines(cons$V1[j],lw1$fitted[j],col="red",lwd=3)
#quantile(lw1$fitted[j],0.95)
#abline(h=quantile(lw1$fitted[j],0.95))

#quantile(lw1$fitted[j],0.05)
#abline(h=quantile(lw1$fitted[j],0.05))

cons$fit = lw1$fitted[j]
# summarize results --------------
cons$ex = paste(cons$chr, ":", cons$V4, "-", cons$V5, sep="")


cons.ag = aggregate(cons$fit, by = list(cons$ex), function(x) c(length(x), mean(x), max(x), min(x), 
	length(x[x<0.72]), length(x[x>0.87])))

cons.ag$GB = rep(ap.nms,nrow(cons.ag))


write.table(cons.ag, header=F,row.names=F)

















############# later ------------------


inc.tab = merge(inc.tab, psi.events, by = "EVENT")

boxplot(log10(inc.tab$N[inc.tab$COORD %in% cons$V1]),inc.tab$N[inc.tab$COORD %in% cons$V1[which(cons$V6>0)]] )


x = inc.tab[inc.tab$COORD %in% cons$V1,]
conser = rep('null', nrow(x))
conser[which(x$COORD %in% cons$V1[which(cons$V6>0)])] = "non_cons"
conser[which(x$COORD %in% cons$V1[which(cons$V7>0)])] = "high_cons"
x$conser = conser
x = x[c(3,4,7,15,13)]
x2 = x[which(x$COMPLEX %in% c("C2", "C3")),]
x2$COMPLEX = NULL
x2 = x2[which(x2$conser!="null"),]
test = melt(x2, measure.vars = c("N", "F"))

tiss.joy = ggplot(test, aes(x = (value), y = conser, fill = variable)) + 
				geom_joy(rel_min_height = 0.1, alpha = 0.8) +

				theme_classic() + 
				theme(#axis.line=element_blank(),
						axis.text.x=element_text(size = 16),
						axis.text.y=element_text(size = 16),
						#axis.ticks.x=element_blank(),
						axis.title.x=element_text(size = 16),
						axis.title.y=element_text(size = 16),
						legend.position="bottom",
						panel.background=element_blank(),
						panel.border=element_blank(),
						panel.grid.major=element_blank(),
						panel.grid.minor=element_blank(),
						plot.background=element_blank()) 


x = aov(test$value~test$conser*tests$variable)


              Df  Sum Sq  Mean Sq F value Pr(>F)    
test$conser    1 0.00643 0.006428   169.3 <2e-16 ***
Residuals   2374 0.09011 0.000038                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> tiss.joy
Picking joint bandwidth of 0.00115
> TUkeyHSD(x)
Error in TUkeyHSD(x) : could not find function "TUkeyHSD"
> TukeyHSD(x)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = test$value ~ test$conser)

$`test$conser`
                          diff         lwr         upr p adj
non_cons-high_cons 0.003693113 0.003136594 0.004249633     0


x = aov(test$value~test$conser*test$variable)
> summary(x)
                            Df  Sum Sq  Mean Sq F value   Pr(>F)    
test$conser                  1 0.00643 0.006428  171.65  < 2e-16 ***
test$variable                1 0.00046 0.000460   12.29 0.000463 ***
test$conser:test$variable    1 0.00083 0.000826   22.05  2.8e-06 ***
Residuals                 2372 0.08882 0.000037                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> TukeyHSD(x)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = test$value ~ test$conser * test$variable)

$`test$conser`
                          diff         lwr        upr p adj
non_cons-high_cons 0.003693113 0.003140347 0.00424588     0

$`test$variable`
            diff          lwr         upr     p adj
F-N 0.0008803923 0.0003880305 0.001372754 0.0004626

$`test$conser:test$variable`
                                 diff           lwr           upr     p adj
non_cons:N-high_cons:N   0.0023692963  0.0013444376  0.0033941550 0.0000000
high_cons:F-high_cons:N  0.0001583102 -0.0005985981  0.0009152185 0.9498262
non_cons:F-high_cons:N   0.0051752407  0.0041503821  0.0062000994 0.0000000
high_cons:F-non_cons:N  -0.0022109861 -0.0032358448 -0.0011861274 0.0000002
non_cons:F-non_cons:N    0.0028059444  0.0015699184  0.0040419705 0.0000000
non_cons:F-high_cons:F   0.0050169306  0.0039920719  0.0060417892 0.0000000



