############################################
###Loading and Separating the FASTA Files###
############################################


#Load Libraries
library(ape)
library(seqinr)

#Creating the Separator File
dict<- c("PB2\\|","PB1\\|", "PA\\|","HA\\|","NP\\|","NA\\|","M1\\|","M2\\|","NS1\\|","NS2\\|")

#Load the file
ali <- read.FASTA("influenza.fasta")

#Remove the special characters
strnames <- names(ali)
spcChar <- c("'","\\(","\\)","\\[","\\]","!","\\.","-"," ")

for(i in 1:length(spcChar)){
  strnames <- gsub(pattern = spcChar[i], replacement = "", x = strnames)  
}
names(ali) <- strnames
rm("strnames","spcChar")

#Separate the Names and Sort
lab <- names(ali)
labsort <- sort(lab)

PB2 <- NULL
PB1 <- NULL
PA <- NULL
HA <- NULL
NP <- NULL
NAT <- NULL
M1 <- NULL
M2 <- NULL
NS1 <- NULL
NS2 <- NULL

for(i in 1:length(labsort)){
  if(grepl(dict[1],names(ali)[i])){
    PB2 <- c(PB2,i)
  }else if(grepl(dict[2],names(ali)[i])){
    PB1 <- c(PB1, i)
  }else if(grepl(dict[3],names(ali)[i])){
    PA <- c(PA,i)
  }else if(grepl(dict[4],names(ali)[i])){
    HA <- c(HA, i)
  }else if(grepl(dict[5],names(ali)[i])){
    NP <- c(NP,i)
  }else if(grepl(dict[6],names(ali)[i])){
    NAT <- c(NAT,i)
  }else if(grepl(dict[7],names(ali)[i])){
    M1 <- c(M1,i)
  }else if(grepl(dict[8],names(ali)[i])){
    M2 <- c(M2,i)
  }else if(grepl(dict[9],names(ali)[i])){
    NS1 <- c(NS1,i)
  }else if(grepl(dict[10],names(ali)[i])){
    NS2 <- c(NS2,i)
  }
}

#Put variables into same list
masterfilt <- list(PB2,PB1,PA,HA,NP,NAT,M1,M2,NS1,NS2)

#Creating Separate Fasta Files
for(i in 1:length(dict)){
  alitmp <- ali[masterfilt[[i]]]
  write.dna(x = alitmp,format = "fasta",file = paste0("Segment",i,".fasta"))
  cat(i, "out of", length(dict), "Done.\n")
}

#Only Use When Debugging and need to remove Clutter
rm(PB2,PB1,PA,HA,NP,NAT,M1,M2,NS1,NS2,lab,labsort,dict,alitmp)

################
####Aligning####
################

file.count <- 1:10

if((Sys.info()["sysname"] == "Windows") == TRUE){
  for (i in 1:length(file.count)){
    shell(paste("muscle -in Segment",i,".fasta -out Output.fasta", sep="")) #Renamed the MUSCLE program to simply muscle.
    file.rename("Output.fasta",paste("Segment",i,".fasta",sep=""))
  }
}else{
  for (i in 1:length(file.count)){
    system(paste("muscle -in Segment",i,".fasta -out Output.fasta", sep=""))
    file.rename("Output.fasta",paste("Segment",i,".fasta",sep=""))
  }
}
readline("Ensure that all the segments are still present (Only applies with incomplete genomes/provided data).\n Press enter when done.")

######################
####Concatenation#####
######################

library(phangorn)

#Reading and pasting the alignment files
ali<-read.aa("Segment1.fasta",format="fasta")
aliname <- names(ali)

#Removing Segment Names
for(j in 1:length(aliname)){
  aliname[j] <- gsub(pattern = "^.*?\\.",replacement = "",x = aliname[j])
}
names(ali) <- aliname


Seq <- as.character(ali)
Seq <- Seq[order(rownames(Seq)),]

#Column labels to be used
seg<- c("PB2","PB1", "PA", "HA","NP","NA","M1","M2","NS1","NS2")

#Adding column names
cnames<-NULL
for(x in 1:ncol(Seq)){
  cnames[x] <- paste0(seg[1],".",x)
}
colnames(Seq)<-cnames

#Automating the first portion for the rest
for(i in 2:length(file.count)){
  ali<-read.aa(paste("Segment",i,".fasta",sep = ""),format="fasta")
  aliname <- names(ali)
  
  #Removing Segment Names
  for(j in 1:length(aliname)){
    aliname[j] <- gsub(pattern = "^.*?\\.",replacement = "",x = aliname[j])
  }
  names(ali) <- aliname
  
  #Extraction
  ali <- as.character(ali)
  ali <- ali[order(rownames(ali)),]
  
  #Setting Column Names
  cnames<-NULL
  for(x in 1:ncol(ali)){
    cnames[x] <- paste(paste(seg[i],".",sep = ""),x, sep = "")
  }
  colnames(ali)<-cnames
  
  #Accounting for differences in Strain number
  tmp <- rownames(Seq) %in% rownames(ali)
  Seq <- Seq[tmp,]
  tmp <- rownames(ali) %in% rownames(Seq)
  ali <- ali[tmp,]
  
  #Checking if the same data is used for each segment
  if(sum(unique(rownames(Seq))!=unique(rownames(ali)))!=0){
    stop("Data Alignment Error at", i-1, "-", i)
  } else {
    cat("Rownames match",paste(i,sep=""),"\n")
  }
  
  #Binding the segments together
  Seq <- cbind(Seq,ali)
}

#####################
####DNA Filtering####
#####################
pol<-c()
for(x in 1:ncol(Seq)){
  if (length(unique(Seq[,x]))==1)
  {pol[x]<-FALSE}
  else 
  {pol[x]<-TRUE}  
}
pol<-colnames(Seq)[pol]

Seqpol<-Seq[,pol]

#############################
####Extracting the Traits####
#############################
seqnames<-rownames(Seqpol)

seqsplit<-strsplit(seqnames,split = "\\.")

#Extracting from Traits
Virulence<-NULL
Transmission<-NULL
Polybasic_Cleavage<-NULL

for(i in 1:length(seqsplit)){
  Virulence[i]<- seqsplit[[i]][5]
  Transmission[i]<- seqsplit[[i]][6]
  Polybasic_Cleavage[i]<- seqsplit[[i]][8]
}

#Factorization
Virulence<-as.factor(Virulence)
Transmission<-as.factor(Transmission)
Polybasic_Cleavage<-as.factor(Polybasic_Cleavage)

PhenotypeCriteria <- data.frame(seqnames,Virulence,Transmission,Polybasic_Cleavage)

write.csv(PhenotypeCriteria,"Phenotype Table.csv")

cat("Traits Extracted \n")

#############################
####The Adaptive Boosting####
#############################

#This was run on as separate Cluster.  Only Virulence at size 125 is shown.  To change, just swap out the trait used
#and the chunk size

library(adabag)
library(TeachingDemos)

SeqV<-cbind.data.frame(Virulence,Seq)

# Prep For Adaptive Boosting####
chunk_size <- 125 # smaller values [e.g., 1000] seem to give more potential features
ali_len <- dim(Seqpol)[2]
n_chunks <- floor(ali_len/chunk_size) + 1
feature_used <- NULL
mfinal_chunk <- 100

#Adaptive Boosting SeqV####
for(i in 1:n_chunks){
  # prep alignment chunk
  lo <- chunk_size * (i - 1) + 1
  hi <- (chunk_size * i)
  if(hi > ali_len){
    hi <- ali_len
  }
  pc <- format(100*i/n_chunks, digits=2, nsmall=2)
  print(paste0("i:",i,"; lo:",lo,"; hi:",hi, " (", pc, "%)"))
  short_ali <- Seqpol[, seq(lo,hi)]
  
  mychunkdf <- data.frame(SeqV[,1], short_ali)
  # to have actual site positions
  colnames(mychunkdf)[1] <- "Virulence"
  
   # run adaptive boosting
  repeat{
    chunck.adaboost <- NULL
    #chunck.adaboost <- tryCatch(boosting(Virulence ~., data = mychunkdf, boos = T, mfinal = mfinal_chunk, coeflearn = "Freund", control=rpart.control(minsplit=2)))
    try(chunck.adaboost <- boosting(Virulence ~., data = mychunkdf, boos = T, mfinal = mfinal_chunk, coeflearn = "Freund", control=rpart.control(minsplit=2)))
    if(length(chunck.adaboost) > 0){
      if(length(chunck.adaboost$trees) != mfinal_chunk){
        rm(chunck.adaboost)
        print("Deleted")
      }else if(length(chunck.adaboost$trees) == mfinal_chunk){
        break
      }
    }
  }
  chunk_importance_sorted <- sort(chunck.adaboost$importance, decreasing=T)
  
  feature_used <- c(feature_used, names(chunk_importance_sorted[chunk_importance_sorted > 1.5]))
}
feature_used


pos_of_interest <- which(colnames(Seqpol) %in% unique(feature_used))
ali_of_interest <- Seqpol[, pos_of_interest]
mypoidf <- data.frame(SeqV[,1], ali_of_interest)
colnames(mypoidf)[1] <- "Virulence"
poi.adaboost.V <- NULL

#poi.adaboost.V <- boosting(Virulence ~., data = mypoidf, mfinal = 50, coeflearn = "Freund", control=rpart.control(minsplit=2))
repeat{
  try(poi.adaboost.V <- boosting(Virulence ~., data = mypoidf, boos = T, mfinal = mfinal_chunk, coeflearn = "Freund", control=rpart.control(minsplit=2)))
  #cat("Test \n")
  if(length(poi.adaboost.V) > 0){
    if(length(poi.adaboost.V$trees) == mfinal_chunk){
      cat("Success \n")
      break
    }
  }
}

importance_sorted.V100 <- sort(poi.adaboost.V$importance, decreasing=T)
importance_sorted.V100

Elapsed <- proc.time() - Begin
save.image(file = "VBoost125.RData")

############################################################
####Plotting the importance values and the Venn Diagrams####
############################################################
library(TeachingDemos)
#Showing the code for a single chunk size in Pathogenicity, Transmission, and Virluence
#Plotting Virulence
pdf("SeqV.adaboost125.pdf", width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,2))
plot(importance_sorted.V125[importance_sorted.V125 > 1], type="h", xlab="Site ranking", ylab="AdaBoost importance", main="Virulence - Chunk Size 125")
text(seq(1:length(importance_sorted.V125[importance_sorted.V125 > 1]))+0.35, importance_sorted.V125[importance_sorted.V125 > 1], names(importance_sorted.V125[importance_sorted.V125 > 1]), cex=.5, srt=90, pos=2)
subplot( 
  plot(importance_sorted.V125, type="h", xlab="", ylab="", cex.axis=0.4), 
  x=grconvertX(c(0.5,1), from="npc"),
  y=grconvertY(c(.5,1), from="npc"),
  type="fig", pars=list( mar=c(1.,1.,0,0)+0.1) )
dev.off()

#Plotting Pathogenicity
pdf("SeqP.adaboost125.pdf", width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,2))
plot(importance_sorted.P125[importance_sorted.P125 > 1], type="h", xlab="Site ranking", ylab="AdaBoost importance", main="Polybasic Cleavage - Chunk Size 125")
text(seq(1:length(importance_sorted.P125[importance_sorted.P125 > 1]))+0.35, importance_sorted.P125[importance_sorted.P125 > 1], names(importance_sorted.P125[importance_sorted.P125 > 1]), cex=.5, srt=90, pos=2)
subplot( 
  plot(importance_sorted.P125, type="h", xlab="", ylab="", cex.axis=0.4), 
  x=grconvertX(c(0.5,1), from="npc"),
  y=grconvertY(c(.5,1), from="npc"),
  type="fig", pars=list( mar=c(1.,1.,0,0)+0.1) )
dev.off()

#Plotting Transmission
pdf("SeqT.adaboost125.pdf", width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,2))
plot(importance_sorted.T125[importance_sorted.T125 > 1], type="h", xlab="Site ranking", ylab="AdaBoost importance", main="Transmission - Chunk Size 125")
text(seq(1:length(importance_sorted.T125[importance_sorted.T125 > 1]))+0.35, importance_sorted.T125[importance_sorted.T125 > 1], names(importance_sorted.T125[importance_sorted.T125 > 1]), cex=.5, srt=90, pos=2)
subplot( 
  plot(importance_sorted.T125, type="h", xlab="", ylab="", cex.axis=0.4), 
  x=grconvertX(c(0.5,1), from="npc"),
  y=grconvertY(c(.5,1), from="npc"),
  type="fig", pars=list( mar=c(1.,1.,0,0)+0.1) )
dev.off()

#Making the Venn Diagrams
#NOTE: Requires several files at different chunk sizes.
library(gplots)

##Virulence
venndi <- vector("list" , 3)

#Chunk 75
load("VBoost75.RData")
venndi[[1]] <- names(head(importance_sorted.V75,25))
rm(list = setdiff(ls(), "venndi"))

#Chunk 125
load("VBoost125.RData")
venndi[[2]] <- names(head(importance_sorted.V125,25))
rm(list = setdiff(ls(), "venndi"))

#Chunk 175
load("VBoost175.RData")
venndi[[3]] <- names(head(importance_sorted.V175,25))
rm(list = setdiff(ls(), "venndi"))

#Naming the lists
names(venndi) <- c("Chunk 75", "Chunk 125", "Chunk175")

#Creating the plot file
VirVenn <- venn(venndi)

##Transmission
venndi <- vector("list" , 3)

#Chunk 75
load("TBoost75.RData")
venndi[[1]] <- names(head(importance_sorted.T75,25))
rm(list = setdiff(ls(), c("VirVenn","venndi")))

#Chunk 125
load("TBoost125.RData")
venndi[[2]] <- names(head(importance_sorted.T125,25))
rm(list = setdiff(ls(), c("VirVenn","venndi")))

#Chunk 175
load("TBoost175.RData")
venndi[[3]] <- names(head(importance_sorted.T175,25))
rm(list = setdiff(ls(), c("VirVenn","venndi")))

#Naming the lists
names(venndi) <- c("Chunk 75", "Chunk 125", "Chunk 175")

#Creating the plot file
TraVenn <- venn(venndi)

##PolyBasic Cleavage
venndi <- vector("list" , 3)

#Chunk 75
load("PBoost75.RData")
venndi[[1]] <- names(head(importance_sorted.P75,25))
rm(list = setdiff(ls(), c("TraVenn","VirVenn","venndi")))

#Chunk 125
load("PBoost130.RData")
venndi[[2]] <- names(head(importance_sorted.P130,25))
rm(list = setdiff(ls(), c("TraVenn","VirVenn","venndi")))

#Chunk 175
load("PBoost170.RData")
venndi[[3]] <- names(head(importance_sorted.P165,25))
rm(list = setdiff(ls(), c("TraVenn","VirVenn","venndi")))

#Naming the lists
names(venndi) <- c("Chunk 75", "Chunk 130", "Chunk 165")

#Creating the plot file
PolVenn <- venn(venndi)

#Saving to PDF
pdf("VennDiagrams.pdf", width=6, height=4)
  par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,2))
  #Virulence/Infectivity
  plot(VirVenn)
  title("Infectivity")
  plot(TraVenn)
  title("Transmission")
  plot(PolVenn)
  title("Pathogenicity")
dev.off()

#Individual Plots - Virulence
pdf("VirVenn.pdf", width=6, height=4)
  plot(VirVenn)
  title("Infectivity")
dev.off()

#Individual Plots - Transmission
pdf("TraVenn.pdf", width=6, height=4)
  plot(TraVenn)
  title("Transmission")
dev.off()

#Individual Plots - Transmission
pdf("PolVenn.pdf", width=6, height=4)
  plot(PolVenn)
  title("Pathogenicity")
dev.off()

######################################################
####Determining the effects of Chunking on runtime####
######################################################

#WD should be where all your Boosting results are located
files <- list.files()
files <- files[grep(x = files, pattern = "*.RData")]
files <- sort(files)
elap <- NULL
for(a in 1:length(files)){
  load(paste(files[a]))
  elap[a] <- Elapsed[3]
  cat(a,"out of", length(files),"\n")
}

#Splitting the Data
pol.y <- elap[grep("P",files)]
pol.y <- log10(pol.y)
pol.x <- files[grep("P",files)]
pol.x <- sapply(strsplit(pol.x,split = "t", fixed = TRUE), function(x) (x[2]))
pol.x <- sapply(strsplit(pol.x,split = ".", fixed = TRUE), function(x) (x[1]))
pol.x <- as.numeric(pol.x)


tra.y <- elap[grep("T",files)]
tra.y <- log10(tra.y)
tra.x <- files[grep("T",files)]
tra.x <- sapply(strsplit(tra.x,split = "t", fixed = TRUE), function(x) (x[2]))
tra.x <- sapply(strsplit(tra.x,split = ".", fixed = TRUE), function(x) (x[1]))
tra.x <- as.numeric(tra.x)

vir.y <- elap[grep("V",files)]
vir.y <- log10(vir.y)
vir.x <- files[grep("V",files)]
vir.x <- sapply(strsplit(vir.x,split = "t", fixed = TRUE), function(x) (x[2]))
vir.x <- sapply(strsplit(vir.x,split = ".", fixed = TRUE), function(x) (x[1]))
vir.x <- as.numeric(vir.x)

#Adapted from http://r-eco-evo.blogspot.ca/2011/08/comparing-two-regression-slopes-by.html
# pol
typ <- rep("pol", length(pol.x))
pol.df <- data.frame(typ, pol.x, pol.y)
colnames(pol.df) <- c("typ", "x", "y")
# tra
typ <- rep("tra", length(tra.x))
tra.df <- data.frame(typ, tra.x, tra.y)
colnames(tra.df) <- c("typ", "x", "y")
# vir
typ <- rep("vir", length(vir.x))
vir.df <- data.frame(typ, vir.x, vir.y)
colnames(vir.df) <- c("typ", "x", "y")

dat.df <- rbind(pol.df, tra.df, vir.df)

# 1. testing the slopes (ancova)
mod1 <- aov(y ~ x * typ, data = dat.df)
summary(mod1)		# no significant interaction bwt runtime and typ
mod2 <- aov(y ~ x + typ, data = dat.df)
summary(mod2)		# typ has a significant effect on runtime (y)
anova(mod1,mod2)	# the most parsimonious model is mod2

# 2. individual regressions, typ by typ
pol <- dat.df[dat.df$typ == "pol", ]
tra <- dat.df[dat.df$typ == "tra", ]
vir <- dat.df[dat.df$typ == "vir", ]

pol.lm <- lm(y ~ x, data = pol)
summary(pol.lm)
tra.lm <- lm(y ~ x, data = tra)
summary(tra.lm)
vir.lm <- lm(y ~ x, data = vir)
summary(vir.lm)

pdf("fig_runtime.pdf", width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1)) 
plot(pol$x, pol$y, pch = 1, col = "orange", xlim=c(70,155), ylim=c(min(c(pol$y, vir$y, tra$y), na.rm=T), max(c(pol$y, vir$y, tra$y), na.rm=T)), xlab = "Chunk Size", ylab = "log10(seconds)")
points(tra$x, tra$y, pch = 2, col = "blue")
points(vir$x, vir$y, pch = 3, col = "red")
abline(pol.lm, col="orange")
abline(tra.lm, col="blue")
abline(vir.lm, col="red")
legend("bottomright", c("pol", "tra", "vir"), pch=c(1,2,3), col=c("orange", "blue", "red"), bty="n",cex=.8)
dev.off()

# fit separate intercepts
nest.lm <- lm(y ~ typ/x - 1, data = dat.df)
summary(nest.lm) # intercepts: pol >> tra > vir