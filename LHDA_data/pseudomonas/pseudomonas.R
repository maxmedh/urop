					######################
					# data preprocessing #
					######################

# load images and files
load("wholeAln.RData")
pheno <- read.csv("Dettman_Phenotypes.csv")
cds_index <- read.csv("alnIndex.csv")


# rename aln positions according to cds_index (takes a few minutes)
colnames(wholeAln)
pos_names <- c()
for(i in 1:length(cds_index[,1])){
	idx <- seq(1, cds_index[i,3])
	pos_names <- c(pos_names, paste0(cds_index[i,2], "_", idx))
}
pos_names
length(pos_names) == length(wholeAln[1,]) # these two vectors should be the same length
colnames(wholeAln) <- pos_names


# eliminate monomorphic sites (takes about a couple of hours)
monositestoberemoved <- NULL
clean_aln <- wholeAln
for(i in 1:length(clean_aln[1,])){
	if(length(unique(clean_aln[,i])) == 1){
		if(!(i %% 10000)){
			print(i)
		}
		monositestoberemoved <- c(monositestoberemoved, i)
	}
}
monositestoberemoved

clean_aln <- clean_aln[,-monositestoberemoved]
length(wholeAln[1,])
length(clean_aln[1,])
save.image("mh_clean_data.RData")

# make sure rownames(clean_aln) and pheno$Isolate are in same order
short_seq_names <- rownames(clean_aln)
short_seq_names <- sub("JD","", short_seq_names)
short_seq_names <- sub("CMsoap.contigs.fa","", short_seq_names)
short_seq_names == as.character(pheno$Isolate) # should all be T -- which it is!




						#######################
						# adaboost on trait 1 #
						#######################

# load libraries
library(tree)
library(ape)
library(adabag)


trait_nb <- 1
trait_nb <- trait_nb + 1 # column number now

cur_trait <- pheno[, trait_nb]

# create data frame for ML analyses
hr_df <- data.frame(cur_trait, clean_aln)

# chuncking preparation
chunk_size <- 5000 
my_mfinal <- 100
ali_len <- dim(clean_aln)[2]
n_chunks <- floor(ali_len/chunk_size) + 1
feature_used_adaboost <- NULL


# adaptive boosting on chuncked data
for(i in 1:n_chunks){
	# prep alignment chunk
	lo <- chunk_size * (i - 1) + 1
	hi <- (chunk_size * i)
	if(hi > ali_len){
		hi <- ali_len
	}
	pc <- format(100*i/n_chunks, digits=2, nsmall=2)
	print(paste0("i:",i,"; lo:",lo,"; hi:",hi, " (", pc, "%)"))
	short_ali <- clean_aln[, seq(lo,hi)]
		
	mychunkdf <- data.frame(cur_trait, short_ali)
	
	# run adabag
	adabag_cur_trait <- boosting(cur_trait ~., data = mychunkdf, boos = T, mfinal = my_mfinal, coeflearn = "Freund", control=rpart.control(minsplit=2))
	chunk_importance_sorted <- sort(adabag_cur_trait$importance, decreasing=T)
	feature_used_adaboost <- c(feature_used_adaboost, names(chunk_importance_sorted[chunk_importance_sorted > 1]))
}
feature_used_adaboost

pos_of_interest_adaboost <- which(colnames(clean_aln) %in% unique(feature_used_adaboost))
ali_of_interest_adaboost <- clean_aln[, pos_of_interest_adaboost]
mypoidf_adaboost <- data.frame(cur_trait, ali_of_interest_adaboost)

adaboost_cur_trait_poi <- boosting(cur_trait ~., data = mypoidf_adaboost, boos = T, mfinal = my_mfinal, coeflearn = "Freund", control=rpart.control(minsplit=2))
trait_1_adaboost <- sort(adaboost_cur_trait_poi$importance, decreasing=T)

# save image
save.image("mh_trait01.RData")

# 10-fold cross-validation
RS.boostcv <- boosting.cv(cur_trait ~., data= mypoidf_adaboost, v=10, mfinal = 5000, coeflearn = "Freund", control=rpart.control(minsplit=2))
RS.boostcv[-1] 

succRate <- 100 - 100*RS.boostcv[-1]$error


# save image
save.image("mh_trait01cv.RData")


					#############
					# post-proc #
					#############


library(RColorBrewer)
library(gplots)
system("rm -fr figures/")
system("mkdir figures")

##################################
# GET TOP-10 SITES AND SITE ID'S #
##################################

mycol <- brewer.pal(3, "Set1")

# repeatability of estimated sites / their importance values by adaboost
top10_trait01 <- list()
all_trait01 <- list()

#trait01
for(ii in 1:4){
	print(paste0("Doing trait01 run0", ii))
	load(paste0("run0", ii, "/mh_trait01.RData"))
	top10_trait01[[ii]] <- trait_1_adaboost[1:10]
	all_trait01[[ii]] <- trait_1_adaboost
}

# save and cleanup a bit...
save(top10_trait01, all_trait01, 
					file="top10_trait_only.RData")
rm(list = ls())
load("top10_trait_only.RData")



#################################################
# importance correlation plots across all sites #
#################################################
mycol <- brewer.pal(3, "Set1")
#trait01
ref <- comp_mat <- as.matrix(all_trait01[[1]])
for(i in 2:4){
	print(paste0("Doing trait01 run0", i))
	compare <- c()
	for(j in 1:length(ref[,1])){
		cur_name <- pos <- NULL
		cur_name <- rownames(ref)[j]
		pos <- grep(cur_name, names(all_trait01[[i]]))
		if(length(pos) > 0){
			compare[j] <- all_trait01[[i]][pos]
		}else{
			compare[j] <- NA
		}
	}
	comp_mat <- cbind(comp_mat, compare)
}
pdf("figures/Fig1_ImportanceComp.pdf", width=12, height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,2))
# trait01
plot(comp_mat[,1], comp_mat[,2], pch=20, col=mycol[1], main="trait01", ylim=c(0, max(comp_mat[1,])), xlab="Importance (reference run)", ylab="Importance (comparative runs)")
points(comp_mat[,1], comp_mat[,3], pch=20, col=mycol[2])
points(comp_mat[,1], comp_mat[,4], pch=20, col=mycol[3])
abline(a = 0, b = 1, lty=2, col="gray")
legend("topleft", c("run A1 vs. A2", "run A1 vs. B1", "run A1 vs. B2"), pch=20, col=mycol, bty = "n")


dev.off()


##################################################
# venn diagrams across sites with importance > 1 #
##################################################

cur_sites_trait01 <- list()
site_list_trait01 <- c()

#trait01
for(i in 1:4){
	print(paste0("Doing trait01 run0", i))
	cur_sites <- NULL
	cur_sites_trait01[[i]] <- names(all_trait01[[i]][all_trait01[[i]] > 1])
}
pdf("figures/Fig2_Importance_1_Venn.pdf", width=12, height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1),mfrow=c(2,2))
venn(list(run_A1=cur_sites_trait01[[1]], run_A2=cur_sites_trait01[[2]], run_B1=cur_sites_trait01[[3]], run_B2=cur_sites_trait01[[4]]))
title("trait01")
site_list_trait01 <- Reduce(intersect, cur_sites_trait01)
site_list_trait <- paste(site_list_trait01, collapse =", ")
text(0,0, paste0("Common sites: ", site_list_trait), cex=.4, pos = 4)


dev.off()



####################
# fetch gene names #
####################

# get gene IDs
# from: http://www.pseudomonas.com/strain/download
if(!file.exists("UCBPP-PA14.csv")){
	download.file("http://www.pseudomonas.com/downloads/pseudomonas/pgd_r_16_2/Pseudomonas/complete/gtf-complete.tar.gz", destfile="gtf-complete.tar.gz")
	system("tar -zxvf gtf-complete.tar.gz")
	system("mv gtf/Pseudomonas_aeruginosa_UCBPP-PA14_109.gtf UCBPP-PA14.csv")
	system("rm -fr gtf/ gtf-complete.tar.gz")
}
PA14_annot <- read.table("UCBPP-PA14.csv")

gene_list_trait01 <- list()
gene_names_trait01 <- c()
gene_PAnames_trait01 <- c()

#trait01
for(i in 1:length(site_list_trait01)){
	tmp_lst1 <- tmp_lst2 <- pos <- NULL
	tmp_lst1 <- site_list_trait01[i]
	tmp_lst2 <- substr(tmp_lst1, 1, 10)
	pos <- grep(tmp_lst2, PA14_annot[,16])
	gene_names_trait01 <- c(gene_names_trait01, as.character(PA14_annot[pos, 19]))
	gene_PAnames_trait01 <- c(gene_PAnames_trait01, tmp_lst2)
}
all_top_gene_names <- gene_names_trait01
gene_names_trait01 <- unique(gene_names_trait01)
gene_PAnames_trait01 <- unique(gene_PAnames_trait01)
load("run01/mh_trait01.RData")
pdf("figures/Fig3_trait01_adaboost.pdf", width=12, height=8)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(2,2)) 
plot(trait_1_adaboost[trait_1_adaboost > 1], type="h", xlab="Site ranking", ylab="AdaBoost importance", main="trait01")
text(seq(1:length(trait_1_adaboost[trait_1_adaboost > 1])), trait_1_adaboost[trait_1_adaboost > 1], all_top_gene_names, cex=.5, srt=90, pos=2)
subplot( 
  plot(trait_1_adaboost, type="h", xlab="", ylab="", cex.axis=0.4), 
  x=grconvertX(c(0.5,1), from="npc"),
  y=grconvertY(c(.5,1), from="npc"),
  type="fig", pars=list( mar=c(1.,1.,0,0)+0.1) )


dev.off()

# report gene names
write.table(gene_names_trait01, file="gene_names_trait01.txt", quote=F, row.names=F, sep = "\t")



#########################################################################################
#########################################################################################
#########################################################################################

##########################
# DIRECT RESULTS FROM ML #
##########################

								###########
								# trait01 #
								###########

# load images and files
load("run01/mh_trait01.RData")
# ababoost
pdf("figures/trait01_adaboost_run01.pdf", width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1)) 
plot(trait_1_adaboost[trait_1_adaboost > 1], type="h", xlab="Site ranking", ylab="AdaBoost importance", main="trait01")
text(seq(1:length(trait_1_adaboost[trait_1_adaboost > 1])), trait_1_adaboost[trait_1_adaboost > 1], names(trait_1_adaboost[trait_1_adaboost > 1]), cex=.5, srt=90, pos=2)
subplot( 
  plot(trait_1_adaboost, type="h", xlab="", ylab="", cex.axis=0.4), 
  x=grconvertX(c(0.5,1), from="npc"),
  y=grconvertY(c(.5,1), from="npc"),
  type="fig", pars=list( mar=c(1.,1.,0,0)+0.1) )
dev.off()

# load images and files
load("run02/mh_trait01.RData")
# ababoost
pdf("figures/trait01_adaboost_run02.pdf", width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1)) 
plot(trait_1_adaboost[trait_1_adaboost > 1], type="h", xlab="Site ranking", ylab="AdaBoost importance", main="trait01")
text(seq(1:length(trait_1_adaboost[trait_1_adaboost > 1])), trait_1_adaboost[trait_1_adaboost > 1], names(trait_1_adaboost[trait_1_adaboost > 1]), cex=.5, srt=90, pos=2)
subplot( 
  plot(trait_1_adaboost, type="h", xlab="", ylab="", cex.axis=0.4), 
  x=grconvertX(c(0.5,1), from="npc"),
  y=grconvertY(c(.5,1), from="npc"),
  type="fig", pars=list( mar=c(1.,1.,0,0)+0.1) )
dev.off()

# load images and files
load("run03/mh_trait01.RData")
# ababoost
pdf("figures/trait01_adaboost_run03.pdf", width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1)) 
plot(trait_1_adaboost[trait_1_adaboost > 1], type="h", xlab="Site ranking", ylab="AdaBoost importance", main="trait01")
text(seq(1:length(trait_1_adaboost[trait_1_adaboost > 1])), trait_1_adaboost[trait_1_adaboost > 1], names(trait_1_adaboost[trait_1_adaboost > 1]), cex=.5, srt=90, pos=2)
subplot( 
  plot(trait_1_adaboost, type="h", xlab="", ylab="", cex.axis=0.4), 
  x=grconvertX(c(0.5,1), from="npc"),
  y=grconvertY(c(.5,1), from="npc"),
  type="fig", pars=list( mar=c(1.,1.,0,0)+0.1) )
dev.off()

# load images and files
# ababoost
pdf("figures/trait01_adaboost_run04.pdf", width=6, height=4)
par(oma = c(0,0,0,0), mar = c(4, 4, 1, 1), mfrow=c(1,1)) 
plot(trait_1_adaboost[trait_1_adaboost > 1], type="h", xlab="Site ranking", ylab="AdaBoost importance", main="trait01")
text(seq(1:length(trait_1_adaboost[trait_1_adaboost > 1])), trait_1_adaboost[trait_1_adaboost > 1], names(trait_1_adaboost[trait_1_adaboost > 1]), cex=.5, srt=90, pos=2)
subplot( 
  plot(trait_1_adaboost, type="h", xlab="", ylab="", cex.axis=0.4), 
  x=grconvertX(c(0.5,1), from="npc"),
  y=grconvertY(c(.5,1), from="npc"),
  type="fig", pars=list( mar=c(1.,1.,0,0)+0.1) )
dev.off()





##################################
save.image("top10_trait.RData")
q(save="no")
##################################

