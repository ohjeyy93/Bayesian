setwd("/Users/adminuser/Desktop/CDC/Bayesian/pyamd")

#### provide genotype datasheet
data = read.csv("Example_haplotype_data_sheet.txt",skip=0,stringsAsFactors = FALSE, sep = "\t")


## provide names of loci as per the genotype datasheet provided
locinames = c("CDC1","CDC2","CDC3","CDC4","X378_PART_A","X378_PART_B","X378_PART_C",
		    "X360i2_PART_A","X360i2_PART_B","X360i2_PART_C","X360i2.PART_D",
		    "Junction","MSR_Left","MSR_Right")
locinames_base = c("CDC1","CDC2","CDC3","CDC4","X378","X360","Junction","MSR")

## provide ploidy of each locus - ordered the same as the locinames variable
ploidy = c(2,2,2,2,2,2,2,2,2,2,2,1,1,1)

### is.na means not missing data
data = data[!is.na(data$Seq_ID) & data$Seq_ID != "",]
ids = data$Seq_ID
nids = length(ids)
nloci = length(locinames)

### transform data 

newdata = c()
datacompleteness = c()
# MaxMOI for each locus 
for (j in 1:nloci) {
	locicolumns = grepl(paste(locinames[j],"",sep=""),colnames(data))
	#print(colnames(data))
	#print(locicolumns)
	raw_alleles = data[,locicolumns]
	print(raw_alleles)
	#print(rowSums(raw_alleles == "X",na.rm=TRUE))
	maxMOI = max(c(2,rowSums(raw_alleles == "X",na.rm=TRUE)))
	#print(maxMOI)
	MOI = rowSums(raw_alleles == "X",na.rm=TRUE)
	nalleles = sum(locicolumns,na.rm=TRUE)
	#print(nalleles)
	newdatatemp = rbind(matrix(NA,nids,maxMOI))
	print(newdatatemp)
	sapply(1:nids,function (x) if (MOI[x] > 0) { newdatatemp[x,1:MOI[x]] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")})
	if (ploidy[j] > 1) {
		sapply(1:nids,function (x) if (MOI[x] == 1) { newdatatemp[x,1:2] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")})
	}
	colnames(newdatatemp) = paste(1:maxMOI,"_Hap_",locinames[j],sep="")
	#print(newdatatemp)
	newdata = cbind(newdata,newdatatemp )
	datacompleteness  = cbind(datacompleteness,MOI)
}

#print(ncol(newdata))
data = cbind(ids,newdata)
data = data.frame(data)

# calculate locus-specific amplification success
		       
datacompleteness_bylocus = sapply(locinames_base,function (x) rowSums(cbind(datacompleteness[,grepl(x,locinames)]))>0)

print("Number of loci sequenced")
table(rowSums(datacompleteness_bylocus))

#print(datacompleteness_bylocus)
#print(rowSums(datacompleteness_bylocus))


# define inclusion criteria below - filters genotype datasheet to exclude specimens that fail to meet these criteria

#print(data)
#print(datacompleteness_bylocus[,c(5,7,8)])
#print(rowSums(datacompleteness_bylocus))
#print(rowSums(datacompleteness_bylocus[,c(5,7,8)]) == 3)
#print(rowSums(datacompleteness_bylocus[,c(5,7,8)]))
#print(rowSums(datacompleteness_bylocus) >= 4)
#print(rowSums(datacompleteness_bylocus[,c(5,7,8)]) == 3 &rowSums(datacompleteness_bylocus) >= 4)
cleandata = data[(rowSums(datacompleteness_bylocus[,c(5,7,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | rowSums(datacompleteness_bylocus) >= 5 | (rowSums(datacompleteness_bylocus[,c(6,7,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | (rowSums(datacompleteness_bylocus[,c(5,6,7)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | (rowSums(datacompleteness_bylocus[,c(5,6,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4),] 
#print(cleandata)
print((rowSums(datacompleteness_bylocus[,c(5,7,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | rowSums(datacompleteness_bylocus) >= 5 | (rowSums(datacompleteness_bylocus[,c(6,7,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | (rowSums(datacompleteness_bylocus[,c(5,6,7)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | (rowSums(datacompleteness_bylocus[,c(5,6,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4))

#write.csv(cleandata,"/Users/adminuser/Desktop/CDC/Bayesian/pyamd/cleandata.csv", row.names = FALSE )
nrow(cleandata)
nrow(data)
ncol(cleandata)

data = cleandata
ids = data$ids
#print(ids)
nids = length(ids)

#print(nids)

####write.csv(data,"transformed_data.csv")
