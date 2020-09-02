### Calculate matrix using Plucinski's Unsupervized Naive Bayes classifier

# calculate allele frequencies

# calculate likelihood of no relation and likelihood of relation

library(parallel)

calculate_loglikelihood = function(v1,v2,p1,p2,ploid){
	# p1 is vector of allele frequencies in sample 1
	# p2 is vector of allele frequencies in sample 2
	n1 = length(p1)
	n2 = length(p2)
	loglikelihood0 = sum(log(p1))+sum(log(p2))
	loglikelihood1 = log(sum(sum(sapply(1:n1, function (i) sapply(1:n2, 
				function (j) 1/n1/n2 * (v1[i] == v2[j])*exp((sum(log(p1[-i]))+sum(log(p2)))))))))
	if (length(v1) > 1 & length(v2) > 1) {
	pairs1 = combn(v1,2,simplify = FALSE)
	pairs2 = combn(v2,2,simplify = FALSE)
	npairs1 = length(pairs1)
	npairs2 = length(pairs2)

	loglikelihood2 = log(sum(sum(sapply(1:npairs1, function (i) sapply(1:npairs2, 
				function (j) 1/npairs1/npairs2 * (sum(sort(pairs1[[i]]) == sort(pairs2[[j]]))==2)*
							exp((sum(log(p1[-match(pairs1[[i]],v1)]))+sum(log(p2)))))))))
	} else { loglikelihood2 = NA}
	if (ploid == 1) {
		loglikelihood2 = loglikelihood1 
	}
	c(loglikelihood0,loglikelihood1,loglikelihood2)
}

##### modify epsilon value here
					    
epsilon = 0.3072  #rate of missed alleles, estimated by Joel May 2 2019

calculate_loglikelihood2 = function(v1,v2,p1,p2,ploid){
	# p1 is vector of allele frequencies in sample 1
	# p2 is vector of allele frequencies in sample 2
	n1 = length(p1)
	n2 = length(p2)
	#print(n1)
	#print(n2)
	loglikelihood0 = sum(log(p1))+sum(log(p2))
	#print(loglikelihood0)
	#print(p1)
	#print(p2)
	#print(n1)
	#print(v1)
	#print(v2)
	loglikelihood1 = log(max(sapply(1:n1, function (i) sapply(1:n2, 
				function (j)  (v1[i] == v2[j])*exp((sum(log(p1[-i]))+sum(log(p2)))))),na.rm=TRUE))
	#print(exp((sum(log(p1[-1]))+sum(log(p2)))))
	#print(sum(log(p1[-1])))
	#print(p1[-1])
	if (length(v1) > 1 & length(v2) > 1) {
	pairs1 = combn(v1,2,simplify = FALSE)
	pairs2 = combn(v2,2,simplify = FALSE)
	npairs1 = length(pairs1)
	npairs2 = length(pairs2)
  #print(p1)
  #print(p2)
	loglikelihood2 = log(max(sapply(1:npairs1, function (i) sapply(1:npairs2, 
				function (j) (sum(sort(pairs1[[i]]) == sort(pairs2[[j]]))==2)*
							exp((sum(log(p1[-match(pairs1[[i]],v1)]))+sum(log(p2)))))),na.rm=TRUE))
	} else { loglikelihood2 = NA}
  #print(loglikelihood1)
	#print(p1)
	#print(v1)
	#print(pairs1)
	#print(-match(pairs1[[1]],v1))
	#print(p1[-match(pairs1[[1]],v1)])
	#print(loglikelihood2)
	#if (length(p1)>2){print(p1[-match(pairs1[[1]],v1)])}
	#if (length(p1)>2){print(p1[-match(pairs1[[2]],v1)])}
	#if (length(p1)>2){print(p1[-match(pairs1[[3]],v1)])}
	#if (length(p1)>2){print(exp((sum(log(p1[-match(pairs1[[1]],v1)]))+sum(log(p2)))))}
	#if (length(p1)>2){print(exp((sum(log(p1[-match(pairs1[[2]],v1)]))+sum(log(p2)))))}
	#if (length(p1)>2){print(exp((sum(log(p1[-match(pairs1[[3]],v1)]))+sum(log(p2)))))}
	#print(sort(pairs1[[1]]))
	#print(sort(pairs2[[1]]))
	#print(sort(pairs1[[1]]) == sort(pairs2[[1]]))
	#print(sum(sort(pairs1[[1]]) == sort(pairs2[[1]])))
	#print(exp(sum(log(p1[-match(pairs1[[1]],v1)]))+sum(log(p2))))
	#print(p1[-match(pairs1[[2]],v1)])
	if (loglikelihood1 == -Inf) {
		loglikelihood1 = log(epsilon * min(c(p1,p2)))
	}
	if (loglikelihood2 == -Inf & ploid > 1) {
		nsharedalleles = max(sapply(1:npairs1, function (i) sapply(1:npairs2, 
						function (j) (sum(sort(pairs1[[i]]) == sort(pairs2[[j]]))))))
		loglikelihood2 = log((epsilon*min(c(p1,p2))) ^ (2-nsharedalleles))
	}
	if (ploid == 1) {
		loglikelihood2 = loglikelihood1 
	}
	c(loglikelihood0,loglikelihood1,loglikelihood2)
}

alleles = list()
frequencies = list()

for (j in 1:nloci) {
	locicolumns = grepl(paste(locinames[j],"",sep=""),colnames(data))
	#print(locinames[j])
	#print(colnames(data))
	#print(locicolumns)
	raw_alleles = c(as.matrix(data[,locicolumns]))
	write.csv(data,"/Users/adminuser/Desktop/CDC/Bayesian/pyamd/data1.csv", row.names = FALSE )
	#print(raw_alleles)
	#print(data)
	#nrow(data)
	#ncol(data)
	#print(data[,locicolumns])
	#print(raw_alleles)
	#write.csv(raw_alleles,"/Users/adminuser/Desktop/CDC/Bayesian/pyamd/raw_allelesn1.csv", row.names = FALSE )
	raw_alleles[raw_alleles == "NA"] = NA
	raw_alleles[raw_alleles == 0] = NA
	print(raw_alleles[raw_alleles == "NA"])
	#print(raw_alleles)
	write.csv(raw_alleles,"/Users/adminuser/Desktop/CDC/Bayesian/pyamd/raw_alleles2.csv", row.names = FALSE )
	alleles[[j]] = unique(raw_alleles[!is.na(raw_alleles)])
	#write.csv(alleles[[j]],paste0("/Users/adminuser/Desktop/CDC/Bayesian/pyamd/allelestest", (j), ".csv"), row.names = FALSE )
	#write.csv(alleles,"/Users/adminuser/Desktop/CDC/Bayesian/pyamd/alleles.csv", row.names = FALSE )
	frequencies[[j]] = sapply(alleles[[j]], function(x) sum(raw_alleles == x,na.rm=TRUE))
	#write.csv(frequencies[[j]],paste0("/Users/adminuser/Desktop/CDC/Bayesian/pyamd/frequencies", (j), ".csv"), row.names = FALSE )
	frequencies[[j]] = frequencies[[j]] / sum(frequencies[[j]])
	#write.csv(frequencies[[j]],paste0("/Users/adminuser/Desktop/CDC/Bayesian/pyamd/frequencies", (j), ".csv"), row.names = FALSE )
	#print(frequencies[[j]])

}


observeddatamatrix = list()
for (j in 1:nloci) {
	locus = locinames[j]
	locicolumns = grepl(paste(locus,"",sep=""),colnames(data))
	oldalleles = as.vector(data[,locicolumns])
	oldalleles [oldalleles == "NA"] = NA
	oldalleles [oldalleles == 0] = NA
	if (length(dim(oldalleles)[2]) == 0) {
		oldalleles = matrix(oldalleles,length(oldalleles),1)
	}
	#print(oldalleles)
	observeddatamatrix[[j]] = oldalleles 
}

pairwisedistance = function(isolate1,isolate2){
	#print(((isolate2-1)*nids+isolate1)/ (nids*nids))
	#print(isolate1)
	#print(isolate2)
	loglik = matrix(NA,nloci,3)
	#print(loglik)
	for (j in 1:nloci) {
		v1 = observeddatamatrix[[j]][isolate1,]
		#print(v1)
		v1 = v1[!is.na(v1)]
		p1 = frequencies[[j]][match(v1,alleles[[j]])]
		v2 = observeddatamatrix[[j]][isolate2,]
		v2 = v2[!is.na(v2)]
		p2 = frequencies[[j]][match(v2,alleles[[j]])]
		#print(isolate1)
		#print(isolate2)
		#print(v1)
		#print(v2)
		#print(p1)
		#print(p2)
		if (length(v1) > 0 & length(v2) > 0) {
			loglik[j,] = calculate_loglikelihood2(v1,v2,p1,p2,ploidy[j])
			print(loglik[j,])
		} else { loglik[j,] = c(NA,NA,NA)}
		#write.csv(loglik[j,],paste0("/Users/adminuser/Desktop/CDC/Bayesian/pyamd/loglik", (j), ".csv"), row.names = FALSE )
	}
	loglik_allloci = (colSums(matrix(loglik,ncol=3),na.rm=TRUE))
	#write.csv(loglik_allloci,("/Users/adminuser/Desktop/CDC/Bayesian/pyamd/loglik_allloci.csv"), row.names = FALSE )
	lik_allloci = exp(loglik_allloci) / sum(exp(loglik_allloci))
	#write.csv(lik_allloci,("/Users/adminuser/Desktop/CDC/Bayesian/pyamd/lik_allloci.csv"), row.names = FALSE )
	sum(lik_allloci*c(0,1,2))
}



allpossiblepairs = expand.grid(1:nids,1:nids)
#print(nids)
#expand the grid from imported data
#print(allpossiblepairs)
#write.csv(allpossiblepairs,"/Users/adminuser/Desktop/CDC/Bayesian/pyamd/allpossiblepairs.csv", row.names = FALSE )
allpossiblepairs = unique(allpossiblepairs[allpossiblepairs[,1] <= allpossiblepairs[,2],])
#print(allpossiblepairs)
#write.csv(allpossiblepairs,"/Users/adminuser/Desktop/CDC/Bayesian/pyamd/allpossiblepairs2.csv", row.names = FALSE )
# pairwisedistancevector = unlist(lapply(1:dim(allpossiblepairs)[1], function (x) pairwisedistance(allpossiblepairs[x,1],allpossiblepairs[x,2]))) # not parallel

				  
				  ###### MODIFY NUMBER OF CORES BELOW - mc.cores=###

#print(allpossiblepairs[11,1])
#print(allpossiblepairs[11,2])

#print(pairwisedistance(allpossiblepairs[1,1],allpossiblepairs[1,2]),mc.cores=12)
pairwisedistancevector = unlist(mclapply(1:dim(allpossiblepairs)[1], function (x) pairwisedistance(allpossiblepairs[x,1],allpossiblepairs[x,2]),mc.cores=12)) # parallel
#print(pairwisedistance(allpossiblepairs[1,1],allpossiblepairs[1,2]),mc.cores=12)

pairwisedistancematrix = matrix(NA,nids,nids)
sapply(1:dim(allpossiblepairs)[1], function (x) pairwisedistancematrix[allpossiblepairs[x,1],allpossiblepairs[x,2]] <<- pairwisedistancevector[x])
sapply(1:dim(allpossiblepairs)[1], function (x) pairwisedistancematrix[allpossiblepairs[x,2],allpossiblepairs[x,1]] <<- pairwisedistancevector[x])

colnames(pairwisedistancematrix) = ids 
rownames(pairwisedistancematrix) = ids
###write.csv(pairwisedistancematrix,"pairwisedistancematrix_Bayesian.csv")

Bayesian_pairwisedistancematrix = pairwisedistancematrix 


pairwisedistancematrix

colv_clustering = rep(rgb(0,0,0), length(ids))


print("Calculation of Bayesian matrix complete")
