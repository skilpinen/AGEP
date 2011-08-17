# AGEP library codes
# 
# Author: Sami Kilpinen & Kalle Ojala
###############################################################################

create.ref.obj<-function(data,annotation,min.data=5,min.bw=-Inf,max.bw=Inf,resolution=512,verbose=TRUE,tos=NULL,return.expr=TRUE){
	# Wrapper function to generate reference data object
	if (ncol(data)!=length(annotation)){stop("Annotation vector needs to be as long as there are columns in the datamatrix!")}
	bws<-calculate.bws(data=data,gr=annotation,min.data=min.data,min.bw=min.bw,max.bw=max.bw,verbose=verbose,tos=tos)
	ref.data.obj<-create.tissue.dist(data=data,gr=annotation,bws=bws,resolution=resolution,verbose=verbose,return.expr=return.expr)
	return(ref.data.obj)
}


do.match<-function(query.sample,query.annotation=NULL,query.id=NULL,ref.data.obj,verbose=TRUE,verbose.res=FALSE){
	# Wrapper function to match one sample against reference data (object)
	if (length(query.id)>1){stop("Function can handle only one query sample at the time!")}
	if (!is.vector(query.sample) | !is.numeric(query.sample)){stop("query.sample needs to be a numeric vector!")}
	match.res<-match.query(query=query.sample,gr.que=query.annotation,que.sample.id=query.id,ref.data.obj=ref.data.obj,verbose.res=verbose.res,verbose=verbose)
	return(match.res)
}

calculate.bws<-function(data,gr,min.data=5,min.bw=-Inf,max.bw=Inf,verbose=FALSE,tos=NULL){
	un.groups<-unique(gr)
	res<-vector("list",length=length(un.groups))
	names(res)<-un.groups
	
	if (is.null(tos)){
		if (verbose){cat("\r Finding maximum values of each variable in each class (can take awhile)\n");flush.console()}
		# Change NA values temporarily to -Inf 
		na.ii<-is.na(data)
		data[na.ii]<--Inf
		if (length(unique(gr))>1){
			maxes<-aggregate(t(data),by=list(gr),FUN=max,na.rm=T)
			# Return those to NA
			data[na.ii]<-NA
			row.names(maxes)<-maxes[,1]
			maxes<-as.matrix(maxes[,-1])
		} else {
			maxes<-as.matrix(apply(data,1,max,na.rm=T))
			row.names(maxes)<-row.names(data)
			# Return those to NA
			data[na.ii]<-NA
		}
	} else {
		maxes<-tos
	}
	
	# Calculate bws
	for (i in 1:length(un.groups)){
		if (verbose){
			cat("\r Calculating bandwidths for class ",i,"/",length(un.groups),"              ")
			flush.console()
		}
		ii<-gr %in% un.groups[i]
		
		res[[i]]<-apply(data[,ii],1,function(x){
					if (length(na.omit(x))>=min.data){
						tmp.bw<-bw.nrd(na.omit(x))
						min.bw<-min.bw
						max.bw<-max.bw
						
						if (tmp.bw<min.bw){return(min.bw)} 
						else if (tmp.bw>max.bw){return(max.bw)}
						else {return(tmp.bw)}
						
					} else {NA}
				})
	}
	
	# Calculate tos
	if (is.null(tos)){
		if (length(unique(gr))>1){
			tmp.frame<-t(data.frame(res))
			row.names(tmp.frame)<-names(res)
		} else {
			tmp.frame<-data.frame(res)
		}
		
		df<-as.matrix(tmp.frame)
		tmp<-(df*3)+maxes[row.names(df),]
		if (length(unique(gr))>1){
			tos<-apply(tmp,2,max,na.rm=T)
		} else {
			tos<-tmp[,1]
		}
	} else {
		tos<-maxes
	}
	if (verbose){
		cat("\n")
		flush.console()
	}
	return(list(bws=res,tos=tos))
}

create.tissue.dist<-function(data,gr,bws,resolution=512,verbose,return.expr=TRUE){	
	tissue.dist<-list()
	un.groups<-unique(gr)
	scaling.factors<-list()
	from<-rep(0,nrow(data))
	to<-bws$tos
	bws<-bws$bws
	
	for (i in 1:length(un.groups)){
		
		ii<-which(gr==un.groups[i])
		i4bw<-which(names(bws)==un.groups[i])
		bws.tmp<-bws[[i4bw]]
		scale.factor<-numeric()
		
		tissue.dist[[i]]<-sapply(1:nrow(data),function(y){
					if (verbose){
						cat("\r Calculating densities for class ",y,"   of ",i," classes   ")
						flush.console()
					}
					if (!is.na(bws.tmp[y]) & bws.tmp[y]>0 & sum(!is.na(data[y,ii]))>=5){
						tmp.dens<-density(data[y,ii],na.rm=TRUE,bw=bws.tmp[y],from=from[y],to=to[y],n=resolution)$y
					} else {
						tmp.dens<-rep(NA,resolution)
					}
					# Scale area to one
					x.fac<-(to[y])/resolution
					scale.factor[y]<<-1/sum(tmp.dens*x.fac)
					tmp.dens<-tmp.dens*scale.factor[y]
					return(tmp.dens)
				})
		names(scale.factor)<-row.names(data)
		scaling.factors[[i]]<-scale.factor
		colnames(tissue.dist[[i]])<-row.names(data)
	}
	names(scaling.factors)<-as.character(un.groups)
	names(tissue.dist)<-as.character(un.groups)
	names(from)<-names(to)
	if (return.expr){
		return(list(data.densities=tissue.dist,froms=from,tos=to,scaling.factors=scaling.factors,resolution=resolution,bws=bws,ref.data=data,annotation=gr))
	} else {
		return(list(data.densities=tissue.dist,froms=from,tos=to,scaling.factors=scaling.factors,resolution=resolution,bws=bws,ref.data=NA,annotation=gr))
	}
}

kd.fast<-function(nmbrs,x,bw){
	# nmbrs is ref.data
	if (length(x)!=nrow(nmbrs)){stop("Query vector not equal in length to ref. data matrix!")}
	na.i<-is.na(x)
	n<-rowSums(!is.na(nmbrs))				
	tmp1<-rowSums((1/sqrt(pi*2)) * ((exp(1)) ** -(0.5 * (((x-nmbrs) / bw) ** 2))),na.rm=TRUE)/(n*bw)
	tmp1[na.i]<-NA
	return(tmp1)
}

match.query<-function(query,gr.que,que.sample.id,ref.data.obj,verbose.res=TRUE,verbose=FALSE){
	
	gr.ref<-ref.data.obj$annotation
	un.groups<-unique(gr.ref)
	
	genes.in.ref.data<-names(ref.data.obj$froms)
	common.genes<-intersect(names(query),genes.in.ref.data)
	if (length(common.genes)<1){stop("Not common variables between query and reference data!")}
	
	res<-matrix(data=NA,nrow=length(genes.in.ref.data),ncol=length(un.groups),dimnames=list(genes.in.ref.data,un.groups))
	
	query<-query[common.genes]	
	froms<-ref.data.obj$froms[common.genes]
	tos<-ref.data.obj$tos[common.genes]
	query.data.i<-!is.na(query)
	res<-res[common.genes,]
	bw<-ref.data.obj$bws
	resolution<-ref.data.obj$resolution
	scaling.factors<-ref.data.obj$scaling.factors
	tissue.dist<-ref.data.obj$data.densities
	ref.data<-ref.data.obj$ref.data[common.genes,]
	
	# Tissue specific part
	for (i in 1:length(un.groups)){
		if(verbose){
			cat("\r Comparing to group ",i,"/",length(un.groups),"                 ")
			flush.console()
		}
		kd4query<-numeric(length=length(common.genes))
		ii<-gr.ref %in% un.groups[i]		
		i4bw<-which(names(bw)==un.groups[i])
		i4scfc<-which(names(scaling.factors)==un.groups[i])
		bws<-bw[[i4bw]]
		bws<-bws[common.genes]
		s.fac<-scaling.factors[[i4scfc]][common.genes]
		tmp.data<-ref.data[,ii]
		
		dens.tmp2<-tissue.dist[[un.groups[i]]][,common.genes]
		kd4query<-kd.fast(nmbrs=tmp.data,x=query,bw=bws)
		kd4query<-kd4query*s.fac
		kd4query<-rowSums(t(dens.tmp2)>=kd4query)/nrow(dens.tmp2)		
		kd4query[(is.na(bws) | bws<0) & !query.data.i]<-NA		
		res[,i]<-kd4query
	}
	
	res.computed<-apply(t(1-res)+0.25,2,function(x){
				tmp<-(1/x) %*% t(x)
				ii<-tmp>=1 & !is.na(tmp)
				tmp[ii]<-(1/tmp[ii])
				tmp<-(1-((tmp-.2)*1.25))
				tmp[ii]<-tmp[ii]*-1
				diag(tmp)<-NA
				rowMeans(tmp,na.rm=T)
			})
	
	row.names(res.computed)<-colnames(res)
	
	final.scores<-sort(apply(res.computed,1,mean,na.rm=T), decreasing = TRUE)
	
	if (verbose.res){		
		return(list(tm.scores=t(res),ts.scores=res.computed,query.id=que.sample.id,query.annotation=gr.que,tissue.scores=final.scores))
	}
	else {
		return(list(query.id=que.sample.id,query.annotation=gr.que,tissue.scores=final.scores))
	}
}

# Function for weighing the AGEP results with gene/class specific weights
weighByGenes <- function (match.object, hit.genes, weight = 1, exponent = 2) {
	
	# Result vector
	results <- c()
	
	# Loop across all tissue of the result
	for (i in 1:nrow(match.object$ts.scores)) {
		
		# Which genes are common to both hit.genes and the AGEP result object
		hit.index  <- which(names(hit.genes) == rownames(match.object$ts.scores)[i])
		
		# Where are these in the 
		hit.these  <- match(colnames(match.object$ts.scores), names(hit.genes[[hit.index]]))
		
		# Calculate multiplier, balanced
		multiplier <- (hit.genes[[hit.index]])[hit.these] ^ exponent / mean((hit.genes[[hit.index]])[hit.these] ^ exponent, na.rm = TRUE) * weight + (1 - weight)
		
		# Calculate results with gene multipliers applied
		results[i] <- mean((match.object$ts.scores[i,]) * multiplier, na.rm = TRUE)
	}
	
	# Give names to categories
	names(results) <- rownames(match.object$ts.scores)
	
	# Resturn results
	return(results)
}

hitGenes <- function (ref.obj, verbose = TRUE) {
	require(bigmemory)
	require(biganalytics)
	
	b <- big.matrix(length(ref.obj$data.densities),prod(dim(ref.obj$data.densities[[1]])))
	for (i in 1:length(ref.obj$data.densities)){
		b[i,]<-as.vector(ref.obj$data.densities[[i]])
	}
	
	hit.genes<-lapply(1:length(ref.obj$data.densities),function(i){
				cat("\r",i)
				b[i,]<-0				
				maxes<-colmax(b,na.rm=T)
				
				bi<-as.vector(ref.obj$data.densities[[i]])
				b[i,]<-bi
				
				tmp.mat<-bi-maxes
				
				dim(tmp.mat)<-dim(ref.obj$data.densities[[i]])
				
				tmp.mat[tmp.mat<0]<-0
				tmp.mat<-sapply(1:ncol(tmp.mat),function(y){tmp.mat[,y]*(ref.obj$tos[y]/512)})
				colnames(tmp.mat)<-names(ref.obj$tos)
				sort(colSums(tmp.mat,na.rm=T))
			})
	names(hit.genes)<-names(ref.obj$data.densities)
	return(hit.genes)
}


