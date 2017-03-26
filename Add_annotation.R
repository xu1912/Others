fn="chr18_WGBS_position.txt"
pos=read.table(fn,header=T,sep="\t")
fn_n=strsplit(fn,"\\.")[[1]][1]
OutputF=paste(fn_n,".annot",sep="")

##Add promoter region annotation.
promoter=read.table("chr18_promoter.bed",header=F)
pos$promoter="N"
sp=1
for(i in 1:nrow(pos)){

	if (pos$start[i]< promoter$V2[sp]){
		next
	}
	
	for(j in sp:nrow(promoter)){
	
		if (pos$start[i]> promoter$V3[j]){
			next
		}
		
		if (pos$start[i] <= promoter$V3[j] && pos$start[i]>= promoter$V2[j]){
			pos$promoter[i]="Y"
			sp=j
			break
		}
		
		if (pos$start[i]< promoter$V2[j]){
			sp=j
			break
		}
	
	}

}


##Add gene region annotation.
pos$gene="N"
gd=read.table("glist-hg19_chr18",header=F)
genel <- gd[order(gd$V2,gd$V3),]
sp=1
for(i in 1:nrow(pos)){

	if (pos$start[i]< genel$V2[sp]){
		next
	}
	
	for(j in sp:nrow(genel)){
	
		if (pos$start[i]> genel$V3[j]){
			next
		}
		
		if (pos$start[i] <= genel$V3[j] && pos$start[i]>= genel$V2[j]){
			pos$gene[i]="Y"
			sp=j
			break
		}
		
		if (pos$start[i]< genel$V2[j]){
			sp=j
			break
		}
	
	}

}


##Add GERP region annotation.
pos$gerp="N"
gerp=read.table("hg19_GERP/hg19_chr18_elems.txt",header=F)
gerp <- gerp[order(gerp$V1,gerp$V2),]
sp=1
for(i in 1:nrow(pos)){

	if (pos$start[i]< gerp$V1[sp]){
		next
	}
	
	for(j in sp:nrow(gerp)){
	
		if (pos$start[i]> gerp$V2[j]){
			next
		}
		
		if (pos$start[i] <= gerp$V2[j] && pos$start[i]>= gerp$V1[j]){
			pos$gerp[i]="Y"
			sp=j
			break
		}
		
		if (pos$start[i]< gerp$V1[j]){
			sp=j
			break
		}
	
	}

}


##Add GC content in percentage annotation.
gccontent=read.table("hg19.gc5Base.chr18.txt",header=F)
sp=1
pos$gc_perct=0
for(i in 1585:nrow(pos)){

	wn=floor((pos$start[i]-1)/5)*5+1
	ln=which(gccontent[,1]==wn)
	if(length(ln)>0){
		pos$gc_perct[i]=gccontent[ln,2]
	}

}


##Add CGI annotation.
pos$cpgisland="N"
cpgisland=read.table("cpgIsland_chr18.pos",header=F)
cpgisland <- cpgisland[order(cpgisland$V1,cpgisland$V2),]
sp=1
for(i in 1:nrow(pos)){

	if (pos$start[i]< cpgisland$V1[sp]){
		next
	}
	
	for(j in sp:nrow(cpgisland)){
	
		if (pos$start[i]> cpgisland$V2[j]){
			next
		}
		
		if (pos$start[i] <= cpgisland$V2[j] && pos$start[i]>= cpgisland$V1[j]){
			pos$cpgisland[i]="Y"
			sp=j
			break
		}
		
		if (pos$start[i]< cpgisland$V1[j]){
			sp=j
			break
		}
	
	}

}


##Generate CGI shore (+-2kb out of island) and CGI shelf (+-2kb out of shore)
cgi_shore=matrix(0,nrow=nrow(cpgisland)*2,ncol=2)
cgi_shelf=matrix(0,nrow=nrow(cpgisland)*2,ncol=2)

for(i in 1:nrow(cpgisland)){

	cgi_shore[(2*i-1),2]=cpgisland[i,1]-1
	cgi_shore[(2*i-1),1]=cgi_shore[(2*i-1),2]-2000
	cgi_shore[(2*i),1]=cpgisland[i,2]+1
	cgi_shore[(2*i),2]=cgi_shore[(2*i),1]+2000
	
	cgi_shelf[(2*i-1),2]=cgi_shore[(2*i-1),1]-1
	cgi_shelf[(2*i-1),1]=cgi_shelf[(2*i-1),2]-2000
	cgi_shelf[(2*i),1]=cgi_shore[(2*i),2]+1
	cgi_shelf[(2*i),2]=cgi_shelf[(2*i),1]+2000

}


##Add CGI shore annotation.
pos$cpg_shore="N"
cgi_shore <- cgi_shore[order(cgi_shore[,1],cgi_shore[,2]),]
sp=1
for(i in 1:nrow(pos)){

	if (pos$start[i]< cgi_shore[sp,1]){
		next
	}
	
	for(j in sp:nrow(cgi_shore)){
	
		if (pos$start[i]> cgi_shore[j,2]){
			next
		}
		
		if (pos$start[i] <= cgi_shore[j,2] && pos$start[i]>= cgi_shore[j,1]){
			pos$cpg_shore[i]="Y"
			sp=j
			break
		}
		
		if (pos$start[i]< cgi_shore[j,1]){
			sp=j
			break
		}
	
	}

}


##Add CGI shelf annotation.
pos$cpg_shelf="N"
cgi_shelf <- cgi_shelf[order(cgi_shelf[,1],cgi_shelf[,2]),]
sp=1
for(i in 1:nrow(pos)){

	if (pos$start[i]< cgi_shelf[sp,1]){
		next
	}
	
	for(j in sp:nrow(cgi_shelf)){
	
		if (pos$start[i]> cgi_shelf[j,2]){
			next
		}
		
		if (pos$start[i] <= cgi_shelf[j,2] && pos$start[i]>= cgi_shelf[j,1]){
			pos$cpg_shelf[i]="Y"
			sp=j
			break
		}
		
		if (pos$start[i]< cgi_shelf[j,1]){
			sp=j
			break
		}
	
	}

}


write.table(pos,OutputF,col.names=TRUE,row.names=FALSE,quote=FALSE,sep=" ")
