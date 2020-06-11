analyseMR <- function(gene, dir, Ngwas=239087, N_eQTLs=32000){
	# This code is lightly adapted from https://github.com/eleporcu/TWMR
  
	out<-c("gene","alpha","SE","P","Nsnps","Ngene")
	
	# Load eqtl matrix
	file<-file.path(dir, paste(gene,"matrix",sep="."))
	print(paste("Loading", file, sep=" "))
	filecluster<-read.table(file,header=T,sep=" ",dec=".")
	beta<-as.matrix(filecluster[,2:(length(filecluster[1,])-1)])
	x<-colSums(abs(beta))
	remove<-which(x==0)
	if(length(remove)>0) {
	beta<-beta[,-remove]}
	beta<-as.matrix(beta)
	gamma<-as.matrix(filecluster[,length(filecluster[1,])])
	LDmatrix<-file.path(dir, paste(gene,"ld",sep="."))
	print(paste("Loading", LDmatrix, sep=" "))
	C<-read.table(LDmatrix,header=F,sep=" ",dec=".")
	C<-as.matrix(C[,1:length(C[,1])])
	S<-t(beta)%*%solve(C)%*%beta
	H<-(1-1/sqrt(3781))*S+(1/sqrt(3781))*diag(length(S[,1]))
	alpha<-solve(H)%*%(t(beta)%*%solve(C)%*%gamma)
	alpha<-as.vector(alpha)
	     C_inv <- solve(C)
	GCG_inv <- t(beta) %*% solve(C) %*% beta
	GCG_inv<-(1-1/sqrt(3781))*GCG_inv+(1/sqrt(3781))*diag(length(GCG_inv[,1]))
	GCG_inv<-solve(GCG_inv)
	     df_dg <- GCG_inv %*% t(beta) %*% C_inv
	     df_dG <- (GCG_inv %x% (t(gamma) %*% C_inv %*% ((beta %*% GCG_inv %*% t(beta)) %*% C_inv + diag(nrow(beta))))) + ((-t(gamma) %*% C_inv %*% beta %*% GCG_inv) %x% (GCG_inv %*% t(beta) %*% C_inv))
	     J <- cbind(df_dG, df_dg)
	 SEs<-c(rep(1/sqrt(32000),length(beta[1,])*length(beta[,1])),rep(1/sqrt(Ngwas),length(gamma[,1])))
	 R<-diag(length(beta[1,])+1)
	 Sigma <- (SEs %*% t(SEs)) * (C %x% R)
	     V <- J %*% Sigma %*% t(J)
	se<- sqrt(V[1,1])
	N=length(beta[,1])
	Ngene=length(beta[1,])
	Z<-alpha[1]/se
	pval<-2*pnorm(abs(Z),lower.tail=FALSE)
	line<-c(colnames(filecluster)[2],alpha[1],se,pval,N,Ngene)
	out<-rbind(out,line)
	out_file <- file.path(dir, paste(gene, "alpha", sep="."))
	print(paste("Writing output to", out_file, sep=" "))
	write.table(out,file=out_file,quote=F,col.names=F,row.names=F)
}
