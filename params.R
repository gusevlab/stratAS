library('VGAM')
arg = commandArgs(trailingOnly=T)
f.COUNTS = arg[1]
f.CNV = arg[2]
f.OUT = arg[3]

cnv = read.table( f.CNV , head=F , as.is=T)
vcf = read.table( f.COUNTS , as.is=T )

# put in phase
al.ref = vcf[,4]
al.alt = vcf[,5]
switch = (vcf[,3] == "1|0")
tmp = al.ref[switch]
al.ref[switch] = al.alt[switch]
al.alt[switch] = tmp

phi = rep(NA,nrow(cnv))
mu = rep(NA,nrow(cnv))
num = rep(NA,nrow(cnv))

for ( c  in 1:nrow(cnv) ) {
	overlap = vcf[,1] == cnv[c,2] & vcf[,2] >= cnv[c,3] & vcf[,2] <= cnv[c,4]
	if ( sum(overlap) > 100 ) {
		fit = vglm(cbind( al.ref[overlap] , al.alt[overlap] ) ~ 1, betabinomialff, trace = FALSE)
		cof = Coef(fit)
		phi[c] = 1/(1+sum(cof))
		mu[c] = cof[1] / sum(cof)
		num[c] = sum(overlap)
		# cat( unlist(cnv[c,]) , phi[c] , mu[c] , '\n' , sep='\t' , file=stderr() )
	}
}

write.table( cbind(cnv,format(cbind(phi,mu),digits=3),num) , quote=F, row.names=F, col.names=F , sep='\t' , file=paste(f.OUT,".cnvphi",sep=''))

fit = vglm(cbind( al.ref , al.alt ) ~ 1, betabinomialff, trace = FALSE)              
cof = Coef(fit)
phi = 1/(1+sum(cof))
mu = cof[1] / sum(cof)
cat( cnv[1,1] , phi , mu , '\n' , sep='\t' , file=paste(f.OUT,".allphi",sep='') )