library('VGAM')
library("optparse")

option_list = list(
  make_option("--inp_counts", action="store", default=NA, type='character',
              help="Path to file containing allelic counts for this individual [required]"),
  make_option("--inp_cnv", action="store", default=NA, type='character',
              help="Path to file containing CNV boundaries for this individual [optional]"),	
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]") 
)
opt = parse_args(OptionParser(option_list=option_list))
    
vcf = read.table( opt$inp_counts , as.is=T )

# put in phase
al.ref = vcf[,4]
al.alt = vcf[,5]
switch = (vcf[,3] == "1|0")
tmp = al.ref[switch]
al.ref[switch] = al.alt[switch]
al.alt[switch] = tmp

# parameters to be estimated
phi = rep(NA,nrow(cnv))
mu = rep(NA,nrow(cnv))
num = rep(NA,nrow(cnv))

if ( !is.na(opt$inp_cnv) ) {
	cnv = read.table( opt$inp_cnv , head=F , as.is=T)

	# read through CNVs 
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
	write.table( cbind(cnv,format(cbind(phi,mu),digits=3),num) , quote=F, row.names=F, col.names=F , sep='\t' , file=paste(opt$out,".local.params",sep=''))
}

# fit all counts
fit = vglm(cbind( al.ref , al.alt ) ~ 1, betabinomialff, trace = FALSE)              
cof = Coef(fit)
phi = 1/(1+sum(cof))
mu = cof[1] / sum(cof)
cat( cnv[1,1] , phi , mu , '\n' , sep='\t' , file=paste(opt$out,".global.params",sep='') )
