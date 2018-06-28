library('VGAM')
library("optparse")

option_list = list(
	make_option("--inp_counts", action="store", default=NA, type='character',
              help="Path to file containing allelic counts for this individual [required]"),
	make_option("--inp_cnv", action="store", default=NA, type='character',
              help="Path to file containing CNV boundaries for this individual [optional]"),	
	make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
	make_option("--min_cov", action="store", default=5 , type='integer',
              help="Minimum number of REF reads and ALT reads to include this site. [default: %default]")	
)
opt = parse_args(OptionParser(option_list=option_list))
    
vcf = read.table( opt$inp_counts , as.is=T )
# filter heterozygous sites with minimum reads
vcf = vcf[ (vcf$HAP == "1|0" | vcf$HAP = "0|1") & vcf$REF.READS >= opt$min_cov & vcf$ALT.READS >= opt$min_cov ,]

# put in phase
al.ref = vcf$REF.READS
al.alt = vcf$ALT.READS
switch = (vcf$HAP == "1|0")
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
	for ( c in 1:nrow(cnv) ) {
		overlap = vcf$CHR == cnv$CHR[c] & vcf$POS >= cnv$P0[c] & vcf$POS <= cnv$P1[c]
		if ( sum(overlap) > 100 ) {
			fit = vglm(cbind( al.ref[overlap] , al.alt[overlap] ) ~ 1, betabinomialff, trace = FALSE)
			cof = Coef(fit)
			phi[c] = 1/(1+sum(cof))
			mu[c] = cof[1] / sum(cof)
			num[c] = sum(overlap)
		}
	}
	write.table( cbind(cnv[,c("CHR","P0","P1")],format(cbind(phi,mu),digits=3),num) , quote=F, row.names=F, col.names=F , sep='\t' , file=paste(opt$out,".local.params",sep=''))
}

# fit all counts
fit = vglm(cbind( al.ref , al.alt ) ~ 1, betabinomialff, trace = FALSE)              
cof = Coef(fit)
phi = 1/(1+sum(cof))
mu = cof[1] / sum(cof)
cat( phi , mu , '\n' , sep='\t' , file=paste(opt$out,".global.params",sep='') )
