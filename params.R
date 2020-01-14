library('VGAM')
library("optparse")

option_list = list(
	make_option("--id", action="store", default=NA, type='character',
              help="Sample identifier [required]"),
	make_option("--inp_counts", action="store", default=NA, type='character',
              help="Path to file containing allelic counts for this individual [required]"),
	make_option("--inp_cnv", action="store", default=NA, type='character',
              help="Path to file containing CNV boundaries for this individual [optional]"),	
	make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
	make_option("--min_cov", action="store", default=5 , type='integer',
              help="Minimum number of REF reads and ALT reads to include this site. [default: %default]"),
	make_option("--min_snps", action="store", default=100 , type='integer',
              help="Minimum number of SNPs in a local CNV region needed to estimate parameters. [default: %default]")              
)
opt = parse_args(OptionParser(option_list=option_list))
    
vcf = read.table( opt$inp_counts , as.is=T , header = TRUE)
# filter heterozygous sites with minimum reads
vcf = vcf[ (vcf$HAP == "1|0" | vcf$HAP == "0|1") & vcf$REF.READS + vcf$ALT.READS >= opt$min_cov ,]

# put in phase
al.ref = vcf$REF.READS
al.alt = vcf$ALT.READS
switch = (vcf$HAP == "1|0")
tmp = al.ref[switch]
al.ref[switch] = al.alt[switch]
al.alt[switch] = tmp

if ( !is.na(opt$inp_cnv) ) {
	cnv = read.table( opt$inp_cnv , header=TRUE , as.is=T)
	
	# local parameters to be estimated
	phi = rep(NA,nrow(cnv))
	mu = rep(NA,nrow(cnv))
	num = rep(NA,nrow(cnv))

	# read through CNVs 
	for ( c in 1:nrow(cnv) ) {
		overlap = vcf$CHR == cnv$CHR[c] & vcf$POS >= cnv$P0[c] & vcf$POS <= cnv$P1[c]
		if ( sum(overlap) > opt$min_snps ) {
			fit = vglm(cbind( al.ref[overlap] , al.alt[overlap] ) ~ 1, betabinomialff, trace = FALSE)
			cof = Coef(fit)
			phi[c] = 1/(1+sum(cof))
			mu[c] = cof[1] / sum(cof)
			num[c] = sum(overlap)
		}
	}
	if ( sum(colnames(cnv) == "CNV") == 1 ) {
	write.table( cbind(opt$id,cnv[,c("CHR","P0","P1","CNV")],format(cbind(phi,mu),digits=3),num) , quote=F, row.names=F, col.names=c("ID","CHR","P0","P1","CNV","PHI","MU","N"), sep='\t' , file=paste(opt$out,".local.params",sep=''))
	} else {
	write.table( cbind(opt$id,cnv[,c("CHR","P0","P1")],format(cbind(phi,mu),digits=3),num) , quote=F, row.names=F, col.names=c("ID","CHR","P0","P1","PHI","MU","N"), sep='\t' , file=paste(opt$out,".local.params",sep=''))
	}
}

# fit all counts
fit = vglm(cbind( al.ref , al.alt ) ~ 1, betabinomialff, trace = FALSE)              
cof = Coef(fit)
phi = 1/(1+sum(cof))
mu = cof[1] / sum(cof)
num = length( al.ref )
write.table( cbind(opt$id,phi , mu, num) , quote=F , row.names=F , sep='\t' , col.names=c("ID","PHI","MU","N") , file=paste(opt$out,".global.params",sep='') )

# --- check overdispersion by quantiles (with enough data)
#if ( opt$debug ) {
#	# split into deciles
#	decs = quantile(al.ref+al.alt, prob = seq(0, 1, length = 3), type = 5)
#	for ( i in 1:4 ) {
#		keep = al.ref+al.alt >= decs[i] & al.ref+al.alt < decs[i+1]
#		fit = vglm(cbind( al.ref[keep] , al.alt[keep] ) ~ 1, betabinomialff, trace = FALSE)
#		cof = Coef(fit)
#		phi = 1/(1+sum(cof))
#		mu = cof[1] / sum(cof)
#		num = length( al.ref )
#		cat( decs[i] , decs[i+1] , phi , mu , num , '\n' ) 
#	}
#}
