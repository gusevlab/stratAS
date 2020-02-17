library('VGAM')
library("optparse")

options(digits=3)

option_list = list(
	make_option("--id", action="store", default=NA, type='character',
              help="Sample identifier [required]"),
	make_option("--inp_counts", action="store", default=NA, type='character',
              help="Path to file containing allelic counts for this individual [required]"),
	make_option("--inp_cnv", action="store", default=NA, type='character',
              help="Path to file containing CNV boundaries for this individual [optional]"),	
	make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),           
	make_option("--group", action="store", default=NA , type='integer',
	            help="Group regions together by fold change (CNV column) and estimate jointly, specify how many groups [default: %default]"), 
	make_option("--group_snp", action="store_true", default=FALSE,
	            help="Group based on SNPs rather than CNVs (useful when there are many short CNVs) [default: %default]"), 
	make_option("--min_cov", action="store", default=5 , type='integer',
              help="Minimum number of REF reads and ALT reads to include this site. [default: %default]"),
	make_option("--min_snps", action="store", default=100 , type='integer',
              help="Minimum number of SNPs in a local CNV region needed to estimate parameters. [default: %default]"),
	make_option("--verbose", action="store_true", default=FALSE,
	            help="Verbose mode [default: %default]")                          
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

# fit all counts
fit = vglm(cbind( al.ref , al.alt ) ~ 1, betabinomialff, trace = FALSE)              
cof = Coef(fit)
phi = 1/(1+sum(cof))
mu = cof[1] / sum(cof)
num = length( al.ref )
write.table( cbind(opt$id,phi , mu, num) , quote=F , row.names=F , sep='\t' , col.names=c("ID","PHI","MU","N") , file=paste(opt$out,".global.params",sep='') )

if ( !is.na(opt$inp_cnv) ) {
	cnv = read.table( opt$inp_cnv , header=TRUE , as.is=T)
	
	# local parameters to be estimated
	phi = rep(NA,nrow(cnv))
	mu = rep(NA,nrow(cnv))
	num = rep(NA,nrow(cnv))

	# estimate grouped
	if ( !is.na(opt$group) ) {
		if ( opt$group < 2 ) {
			stop("ERROR: --group parameter must be at least 2\n")		
		}
		if ( sum(colnames(cnv) == "CNV") == 0 ) {
			stop("ERROR: --group requires a CNV column in --inp_cnv file on which to group\n")
		}

		if ( opt$group_snp ) {
			# assign every snp a CNV
			snp.cnv = rep(NA , nrow(vcf))
			for ( c in 1:nrow(cnv) ) {
				overlap = vcf$CHR == cnv$CHR[c] & vcf$POS >= cnv$P0[c] & vcf$POS <= cnv$P1[c]
				snp.cnv[ overlap ] = cnv$CNV[c]
			}	
			# split CNVs into deciles
			decs = quantile( snp.cnv[ !is.na(snp.cnv) ] , prob = seq(0, 1, length = (opt$group+1) ), type = 5)
		} else {			
			# split CNVs into deciles
			decs = quantile( cnv$CNV , prob = seq(0, 1, length = (opt$group+1) ), type = 5)
		}
		
		# for each decile
		for ( i in 1:(opt$group) ) {
			if( opt$verbose ) cat( "Grouped mode: Testing quantile ",i,": ",decs[i]," to ",decs[i+1],sep='',file=stderr() )
			
			# grab the CNVs
			keep = cnv$CNV >= decs[i] & cnv$CNV < decs[i+1]
			
			# merge the respective variants
			cur.al.ref = vector()
			cur.al.alt = vector()
			cur.tot = 0
			cnv.keep = rep( FALSE , nrow(cnv) )
			
			for ( c in which(keep) ) {
				overlap = vcf$CHR == cnv$CHR[c] & vcf$POS >= cnv$P0[c] & vcf$POS <= cnv$P1[c]
				cnv.keep[ c ] = TRUE
				cur.al.ref = c(cur.al.ref , al.ref[overlap])
				cur.al.alt = c(cur.al.alt , al.alt[overlap])
				cur.tot = cur.tot + sum(overlap)
			}
			
			if( opt$verbose ) cat( ", ",cur.tot," total SNPs",sep='',file=stderr() )
			
			if ( cur.tot > opt$min_snps ) {
				fit = vglm(cbind( cur.al.ref , cur.al.alt ) ~ 1, betabinomialff, trace = FALSE)
				cof = Coef(fit)
				cur.phi = 1/(1+sum(cof))
				cur.mu = cof[1] / sum(cof)
				
				phi[ cnv.keep ] = cur.phi
				mu[ cnv.keep ] = cur.mu
				num[ cnv.keep ] = sum(overlap)

				if( opt$verbose ) cat( ":\tPHI=", cur.phi , ",\tMU=", cur.mu , "\n",sep='',file=stderr() )
			} else {
				if( opt$verbose ) cat( "\n",sep='',file=stderr() )			
			}
		}
	} else {
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
	}
	
	if ( sum(colnames(cnv) == "CNV") == 1 ) {
	write.table( cbind(opt$id,cnv[,c("CHR","P0","P1","CNV")],format(cbind(phi,mu),digits=3),num) , quote=F, row.names=F, col.names=c("ID","CHR","P0","P1","CNV","PHI","MU","N"), sep='\t' , file=paste(opt$out,".local.params",sep=''))
	} else {
	write.table( cbind(opt$id,cnv[,c("CHR","P0","P1")],format(cbind(phi,mu),digits=3),num) , quote=F, row.names=F, col.names=c("ID","CHR","P0","P1","PHI","MU","N"), sep='\t' , file=paste(opt$out,".local.params",sep=''))
	}
}