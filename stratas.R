library(VGAM)
library(optparse)

option_list = list(
	make_option("--input", action="store", default=NA, type='character',
              help="Path to input file [required]"),
	make_option("--samples", action="store", default=NA, type='character',
              help="Path to sample identifier file, must have ID and CONDITION columns [required]"),	
	make_option("--peaks", action="store", default=NA, type='character',
              help="Path to file containing peak/gene boundaries, must contain CHR P0 P1 NAME CENTER columns [required]"),	
	make_option("--global_param", action="store", default=NA, type='character',
              help="Path to global parameter file [required]"),
	make_option("--local_param", action="store", default=NA, type='character',
              help="Path to local parameter file"),
	make_option("--window", action="store", default=100e3 , type='integer',
              help="Window (in bp) for SNPs to test around the peak boundary. [default: %default]"),
	make_option("--perm", action="store", default=0 , type='integer',
              help="# of permutations to shuffle the allele labels (0=off). [default: %default]"),
	make_option("--perm_cond", action="store", default=0 , type='integer',
	            help="# of permutations to shuffle the condition labels (0=off). [default: %default]"),	
	make_option("--min_cov", action="store", default=1 , type='integer',
              help="Individuals must have at least this many reads (for both alleles) to be tested. [default: %default]"),
	make_option("--min_maf", action="store", default=0.01 , type='double',
              help="Minimum minor allele frequency for test SNP. [default: %default]"),
	make_option("--min_het", action="store", default=0.01 , type='double',
              help="Minimum minor heterozygous frequency for test SNP. [default: %default]"),
	make_option("--max_rho", action="store", default=0.10 , type='double',
              help="Maximum local/global over-dispersion parameter for which to include individual in test. [default: %default]"),
	make_option("--binom", action="store_true", default=FALSE,
              help="Also perform a standard binomial test. [default: %default]"),
	make_option("--bbreg", action="store_true", default=FALSE,
	            help="Also perform a beta binomial regression, requires library(aod). [default: %default]"),	
	make_option("--indiv", action="store_true", default=FALSE,
	            help="Also report the per-individual allele fractions. [default: %default]"),	
	make_option("--exclude", action="store_true", default=75 , type='integer',
              help="The mimium distance between SNPs allowed in the haplotype. [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

PAR.WIN = opt$window

NUM.PERM = opt$perm
NUM.PERM_COND = opt$perm_cond

MIN.MAF = opt$min_maf
MIN.HET = opt$min_het
MIN.COV = opt$min_cov
MAX.RHO = opt$max_rho
DO.BINOM = opt$binom
DO.INDIV = opt$indiv
DO.BBREG = opt$bbreg

message("-----------------------------------")
message("stratAS")
message("https://github.com/gusevlab/stratAS")
message("-----------------------------------")

# --- error checks:
if ( DO.BBREG ) { 
  library("aod")
  if ( is.na(opt$local_param) ) {
    stop("ERROR: --bbreg requires --local_param file\n")
  }
}

if ( NUM.PERM > 0 && NUM.PERM_COND > 0 ) {
  stop("ERROR: --perm and --perm_cond cannot both be set\n")
}

if ( NUM.PERM < 0 || NUM.PERM_COND < 0 ) {
  stop("ERROR: --perm or --perm_cond cannot be negative\n")
}
# ---

message("Files being read in")
peaks = read.table( opt$peaks , head=T , as.is=T)
mat = read.table( opt$input , as.is=T )
phe = read.table( opt$samples , head=T , as.is=T)
cnv.all = read.table( opt$global_param , head=T ,as.is=T)
message("Files have been read in")

PERM.PVTHRESH = 0.05 / nrow(mat)

bb.loglike = function( mu , ref , alt, rho ) {
        keep = !is.na(ref + alt) & ref + alt > 0
        -1 * sum( dbetabinom(alt[keep], (ref+alt)[keep], mu, rho = rho[keep], log = T) )
}

bbinom.test = function( ref , alt , rho ) {
	if ( length(ref) > 0 && length(alt) > 0 && length(rho) > 0 ) {
		opt = optimize( bb.loglike , interval=c(0,1) , ref , alt, rho )
		opt$lrt = 2 * (opt$objective - bb.loglike( 0.5 , ref , alt , rho) )
		opt$pv = pchisq( abs(opt$lrt) , df=1 , lower.tail=F )
	} else {
		opt = list( "lrt" = NA , "pv" = NA , "min" = NA , "objective" = NA )
	}
	return( opt )
}

bbreg.test = function( ind , ref , alt , rho , cond , covar ) {
  if ( length(unique(cond)) > 1 && length(ref) > 0 && length(alt) > 0 && length(rho) > 0 ) {
    df = data.frame( y=alt , n=ref+alt , cond = cond , covar = covar , phi.group = as.factor( ind ) )
    n = length(rho)
    # test with a fixed overd parameter for each individual (for some reason this is very SLOW!)
    # reg = betabin( cbind( y , n - y  ) ~ 1 + cond , ~ phi.group , df , fixpar = list( 3:(n+2) , rho ) )
    
    # for efficiency, discretize the overd parameters into five groups
    nq = 5
    phi.q = rep(NA,length(ind))
    rho.q = rep(NA,nq)
    qq = quantile(unique(rho), probs = seq(0, 1, .2))
    for ( i in 1:nq ) {
      keep = rho >= qq[i] & rho <= qq[i+1]
      phi.q[ keep ] = i
      rho.q[ i ] = mean( rho[keep] )
    }
    df$phi.q = as.factor(phi.q)
    
    reg = betabin( cbind( y , n - y  ) ~ 1 + cond + covar , ~ phi.q , df , fixpar = list( 4:(nq+3) , rho.q ) )
    
    zscores = coef(reg) / sqrt(diag(vcov( reg )))
    pvals = 2*(pnorm( abs(zscores) , lower.tail=F))
    opt = list( "pv" = pvals )
  } else {
    opt = list( "pv" = NA )
  }
  return( opt )
}

cur.chr = unique(mat[,1])
m = match(peaks$CHR , cur.chr )
peaks = peaks[!is.na(m),]

if ( !is.na(opt$local_param) )  {
	cnv.local = read.table( opt$local_param , head=T ,as.is=T)
	m = match(cnv.local$CHR , cur.chr )
	cnv.local = cnv.local[!is.na(m),]
}

N = (ncol(mat) - 5)/4
M = nrow(mat)

HAPS = list()
GENO.H1 = matrix( 0 , nrow=M , ncol=N )
GENO.H2 = matrix( 0 , nrow=M , ncol=N )

PHENO = phe$CONDITION
m = match( phe$ID , cnv.all$ID )
cnv.all = cnv.all[m,]
RHO.ALL = cnv.all$PHI

for ( h in 1:2 ) {
	HAPS[[h]] = matrix( 0 , nrow=M , ncol=N )
}

# standardize matrix to the same haplotypes
for ( i in 1:N ) {
	GENO.H1[,i] = mat[ , 6 + 4*(i-1) ]
	GENO.H2[,i] = mat[ , 6 + 4*(i-1) + 1 ]
	HET = GENO.H1[,i] != GENO.H2[,i]
	cur.ALT = mat[ , 6 + 4*(i-1) ] == 0
	HAPS[[1]][HET & cur.ALT,i] = mat[ HET & cur.ALT , 6 + 4*(i-1) + 2 ]
	HAPS[[1]][HET & !cur.ALT,i] = mat[ HET & !cur.ALT , 6 + 4*(i-1) + 3 ]
	HAPS[[2]][HET & cur.ALT,i] = mat[ HET & cur.ALT , 6 + 4*(i-1) + 3 ]
	HAPS[[2]][HET & !cur.ALT,i] = mat[ HET & !cur.ALT , 6 + 4*(i-1) + 2 ]	
}

options( digits = 4 )

RHO = RHO.ALL
RHO[ RHO > MAX.RHO ] = NA

COL.HEADER = c("CHR","POS","RSID","P0","P1","NAME","CENTER","N.HET","N.READS","ALL.AF","ALL.BBINOM.P","C0.AF","C0.BBINOM.P","C1.AF","C1.BBINOM.P","DIFF.BBINOM.P")
COL.HEADER.BINOM = c("ALL.BINOM.P","ALL.C0.BINOM.P","ALL.C1.BINOM.P","FISHER.OR","FISHER.DIFF.P")
COL.HEADER.INDIV = c("IND.C0","IND.C0.COUNT.REF","IND.C0.COUNT.ALT","IND.C1","IND.C1.COUNT.REF","IND.C1.COUNT.ALT")
COL.HEADER.BBREG = c("ALL.BBREG.P","DIFF.BBREG.P","CNV.BBREG.P")

HEAD = COL.HEADER
if ( DO.BINOM ) HEAD = c(HEAD,COL.HEADER.BINOM)
if ( DO.BBREG ) HEAD = c(HEAD,COL.HEADER.BBREG)
if ( DO.INDIV ) HEAD = c(HEAD,COL.HEADER.INDIV)

cat( HEAD , sep='\t')
cat('\n')

for ( p in 1:nrow(peaks) ) {
	cur = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$P0[p] & mat[,2] <= peaks$P1[p]
	if(sum(cur) == 0 ) next()

	if ( !is.na(opt$local_param) )  {
		cur.cnv = cnv.local[ cnv.local$CHR == peaks$CHR[p] & cnv.local$P0 < peaks$P0[p] & cnv.local$P1 > peaks$P1[p] , ]
		m = match( phe$ID , cur.cnv$ID )
		cur.cnv = cur.cnv[m,]
		RHO = cur.cnv$PHI
		RHO[ RHO > MAX.RHO ] = NA
		COVAR = cur.cnv$CNV
	}

	# collapse reads at this peak
	cur.h1 = vector()
	cur.h2 = vector()
	cur.i = vector()
	
	for ( i in 1:N ) {
		reads.keep = GENO.H1[cur,i] != GENO.H2[cur,i] & HAPS[[1]][cur,i] >= MIN.COV & HAPS[[2]][cur,i] >= MIN.COV
		reads.keep[reads.keep] <- !c(FALSE, diff(mat[cur, 2][reads.keep]) < opt$exclude)
		cur.h1 = c( cur.h1 , (HAPS[[1]][cur,i])[reads.keep] )
		cur.h2 = c( cur.h2 , (HAPS[[2]][cur,i])[reads.keep] )
		cur.i = c( cur.i , rep( i , sum(reads.keep)) )
	}
	
	if ( length(unique(cur.i)) > MIN.MAF*N && sum(cur.h1) + sum(cur.h2) > 0 ) {
		# test all nearby SNPs

		if( PAR.WIN == -1 ) {
			cur.snp = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$P0[p] & mat[,2] <= peaks$P1[p]
		} else {
			cur.snp = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$CENTER[p] - PAR.WIN & mat[,2] <= peaks$CENTER[p] + PAR.WIN
		}

		for ( s in which(cur.snp) ) {
			# restrict to hets
			HET = GENO.H1[s,] != GENO.H2[s,] & !is.na(RHO)
			if ( sum(HET) > MIN.HET*N ) {
				# collect REF/ALT heterozygous haplotypes
				m1 = match( cur.i , which(HET & GENO.H1[s,] == 0) )
				m2 = match( cur.i , which(HET & GENO.H2[s,] == 0) )
				CUR.REF = c( cur.h1[ !is.na(m1) ] , cur.h2[ !is.na(m2) ])
				CUR.ALT = c( cur.h2[ !is.na(m1) ] , cur.h1[ !is.na(m2) ])
				CUR.IND = c( cur.i[ !is.na(m1) ] , cur.i[ !is.na(m2)] )

				if ( sum(CUR.REF) + sum(CUR.ALT) > 0 ) {
					for ( perm in 0:max(NUM.PERM,NUM.PERM_COND) ) {
					  
					  CUR.PHENO = PHENO
					  
						if ( perm > 0 ) {
							# randomly swap REF/ALT alleles
						  if ( NUM.PERM > 0 ) {
  							cur.swap = unique(CUR.IND)[ as.logical(rbinom( length(unique(CUR.IND)) , 1 , 0.5 )) ]
  							cur.swap = !is.na(match( CUR.IND , cur.swap ))
  							tmp = CUR.REF[ cur.swap ]
  							CUR.REF[ cur.swap ] = CUR.ALT[ cur.swap ]
  							CUR.ALT[ cur.swap ] = tmp
  						# shuffle the phenotype label
  						} else if (NUM.PERM_COND > 0 ) {
						    CUR.PHENO = sample(PHENO)
						  }
						}

						m = !is.na(match( CUR.IND , which(CUR.PHENO==0) ))
						CUR.REF.C0 = CUR.REF[m]
						CUR.ALT.C0 = CUR.ALT[m]
						CUR.IND.C0 = CUR.IND[m]

						m = !is.na(match( CUR.IND , which(CUR.PHENO==1) ))
						CUR.REF.C1 = CUR.REF[m]
						CUR.ALT.C1 = CUR.ALT[m]
						CUR.IND.C1 = CUR.IND[m]		
						
						# --- perform beta-binomial test
						tst.bbinom.C0 = bbinom.test( CUR.REF.C0 , CUR.ALT.C0 , RHO[CUR.IND.C0] )
						tst.bbinom.C1 = bbinom.test( CUR.REF.C1 , CUR.ALT.C1 , RHO[CUR.IND.C1] )
						tst.bbinom.ALL = bbinom.test( CUR.REF , CUR.ALT , RHO[CUR.IND] )
						lrt.BOTH = 2 * (tst.bbinom.C0$objective + tst.bbinom.C1$objective - tst.bbinom.ALL$objective)
						pv.BOTH = pchisq( abs(lrt.BOTH) , df=1 , lower.tail=F )			

						# --- print main output
						cat( unlist(mat[s,1:3]) , unlist(peaks[p,c("P0","P1","NAME","CENTER")]) , sum(HET) , sum(CUR.REF) + sum(CUR.ALT) , tst.bbinom.ALL$min , tst.bbinom.ALL$pv , tst.bbinom.C0$min , tst.bbinom.C0$pv , tst.bbinom.C1$min , tst.bbinom.C1$pv , pv.BOTH , sep='\t' )
						
						# --- perform binomial test
						if ( DO.BINOM ) {
							tst.binom = binom.test( sum(CUR.ALT),sum(CUR.REF)+sum(CUR.ALT) )
							if ( sum(CUR.REF.C0)+sum(CUR.ALT.C0) > 0 ) tst.binom.C0 = binom.test( sum(CUR.ALT.C0),sum(CUR.REF.C0)+sum(CUR.ALT.C0) ) else tst.binom.C0 = list("p.value"=NA)
							if ( sum(CUR.REF.C1)+sum(CUR.ALT.C1) > 0 ) tst.binom.C1 = binom.test( sum(CUR.ALT.C1),sum(CUR.REF.C1)+sum(CUR.ALT.C1) ) else tst.binom.C1 = list("p.value"=NA)
							# --- perform fisher's exact test between conditions					
							tst.fisher = fisher.test( cbind( c(sum(CUR.REF.C0),sum(CUR.ALT.C0)) , c(sum(CUR.REF.C1),sum(CUR.ALT.C1)) ) )

							cat( "" , tst.binom$p.value ,  tst.binom.C0$p.value , tst.binom.C1$p.value , tst.fisher$est , tst.fisher$p.value , sep='\t' )
						}
						
						# --- perform beta-binom regression
						if ( DO.BBREG ) {
						  tst.bbreg = bbreg.test( CUR.IND , CUR.REF , CUR.ALT , RHO[ CUR.IND ] , CUR.PHENO[CUR.IND] , COVAR[CUR.IND] )
						  cat( "" , tst.bbreg$pv , sep='\t' )
						}
						
						# --- print individual counts
						if ( DO.INDIV ) {
						  cat( "" , paste(CUR.IND.C0,collapse=',') , paste(CUR.REF.C0,collapse=',') , paste(CUR.ALT.C0,collapse=',') , paste(CUR.IND.C1,collapse=',') , paste(CUR.REF.C1,collapse=',') , paste(CUR.ALT.C1,collapse=',') , sep='\t' )
						}
						
						cat('\n')
					}
				}
			}
		}
	}
}