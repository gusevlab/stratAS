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
	make_option("--keep", action="store", default=NA, type='character',
              help="Set of individuals to keep for prediction"),
	make_option("--predict_snps", action="store", default=NA, type='character',
              help="Set of SNPs to use for prediction"),
	make_option("--update_snps", action="store", default=NA, type='character',
              help="Update the SNP identifiers (col1 = chr:pos, col2=updated rsid)"),
	make_option("--predict", action="store_true", default=FALSE,
              help="Build TWAS/predictive models"),
	make_option("--predict_only", action="store_true", default=FALSE,
              help="Skip all testing steps"),
	make_option("--total_matrix", action="store", default=NA, type='character',
              help="Path to matrix of total activity, enables the linear model. Note: if --local_param is on, then CNV is included as covariate."),
	make_option("--total_rn", action="store_true", default=FALSE,
              help="Rank normalize the total expression phenotype."),
	make_option("--covar", action="store", default=NA, type='character',
              help="Path to covariates for total activity"),
	make_option("--combine", action="store", default="store_true", type=FALSE,
              help="Output combined BBINOM and TOTAL statistics by Stouffer's method, must have all relevant flags for basic and total input"),
	make_option("--window", action="store", default=100e3 , type='integer',
              help="Window (in bp) for SNPs to test around the peak CENTER value. Set to -1 to only test SNPs inside the peak boundary. [default: %default]"),
	make_option("--perm", action="store", default=0 , type='integer',
              help="# of permutations to shuffle the allele labels (0=off). [default: %default]"),
	make_option("--seed", action="store", default=NA , type='integer',
              help="Random seed. [default: %default]"),
    make_option("--perm_cond", action="store", default=0 , type='integer',
	            help="# of permutations to shuffle the condition labels (0=off). [default: %default]"),	
	make_option("--min_cov", action="store", default=1 , type='integer',
              help="Individuals must have at least this many reads (for both alleles) to be tested. [default: %default]"),
	make_option("--min_maf", action="store", default=0.01 , type='double',
              help="Minimum minor allele frequency for test SNP. [default: %default]"),
	make_option("--min_het", action="store", default=0.01 , type='double',
              help="Minimum minor heterozygous frequency for test SNP. [default: %default]"),
	make_option("--min_n_pred", action="store", default=10 , type='integer',
              help="Minimum heterozygous count for prediction. [default: %default]"),
	make_option("--max_rho", action="store", default=0.10 , type='double',
              help="Maximum local/global over-dispersion parameter for which to include individual in test. [default: %default]"),
	make_option("--binom", action="store_true", default=FALSE,
              help="Also perform a standard binomial test. [default: %default]"),
	make_option("--mbased", action="store_true", default=FALSE,
              help="Also perform the MBASED test for differences (requires MBASED libraries). [default: %default]"),
    make_option("--bbreg", action="store_true", default=FALSE,
	            help="Also perform a beta binomial regression, requires library(aod). [default: %default]"),	
	make_option("--fill_cnv", action="store_true", default=FALSE,
	            help="Set individuals with missing CNV calls to diploid and \rho=0.01 [default: %default]"),	
	make_option("--indiv", action="store_true", default=FALSE,
	            help="Also report the per-individual allele fractions. [default: %default]"),
	make_option("--collapse_reads", action="store_true", default=FALSE,
	            help="Merge all sites for each individual prior to analysis. [default: %default]"),
	make_option("--sim", action="store_true", default=FALSE,
	            help="Simulate imbalance and test. Allelic-fraction specified by --sim_af. [default: %default]"),	            
	make_option("--sim_cnv", action="store_true", default=FALSE,
	            help="Add local CNVs to simulations. [default: %default]"),
	make_option("--sim_cnv_allelic", action="store_true", default=FALSE,
	            help="CNVs always impact one allele. [default: %default]"),
	make_option("--sim_af1", action="store", default=0.50 , type='double',
              help="Specify the allelic fraction for CONDITION==1. [default: %default]"),              
	make_option("--sim_af0", action="store", default=0.50 , type='double',
              help="Specify the allelic fraction for CONDITION==0. [default: %default]"), 
	make_option("--cond_cnv_mean", action="store_true", default=FALSE,
	            help="Use MU as the covariate for BBREG, otherwise uses CNV as the covariate (requires --bbreg). [default: %default]"),	            
	make_option("--min_cnv_round", action="store", default=NA , type='double',
              help="Absolute CNV values below this cutoff get set to zero. [default: %default]"),
	make_option("--mask_cnv_cutoff", action="store", default=NA , type='double',
              help="Mask out any sites that have an absolute CNV value above this cutoff (NA = no masking). [default: %default]"),
	make_option("--exclude", action="store", default=75 , type='integer',
              help="The mimium distance between SNPs allowed in the haplotype. [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

if ( !is.na(opt$seed) ) set.seed( opt$seed )

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
DO.MBASED = opt$mbased
DO.TOTAL = FALSE
DO.COMBINE = opt$combine

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

if ( DO.MBASED ) {
  library("MBASED")
}

if ( NUM.PERM > 0 && NUM.PERM_COND > 0 ) {
  stop("ERROR: --perm and --perm_cond cannot both be set\n")
}

if ( NUM.PERM < 0 || NUM.PERM_COND < 0 ) {
  stop("ERROR: --perm or --perm_cond cannot be negative\n")
}
# ---

if ( opt$predict ) {
	library("glmnet")
	
	pred.lasso.binom = function( x , y1 , y2  ) {
		keep = !is.na( apply(x,2,sd) )
		coef = rep(0,length(keep))
		if ( sum(keep) >= 5 ) {		
			glm = cv.glmnet( y = cbind(y1,y2) , x = x[,keep] , alpha=1 , family="binomial" )
			coef[ keep ] = coef( glm , s="lambda.min" )[-1]
		}
		return( coef )
	}			
	pred.lasso = function( x , y , w = NULL ) {
		keep = !is.na( apply(x,2,sd) )
		coef = rep(0,length(keep))
		if ( sum(keep) >= 5 && sd(y) != 0 ) {		
			if ( is.null(w) ) glm = cv.glmnet( y = y , x = x[,keep] , alpha=1 )
			else glm = cv.glmnet( y = y , x = x[,keep] , alpha=1 , weights=w )
			coef[ keep ] = coef( glm , s="lambda.min" )[-1]
		}
		return( coef )
	}
	pred.combined = function( x1 , y1 , x2 , y2 ) {
		x = rbind( x1 , x2 )
		y = c( y1 , y2 )
		return( pred.lasso( x , y ) )
	}
	pred.interaction = function( x1 , y1 , x2 , y2 ) {
		M = ncol(x1)
		x.a = rbind( x1 , x2 )
		x.b = rbind( x1 , matrix(0,nrow=nrow(x2),ncol=ncol(x2)) )
		x.i = c( rep(1,nrow(x1)) , rep(0,nrow(x2)) )
		x = cbind(x.a,x.b,x.i)
		y = c( y1 , y2 )
		
		coef = pred.lasso( x , y )
		coef.sum = coef[ 1:(M) ] + coef[ (M+1):(2*M) ]
		return( coef.sum )
	}
	
	pred.marginal = function( x , y , top = TRUE ) {
		eff.wgt = t( x ) %*% (y) / sqrt( length(y) - 1 )
		eff.wgt[ is.na(eff.wgt) ] = 0
		if ( top ) eff.wgt[ - which.max( eff.wgt^2 ) ] = 0
		return( eff.wgt )
	}

	hereg = function( x , y ) {
		sx = scale(x)
		sx = sx[ , apply( is.na(sx) , 2 , sum ) == 0 , drop=F ]

		K = sx %*% t( sx ) / ncol(sx)
		phe.prod  = scale(y) %*% t(scale(y))
		diag(phe.prod) = NA
		vec.phe = c(phe.prod[lower.tri(phe.prod)])
		vec.grm = c(K[lower.tri(K)])
		keep = !is.na(vec.phe)
		vec.phe = vec.phe[keep]
		vec.grm = vec.grm[keep]
		hsq = lm(vec.phe ~ vec.grm)		
		return( hsq$coef[2] )
	}
	
	pred.train = function( x.tot , x.hap , y.tot , y.h1 , y.h2 , output , snps , folds=5 ) {
		options(warn=-1)
		
		hap.wgt = y.h1 + y.h2		
		y.hap = log( y.h1 / y.h2 )
		
		# remove NA's
		tot.keep = !is.na(y.tot)
		x.tot = x.tot[tot.keep,]
		y.tot = y.tot[tot.keep]
		hap.keep = !is.na(y.hap)
		x.hap = x.hap[hap.keep,]
		y.hap = y.hap[hap.keep]
		hap.wgt = hap.wgt[hap.keep]
		y.h1 = y.h1[ hap.keep ]
		y.h2 = y.h2[ hap.keep ]
		
		# cross validation
		tot.ord = sample(1:nrow(x.tot))
		tot.cut = floor(seq(0,length(tot.ord),length.out=(folds+1)))

		hap.ord = sample(1:nrow(x.hap))
		hap.cut = floor(seq(0,length(hap.ord),length.out=(folds+1)))

		models = c("lasso","lasso.as","lasso.combined","top1.as","top1","top1.combined")
		n.models = length(models)
						
		tot.pred = matrix(NA,nrow=nrow(x.tot),ncol=n.models)
		hap.pred = matrix(NA,nrow=nrow(x.hap),ncol=n.models)
		
		x.tot.scaled = scale(x.tot)
		x.hap.scaled = scale(x.hap)
		y.tot.scaled = scale(y.tot)
		y.hap.scaled = scale(y.hap)
		
		x.hap.scaled[ is.na(x.hap.scaled) ] = 0

		for ( i in 2:(folds+1) ) {
			batch = (tot.cut[i-1]+1):tot.cut[i]
			tot.heldout = tot.ord[ batch ]
			
			if ( sum(tot.keep) >= opt$min_n_pred ) {
				# lasso
				c.tot = pred.lasso( x = x.tot[-tot.heldout,] , y = y.tot.scaled[-tot.heldout] )
				tot.pred[ tot.heldout , 1 ] = x.tot[tot.heldout,] %*% c.tot
				c.tot.top = pred.marginal( x = x.tot.scaled[-tot.heldout,] , y = y.tot.scaled[-tot.heldout] , top=FALSE )
				
				# top1
				c.tot.top.single = c.tot.top
				c.tot.top.single[ - which.max( c.tot.top^2 ) ] = 0
				tot.pred[ tot.heldout , 5 ] = x.tot.scaled[tot.heldout,] %*% c.tot.top.single
			}
			
			if ( sum(hap.keep) >= opt$min_n_pred ) {
				batch = (hap.cut[i-1]+1):hap.cut[i]
				hap.heldout = hap.ord[ batch ]
				
				c.hap = pred.lasso( x = x.hap[-hap.heldout,] , y = y.hap[-hap.heldout] , w = sqrt(hap.wgt[-hap.heldout]) )
				c.hap.top = pred.marginal( x = x.hap.scaled[-hap.heldout,] , y = y.hap.scaled[-hap.heldout] , top=FALSE )

				c.combined = pred.combined( x1 = x.tot.scaled[-tot.heldout,] , x2 = x.hap.scaled[-hap.heldout,] , y1 = y.tot.scaled[-tot.heldout] , y2 = y.hap.scaled[-hap.heldout] )

				# combined statistic
				c.combined.top = ( c.tot.top + c.hap.top ) / sqrt(2)
				# remove cases where either test failed
				c.combined.top[ c.tot.top == 0 | c.hap.top == 0 ] = 0
				# get the max QTL
				c.combined.top[ - which.max( c.combined.top^2 ) ] = 0
				c.hap.top[ - which.max( c.hap.top^2 ) ] = 0
				
				# lasso.as
				tot.pred[ tot.heldout , 2 ] = (x.tot[tot.heldout,] - 1) %*% c.hap
				# lasso.combined
				tot.pred[ tot.heldout , 3 ] = x.tot.scaled[tot.heldout,] %*% c.combined
				# top1.as
				tot.pred[ tot.heldout , 4 ] = x.tot.scaled[tot.heldout,] %*% c.hap.top

				# lasso
				hap.pred[ hap.heldout , 1 ] = x.hap[hap.heldout,] %*% c.tot
				# lasso.as
				hap.pred[ hap.heldout , 2 ] = x.hap[hap.heldout,] %*% c.hap
				# lasso.combined
				hap.pred[ hap.heldout , 3 ] = x.hap.scaled[hap.heldout,] %*% c.combined	
				# top1.as
				hap.pred[ hap.heldout , 4 ] = x.hap.scaled[hap.heldout,] %*% c.hap.top
				
				# top1.combined
				hap.pred[ hap.heldout , 6 ] = x.hap.scaled[hap.heldout,] %*% c.combined.top
				# top1.combined
				tot.pred[ tot.heldout , 6 ] = x.tot.scaled[tot.heldout,] %*% c.combined.top
				
				if ( sum(tot.keep) >= opt$min_n_pred ) {
					hap.pred[ hap.heldout , 5 ] = x.hap[hap.heldout,] %*% c.tot.top.single
				}
			}
		}
		
		cv.performance = matrix( NA , nrow=4 , ncol=n.models )
		rownames(cv.performance) = c("hap.rsq","hap.pval","rsq","pval")
		colnames(cv.performance) = models
		for ( i in 1:n.models ) {
			try( { tst = cor.test(hap.pred[,i],y.hap); cv.performance[1,i] = tst$est^2 - 1/sum(hap.keep); cv.performance[2,i] = tst$p.value } , silent=T )
			try( { tst = cor.test(tot.pred[,i],y.tot); cv.performance[3,i] = tst$est^2 - 1/sum(tot.keep); cv.performance[4,i] = tst$p.value } , silent=T )
		}

		# --- final training
		wgt.matrix = matrix( NA , nrow=ncol(x.tot) , ncol=n.models )
		rownames( wgt.matrix ) = snps[,2]
		colnames( wgt.matrix ) = models
		# lasso
		wgt.matrix[,1] = pred.lasso( x = x.tot , y = y.tot.scaled )
		# top1
		wgt.matrix[,5] = pred.marginal( x = x.tot.scaled , y = y.tot.scaled , top=FALSE )
		# hap models
		if ( sum(hap.keep) >= opt$min_n_pred ) {
			# lasso as
			wgt.matrix[,2] = pred.lasso( x = x.hap , y = y.hap , w = sqrt(hap.wgt) )
			# lasso combined
			wgt.matrix[,3] = pred.combined( x1 = x.tot.scaled , x2 = x.hap.scaled , y1 = y.tot.scaled , y2 = y.hap.scaled )
			# top1 as
			wgt.matrix[,4] = pred.marginal( x = x.hap.scaled , y = y.hap.scaled , top=FALSE )
			# top1 combined
			wgt.matrix[,6] = (wgt.matrix[,4] + wgt.matrix[,5]) / sqrt(2)
		}

		N.as = nrow(x.hap)
		N.qt = nrow(x.tot)
		N.tot = N.as + N.qt
		
		# heritability
		hsq.as = NA
		hsq.qtl = NA
		hsq = 0		
		hsq.pv = 0
		
		try( {hsq.as = hereg( x.hap , y.hap )} , silent=TRUE )
		try( {hsq.qtl = hereg( x.tot , y.tot )} , silent=TRUE )
		# ---
		
		# check if any models have non-zero weights
		if( sum(apply(wgt.matrix!=0,2,sum,na.rm=T)!=0,na.rm=T) > 0 ) {
			save( wgt.matrix , snps , cv.performance , hsq, hsq.pv , hsq.as , hsq.qtl ,  N.tot , N.as , N.qt , file = paste( "WEIGHTS/" , output , ".wgt.RDat" , sep='' ) )
			cat( output , hsq.as , hsq.qtl , c(cv.performance) , '\n' , sep='\t' )
		}
		options(warn=0)
	}	
}

message("Files being read in")
peaks = read.table( opt$peaks , head=T , as.is=T)
mat = read.table( opt$input , as.is=T )
phe = read.table( opt$samples , head=T , as.is=T)

cnv.all = read.table( opt$global_param , head=T ,as.is=T)
if ( !is.na(opt$total_matrix) ) {
	total.mat = as.matrix( read.table( opt$total_matrix , head=T , check.names=FALSE , row.names=1 ) )
	DO.TOTAL = TRUE
}
if ( !is.na(opt$covar) ) {
	covar.mat = read.table( opt$covar , head=T , check.names=FALSE , row.names=1 )
}

if ( DO.COMBINE & !DO.TOTAL ) {
	stop("ERROR: Cannot have --combine without proper --total_matrix input")
}

message("Files have been read in")

PERM.PVTHRESH = 0.05 / nrow(mat)

bb.simulate = function( tot , rho , mu ) {
        keep = !is.na(tot) & tot > 0

        ret = list()
        ret[["ALT"]] = rep(0,length(keep))
        ret[["ALT"]][ keep ] = rbetabinom(n = sum(keep), size = tot[keep], prob = mu[keep] , rho = rho[keep] )
        ret[["REF"]] = as.vector(tot - ret[["ALT"]])       
        return( ret )
}

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

# ref = CUR.REF.C0; alt = CUR.ALT.C0; rho = RHO[ CUR.IND.C0 ]; covar = COVAR.C0 

bbreg.test = function( ref , alt , rho , covar , cond=NULL ) {
    if ( !is.null(cond) ) {
	ret = list( "pv" = c(NA,NA,NA) )
    	if ( sum(!is.na(covar)) == 0 ) return( ret )
      df = data.frame( y=alt , n=ref+alt , cond = cond , covar = covar , rho=rho )      
    } else {
	ret = list( "pv" = c(NA,NA) )
    	if ( sum(!is.na(covar)) == 0 ) return( ret )
      df = data.frame( y=alt , n=ref+alt , covar = covar , rho=rho )
    }

    # test with a fixed overd parameter for each individual (for some reason this is very SLOW!)
    # reg = betabin( cbind( y , n - y  ) ~ 1 + cond , ~ phi.group , df , fixpar = list( 3:(n+2) , rho ) )
    df = df[ !is.na(df$covar) & df$n > 0, ]
    
    # for efficiency, discretize the overd parameters into five groups
    nq = min( 5 , length(unique(rho)) )
    phi.q = rep(NA,nrow(df))
    rho.q = rep(NA,nq)
    qq = quantile(unique(df$rho), probs = seq(0, 1, 1/nq))
    for ( i in 1:nq ) {
      keep = df$rho >= qq[i] & df$rho <= qq[i+1]
      phi.q[ keep ] = i
      rho.q[ i ] = mean( df$rho[keep] )
    }
    df$phi.q = as.factor(phi.q)

    ## wrapping the betabin call in a silent try to supress convergance issues
    ret = tryCatch( {
    if ( !is.null(cond) ) {
#      if ( sd(df$covar) == 0 || sd(df$cond) == 0 || cor(df$covar,df$cond) == 1 ) return( ret.null )

	  # test each marginal term to ensure convergence
      reg.2 = betabin( cbind( y , n - y ) ~ 1 + cond , ~ phi.q , df , fixpar = list( 3:(nq+2) , rho.q ) )
      reg.3 = betabin( cbind( y , n - y ) ~ 1 + covar , ~ phi.q , df , fixpar = list( 3:(nq+2) , rho.q ) )
      reg = betabin( cbind( y , n - y ) ~ 1 + cond + covar , ~ phi.q , df , fixpar = list( 4:(nq+3) , rho.q ) )
    } else {
#      if ( sd(df$covar) == 0 ) return( ret.null )
      reg = betabin( cbind( y , n - y  ) ~ 1 + covar , ~ phi.q , df , fixpar = list( 3:(nq+2) , rho.q ) )
    }

#    if ( sum(!is.na(reg@varparam)) == 0 ) return( ret.null )   
    zscores = coef(reg) / sqrt(diag(vcov( reg )))
    pvals = 2*(pnorm( abs(zscores) , lower.tail=F))
    opt = list( "pv" = pvals )
    opt } , silent = TRUE , warning = function(w){ return(ret) } , error = function(e) { return(ret) } )
    return( ret )
}

mbased.test = function( CUR.REF.C0 , CUR.ALT.C0 , CUR.REF.C1 , CUR.ALT.C1 , RHO.C0 , RHO.C1 , MAX.SITES = NULL ) {

	MB.N = max(length(CUR.REF.C0),length(CUR.REF.C1))

	MB.REF.C0 = rep(0,MB.N)
	MB.REF.C1 = MB.REF.C0
	MB.ALT.C0 = MB.REF.C0
	MB.ALT.C1 = MB.REF.C0
	MB.RHO.C0 = MB.REF.C0
	MB.RHO.C1 = MB.REF.C0

	MB.REF.C0[ 1:length(CUR.REF.C0) ] = CUR.REF.C0
	MB.REF.C1[ 1:length(CUR.REF.C1) ] = CUR.REF.C1
	MB.ALT.C0[ 1:length(CUR.ALT.C0) ] = CUR.ALT.C0
	MB.ALT.C1[ 1:length(CUR.ALT.C1) ] = CUR.ALT.C1
	MB.RHO.C0[ 1:length(RHO.C0) ] = RHO.C0
	MB.RHO.C1[ 1:length(RHO.C1) ] = RHO.C1

	if ( !is.null( MAX.SITES ) ) {
	MB.N = min( MAX.SITES , MB.N )
	MB.REF.C0 = MB.REF.C0[1:MB.N]
	MB.REF.C1 = MB.REF.C1[1:MB.N]
	MB.ALT.C0 = MB.ALT.C0[1:MB.N]
	MB.ALT.C1 = MB.ALT.C1[1:MB.N]
	MB.RHO.C0 = MB.RHO.C0[1:MB.N]
	MB.RHO.C1 = MB.RHO.C1[1:MB.N]
	}

	mySNVs_2s <- GRanges(
	  seqnames=rep("chr1",MB.N),
	  ranges=IRanges(start=1:MB.N, width=1),
	  aseID=rep("gene1",MB.N),
	  allele1=rep("A",MB.N),
	  allele2=rep("T",MB.N)
	)
	names(mySNVs_2s) <- paste( rep("gene1_SNV",MB.N) , 1:MB.N , sep='')

	## create input RangedSummarizedExperiment object
	myTumorNormalExample <- SummarizedExperiment(
	  assays=list(
		lociAllele1Counts=matrix(
		  c(
			MB.REF.C0,
			MB.REF.C1
		  ),
		  ncol=2,
		  dimnames=list(
		   names(mySNVs_2s), 
		   c('tumor', 'normal')
		  )
		),
		lociAllele2Counts=matrix(
		  c(
			MB.ALT.C0,
			MB.ALT.C1
		  ),
		  ncol=2, 
		  dimnames=list(
		   names(mySNVs_2s), 
		   c('tumor', 'normal')
		  )
		),
		lociCountsDispersions=matrix(
		  c( MB.RHO.C0 , MB.RHO.C1 ),
		  ncol=2, 
		  dimnames=list(
		   names(mySNVs_2s), 
		   c('tumor', 'normal')
		  )
		)   
	  ),
	  rowRanges=mySNVs_2s
	)

	ASEresults_2s <- runMBASED(
	  ASESummarizedExperiment=myTumorNormalExample,
	  isPhased=TRUE,
	  numSim=0
	)
	ret = list( "pv1" = assays(ASEresults_2s)$pValueASE , "min1" = assays(ASEresults_2s)$majorAlleleFrequencyDifference )
	
	# ---
	myTumorNormalExample <- SummarizedExperiment(
	  assays=list(
		lociAllele1Counts=matrix(
		  c(
			MB.REF.C1,
			MB.REF.C0
		  ),
		  ncol=2,
		  dimnames=list(
		   names(mySNVs_2s), 
		   c('tumor', 'normal')
		  )
		),
		lociAllele2Counts=matrix(
		  c(
			MB.ALT.C1,
			MB.ALT.C0
		  ),
		  ncol=2, 
		  dimnames=list(
		   names(mySNVs_2s), 
		   c('tumor', 'normal')
		  )
		),
		lociCountsDispersions=matrix(
		  c( MB.RHO.C1 , MB.RHO.C0 ),
		  ncol=2, 
		  dimnames=list(
		   names(mySNVs_2s), 
		   c('tumor', 'normal')
		  )
		)   
	  ),
	  rowRanges=mySNVs_2s
	)

	ASEresults_2s <- runMBASED(
	  ASESummarizedExperiment=myTumorNormalExample,
	  isPhased=TRUE,
	  numSim=0
	)
	
	ret$pv2 = assays(ASEresults_2s)$pValueASE
	ret$min2 = assays(ASEresults_2s)$majorAlleleFrequencyDifference
	
	return( ret )
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

# order the total matrix
if ( DO.TOTAL ) {
	total.mat = total.mat[ , match( phe$ID , colnames(total.mat)) ,drop=F]
	total.mat = total.mat[ match( peaks$NAME , rownames(total.mat)) , ,drop=F]
}

if ( !is.na(opt$covar) ) {
	covar.mat = t( as.matrix(covar.mat)[ , match( phe$ID , colnames(covar.mat)) ] )
}

options( digits = 4 )

COL.HEADER = c("CHR","POS","RSID","P0","P1","NAME","CENTER","N.HET","N.READS","ALL.AF","ALL.BBINOM.P","C0.AF","C0.BBINOM.P","C1.AF","C1.BBINOM.P","DIFF.BBINOM.P")
COL.HEADER.BINOM = c("ALL.BINOM.P","C0.BINOM.P","C1.BINOM.P","FISHER.OR","FISHER.DIFF.P")
COL.HEADER.INDIV = c("IND.C0","IND.C0.COUNT.REF","IND.C0.COUNT.ALT","IND.C1","IND.C1.COUNT.REF","IND.C1.COUNT.ALT")
COL.HEADER.BBREG = c("C0.BBREG.P","C0.CNV.BBREG.P","C1.BBREG.P","C1.CNV.BBREG.P","ALL.BBREG.P","DIFF.BBREG.P","DIFF.CNV.BBREG.P")
COL.HEADER.TOTAL = c("ALL.TOTAL.Z","ALL.TOTAL.P","C0.TOTAL.Z","C0.TOTAL.P","C1.TOTAL.Z","C1.TOTAL.P","DIFF.TOTAL.Z","DIFF.TOTAL.P")
COL.HEADER.COMBINE = c("ALL.COMBINE.Z","ALL.COMBINE.P","C0.COMBINE.Z","C0.COMBINE.P","C1.COMBINE.Z","C1.COMBINE.P","DIFF.COMBINE.Z","DIFF.COMBINE.P")
COL.HEADER.MBASED = c("DIFF.MBASED1.P","DIFF.MBASED1.AF","DIFF.MBASED2.P","DIFF.MBASED2.AF")

HEAD = COL.HEADER
if ( DO.BINOM ) HEAD = c(HEAD,COL.HEADER.BINOM)
if ( DO.MBASED ) HEAD = c(HEAD,COL.HEADER.MBASED)
if ( DO.BBREG ) HEAD = c(HEAD,COL.HEADER.BBREG)
if ( DO.INDIV ) HEAD = c(HEAD,COL.HEADER.INDIV)
if ( DO.TOTAL ) HEAD = c(HEAD,COL.HEADER.TOTAL)
if ( DO.COMBINE ) HEAD = c(HEAD,COL.HEADER.COMBINE)

if ( !is.na(opt$predict_snps) ) {
	predict_snps = rep(F,nrow(mat))
	predict_snps.lst = read.table(opt$predict_snps,as.is=T)[,1]
	# search by SNP id
	predict_snps[ !is.na( match( mat[,3] , predict_snps.lst ) ) ] = T
	# or by position
	predict_snps[ !is.na( match( paste(mat[,1],mat[,2],sep=':') , predict_snps.lst ) ) ] = T
} else {
	predict_snps = rep(T,nrow(mat))
}

if ( !is.na(opt$update_snps) ) {
	update_snps.lst = read.table(opt$update_snps,as.is=T)
	update_snps.lst = update_snps.lst[ match( paste(mat[,1],mat[,2],sep=':') , update_snps.lst[,1] ) , 2 ]
	update_snps.keep = !is.na(update_snps.lst)
	mat[update_snps.keep,3] = update_snps.lst[ update_snps.keep ]
}

if ( !is.na(opt$keep) ) {
	ind.keep = rep(F,nrow(phe))
	ind.keep.lst = read.table(opt$keep,as.is=T)[,1]
	ind.keep[ !is.na( match( phe$ID , ind.keep.lst ) ) ] = T
} else {
	ind.keep = rep(T,nrow(phe))
}
	
if ( ! opt$predict_only ) {
	cat( HEAD , sep='\t')
	cat('\n')
}

ALL.MAF = apply( (GENO.H1 + GENO.H2)/2,1,mean)

for ( p in 1:nrow(peaks) ) {
	cur = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$P0[p] & mat[,2] <= peaks$P1[p]
	if(sum(cur,na.rm=T) == 0 ) next

	if ( !is.na(opt$local_param) )  {
		cnv.region = cnv.local[ cnv.local$CHR == peaks$CHR[p] & cnv.local$P0 < peaks$P0[p] & cnv.local$P1 > peaks$P1[p] , ]
		m = match( phe$ID , cnv.region$ID )
		cnv.region = cnv.region[m,]
		CUR.RHO = cnv.region$PHI
		CUR.CNV = cnv.region$CNV
		CUR.CNV_MU = cnv.region$MU
		
		if ( opt$fill_cnv ) {
			CUR.RHO[ is.na(CUR.RHO) ] = 0.01
			CUR.CNV[ is.na(CUR.CNV) ] = 0
			CUR.CNV_MU[ is.na(CUR.CNV) ] = 0.5
		}
		if ( !is.na(opt$mask_cnv_cutoff) ) {
			mask = abs(CUR.CNV) > opt$mask_cnv_cutoff
			CUR.CNV[ mask ] = NA
			CUR.RHO[ mask ] = NA
			CUR.CNV_MU[ mask ] = NA
		}
		if ( !is.na(opt$min_cnv_round) ) {
			mask = abs(CUR.CNV) < opt$min_cnv_round
			CUR.CNV[ mask ] = 0
		}
	} else {
		CUR.RHO = RHO.ALL
	}
	CUR.RHO[ CUR.RHO > MAX.RHO ] = NA

	# --- collapse reads at this peak
	cur.h1 = vector()
	cur.h2 = vector()
	cur.i = vector()
	
	cur.h1.tot = rep(0,N)
	cur.h2.tot = rep(0,N)
	
	# !is.na(CUR.RHO) : skip this sample if they have been masked out / have no overdispersion estimate
	for ( i in (1:N)[ ind.keep & !is.na(CUR.RHO) ] ) {
		
		reads.keep = GENO.H1[cur,i] != GENO.H2[cur,i] & HAPS[[1]][cur,i] >= MIN.COV & HAPS[[2]][cur,i] >= MIN.COV & HAPS[[1]][cur,i] + HAPS[[2]][cur,i] > 0
		reads.keep[reads.keep] <- !c(FALSE, diff(mat[cur, 2][reads.keep]) < opt$exclude)
		cur.h1 = c( cur.h1 , (HAPS[[1]][cur,i])[reads.keep] )
		cur.h2 = c( cur.h2 , (HAPS[[2]][cur,i])[reads.keep] )
		cur.i = c( cur.i , rep( i , sum(reads.keep)) )
		
		cur.h1.tot[i] = sum( (HAPS[[1]][cur,i])[reads.keep] )
		cur.h2.tot[i] = sum( (HAPS[[2]][cur,i])[reads.keep] )
	}

	# --- construct the total activity phenotype and covariates
	if ( DO.TOTAL || opt$predict ) {
		TOT.Y = total.mat[ p , ]
		TOT.Y[ !ind.keep ] = NA
		
		# rank normalize the phenotype
		if ( opt$total_rn ) TOT.Y = scale( rank(TOT.Y) / length(TOT.Y) )

		if ( !is.na(opt$covar) && !is.na(opt$local_param) ) {
			LM.COVAR = cbind( covar.mat , CUR.CNV )
		} else if ( !is.na(opt$covar) ) {
			LM.COVAR = covar.mat
		} else if ( !is.na(opt$local_param) ) {
			LM.COVAR = CUR.CNV
		} else {
			LM.COVAR = null
		}
		
		# test that covariates have variance
		covar.sds = apply( LM.COVAR,2,sd,na.rm=T )
		if ( sum(is.na(covar.sds) | covar.sds == 0) > 0 ) {
			cat( "WARNING: Covariates " , colnames( LM.COVAR ) [ is.na(covar.sds) | covar.sds == 0 ] , "had zero variance for peak",p,", skipping this peak\n" , file = stderr())
			next
		}
	}

	# --- define nearby SNPs
	if( PAR.WIN == -1 ) {
		cur.snp = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$P0[p] & mat[,2] <= peaks$P1[p]
	} else {
		cur.snp = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$CENTER[p] - PAR.WIN & mat[,2] <= peaks$CENTER[p] + PAR.WIN
	}
	cur.snp = cur.snp & ALL.MAF > MIN.MAF & ALL.MAF < (1-MIN.MAF)
	if ( sum(cur.snp) == 0 ) next
	
	# --- build predictive models
	if ( opt$predict ) {
		cur.pred.snp = cur.snp & predict_snps
		PRED.GEN = t(GENO.H1[cur.pred.snp,] + GENO.H2[cur.pred.snp,])
		PRED.HAP = t(GENO.H1[cur.pred.snp,] - GENO.H2[cur.pred.snp,])
		#AS.Y = log( cur.h1.tot / cur.h2.tot )

		if ( !is.null(LM.COVAR) & sum(!is.na(TOT.Y)) > 2 ) {
			CUR.TOT.Y = resid( lm( TOT.Y ~ LM.COVAR , na.action="na.exclude" ) )
		} else {
			CUR.TOT.Y = TOT.Y
		}

		# x.tot = PRED.GEN ; x.hap = PRED.HAP ; y.tot = CUR.TOT.Y ; y.h1 = cur.h1.tot ; y.h2 = cur.h2.tot ; output = peaks$NAME[p] ; snps = mat[cur.pred.snp,c(1,3,2,2,4,5)]
		try( { pred.train( x.tot = PRED.GEN , x.hap = PRED.HAP , y.tot = CUR.TOT.Y , y.h1 = cur.h1.tot , y.h2 = cur.h2.tot , output = peaks$NAME[p] , snps = mat[cur.pred.snp,c(1,3,2,2,4,5)] ) } , silent=T )
	}
	if ( opt$predict_only ) next	
	
	# --- test each nearby SNP
	if ( length(unique(cur.i)) > MIN.MAF*N && sum(cur.h1) + sum(cur.h2) > 0 ) {
		for ( s in which(cur.snp) ) {
			
			tst.bbinom.C0 = c(NA,NA)
			tst.bbinom.C1 = c(NA,NA)
			tst.bbinom.ALL = c(NA,NA)
			pv.BOTH = c(NA)
			
			# restrict to hets
			HET = GENO.H1[s,] != GENO.H2[s,] & !is.na(CUR.RHO)

			if ( sum(HET) > MIN.HET*N ) {
				# collect REF/ALT heterozygous haplotypes
				m1 = match( cur.i , which(HET & GENO.H1[s,] == 0) )
				m2 = match( cur.i , which(HET & GENO.H2[s,] == 0) )
				CUR.REF = c( cur.h1[ !is.na(m1) ] , cur.h2[ !is.na(m2) ])
				CUR.ALT = c( cur.h2[ !is.na(m1) ] , cur.h1[ !is.na(m2) ])
				CUR.IND = c( cur.i[ !is.na(m1) ] , cur.i[ !is.na(m2)] )
						
				# re-orient the CNV mean
				if ( !is.na(opt$local_param) ) {
					CUR.ALLELE.CNV_MU = c( CUR.CNV_MU[ cur.i[ !is.na(m1) ] ] - 0.5 , 0.5 - CUR.CNV_MU[ cur.i[ !is.na(m2) ] ] )
				}
				
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

						m.C0 = !is.na(match( CUR.IND , which(CUR.PHENO==0) ))
						m.C1 = !is.na(match( CUR.IND , which(CUR.PHENO==1) ))

						# --- simulation
						if ( opt$sim ) {
							sim.mu = rep( 0.5 , length(CUR.IND) )
							sim.mu[ m.C0 ] = opt$sim_af0
							sim.mu[ m.C1 ] = opt$sim_af1

							if ( opt$sim_cnv ) {
								# alternative read count (relative to 1)
								sim.alt.ct = sim.mu / (1 - sim.mu)
								# reference read count (1 at baseline)
								sim.ref.ct = rep(1,length(sim.alt.ct))
	
								# cnv multiplier
								if ( opt$sim_cnv_allelic ) {
									sim.cnv.ct = 2 * ( 2^( abs(CUR.CNV[CUR.IND]) ) ) - 1
								} else {
									sim.cnv.ct = 2 * ( 2^( CUR.CNV[CUR.IND] ) ) - 1
								}							
	
								# total fraction = expected alt count * cnv multiplier
								sim.mu = (sim.alt.ct * sim.cnv.ct) / (sim.ref.ct + sim.alt.ct*sim.cnv.ct)
							}
							sim.mu[ sim.mu > 1 ] = 1
							sim.mu[ sim.mu < 0 ] = 0						
							cur.sim = bb.simulate( CUR.REF + CUR.ALT , CUR.RHO[CUR.IND] , sim.mu )
							CUR.REF = cur.sim$REF
							CUR.ALT = cur.sim$ALT
						}

						# merge calls from the same individual
						if ( opt$collapse_reads ) {
							NEW.REF = vector()
							NEW.ALT = vector()
							NEW.IND = vector()
							for ( i in unique(CUR.IND) ) {
								NEW.IND = c(NEW.IND,i)
								NEW.REF = c(NEW.REF,sum(CUR.REF[CUR.IND==i]))
								NEW.ALT = c(NEW.ALT,sum(CUR.ALT[CUR.IND==i]))
							}
							CUR.REF = NEW.REF
							CUR.ALT = NEW.ALT
							CUR.IND = NEW.IND
							
							m.C0 = !is.na(match( CUR.IND , which(CUR.PHENO==0) ))
							m.C1 = !is.na(match( CUR.IND , which(CUR.PHENO==1) ))
						}
				
						CUR.REF.C0 = CUR.REF[m.C0]
						CUR.ALT.C0 = CUR.ALT[m.C0]
						CUR.IND.C0 = CUR.IND[m.C0]

						CUR.REF.C1 = CUR.REF[m.C1]
						CUR.ALT.C1 = CUR.ALT[m.C1]
						CUR.IND.C1 = CUR.IND[m.C1]
						
						# --- perform beta-binomial test
						tst.bbinom.C0 = bbinom.test( CUR.REF.C0 , CUR.ALT.C0 , CUR.RHO[CUR.IND.C0] )
						tst.bbinom.C1 = bbinom.test( CUR.REF.C1 , CUR.ALT.C1 , CUR.RHO[CUR.IND.C1] )
						tst.bbinom.ALL = bbinom.test( CUR.REF , CUR.ALT , CUR.RHO[CUR.IND] )
						lrt.BOTH = 2 * (tst.bbinom.C0$objective + tst.bbinom.C1$objective - tst.bbinom.ALL$objective)
						pv.BOTH = pchisq( abs(lrt.BOTH) , df=1 , lower.tail=F )

						# --- print main output
						if ( perm > 0 ) cat( "PERM_" , perm , ":" , sep='' )
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
						
						# --- perform MBASED test
						if ( DO.MBASED ) {
							tst.mbased = list( "pv1" = NA , "pv2" = NA , "min1" = NA , "min2" = NA )
							try( { tst.mbased = mbased.test( CUR.REF.C0 , CUR.ALT.C0 , CUR.REF.C1 , CUR.ALT.C1 , CUR.RHO[CUR.IND.C0] , CUR.RHO[CUR.IND.C1] ) } , silent=T )
							cat( "" , tst.mbased$min1 , tst.mbased$pv1 , tst.mbased$min2 , tst.mbased$pv2 , sep='\t' )
						}
							
						# --- perform beta-binom regression							
						if ( DO.BBREG ) {
							# --- set up covariates for BBREG		
							if ( opt$cond_cnv_mean ) {
								# mean is already computed per allele
								BBREG.COVAR.C0 = CUR.ALLELE.CNV_MU[m.C0]
								BBREG.COVAR.C1 = CUR.ALLELE.CNV_MU[m.C1]
								BBREG.COVAR = CUR.ALLELE.CNV_MU
							} else {	
								# CNV is already computed per individual						
								BBREG.COVAR.C0 = CUR.CNV[CUR.IND.C0]
								BBREG.COVAR.C1 = CUR.CNV[CUR.IND.C1]
								BBREG.COVAR = CUR.CNV[CUR.IND]					
							}
							
							tst.bbreg.c0 = list("pv" = c(NA,NA) )
							tst.bbreg.c1 = list("pv" = c(NA,NA) )
							tst.bbreg = list("pv" = c(NA,NA,NA) )
							
							try( { tst.bbreg.c0 = bbreg.test( CUR.REF.C0 , CUR.ALT.C0 , CUR.RHO[ CUR.IND.C0 ] , BBREG.COVAR.C0 ) } , silent = F )
							try( { tst.bbreg.c1 = bbreg.test( CUR.REF.C1 , CUR.ALT.C1 , CUR.RHO[ CUR.IND.C1 ] , BBREG.COVAR.C1 ) } , silent = F )
							if( !is.na(min(tst.bbreg.c0$pv)) && !is.na(min(tst.bbreg.c1$pv)) ) {
								try( { tst.bbreg = bbreg.test( CUR.REF , CUR.ALT , CUR.RHO[ CUR.IND ] , BBREG.COVAR , CUR.PHENO[CUR.IND] ) } , silent = F )
							}
							cat( "" , tst.bbreg.c0$pv , tst.bbreg.c1$pv , tst.bbreg$pv , sep='\t' )
						}
						
						if ( DO.TOTAL ) {
							GEN = GENO.H1[s,] + GENO.H2[s,]
							if ( perm > 0 ) GEN = sample(GEN)
							
							if ( !is.null( LM.COVAR ) ) {
								reg.out.all = c(NA,NA)
								try( { reg.out.all = summary(lm( TOT.Y ~ GEN + LM.COVAR ))$coef[2,c(3,4)] } , silent=T )
								cat( "" , reg.out.all , sep='\t' )
								
								reg.out.c0 = c(NA,NA)
								try( { reg.out.c0 = summary(lm( TOT.Y[CUR.PHENO==0] ~ GEN[CUR.PHENO==0] + LM.COVAR[CUR.PHENO==0,] ))$coef[2,c(3,4)] } , silent=T )
								cat( "" , reg.out.c0 , sep='\t' )
		
								reg.out.c1 = c(NA,NA)
								try( { reg.out.c1 = summary(lm( TOT.Y[CUR.PHENO==1] ~ GEN[CUR.PHENO==1] + LM.COVAR[CUR.PHENO==1,] ))$coef[2,c(3,4)] } , silent=T )
								cat( "" , reg.out.c1 , sep='\t' )

								reg.out.d = c(NA,NA)
								# check if the interaction term has variance
								if ( !is.na(sd(GEN*CUR.PHENO,na.rm=T)) & sd(GEN*CUR.PHENO,na.rm=T) > 0 ) {
									try( { reg.out.d = summary(lm( TOT.Y ~ GEN*CUR.PHENO + GEN + CUR.PHENO + LM.COVAR ))$coef[2,c(3,4)] } ,silent=T )
								}
								cat( "" , reg.out.d , sep='\t' )
							} else {
								reg.out.all = c(NA,NA)
								try( { reg.out.all = summary(lm( TOT.Y ~ GEN ))$coef[2,c(3,4)] } , silent=T )
								cat( "" , reg.out.all , sep='\t' )
								
								reg.out.c0 = c(NA,NA)
								try( { reg.out.c0 = summary(lm( TOT.Y[CUR.PHENO==0] ~ GEN[CUR.PHENO==0] ))$coef[2,c(3,4)] } , silent=T )
								cat( "" , reg.out.c0 , sep='\t' )
		
								reg.out.c1 = c(NA,NA)
								try( { reg.out.c1 = summary(lm( TOT.Y[CUR.PHENO==1] ~ GEN[CUR.PHENO==1] ))$coef[2,c(3,4)] } , silent=T )
								cat( "" , reg.out.c1 , sep='\t' )

								reg.out.d = c(NA,NA)
								# check if the interaction term has variance
								if ( !is.na(sd(GEN*CUR.PHENO,na.rm=T)) & sd(GEN*CUR.PHENO,na.rm=T) > 0 ) {
									try( { reg.out.d = summary(lm( TOT.Y ~ GEN*CUR.PHENO + GEN + CUR.PHENO ))$coef[2,c(3,4)] } ,silent=T )
								}
								cat( "" , reg.out.d , sep='\t' )
							}
							
							# Combine BBINOM and TOTAL tests
							if ( DO.COMBINE ) {
								comb.out = c(NA,NA)
								try ( {
									bbinom.z = qnorm(tst.bbinom.ALL$pv/2,lower.tail=F)
									# sign the Z-score based on the AF
									if ( tst.bbinom.ALL$min < 0.5 ) bbinom.z = -1 * bbinom.z
									# stouffers method
									comb.out[1] = (bbinom.z + reg.out.all[1])/sqrt(2)
									# convert to p-value
									comb.out[2] = 2*(pnorm( abs(comb.out[1]) , lower.tail=F))
									},silent=T )
								cat( "" , comb.out , sep='\t' )
								
								comb.out = c(NA,NA)
								try ( {
									bbinom.z = qnorm(tst.bbinom.C0$pv/2,lower.tail=F)
									if ( tst.bbinom.C0$min < 0.5 ) bbinom.z = -1 * bbinom.z
									comb.out[1] = (bbinom.z + reg.out.c0[1])/sqrt(2)
									comb.out[2] = 2*(pnorm( abs(comb.out[1]) , lower.tail=F))
									},silent=T )
								cat( "" , comb.out , sep='\t' )
								
								comb.out = c(NA,NA)
								try ( {
									bbinom.z = qnorm(tst.bbinom.C1$pv/2,lower.tail=F)
									if ( tst.bbinom.C1$min < 0.5 ) bbinom.z = -1 * bbinom.z
									comb.out[1] = (bbinom.z + reg.out.c1[1])/sqrt(2)
									comb.out[2] = 2*(pnorm( abs(comb.out[1]) , lower.tail=F))
									},silent=T )
								cat( "" , comb.out , sep='\t' )

								comb.out = c(NA,NA)
								try ( {
									bbinom.z = qnorm(pv.BOTH/2,lower.tail=F)
									# sign the z-score based on difference in AFs
									if ( tst.bbinom.C1$min < tst.bbinom.C0$min ) bbinom.z = -1 * bbinom.z
									comb.out[1] = (bbinom.z + reg.out.d[1])/sqrt(2)
									comb.out[2] = 2*(pnorm( abs(comb.out[1]) , lower.tail=F))
									},silent=T )
								cat( "" , comb.out , sep='\t' )
							}
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
