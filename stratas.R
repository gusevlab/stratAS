library('VGAM')

args = commandArgs(trailingOnly=T)
f.mat = args[1]

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
		opt = list( "lrt" = NA , "pv" = NA , "min" = NA )
	}
	return( opt )
}
					
peaks = read.table("gencode.protein_coding.transcripts.bed",as.is=T)
mat = read.table(f.mat,as.is=T)
phe = read.table("WASP/ASVCF/KIRC.ALL.AS.PHE",as.is=T)
map = read.table(paste(f.mat,".HIGHMAP.bed",sep=''),as.is=T)

m = match( paste(mat[,1],mat[,2]) , paste(map[,1],map[,3]) )
highmap = !is.na(m)
cat( mean(highmap) , "high mappability sites\n" , file=stderr() )

cnv.all = read.table("WASP/ASVCF/KIRC.ALL.AS.CNV",as.is=T)
cnv.local = read.table("WASP/ASVCF/KIRC.ALL.AS.CNVLOCAL",as.is=T)
#cnv.starts = c(1,cumsum(rle(tbl[,1])$len))
#cnv.starts = cnv.ctr[-length(cnv.ctr)]
cnv.local[,2] = paste("chr",cnv.local[,2],sep='')

cur.chr = unique(mat[,1])
m = match(peaks[,1] , cur.chr )
peaks = peaks[!is.na(m),]

m = match(cnv.local[,2] , cur.chr )
cnv.local = cnv.local[!is.na(m),]

PAR.WIN = 100e3
NUM.PERM = 0
MIN.MAF = 0.01
MIN.HET = 0.01
MAX.RHO = 0.10
DO.BINOM = FALSE
PERM.PVTHRESH = 0.05 / nrow(mat)

N = (ncol(mat) - 5)/4
M = nrow(mat)

HAPS = list()
GENO.H1 = matrix( 0 , nrow=M , ncol=N )
GENO.H2 = matrix( 0 , nrow=M , ncol=N )

PHENO = phe[,2]
m = match( phe[,1] , cnv.all[,1] )
cnv.all = cnv.all[m,]
RHO.ALL = cnv.all[,2]

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

options( digits = 2 )

for ( p in 1:nrow(peaks) ) {
	cur = mat[,1] == peaks[p,1] & mat[,2] >= peaks[p,2] & mat[,2] <= peaks[p,3] & highmap
	cur.cnv = cnv.local[ cnv.local[,2] == peaks[p,1] & cnv.local[,3] < peaks[p,2] & cnv.local[,4] > peaks[p,3] , ]
	m = match( phe[,1] , cur.cnv[,1] )
	cur.cnv = cur.cnv[m,]
	RHO = cur.cnv[,5]
	RHO[ RHO > MAX.RHO ] = NA
	
	# collapse reads at this peak
	cur.h1 = vector()
	cur.h2 = vector()
	cur.i = vector()
	# Remove any sites that are too close together
	snps.overlap = ( c(mat[cur,2],NA) - c(NA,mat[cur,2]) < 100 )[ 1:sum(cur) ]
	snps.overlap[1] = FALSE
	cur[ cur ] = !snps.overlap
	
	for ( i in 1:N ) {
		reads.keep = GENO.H1[cur,i] != GENO.H2[cur,i] & HAPS[[1]][cur,i] >= 5 & HAPS[[2]][cur,i] >= 5
		cur.h1 = c( cur.h1 , (HAPS[[1]][cur,i])[reads.keep] )
		cur.h2 = c( cur.h2 , (HAPS[[2]][cur,i])[reads.keep] )
		cur.i = c( cur.i , rep( i , sum(reads.keep)) )
	}
	
	if ( length(unique(cur.i)) > MIN.MAF*N && sum(cur.h1) + sum(cur.h2) > 0 ) {
		cat( p , p / nrow(peaks) , unlist(peaks[p,]) , '\n' , file=stderr() )
		# test all nearby SNPs
		cur.snp = mat[,1] == peaks[p,1] & mat[,2] >= peaks[p,5] - PAR.WIN & mat[,2] <= peaks[p,5] + PAR.WIN
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
					
					for ( perm in 0:NUM.PERM ) {
						if ( perm > 0 ) {
							# randomly swap REF/ALT alleles
							cur.swap = unique(CUR.IND)[ as.logical(rbinom( length(unique(CUR.IND)) , 1 , 0.5 )) ]
							cur.swap = !is.na(match( CUR.IND , cur.swap ))
							tmp = CUR.REF[ cur.swap ]
							CUR.REF[ cur.swap ] = CUR.ALT[ cur.swap ]
							CUR.ALT[ cur.swap ] = tmp
						}

						m = !is.na(match( CUR.IND , which(PHENO==0) ))
						CUR.REF.C0 = CUR.REF[m]
						CUR.ALT.C0 = CUR.ALT[m]
						CUR.IND.C0 = CUR.IND[m]

						m = !is.na(match( CUR.IND , which(PHENO==1) ))
						CUR.REF.C1 = CUR.REF[m]
						CUR.ALT.C1 = CUR.ALT[m]
						CUR.IND.C1 = CUR.IND[m]		
						
						# --- perform beta-binomial test
						tst.bbinom.C0 = bbinom.test( CUR.REF.C0 , CUR.ALT.C0 , RHO[CUR.IND.C0] )
						tst.bbinom.C1 = bbinom.test( CUR.REF.C1 , CUR.ALT.C1 , RHO[CUR.IND.C1] )
						tst.bbinom.BOTH = bbinom.test( CUR.REF , CUR.ALT , RHO[CUR.IND] )
						lrt.BOTH = 2 * (tst.bbinom.C0$objective + tst.bbinom.C1$objective - tst.bbinom.BOTH$objective)
						pv.BOTH = pchisq( abs(lrt.BOTH) , df=1 , lower.tail=F )			

						# --- perform binomial test
						if ( DO.BINOM ) {
							tst.binom = binom.test( sum(CUR.ALT),sum(CUR.REF)+sum(CUR.ALT) )
							if ( sum(CUR.REF.C0)+sum(CUR.ALT.C0) > 0 ) tst.binom.C0 = binom.test( sum(CUR.ALT.C0),sum(CUR.REF.C0)+sum(CUR.ALT.C0) ) else tst.binom.C0 = list("p.value"=NA)
							if ( sum(CUR.REF.C1)+sum(CUR.ALT.C1) > 0 ) tst.binom.C1 = binom.test( sum(CUR.ALT.C1),sum(CUR.REF.C1)+sum(CUR.ALT.C1) ) else tst.binom.C1 = list("p.value"=NA)
							# --- perform fisher's exact test between conditions					
							tst.fisher = fisher.test( cbind( c(sum(CUR.REF.C0),sum(CUR.ALT.C0)) , c(sum(CUR.REF.C1),sum(CUR.ALT.C1)) ) )
							
							cat( unlist(mat[s,1:3]) , unlist(peaks[p,-1]) , sum(HET) , sum(CUR.REF) + sum(CUR.ALT) , tst.binom$est , tst.binom$p.value , tst.bbinom.BOTH$min , tst.bbinom.BOTH$pv , tst.fisher$est , tst.binom.C0$p.value , tst.bbinom.C0$pv , tst.binom.C1$p.value , tst.bbinom.C1$pv , tst.fisher$p.value , pv.BOTH , '\n' , sep='\t' )
						} else {
							cat( unlist(mat[s,1:3]) , unlist(peaks[p,-1]) , sum(HET) , sum(CUR.REF) + sum(CUR.ALT) , tst.bbinom.C0$min , tst.bbinom.C0$pv , tst.bbinom.C1$min , tst.bbinom.C1$pv , pv.BOTH , '\n' , sep='\t' )
						}
						
						#cat( perm , p , s , sum(HET) , sum(CUR.REF) + sum(CUR.ALT) , unlist(mat[s,1:3]) , unlist(peaks[p,-1]) , tst.binom$est , tst.binom$p.value , tst.bbinom.BOTH$min , tst.bbinom.BOTH$pv , tst.fisher$est , tst.binom.C0$p.value , tst.bbinom.C0$pv , tst.binom.C1$p.value , tst.bbinom.C1$pv , tst.fisher$p.value , pv.BOTH , '\n' , sep='\t' )
						# if ( perm == 0 && min( tst.binom$p.value , tst.bbinom.BOTH$pv , tst.fisher$p.value , tst.bbinom.C0$pv , tst.bbinom.C1$pv , pv.BOTH , na.rm=T) > PERM.PVTHRESH ) break()
					}
				}
			}
		}
	}
}
