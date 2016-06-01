
source('core.lib.R')
library('binom')
library('VGAM')

read.data <- function(in.file){
  if(length(grep('\\.rds$', in.file)) > 0){
    vars.filt = readRDS(in.file)
  }else{
    colClasses = c(rep('character', 2), rep('integer', 2))
    vars.filt = read.table(in.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, colClasses = colClasses)
    colnames(vars.filt) = c('gene', 'snp', 'alt.dp', 'ref.dp')
    if(ncol(vars.filt) == 5){
      colnames(vars.filt)[5] = 'p.sample.null'
      vars.filt[, 5] = as.numeric(vars.filt[, 5])
    }
    if(ncol(vars.filt) == 6){ #simulated data where class is known
      colnames(vars.filt)[5:6] = c('p.sample.null', 'class')
      vars.filt[, 5] = as.numeric(vars.filt[, 5])
    }

  }
  
  return(vars.filt)
}

get.main.pval <- function(obs.df, method = 'wilson', n.samples, nmax.vars, betabin.rho = 0, p.sample.null = 0.5, testval.col = 'z', gene.null.nvar.list = NULL){
  
  ###
  #Obs dist
  ###

  #Adds dp col (dp = alt.dp + ref.dp) if a dp col doesnt already exist
  obs.df = add.dp.col(obs.df)
  
  #Get SNV test-stat (z)
  obs.df = get.snv.testval(obs.df, method)

  #Add column: n.vars
  obs.df = add.nvar.column(obs.df)
    
  #Calc gene based measures for each gene based on the snv.testval (median, mean, cv, Stouffer)
  obs.gene.df = get.gene.testval(obs.df, testval.col = testval.col) 
 
  #Add column 'nvars.bin'
  obs.gene.df = add.nvars.bin(obs.gene.df, nmax.vars)

  
  ###
  #Null dist
  ###
  #Generate null by binomial sampling with same depths as the snvs within each gene, that is, one null per snv combined by stouffer into one null per gene
  #(Alt: compare the frac:s directly (estimated vs sampled), then Stouffer on those p-values)

  #get p.sample.null
  cols = colnames(obs.df)
  if(length(which(cols == 'p.sample.null')) > 0){
    p.sample.null = 'p.sample.null'
  }
  
  if(is.null(gene.null.nvar.list)){
    liptak.s.p.nom = get.binresamp.pval(obs.df, obs.gene.df, n.samples, method, betabin.rho, p.sample.null)
    obs.gene.df = cbind(obs.gene.df, liptak.s.p.nom)
  }
  else{ #pre-calculated null dists

    #split and loop on nvars.bin
    obs.gene.df.nvar.list = split(obs.gene.df, obs.gene.df[, 'nvars.bin'])
    nvars.bins = names(obs.gene.df.nvar.list)
    for(j.bin in nvars.bins){
      gene.null = gene.null.nvar.list[[j.bin]]
      obs.gene.df.nvar = obs.gene.df.nvar.list[[j.bin]]
      
      #Get p-val by comparing to null
      liptak.s.p.nom = sapply(obs.gene.df.nvar[, 'liptak.s'], function(obs.val, gene.null){p.nom = length(which(gene.null >= obs.val)) / length(gene.null)}, gene.null = gene.null)
      obs.gene.df.nvar.list[[j.bin]] = cbind(obs.gene.df.nvar, liptak.s.p.nom)
    }

    #rbind
    obs.gene.df = as.data.frame(matrix(nrow = 0, ncol = ncol(obs.gene.df.nvar.list[[1]])))
    for(j.bin in nvars.bins){
      obs.gene.df = rbind(obs.gene.df, obs.gene.df.nvar.list[[j.bin]])
    }
  }


  #Multi-test Adjust
  p.nom.cols = colnames(obs.gene.df)[grep('p.nom', colnames(obs.gene.df))]
  for(jcol in p.nom.cols){
    p.bh = p.adjust(obs.gene.df[, jcol], method = 'BH');
    obs.gene.df = cbind(obs.gene.df, p.bh);
    colnames(obs.gene.df)[ncol(obs.gene.df)] = paste(jcol, 'bh', sep = '.');
  }                                                                                                

  
  #Merge with gene class annot (neg/pos)
  #Only applicabable in the case of simulated data (not real data)
  cols = colnames(obs.df)
  if(length(grep('^class$', cols) > 0)){
    gene.annot = unique(obs.df[, c('gene', 'class')])
    obs.gene.df = merge(obs.gene.df, gene.annot, by = 'gene')
  }
  
  return(obs.gene.df)
}

get.snv.testval <- function(data.df, method = 'wilson', eps = 1e-16, alt.dp.col = 'alt.dp', site.dp.col = 'dp', p.binom.null = 0.5){

  #make ref.dp col (doesn't always exist)
  ref.dp = data.df[, site.dp.col] - data.df[, alt.dp.col]
  
  #add pseudo-count if alt.dp or ref.dp == 0
  pseudo.data.df = cbind(data.df[, site.dp.col], data.df[, alt.dp.col], ref.dp)
  colnames(pseudo.data.df) = c('dp', 'alt.dp', 'ref.dp')
  pseudo.data.df = add.pseudo(pseudo.data.df)

  #get dp cols
  dp = pseudo.data.df[, 'dp']
  alt.dp = pseudo.data.df[, 'alt.dp']
  ref.dp = pseudo.data.df[, 'ref.dp']
  
  if(method == 'p.binom'){ #use p-value from binom.test to calculate Z
    
    #get p-value from binom.test
    x = cbind(alt.dp, dp)
    bt.res = t(apply(x, 1, function(x, p){res = binom.test(x['alt.dp'], x['dp'], p = p, alternative = 'two.sided'); return(c(res$p.value, res$estimate, res$conf.int[1:2]));}, p = p.binom.null))
    colnames(bt.res) = c('pval', 'p', 'p.ci.lower', 'p.ci.upper')
    p = bt.res[, 'p']
    p.ci.upper = bt.res[, 'p.ci.upper']
    p.ci.lower = bt.res[, 'p.ci.lower']
    
    #log.odds. NB: abs needed since unphased data
    abs.log.odds = abs(log(p / (1 - p)))

    #log.odds.se
    log.odds.se = (log(p.ci.upper) - log(p.ci.lower)) / (2 * 1.96)

    #Z score from inv.norm(1-pval)    
    z = pval2zscore(bt.res[, 'pval'], eps)
  }
  else{ #binom.confint
    
    #p
    p = alt.dp / dp

    #p.CI
    ci = binom.confint(alt.dp, dp, methods = method) #Agresti-Coull no good since some estimates < 0 and >1!
    p.ci.upper = ci[, 'upper']
    p.ci.lower = ci[, 'lower']

    #log.odds. NB: abs needed since unphased data
    abs.log.odds = abs(log(p / (1 - p)))

    #log.odds.se
    log.odds.se = (log(p.ci.upper) - log(p.ci.lower)) / (2 * 1.96)
    
    #half z (since absolute value)
    z = abs.log.odds / log.odds.se
  }

  #bind all
  data.df = cbind(data.df, z, abs.log.odds, log.odds.se, p, p.ci.upper, p.ci.lower)

  return(data.df)
}

get.binresamp.pval <- function(obs.df, obs.gene.df, n.samples, method = 'wilson', betabin.rho = 0, p.sample.null = 0.5){
  
  obs.df.gene.list = split(obs.df, obs.df[, 'gene'])
  genes = names(obs.df.gene.list)
  n.genes = length(genes)

  p.nom = vector(length = n.genes, mode = 'numeric')
  names(p.nom) = genes
  j.gene.it = 0
  for(jgene in genes){
    j.gene.it = j.gene.it + 1

    if((j.gene.it %% 100) == 0){
      cat(paste('j.gene: ', j.gene.it, '\n', sep = ''))
    }
    
    #Get null gene testval (liptak.stouffer of snv testvals only, now not eg cv of snv testvals returned, due to speed and size)
    obs.df.gene = obs.df.gene.list[[jgene]]
    if(is.character(p.sample.null)){
      gene.null = get.gene.null(obs.df.gene[, 'dp'], n.samples, method, obs.df.gene[, p.sample.null], betabin.rho) #vec: len: n.samples
    }else{
      gene.null = get.gene.null(obs.df.gene[, 'dp'], n.samples, method, p.sample.null, betabin.rho) #vec: len: n.samples
    }
      
    #Get p-val by comparing to null
    obs.val = obs.gene.df[which(obs.gene.df[, 'gene'] == jgene), 'liptak.s']
    p.nom[jgene] = length(which(gene.null >= obs.val)) / length(gene.null)
  }
  return(p.nom)
}

get.gene.null <- function(gene.snv.dps, n.samples, method = 'wilson', p.sample.null = 0.5, betabin.rho = 0){

  #get null for each snv
  n.snvs = length(gene.snv.dps)
  snv.null.dist = matrix(nrow = n.samples, ncol = n.snvs)
  if(length(p.sample.null) > 1){ #vector with varying p.sample.null
    for(j.snv in 1:n.snvs){    
      null.df = get.snv.null(gene.snv.dps[j.snv], n.samples, p.sample.null[j.snv], betabin.rho) #out: n.samples alt.dp and dp
      null.df = get.snv.testval(null.df, method)
      snv.null.dist[, j.snv] = null.df[, 'z']
    }

  }else{ #single p.sample.null
    for(j.snv in 1:n.snvs){    
      null.df = get.snv.null(gene.snv.dps[j.snv], n.samples, p.sample.null, betabin.rho) #out: n.samples alt.dp and dp
      null.df = get.snv.testval(null.df, method)
      snv.null.dist[, j.snv] = null.df[, 'z']
    }
  }
  #get gene testval by combining snv null testvals by liptak.stouffer
  gene.null = apply(snv.null.dist, 1, liptak.stouffer) #vec: len: n.samples

  return(gene.null)
}

get.snv.null <- function(j.dp, n.samples, p.sample.null = 0.5, betabin.rho = 0){

  #sample binomial, according to p = p.sample.null  
  #alt.dp = rbinom(n.samples, j.dp, p.sample.null)
  alt.dp = rbetabinom(n.samples, j.dp, p.sample.null, betabin.rho)
  dp = rep(j.dp, n.samples)

  data.df = cbind(alt.dp, dp)

  return(data.df)
}

add.dp.col <- function(data.df, alt.dp.col = 'alt.dp', ref.dp.col = 'ref.dp'){

  cols = colnames(data.df)
  if(length(which(cols == 'dp')) == 0){
    dp = data.df[, alt.dp.col] + data.df[, ref.dp.col]
    data.df = cbind(data.df, dp)
  }
  
  return(data.df)
}

add.pseudo <- function(data.df, alt.dp.col = 'alt.dp', site.dp.col = 'dp'){

  zero.ind = which(data.df[, alt.dp.col] == 0)
  data.df[zero.ind, alt.dp.col] = 1
  data.df[zero.ind, site.dp.col] = data.df[zero.ind, site.dp.col] + 1
  
  zero.ind = which(data.df[, 'ref.dp'] == 0)
  data.df[zero.ind, 'ref.dp'] = 1
  data.df[zero.ind, site.dp.col] = data.df[zero.ind, site.dp.col] + 1

  return(data.df)
}
