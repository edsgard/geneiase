
source('core.lib.R')
library('VGAM')

read.data <- function(in.file){
  if(length(grep('\\.rds$', in.file)) > 0){
    vars.filt = readRDS(in.file)
  }else{
    colClasses = c(rep('character', 2), rep('integer', 4))
    vars.filt = read.table(in.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, colClasses = colClasses)
    colnames(vars.filt) = c('gene', 'snp', 'U.alt.dp', 'U.ref.dp', 'T.alt.dp', 'T.ref.dp')
    if(ncol(vars.filt) == 7){
      colnames(vars.filt)[7] = 'class'
    }
  }

  return(vars.filt)
}

get.main.pval <- function(obs.df, method, n.samples, nmax.vars, betabin.rho, gene.null.nvar.list = NULL){

  ###
  #Obs dist
  ###

  #Adds dp cols (dp = alt.dp + ref.dp) if dp cols doesnt already exist (U.dp, T.dp)
  obs.df = add.dp.col(obs.df)

  #Add p.sample.null column
  obs.df = add.pnull.col(obs.df, p.sample.null = 0.5)
  
  #Get snv testval (z)
  obs.df = get.snv.testval(obs.df, method)

  #Add column: n.vars
  obs.df = add.nvar.column(obs.df)
  
  #Calc gene based measures for each gene based on the snv.testval (median, mean, cv, Stouffer)
  obs.gene.df = get.gene.testval(obs.df) 
 
  #Add column 'nvars.bin'
  obs.gene.df = add.nvars.bin(obs.gene.df, nmax.vars)
  
  
  ###
  #Null dist by sampling
  ###  
  if(is.null(gene.null.nvar.list)){
    liptak.s.p.nom = get.binresamp.pval(obs.df, obs.gene.df, n.samples, method, betabin.rho)
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

get.snv.testval <- function(data.df, method, u.alt.dp.col = 'U.alt.dp', u.ref.dp.col = 'U.ref.dp', t.alt.dp.col = 'T.alt.dp', t.ref.dp.col = 'T.ref.dp', eps = 1e-16){

  dp.cols = c(u.alt.dp.col, u.ref.dp.col, t.alt.dp.col, t.ref.dp.col)
    
  #generic col names
  gen.dp.cols = c('U.alt.dp', 'U.ref.dp', 'T.alt.dp', 'T.ref.dp')
  
  #add pseudo-count if alt.dp or ref.dp == 0
  pseudo.data.df = t(apply(data.df[, dp.cols], 1, function(var.counts){if(any(var.counts == 0)){ var.counts = var.counts + 1;}; return(var.counts);}))
  colnames(pseudo.data.df) = gen.dp.cols
  
  if(method == 'or.dp'){
    or = (pseudo.data.df[, 'U.alt.dp'] / pseudo.data.df[, 'U.ref.dp']) / (pseudo.data.df[, 'T.alt.dp'] / pseudo.data.df[, 'T.ref.dp']) #since fc = p/1-p = odds

    #abs and log of OR. #NB: abs, since phase unknown
    abs.log.or = abs(log(or))  

    #SE from total count
    log.or.se = sqrt(apply(1 / pseudo.data.df[, gen.dp.cols], 1, sum));
    
    #half Z (half-normal, since abs(log(OR)))
    z = abs.log.or / log.or.se
  }
  if(method == 'or.fisher'){
    ft.res = calc.or(pseudo.data.df, gen.dp.cols)
    or = ft.res[, 'or']
    or.ci.lower = ft.res[, 'or.ci.lower']
    or.ci.upper = ft.res[, 'or.ci.upper']

    #abs and log of OR. #NB: abs, since phase unknown
    abs.log.or = abs(log(or))
  
    #SE (standard error) of log(OR)
    log.or.se = (log(or.ci.upper) - log(or.ci.lower)) / (2 * 1.96)

    #half Z (half-normal, since abs(log(OR)))
    z = abs.log.or / log.or.se
  }
  if(method == 'p.fisher'){
    ft.res = calc.or(pseudo.data.df, gen.dp.cols)
    pval = ft.res[, 'pval']
    or = ft.res[, 'or']
    or.ci.lower = ft.res[, 'or.ci.lower']
    or.ci.upper = ft.res[, 'or.ci.upper']

    #abs and log of OR. #NB: abs, since phase unknown
    abs.log.or = abs(log(or))
  
    #SE (standard error) of log(OR)
    log.or.se = (log(or.ci.upper) - log(or.ci.lower)) / (2 * 1.96)

    #Z score from inv.norm(1-pval)
    z = pval2zscore(pval, eps)
  }
  
  #bind all
  data.df = cbind(data.df, z, abs.log.or, log.or.se, or)

  return(data.df)
}

get.binresamp.pval <- function(obs.df, obs.gene.df, n.samples, method = 'or.fisher', betabin.rho = 0){
  
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
    gene.null = get.gene.null(obs.df.gene[, c('U.dp', 'T.dp')], n.samples, method, obs.df.gene[, 'p.sample.null'], betabin.rho) #vec: len: n.samples
      
    #Get p-val by comparing to null
    obs.val = obs.gene.df[which(obs.gene.df[, 'gene'] == jgene), 'liptak.s']
    p.nom[jgene] = length(which(gene.null >= obs.val)) / length(gene.null)
  }
  return(p.nom)
}

get.gene.null <- function(gene.snv.dps, n.samples, method = 'or.fisher', p.sample.null, betabin.rho = 0){
  
  #get null for each snv
  n.snvs = nrow(gene.snv.dps)
  snv.null.dist = matrix(nrow = n.samples, ncol = n.snvs)
  for(j.snv in 1:n.snvs){
    null.df = get.snv.null(gene.snv.dps[j.snv, ], n.samples, p.sample.null[j.snv], betabin.rho) #out: n.samples x (U.alt.dp, U.ref.dp, U.dp, T.alt.dp, T.ref.dp, T.dp)
    null.df = get.snv.testval(null.df, method)
    snv.null.dist[, j.snv] = null.df[, 'z']
  }

  #get gene testval by combining snv null testvals by liptak.stouffer
  gene.null = apply(snv.null.dist, 1, liptak.stouffer) #vec: len: n.samples

  return(gene.null)
}

get.snv.null <- function(j.dps, n.samples, p.sample.null = 0.5, betabin.rho = 0){
#sample betabinomial  

  if(is.character(betabin.rho)){
    betabin.rho = as.numeric(strsplit(betabin.rho, ',')[[1]])
    if(length(betabin.rho) > 1){
      betabin.rho.u = betabin.rho[1]
      betabin.rho.t = betabin.rho[2]
    }
    else{
      betabin.rho.u = betabin.rho
      betabin.rho.t = betabin.rho      
    }
  }
  else{
    betabin.rho.u = betabin.rho
    betabin.rho.t = betabin.rho
  }
  
  U.j.dp = as.integer(j.dps[1])
  T.j.dp = as.integer(j.dps[2])
  
  U.alt.dp = rbetabinom(n.samples, U.j.dp, p.sample.null, betabin.rho.u)
  U.dp = rep(U.j.dp, n.samples)
  U.ref.dp = U.dp - U.alt.dp
  
  T.alt.dp = rbetabinom(n.samples, T.j.dp, p.sample.null, betabin.rho.t)
  T.dp = rep(T.j.dp, n.samples)
  T.ref.dp = T.dp - T.alt.dp
  
  data.df = cbind(U.alt.dp, U.ref.dp, U.dp, T.alt.dp, T.ref.dp, T.dp)

  return(data.df)
}

calc.or <- function(data.df, dp.cols){
  ft.res = t(apply(data.df[, dp.cols], 1, function(var.counts){mat = matrix(var.counts, byrow=TRUE, nrow=2); res = fisher.test(mat); return(c(res$p.value, res$estimate, res$conf.int[1:2]));}))
  colnames(ft.res) = c('pval', 'or', 'or.ci.lower', 'or.ci.upper')
  return(ft.res)  
}

add.pnull.col <- function(data.df, p.sample.null = 0.5, u.alt.dp.col = 'U.alt.dp', u.site.dp.col = 'U.dp', t.alt.dp.col = 'T.alt.dp', t.site.dp.col = 'T.dp'){

  if(is.null(p.sample.null)){
    #add pseudo-count if alt.dp or ref.dp == 0
    dp.cols = c(u.alt.dp.col, u.site.dp.col, t.alt.dp.col, t.site.dp.col)
    pseudo.data.df = t(apply(data.df[, dp.cols], 1, function(var.counts){if(any(var.counts == 0)){ var.counts = var.counts + 1;}; return(var.counts);}))

    #estimate p.sample.null from observed counts
    p.sample.null = (pseudo.data.df[, u.alt.dp.col] + pseudo.data.df[, t.alt.dp.col]) / (pseudo.data.df[, u.site.dp.col] + pseudo.data.df[, t.site.dp.col])
  }else{
    p.sample.null = rep(p.sample.null, nrow(data.df))
  }

  #bind
  data.df = cbind(data.df, p.sample.null)

  return(data.df)  
}

add.dp.col <- function(data.df){

  cols = colnames(data.df)
  if(length(which(cols == 'U.dp')) == 0){
    U.dp = data.df[, 'U.alt.dp'] + data.df[, 'U.ref.dp']
    data.df = cbind(data.df, U.dp)
  }
  if(length(which(cols == 'T.dp')) == 0){
    T.dp = data.df[, 'T.alt.dp'] + data.df[, 'T.ref.dp']
    data.df = cbind(data.df, T.dp)
  }

  return(data.df)
}
