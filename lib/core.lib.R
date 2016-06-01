
get.gene.testval <- function(data.df, testval.col = 'z'){
#Calc gene based measures for each gene based on the snv.testval (median, mean, cv, Stouffer)
  
  mean.s = tapply(data.df[, testval.col], data.df[, 'gene'], mean);
  median.s = tapply(data.df[, testval.col], data.df[, 'gene'], median);
  sd.s = tapply(data.df[, testval.col], data.df[, 'gene'], sd);
  cv.s = sd.s / mean.s;
  cv.s[which(sd.s == 0)] = 0 #To avoid NAs if mean also == 0 (0/0)
  
  #Liptak-Stouffer
  gene.list = split(data.df, data.df[, 'gene'])  
  liptak.s = unlist(lapply(gene.list, function(jgene.data){liptak.stouffer(jgene.data[, testval.col])})) #equal weighing 
    
  #bind all values into matrix
  gene.vals = cbind(mean.s, median.s, sd.s, cv.s, liptak.s);

  #add col: n.vars
  gene2nvars = table(data.df[, 'gene'])
  gene = names(gene2nvars)
  gene2nvars = cbind(as.data.frame(as.matrix(gene2nvars)), gene, stringsAsFactors = FALSE)
  colnames(gene2nvars)[1] = 'n.vars'
  gene.vals = merge(gene2nvars, gene.vals, by.x = 'gene', by.y = 'row.names')
  rownames(gene.vals) = gene.vals[, 'gene']
  
  return(gene.vals)  
}

liptak.stouffer <- function(beta, w = NA, na.rm = FALSE){

  if(na.rm){
    beta = beta[which(!is.na(beta))]
  }
  if(is.na(w)){
    w = rep(1, length(beta))
  }
  z = sum(w * beta) / sqrt(sum(w^2))
  return(z)
}

pval2zscore <- function(p, eps = 1e-16){
#qnorm is based on Wichura's algorithm AS 241 which provides precise results up to about 16 digits
  p[which(p > (1-eps))] = 1 - eps; 
  p[which(p < eps)] = eps;
  z = qnorm(1 - p)

  return(z)
}

fisher.combo <- function(p){
  Xsq = -2*sum(log(p))
  combo.p.nom = 1 - pchisq(Xsq, df = 2*length(p))
  return(combo.p.nom)
}

add.nvars.bin <- function(vars.data, nmax.vars = 11){

  nvars.bin = vars.data[, 'n.vars'];
  nvars.bin[which(nvars.bin >= nmax.vars)] = nmax.vars;
  vars.data = cbind(vars.data, nvars.bin)

  return(vars.data)
}

add.nvar.column <- function(data.df){

  #n.vars
  cols = colnames(data.df)
  if(length(grep('^n.vars$', cols)) == 0){
    gene2nvars = table(data.df[, 'gene'])
    gene = names(gene2nvars)
    gene2nvars = cbind(as.data.frame(as.matrix(gene2nvars)), gene, stringsAsFactors = FALSE)
    colnames(gene2nvars)[1] = 'n.vars'
    data.df = merge(data.df, gene2nvars, by = 'gene')    
  }

  return(data.df)
}

min.vars.filter <- function(vars.filt, min.het.vars){
  gene2nvars = table(vars.filt[, 'gene']);
  genes.pass = names(gene2nvars)[which(gene2nvars >= min.het.vars)];
  vars.filt = merge(vars.filt, as.matrix(genes.pass), by.x = 'gene', by.y = 1)

  return(vars.filt)
}

dump.tab <- function(pval.df, pval.file, sim = FALSE){
  if(length(grep('\\.rds$', pval.file)) > 0){
    saveRDS(pval.df, file = pval.file)
  }else{
    colnames(pval.df) = sub('^gene$', 'feat', colnames(pval.df))
    colnames(pval.df) = sub('^liptak.s.p.nom$', 'p.nom', colnames(pval.df))
    colnames(pval.df) = sub('^liptak.s.p.nom.bh$', 'fdr', colnames(pval.df))
    dump.cols = c('feat', 'n.vars', 'mean.s', 'median.s', 'sd.s', 'cv.s', 'liptak.s', 'p.nom', 'fdr')
    if(sim){
      dump.cols = c(dump.cols, 'class')
    }
    write.table(pval.df[, dump.cols], quote = FALSE, row.names = FALSE, sep = '\t', file = pval.file)
  }
}

get.tpr.fpr <- function(obs.mat, alpha.vec, measure, feat.col = 'gene'){

  #classify based on measure
  gene2sig.alpha = matrix(nrow = nrow(obs.mat), ncol = length(alpha.vec), dimnames = list(obs.mat[, feat.col], alpha.vec))
  for(jalpha in alpha.vec){
    gene2sig.alpha[, as.character(jalpha)] = get.sig(obs.mat, jalpha, measure)    
  }
  
  #Merge with gene class annot
  gene.annot = unique(obs.mat[, c(feat.col, 'class')])
  gene2sig.alpha = merge(gene.annot, gene2sig.alpha, by.x = feat.col, by.y = 'row.names')  

  #Get fpr and sens
  fpr.sens = t(apply(gene2sig.alpha[, as.character(alpha.vec)], 2, function(call, truth){fpr = length(which(call & (truth == 'neg'))) / length(which(truth == 'neg')); sens = length(which(call & (truth == 'pos'))) / length(which(truth == 'pos')); return(c(fpr, sens));}, truth = gene2sig.alpha[, 'class']))
  colnames(fpr.sens) = c('fpr', 'tpr')
  
  return(fpr.sens)
}

get.mcc <- function(obs.mat, alpha.vec, measure, feat.col = 'gene'){

  #classify based on measure
  gene2sig.alpha = matrix(nrow = nrow(obs.mat), ncol = length(alpha.vec), dimnames = list(obs.mat[, feat.col], alpha.vec))
  for(jalpha in alpha.vec){
    gene2sig.alpha[, as.character(jalpha)] = get.sig(obs.mat, jalpha, measure)    
  }
  
  #Merge with gene class annot
  gene.annot = unique(obs.mat[, c(feat.col, 'class')])
  gene2sig.alpha = merge(gene.annot, gene2sig.alpha, by.x = feat.col, by.y = 'row.names')  

  #Get Matthews correlation coef
  mcc = t(apply(gene2sig.alpha[, as.character(alpha.vec)], 2, function(call, truth){fp = length(which(call & (truth == 'neg'))); tp = length(which(call & (truth == 'pos'))); fn = length(which(!call & (truth == 'pos'))); tn = length(which(!call & (truth == 'neg'))); denom = (as.numeric(tp) + as.numeric(fp)) * (as.numeric(tp) + as.numeric(fn)) * (as.numeric(tn) + as.numeric(fp)) * (as.numeric(tn) + as.numeric(fn)); if(denom == 0){denom = 1}; matt = ((as.numeric(tp) * as.numeric(tn)) - (as.numeric(fp) * as.numeric(fn))) / sqrt(denom); return(c(matt, fp, tp, fn, tn));}, truth = gene2sig.alpha[, 'class']))
  colnames(mcc) = c('mcc', 'fp', 'tp', 'fn', 'tn')
  
  return(mcc)
}


get.sig <- function(obs.mat, alpha, pval.cols){
  n.tests = nrow(obs.mat)
  sig.res = matrix(nrow = n.tests, ncol = length(pval.cols))
  colnames(sig.res) = pval.cols
  for(jcol in pval.cols){
    sig.col = rep(FALSE, n.tests)
    sig.col[which(obs.mat[, jcol] <= alpha)] = TRUE
    sig.res[, jcol] = sig.col
  }                                                                                                  
  
  #get genes which are sig in both mean and sd
  multisig = apply(sig.res[, pval.cols, drop = FALSE], 1, function(sig.vals){all(as.logical(sig.vals))});

  return(multisig)
}
