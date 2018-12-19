################functions###################

## prognosis related functions
progPathway <- function(os.surv,rs.surv,ctf1=10e-3,ctf2=10e-3){
  if(is.null(os.surv) && is.null(rs.surv)) stop("No survival information")
  else if(is.null(os.surv)){
    sig_overall <- get_pz(os.surv,kegg,ctf1)
    sig_path <- sig_overall[[1]]
    direct_path <- as.numeric(sig_overall[[2]]>0)*2-1 
  }
  else if(is.null(rs.surv)){
    sig_recu <- get_pz(os.surv,kegg,ctf2)
    sig_path <- sig_recu[[1]]
    direct_path <- as.numeric(sig_recu[[2]]>0)*2-1
  }
  else{
    sig_overall <- get_pz(os.surv,kegg,ctf1)
    sig_overall_p <- sig_overall[[1]]
    sig_overall_w <- as.numeric(sig_overall[[2]]>0)*2-1
    
    sig_recu <- get_pz(os.surv,kegg,ctf2)
    sig_recu_p <- sig_recu[[1]]
    sig_recu_w <- as.numeric(sig_recu[[2]]>0)*2-1
    #intersection
    sig_path <- intersect(names(sig_overall_p),names(sig_recu_p)) 
    a <- sig_overall_w[which(names(sig_overall_p) %in% sig_path)]
    b <- sig_recu_w[which(names(sig_recu_p) %in% sig_path)]
    #pick direction consistent pathways
    direct_path <- as.numeric(a[which(a==b)]>0)*2-1
    sig_path <- sig_path[which(a==b)]
  }
  return(list(sig_path,direct_path))
}

progGene <- function(os.surv,rs.surv,ctf1=10e-7,ctf2=10e-7){
  if(is.null(os.surv) && is.null(rs.surv)){
    overall_sig <- geneValue(ex,os.surv)
    ov_p <- overall_sig[[1]]
    ov_z <- overall_sig[[2]]
    ov_p2 <- ov_p[!duplicated(names(ov_p))]
    sig_gene <- ov_p2[which(p.adjust(ov_p2)<ctf1)]
    gene_direct <- ov_z[!duplicated(names(ov_p))][which(p.adjust(ov_p2)<ctf1)]
  }
  else if(is.null(rs.surv)){
    recu_sig <- geneValue(ex,os.surv)
    rc_p <- recu_sig[[1]]
    rc-z <- recu_sig[[2]]
    rc_p2 <- rc_p[!duplicated(names(rc_p))]
    sig_gene <- rc_p2[which(p.adjust(rc_p2)<ctf2)] 
    gene_direct <- rc_z[!duplicated(names(rc_p))][which(p.adjust(rc_p2)<ctf1)]
  }
  else{
    overall_sig <- geneValue(ex,os.surv)
    ov_p <- overall_sig[[1]]
    ov_z <- overall_sig[[2]]
    ov_p2 <- ov_p[!duplicated(names(ov_p))]
    sig_op <- ov_p2[which(p.adjust(ov_p2)<ctf1)]
    sig_oz <- ov_z[!duplicated(names(ov_p))][which(p.adjust(ov_p2)<ctf1)]
    
    recu_sig <- geneValue(ex,os.surv)
    rc_p <- recu_sig[[1]]
    rc_z <- recu_sig[[2]]
    rc_p2 <- rc_p[!duplicated(names(rc_p))]
    sig_rp <- rc_p2[which(p.adjust(rc_p2)<ctf2)] 
    sig_rz <- rc_z[!duplicated(names(rc_p))][which(p.adjust(rc_p2)<ctf1)]
    #intersection
    sig_gene <- intersect(names(sig_op),names(sig_rp))
    #pick direction consistent genes
    direct_gene <- as.numeric(sig_oz[which(sig_oz==sig_rz)]>0)*2-1
    sig_gene <- sig_gene[which(sig_oz==sig_rz)]
  }
  return(list(sig_gene,direct_gene))
}

pathwayValue <- function(sur_info,ex_info,ctf){
  res <- apply(ex_info,1,function(x)summary(coxph(sur_info~x)))
  pvalue <- lapply(res,function(x)x[[10]][3])
  zscore <- lapply(res,function(x)x[[7]][4])
  p <- as.numeric(pvalue)
  names(p)<-names(pvalue)
  z <- as.numeric(zscore)
  names(z) <- names(zscore)
  adj_p <- p.adjust(p)
  sig_p <- adj_p[which(adj_p<ctf)]
  sig_z <- z[which(adj_p<ctf)]
  return(list(sig_p,sig_z))
}

geneValue <- function(gene_info,dirc_info){
  all_p <- numeric()
  all_z <- numeric()
  for(i in 1:length(gene_info)){
    p2symbol <- kegg_path$symbol[which(kegg_path$path_id %in% gene_info[i])]
    g2probe <- out$PROBEID[which(out$SYMBOL %in% p2symbol)]
    pos <- which(row.names(ex) %in% g2probe)
    e <- ex[pos,]
    
    sig_gene <- numeric()
    for(j in 1:length(p2symbol)){
      tmp <- out$PROBEID[which(out$SYMBOL %in% p2symbol[j])]
      if(length(tmp)==1) sig_gene <- rbind(sig_gene,e[tmp,])
      else sig_gene <- rbind(sig_gene,apply(e[tmp,],2,function(x)exp(mean(log(x)))))
    }
    row.names(sig_gene) <- p2symbol
    sig_gene2 <- na.omit(sig_gene)
    res <- apply(sig_gene2,1,function(x)summary(coxph(surtmp~x)))
    pvalue <- lapply(res,function(x)x[[10]][3])
    zscore <- lapply(res,function(x)x[[7]][4])
    p <- as.numeric(pvalue) 
    names(p) <- row.names(sig_gene2)
    z <- as.numeric(zscore) 
    pos2 = which(((as.numeric(z>0))*2-1)==dirc_info[i])
    all_p <- c(all_p,p[pos2])
    all_z <- c(all_z,z[pos2])
  }
  return(list(all_p,all_z))
}
checkOrder <- function(ex,pheno){
  if(colnames(ex)!=pheno$acc)
    stop("The order of sample information must be consistent,function orderDeal is offered to deal with one of these conditions.") 
}
orderDeal <- function(ex,ex.kegg,pheno){ 
  truncate.sample.id <- function(s)
  {
    a = unlist(strsplit(s, "_"))[1]
    a
  } 
  accession = apply(as.matrix(unlist(strsplit(colnames(ex), ".CEL"))), 1, truncate.sample.id)
  pheno_acc = pheno$acc
  Posi <- match(pheno_acc,accession)
  colnames(ex) <- accession
  ex <- ex[,Posi]
  colnames(kegg) <- accession
  kegg <- kegg[,Posi]
}

## diagnosis related functions
diagPathway(adc_pos,ctf){ 
  dia_p <- sapply(1:nrow(kegg),function(x)t.test(kegg[x,adc_pos],kegg[x,-adc_pos],paired = F))
  pt <- as.numeric(dia_p[1,])
  pp <- as.numeric(dia_p[3,])
  adj_pp <- p.adjust(pp)
  sig_path <- row.names(kegg)[which(adj_pp<ctf)]
  direct_path  <- as.numeric(pt[which(adj_pp<ctf)]>0)  
  return(sig_path,direct_path)
}


checkDiagData <- function(ex){ 
  ni2 <- numeric()
  for(i in 1:nrow(ex[ni,])){
    if(length(unique(as.numeric(ex[ni[i],])))==1) ni2 <- c(ni2,ni[i])
  }
  ex2 <- ex[-ni2,]
}
diagGene <- function(adc_pos,ctf){
  p2g <- numeric()
  overall_p <- numeric()
  overall_t <- numeric()
  for(i in 1:length(sig_path)){
    sig_path <- diagPathway[[1]]
    direct_path <- diagPathway[[2]]
    pvsymbol <- kegg_path$symbol[which(kegg_path$path_id %in% sig_path[i])]
    g2probe <- out$PROBEID[which(out$SYMBOL %in% pvsymbol)]
    pos <- which(row.names(ex2) %in% g2probe)
    e <- ex2[pos,]
    row.names(e)<-row.names(ex2)[pos]
    sig_gene <- numeric()
    for(j in 1:length(pvsymbol)){
      tmp <- out$PROBEID[which(out$SYMBOL %in% pvsymbol[j])]
      if(length(tmp)==1) sig_gene <- rbind(sig_gene,e[tmp,])
      else sig_gene <- rbind(sig_gene,apply(e[tmp,],2,function(x)exp(mean(log(x)))))
    }
    row.names(sig_gene) <- pvsymbol
    sig_gene2 <- na.omit(sig_gene)
    res <- sapply(c(1:nrow(sig_gene2)),function(x)t.test(sig_gene2[x,adc_pos],sig_gene2[x,-adc_pos],paired = F))
    pvalue <- as.numeric(res[3,])
    names(pvalue) <- row.names(sig_gene2)
    tscore <- as.numeric(res[1,])
    
    pos2 = which(((as.numeric(tscore>0))*2-1)==direct_path[i])
    pp <- pvalue[pos2]
    tt <- tscore[pos2]
    overall_p <- c(overall_p,pp)
    overall_t <- c(overall_t,tt)
    p2g <- c(p2g,list(pp,tt))
  }
  
  pp2 <- overall_p[!duplicated(names(overall_p))]
  pt2 <- overall_t[!duplicated(names(overall_p))]
  sig_diag <- names(pp2[which(p.adjust(pp2)<ctf)]) 
  diag_w <- as.numeric(pt2[which(p.adjust(pp2)<ctf)]>0)*2-1
  
}
