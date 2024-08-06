###########################
#-- Auxiliary functions --#
###########################

split_elbo = function(y,Xsin,Xcos,L,cp,sigma2=0){
  
  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  elbo = rep(NA,n_cp+1)
  
  for (i in 1:(n_cp+1)) {
    segm = (pts[i]+1):pts[i+1]
    yW = y[segm]
    XsinW = matrix(Xsin[segm,],nrow=length(yW))
    XcosW = matrix(Xcos[segm,],nrow=length(yW))
    
    fitW = SuSiE_group(yW,XsinW,XcosW,L,sigma2=sigma2) 
    elbo[i] = fitW$elbo
  }
  
  sum(elbo)
}

split_elbo_MV = function(Y,Xsin,Xcos,L,cp,sigma2=0){
  n = nrow(Y)
  d = ncol(Y)
  
  elbo = 0
  for (j in 1:d) elbo = elbo + split_elbo(Y[,j],Xsin[,,j],Xcos[,,j],L,cp,sigma2=0)
  
  elbo
}

split_IC = function(y,Xsin,Xcos,L,cp,IC_name,sigma2=0){
  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  if (n_cp > 0) {
    pi = llik = rep(NA,n_cp+1)
    
    for (i in 1:(n_cp+1)) {
      segm = (pts[i]+1):pts[i+1]
      yW = y[segm]
      XsinW = matrix(Xsin[segm,],nrow=length(yW))
      XcosW = matrix(Xcos[segm,],nrow=length(yW))
      
      fitW = SuSiE_group(yW,XsinW,XcosW,L,sigma2=sigma2) 
      if (fitW$sigma2 <= 0) fitW$sigma2 = var(y)
      
      pi[i] =  2*sum(1-apply(fitW$alpha,1,function(x) prod(1-x))) + 1
      llik[i] = sum(dnorm(yW,fitW$yhat,sqrt(fitW$sigma2),log=T))
    }
    p = sum(pi)
    
  } else {
    segm = 1:n
    yW = y[segm]
    XsinW = matrix(Xsin[segm,],nrow=length(yW))
    XcosW = matrix(Xcos[segm,],nrow=length(yW))
    
    fitW = SuSiE_group(yW,XsinW,XcosW,L,sigma2=sigma2) 
    if (fitW$sigma2 <= 0) fitW$sigma2 = var(y)
    
    pi = 2*sum(1-apply(fitW$alpha,1,function(x) prod(1-x))) + 1
    p = pi
    llik = sum(dnorm(yW,fitW$yhat,sqrt(fitW$sigma2),log=T))
  }
  
  if (IC_name == 'MDL') IC = -sum(llik) + log(n_cp) + (n_cp+1)*log(n) + sum(log(pi)) + sum((pi/2+1)*log(diff(pts)))
  if (IC_name == 'BIC') IC = -sum(llik) + p/2*log(n) 
  if (IC_name == 'mBIC') IC = -sum(llik) + p/2*log(n) + 1/2*sum(pi*log(diff(pts/n)))
  
  IC
}

split_IC_MV = function(Y,Xsin,Xcos,L,cp,IC_name,sigma2=0){
  n = nrow(Y)
  d = ncol(Y)
  
  IC = 0
  for (j in 1:d) IC = IC + split_IC(Y[,j],Xsin[,,j],Xcos[,,j],L,cp,IC_name,sigma2=0)
  
  IC
}

#########################
#-- Optimistic search --#
#########################

nOS_internal = function(y,Xsin,Xcos,L,lt,t,rt,l,r,nu=1/2,sigma2=0){
  
  n = nrow(Xsin)
  
  if ((rt-lt) <= 5){
    elbo = c()
    if (lt == l) lt = l+1
    if (rt == r) rt = r-1
    for (i in lt:rt){
      elbo = c(elbo,split_elbo(y,Xsin,Xcos,L,cp=i,sigma2=sigma2))
    }
    cp = which.max(elbo)+lt-1
    out = list(cp=cp,
               elbo=max(elbo,na.rm=T))
    return(out)
  } else {
    
    if ((rt-t) > (t-lt)) {
      w = floor(rt - (rt-t)*nu)
      elbo_w = split_elbo(y,Xsin,Xcos,L,cp=w,sigma2=sigma2)
      elbo_t = split_elbo(y,Xsin,Xcos,L,cp=t,sigma2=sigma2)
      if (elbo_w > elbo_t) { nOS_internal(y,Xsin,Xcos,L,t,w,rt,l,r,nu,sigma2=sigma2) } else { nOS_internal(y,Xsin,Xcos,L,lt,t,w,l,r,nu,sigma2=sigma2) }
    } else {
      w = floor(lt + (t-lt)*nu)
      elbo_w = split_elbo(y,Xsin,Xcos,L,cp=w,sigma2=sigma2)
      elbo_t = split_elbo(y,Xsin,Xcos,L,cp=t,sigma2=sigma2)
      if (elbo_w > elbo_t) { nOS_internal(y,Xsin,Xcos,L,lt,w,t,l,r,nu,sigma2=sigma2) } else { nOS_internal(y,Xsin,Xcos,L,w,t,rt,l,r,nu,sigma2=sigma2) }
    }
    
  }
}

nOS = function(y,Xsin,Xcos,L,l,r,nu=1/2,thresh=0,sigma2=0){
  
  n = length(y)
  t = floor((l+nu*r/(1+nu)))
  out_nOS = nOS_internal(y,Xsin,Xcos,L=L,
                         lt=l,t=t,rt=r,l=l,r=r,nu=nu,sigma2=sigma2)
  
  cp_hat = NA
  
  fit0 = SuSiE_group(y,Xsin,Xcos,L,sigma2=sigma2)
  
  if (out_nOS$elbo - fit0$elbo > thresh)  cp_hat = out_nOS$cp
  
  c(cp_hat,out_nOS$elbo-fit0$elbo)
}


nOS_MV_internal = function(Y,Xsin,Xcos,L,lt,t,rt,l,r,nu=1/2,sigma2=0){
  
  n = dim(Xsin)[1]
  d = dim(Xsin)[3]
  
  if ((rt-lt) <= 5){
    elbo = c()
    if (lt == l) lt = l+1
    if (rt == r) rt = r-1
    for (i in lt:rt){
      elbo = c(elbo,split_elbo_MV(Y,Xsin,Xcos,L,cp=i,sigma2=sigma2))
    }
    cp = which.max(elbo)+lt-1
    out = list(cp=cp,
               elbo=max(elbo,na.rm=T))
    return(out)
  } else {
    
    if ((rt-t) > (t-lt)) {
      w = floor(rt - (rt-t)*nu)
      elbo_w = split_elbo_MV(Y,Xsin,Xcos,L,cp=w,sigma2=sigma2)
      elbo_t = split_elbo_MV(Y,Xsin,Xcos,L,cp=t,sigma2=sigma2)
      if (elbo_w > elbo_t) { nOS_MV_internal(Y,Xsin,Xcos,L,t,w,rt,l,r,nu,sigma2=sigma2) } else { nOS_MV_internal(Y,Xsin,Xcos,L,lt,t,w,l,r,nu,sigma2=sigma2) }
    } else {
      w = floor(lt + (t-lt)*nu)
      elbo_w = split_elbo_MV(Y,Xsin,Xcos,L,cp=w,sigma2=sigma2)
      elbo_t = split_elbo_MV(Y,Xsin,Xcos,L,cp=t,sigma2=sigma2)
      if (elbo_w > elbo_t) { nOS_MV_internal(Y,Xsin,Xcos,L,lt,w,t,l,r,nu,sigma2=sigma2) } else { nOS_MV_internal(Y,Xsin,Xcos,L,w,t,rt,l,r,nu,sigma2=sigma2) }
    }
    
  }
}


nOS_MV = function(Y,Xsin,Xcos,L,l,r,nu=1/2,thresh=0,sigma2=0){
  
  n = nrow(Y)
  d = ncol(Y)
  t = floor((l+nu*r/(1+nu)))
  out_nOS = nOS_MV_internal(Y,Xsin,Xcos,L=L,
                         lt=l,t=t,rt=r,l=l,r=r,nu=nu,sigma2=sigma2)
  cp_hat = NA
  
  elbo0 = rep(NA,d)
  for (j in 1:d) {
    fit0 = SuSiE_group(Y[,j],Xsin[,,j],Xcos[,,j],L,sigma2=sigma2)
    elbo0[j] = fit0$elbo
  }
  elbo0 = sum(elbo0)
  
  if (out_nOS$elbo - elbo0 > thresh)  cp_hat = out_nOS$cp
  
  c(cp_hat,out_nOS$elbo-elbo0)
}


getCandidateSetU = function(y,Xsin,Xcos,Ne,sigma2=0,nu=1/2,thresh=0,rate_c=0.1,Trace=0) {
  L = Ne
  n = length(y)
  w = rate_c*n
  
  it = 1
  if (Trace == 1) cat(paste0("* Searching candidate CP...\n"))
  eval_sgm = matrix(c(1,n),1,2)
  
  os = nOS(y,Xsin,Xcos,L,l=1,r=n,nu,thresh,sigma2=sigma2)
  cp_hat = os[1]
  elbo_hat = os[2]
  
  cp_set = NULL
  cp_list = NULL
  elbo_list = NULL
  if (!is.na(cp_hat)) {
    cp_set = cp_hat
    cp_list = list(cp_hat)
    elbo_list = list(elbo_hat)
    s_set = 1
  } else {
    s_set = 0
  }
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      it = it + 1
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      elbs = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        
        if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
          eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
          
          os = nOS(y[x],Xsin[x,],Xcos[x,],L,1,length(x),nu,thresh,sigma2=sigma2)
          cp_hat = os[1] + pts[k]
          elbo_hat = os[2]
          
          if (!is.na(cp_hat)) if(min(abs(cp_set-cp_hat)) > w) if (cp_hat-w > 1) if (cp_hat+w < n) {
            cps[k] = cp_hat 
            elbs[k] = elbo_hat
          } 
        }
      }
      
      elbs = elbs[!is.na(cps)]
      cps = cps[!is.na(cps)]
      cp_set = unique(sort(c(cp_set,cps)))
      if (length(cps) > 0) {
        cp_list[[it]] = cps[sort.int(elbs,decreasing=T,index.return=T)$ix]
        elbo_list[[it]] = elbs[sort.int(elbs,decreasing=T,index.return=T)$ix]
      }
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
      
    }
  }
  
  if (Trace == 1) cat(paste0("Candidate CP set: {",paste(cp_set,collapse=','),"}\n"))
  list(cp=cp_list,delbo=elbo_list)
}

pruneSetU = function(cp_set,y,Xsin,Xcos,Ne,sigma2=0,IC="MDL",option="joint",Trace=0,pl=FALSE,re_test=3,fixCP=0,maxEval=8) {
  L = Ne
  
  if (length(cp_set) == 0) out = NULL
  
  if (length(cp_set) > 0) {
    
    if (option == "joint") {
      n = length(y)
      
      ind_elbo = sort.int(unlist(cp_set$delbo),decreasing = TRUE,index.return = TRUE)$ix
      cp_set = unlist(cp_set$cp)[ind_elbo]
      
      if (length(cp_set) > maxEval+fixCP) cp_set = cp_set[1:(maxEval+fixCP)]
      
      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      if (length(cp_set) > 1) {
        cp_mat = as(expand.grid(rep(list(0:1),length(cp_set))),"matrix")
        cp_mat = cp_mat[-1,]
        if (fixCP > 0) cp_mat = cp_mat[apply(cp_mat[,1:fixCP],1,function(x) all(x==1)),]
        
        IC_star = rep(NA,nrow(cp_mat))
        for (i in 1:nrow(cp_mat)) {
          IC_star[i] = split_IC(y,Xsin,Xcos,L,cp=sort(cp_set[cp_mat[i,]==1]),IC,sigma2=sigma2)
        }
        
        out = sort(cp_set[cp_mat[which.min(IC_star),]==1])
      } else {
        out = sort(cp_set)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    if (option == 'conditional split') {
      n = length(y)
      
      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      cp_set = cp_set$cp
      
      lev = length(cp_set)
      cp_hat = cp_set[[1]]
      IC_hat = split_IC(y,Xsin,Xcos,L,cp=cp_hat,IC,sigma2=sigma2)
      
      if (lev > 1) {
        conv = 0
        lv = 1
        while (conv == 0) {
          lv = lv + 1
          cps = cp_set[[lv]]
          for (i in 1:length(cps)) {
            cp_star = unique(sort(c(cp_hat,cps[1:i])))
            IC_star = split_IC(y,Xsin,Xcos,L,cp=cp_star,IC,sigma2=sigma2)
            if (IC_star < IC_hat) {
              cp_hat = cp_star
              IC_hat = IC_star
            } else {
              conv = 1
              break
            }
          }
          if (lv == lev) conv = 1
        }
        
        out = sort(cp_hat)
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    
    if (option == 'conditional elbo') {
      n = length(y)
      
      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      ind_elbo = sort.int(unlist(cp_set$delbo),decreasing = TRUE,index.return = TRUE)$ix
      cp_set = unlist(cp_set$cp)[ind_elbo]
      
      ncp = length(cp_set)
      cp_hat = cp_set[1]
      IC_hat = split_IC(y,Xsin,Xcos,L,cp=cp_hat,IC,sigma2=sigma2)
      
      if (ncp > 1) {
        conv = 0
        lv = 1
        while (conv == 0) {
          lv = lv + 1
          cps = cp_set[lv]
          cp_star = unique(sort(c(cp_hat,cps)))
          IC_star = split_IC(y,Xsin,Xcos,L,cp=cp_star,IC,sigma2=sigma2)
          if (IC_star < IC_hat) {
            cp_hat = cp_star
            IC_hat = IC_star
          } else {
            conv = 1
            break
          }
          if (lv == ncp) conv = 1
        }
        
        out = sort(cp_hat)
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    
    if (option == 'best IC') {
      n = length(y)
      
      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      ind_elbo = sort.int(unlist(cp_set$delbo),decreasing = TRUE,index.return = TRUE)$ix
      cp_set = unlist(cp_set$cp)[ind_elbo]
      
      ncp = length(cp_set)
      IC_story = rep(NA,ncp)
      
      if (ncp > 1) {
        for (lv in 1:ncp) {
          cp_star = unique(sort(cp_set[1:lv]))
          IC_star = split_IC(y,Xsin,Xcos,L,cp=cp_star,IC,sigma2=sigma2)
          IC_story[lv] = IC_star
        }
        
        if (pl == TRUE) {
          library(ggplot2)
          pl = ggplot(data=data.frame(x=1:ncp,y=IC_story),aes(x=x,y=y)) + theme_bw() +
            geom_line(size=1) + geom_point(pch=21,fill="white",col="black",size=4) +
            ylab("IC") + xlab("CP") + scale_x_continuous(breaks=1:ncp,labels=cp_set) +
            theme(text=element_text(size=15))
          plot(pl)
        }
        cp_hat = cp_set[1:which.min(IC_story)]
        
        out = sort(cp_hat)
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    if (option == 'best IC refine') {
      n = length(y)
      
      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      ind_elbo = sort.int(unlist(cp_set$delbo),decreasing = TRUE,index.return = TRUE)$ix
      cp_set = unlist(cp_set$cp)[ind_elbo]
      
      ncp = length(cp_set)
      IC_story = rep(NA,ncp)
      
      if (ncp > 1) {
        for (lv in 1:ncp) {
          cp_star = unique(sort(cp_set[1:lv]))
          IC_star = split_IC(y,Xsin,Xcos,L,cp=cp_star,IC,sigma2=sigma2)
          IC_story[lv] = IC_star
        }
        
        if (pl == TRUE) {
          library(ggplot2)
          pl = ggplot(data=data.frame(x=1:ncp,y=IC_story),aes(x=x,y=y)) + theme_bw() +
            geom_line(size=1) + geom_point(pch=21,fill="white",col="black",size=4) +
            ylab("IC") + xlab("CP") + scale_x_continuous(breaks=1:ncp,labels=cp_set) +
            theme(text=element_text(size=15))
          plot(pl)
        }
        
        if (length(cp_set) > 5) {
          to_test = sort(cp_set[c(1,which(diff(IC_story) < 0) + 1)])
          re_test = min(re_test,length(to_test)-1)
          id_eval = sort.int(diff(to_test),index.return=TRUE)$ix[1:re_test]
          re_eval = unique(sort(to_test[c(id_eval,id_eval+1)]))
          cp_fix = setdiff(to_test,re_eval)
          cp_set = c(cp_fix,re_eval)
          
          cp_mat = as(expand.grid(rep(list(0:1),length(cp_set))),"matrix")
          cp_mat = cp_mat[-1,]
          if (length(cp_fix) > 0) cp_mat = cp_mat[apply(cp_mat[,1:length(cp_fix)],1,function(x) all(x==1)),]
          
          IC_star = rep(NA,nrow(cp_mat))
          for (i in 1:nrow(cp_mat)) {
            IC_star[i] = split_IC(y,Xsin,Xcos,L,cp=sort(cp_set[cp_mat[i,]==1]),IC,sigma2=sigma2)
          }
          
          out = sort(cp_set[cp_mat[which.min(IC_star),]==1])
        } else {
          out = cp_set[1:which.min(IC_story)]
        }
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
  }
  
  out
}

getCandidateSetMV = function(Y,Xsin,Xcos,Ne,sigma2=0,nu=1/2,thresh=0,rate_c=0.1,Trace=0) {
  L = Ne
  
  n = nrow(Y)
  d = ncol(Y)
  
  w = rate_c*n
  
  it = 1
  if (Trace == 1) cat(paste0("* Searching candidate CP...\n"))
  eval_sgm = matrix(c(1,n),1,2)
  
  os = nOS_MV(Y,Xsin,Xcos,L,l=1,r=n,nu,thresh,sigma2=sigma2)
  cp_hat = os[1]
  elbo_hat = os[2]
  
  cp_set = NULL
  cp_list = NULL
  elbo_list = NULL
  if (!is.na(cp_hat)) {
    cp_set = cp_hat
    cp_list = list(cp_hat)
    elbo_list = list(elbo_hat)
    s_set = 1
  } else {
    s_set = 0
  }
  
  if (length(cp_set) > 0) {
    conv = 0
    while (conv == 0) {
      it = it + 1
      
      pts = c(1,cp_set,n)
      n_segms = length(pts)-1
      cps = rep(NA,n_segms)
      elbs = rep(NA,n_segms)
      
      for (k in 1:n_segms) {
        x = pts[k]:pts[k+1]
        
        if (all(!(apply(eval_sgm,1,function(x) all(x==c(pts[k],pts[k+1])))))) {
          eval_sgm = rbind(eval_sgm,c(pts[k],pts[k+1]))
          
          os = nOS_MV(Y[x,],Xsin[x,,],Xcos[x,,],L,1,length(x),nu,thresh,sigma2=sigma2)
          cp_hat = os[1] + pts[k]
          elbo_hat = os[2]
          
          if (!is.na(cp_hat)) if(min(abs(cp_set-cp_hat)) > w) if (cp_hat-w > 1) if (cp_hat+w < n) {
            cps[k] = cp_hat 
            elbs[k] = elbo_hat
          } 
        }
      }
      
      elbs = elbs[!is.na(cps)]
      cps = cps[!is.na(cps)]
      cp_set = unique(sort(c(cp_set,cps)))
      if (length(cps) > 0) {
        cp_list[[it]] = cps[sort.int(elbs,decreasing=T,index.return=T)$ix]
        elbo_list[[it]] = elbs[sort.int(elbs,decreasing=T,index.return=T)$ix]
      }
      
      if (length(cp_set) > s_set) {
        s_set = length(cp_set)
      } else {
        conv = 1
      }
      
    }
  }
  
  if (Trace == 1) cat(paste0("Candidate CP set: {",paste(cp_set,collapse=','),"}\n"))
  list(cp=cp_list,delbo=elbo_list)
}

pruneSetMV = function(cp_set,Y,Xsin,Xcos,Ne,sigma2=0,IC="MDL",option="joint",Trace=0,pl=FALSE,re_test=3,fixCP=0,maxEval=8) {
  L = Ne
  
  n = nrow(Y)
  d = ncol(Y)
  
  if (length(cp_set) == 0) out = NULL
  
  if (length(cp_set) > 0) {
    
    if (option == "joint") {
      ind_elbo = sort.int(unlist(cp_set$delbo),decreasing = TRUE,index.return = TRUE)$ix
      cp_set = unlist(cp_set$cp)[ind_elbo]
      
      if (length(cp_set) > maxEval+fixCP) cp_set = cp_set[1:(maxEval+fixCP)]
      
      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      if (length(cp_set) > 1) {
        cp_mat = as(expand.grid(rep(list(0:1),length(cp_set))),"matrix")
        cp_mat = cp_mat[-1,]
        if (fixCP > 0) cp_mat = cp_mat[apply(cp_mat[,1:fixCP],1,function(x) all(x==1)),]
        
        IC_star = rep(NA,nrow(cp_mat))
        for (i in 1:nrow(cp_mat)) {
          IC_star[i] = split_IC_MV(Y,Xsin,Xcos,L,cp=sort(cp_set[cp_mat[i,]==1]),IC,sigma2=sigma2)
        }
        
        out = sort(cp_set[cp_mat[which.min(IC_star),]==1])
      } else {
        out = sort(cp_set)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    if (option == 'conditional split') {

      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      cp_set = cp_set$cp
      
      lev = length(cp_set)
      cp_hat = cp_set[[1]]
      IC_hat = split_IC_MV(Y,Xsin,Xcos,L,cp=cp_hat,IC,sigma2=sigma2)
      
      if (lev > 1) {
        conv = 0
        lv = 1
        while (conv == 0) {
          lv = lv + 1
          cps = cp_set[[lv]]
          for (i in 1:length(cps)) {
            cp_star = unique(sort(c(cp_hat,cps[1:i])))
            IC_star = split_IC_MV(Y,Xsin,Xcos,L,cp=cp_star,IC,sigma2=sigma2)
            if (IC_star < IC_hat) {
              cp_hat = cp_star
              IC_hat = IC_star
            } else {
              conv = 1
              break
            }
          }
          if (lv == lev) conv = 1
        }
        
        out = sort(cp_hat)
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    
    if (option == 'conditional elbo') {

      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      ind_elbo = sort.int(unlist(cp_set$delbo),decreasing = TRUE,index.return = TRUE)$ix
      cp_set = unlist(cp_set$cp)[ind_elbo]
      
      ncp = length(cp_set)
      cp_hat = cp_set[1]
      IC_hat = split_IC_MV(Y,Xsin,Xcos,L,cp=cp_hat,IC,sigma2=sigma2)
      
      if (ncp > 1) {
        conv = 0
        lv = 1
        while (conv == 0) {
          lv = lv + 1
          cps = cp_set[lv]
          cp_star = unique(sort(c(cp_hat,cps)))
          IC_star = split_IC_MV(Y,Xsin,Xcos,L,cp=cp_star,IC,sigma2=sigma2)
          if (IC_star < IC_hat) {
            cp_hat = cp_star
            IC_hat = IC_star
          } else {
            conv = 1
            break
          }
          if (lv == ncp) conv = 1
        }
        
        out = sort(cp_hat)
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    
    if (option == 'best IC') {

      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      ind_elbo = sort.int(unlist(cp_set$delbo),decreasing = TRUE,index.return = TRUE)$ix
      cp_set = unlist(cp_set$cp)[ind_elbo]
      
      ncp = length(cp_set)
      IC_story = rep(NA,ncp)
      
      if (ncp > 1) {
        for (lv in 1:ncp) {
          cp_star = unique(sort(cp_set[1:lv]))
          IC_star = split_IC_MV(Y,Xsin,Xcos,L,cp=cp_star,IC,sigma2=sigma2)
          IC_story[lv] = IC_star
        }
        
        if (pl == TRUE) {
          library(ggplot2)
          pl = ggplot(data=data.frame(x=1:ncp,y=IC_story),aes(x=x,y=y)) + theme_bw() +
            geom_line(size=1) + geom_point(pch=21,fill="white",col="black",size=4) +
            ylab("IC") + xlab("CP") + scale_x_continuous(breaks=1:ncp,labels=cp_set) +
            theme(text=element_text(size=15))
          plot(pl)
        }
        cp_hat = cp_set[1:which.min(IC_story)]
        
        out = sort(cp_hat)
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
    if (option == 'best IC refine') {

      if (Trace == 1) cat(paste0("* Pruning the CP set based on IC...\n"))
      
      ind_elbo = sort.int(unlist(cp_set$delbo),decreasing = TRUE,index.return = TRUE)$ix
      cp_set = unlist(cp_set$cp)[ind_elbo]
      
      ncp = length(cp_set)
      IC_story = rep(NA,ncp)
      
      if (ncp > 1) {
        for (lv in 1:ncp) {
          cp_star = unique(sort(cp_set[1:lv]))
          IC_star = split_IC_MV(Y,Xsin,Xcos,L,cp=cp_star,IC,sigma2=sigma2)
          IC_story[lv] = IC_star
        }
        
        if (pl == TRUE) {
          library(ggplot2)
          pl = ggplot(data=data.frame(x=1:ncp,y=IC_story),aes(x=x,y=y)) + theme_bw() +
            geom_line(size=1) + geom_point(pch=21,fill="white",col="black",size=4) +
            ylab("IC") + xlab("CP") + scale_x_continuous(breaks=1:ncp,labels=cp_set) +
            theme(text=element_text(size=15))
          plot(pl)
        }
        
        if (length(cp_set) > 5) {
          to_test = sort(cp_set[c(1,which(diff(IC_story) < 0) + 1)])
          re_test = min(re_test,length(to_test)-1)
          id_eval = sort.int(diff(to_test),index.return=TRUE)$ix[1:re_test]
          re_eval = unique(sort(to_test[c(id_eval,id_eval+1)]))
          cp_fix = setdiff(to_test,re_eval)
          cp_set = c(cp_fix,re_eval)
          
          cp_mat = as(expand.grid(rep(list(0:1),length(cp_set))),"matrix")
          cp_mat = cp_mat[-1,]
          if (length(cp_fix) > 0) cp_mat = cp_mat[apply(cp_mat[,1:length(cp_fix)],1,function(x) all(x==1)),]
          
          IC_star = rep(NA,nrow(cp_mat))
          for (i in 1:nrow(cp_mat)) {
            IC_star[i] = split_IC_MV(Y,Xsin,Xcos,L,cp=sort(cp_set[cp_mat[i,]==1]),IC,sigma2=sigma2)
          }
          
          out = sort(cp_set[cp_mat[which.min(IC_star),]==1])
        } else {
          out = cp_set[1:which.min(IC_story)]
        }
      } else {
        out = sort(cp_hat)
      }
      if (Trace == 1) cat(paste0("CP set: {",paste(out,collapse=','),"}\n"))
    }
    
  }
  
  out
}


#########################
#-- Plotting function --#
#########################

plotSearch = function(cp_set,n) {
  library(ggplot2)
  nlev = length(cp_set$cp)
  w = n/14
  h = nlev/15
  
  if (nlev > 1) {
    xt = cp_set$cp[[1]]
    yt = nlev
    tx = paste(cp_set$cp[[1]],"\n elbo:",round(cp_set$delbo[[1]],2))
    for (lv in 2:nlev) {
      xt = c(xt,cp_set$cp[[lv]])
      yt = c(yt,rep(nlev-lv+1,length(cp_set$cp[[lv]])))
      tx = c(tx,paste(cp_set$cp[[lv]],"\n elbo:",round(cp_set$delbo[[lv]],2)))
    }
    
    
    pl = ggplot() + theme_bw() +
      geom_point(data=data.frame(x=c(1,n),y=c(0,nlev)),aes(x=x,y=y),col="transparent") +
      geom_segment(aes(x = xt, xend = xt,
                       y = rep(0,length(xt)), yend = yt-h), linetype='dashed') +
      geom_rect(aes(xmin = xt - w, xmax = xt + w,
                    ymin = yt - h, ymax = yt + h), col = "black", fill = "white") +
      geom_text(aes(x=xt,y=yt,label=tx)) +
      geom_point(data=data.frame(x=unlist(cp_set$cp),y=c(rep(0,length(unlist(cp_set$cp))))),
                 aes(x=x,y=y),col="black",fill='red',pch=21,size=3) +
      geom_hline(yintercept=seq(0.5,nlev),col="grey",linetype='dashed') +
      xlab("") + ylab("") +
      scale_y_continuous(breaks=1:nlev,labels=nlev:1) +
      theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5)) +
      ggtitle("Candidate set search")
  }
  
  if (nlev == 1) {
    xt = cp_set$cp[[1]]
    yt = nlev
    tx = paste(cp_set$cp[[1]],"\n elbo:",round(cp_set$delbo[[1]],2))
    
    pl = ggplot() + theme_bw() +
      geom_point(data=data.frame(x=c(1,n),y=c(0,nlev)),aes(x=x,y=y),col="transparent") +
      geom_segment(aes(x = xt, xend = xt,
                       y = rep(0,length(xt)), yend = yt-h), linetype='dashed') +
      geom_rect(aes(xmin = xt - w, xmax = xt + w,
                    ymin = yt - h, ymax = yt + h), col = "black", fill = "white") +
      geom_text(aes(x=xt,y=yt,label=tx)) +
      geom_point(data=data.frame(x=unlist(cp_set$cp),y=c(rep(0,length(unlist(cp_set$cp))))),
                 aes(x=x,y=y),col="black",fill='red',pch=21,size=3) +
      geom_hline(yintercept=seq(0.5,nlev),col="grey",linetype='dashed') +
      xlab("") + ylab("") +
      scale_y_continuous(breaks=1:nlev,labels=nlev:1) +
      theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5)) +
      ggtitle("Candidate set search")
  }
  
  pl
}

plotSpectra = function(y,L=1,cp=NULL,time=NULL,n_om=100,maxOm=0.5) {
  n = length(y)
  if (is.null(time)) time = 1:n
  knots = c(time[1],time[cp],time[n])
  
  om_grid = seq(0,maxOm,length.out=n_om)
  x = seq(1,n)
  Xout = getX(x,om_grid)
  Xsin = Xout$Xsin
  Xcos = Xout$Xcos
  
  ft = fit_segments(y,Xsin,Xcos,L=L,cp=cp,pl=FALSE)
  
  require(ggplot2)
  require(reshape)
  
  pl = ggplot() + theme_bw()
  
  for (sgm in 1:length(ft$models)) {
    xd = seq(knots[sgm],knots[sgm+1],length.out=100)
    zd = matrix(0,n_om,100)
    for (i in 1:ncol(zd)) zd[,i] = abs(ft$models[[sgm]]$B1)
    
    Zd = melt(t(zd))
    Zd$X1 = rep(xd,n_om)
    Zd$X2 = rep(om_grid,each=100)
    
    pl = pl + stat_contour_filled(data=Zd,aes(x=X1,y=X2,z=value))
    
  }
  
  pl = pl +
    geom_vline(xintercept=time[cp],linewidth=1,col="blue",linetype='dashed') +
    labs(colour="Amplitude") + xlab("Time") + ylab(expression(omega)) +
    theme(text=element_text(size=15), 
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position='none') + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(limits=c(-0.01,maxOm),expand=c(0,0))
  
  pl
}

