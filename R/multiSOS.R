###############################################################
# function to perform 2D & 3D analysis
###############################################################
multiSOS <- function(X, Y, family, fil.num, min.level, sizes, lam, maxiter, tol){
  if(missing(family)){family="DaubLeAsymm"}
  if(missing(fil.num)){fil.num=4}
  if(missing(min.level)){min.level=0}
  if(missing(sizes)){cat("Error: missing group sizes.","\n")}
  if(missing(lam)){cat("Error: missing tuning parameter lambda.","\n")}
  if(missing(maxiter)){maxiter=200}
  if(missing(tol)){tol=1e-8}

  ngroup = length(sizes)
  cumsize = c(0,cumsum(sizes))

  # wavelet transform on the design matrix X
  dims = dim(X[[1]])
  if(length(dims)==2){
    C = decomp.2D(X,family=family,fil.num=fil.num,min.level=min.level)
  }else if(length(dims)==3){
    C = decomp.3D(X,family=family,fil.num=fil.num,min.level=min.level)
  }
  c.norm = apply(C,2,function(a) sqrt(sum(a^2)))
  ind = which(c.norm>0)

  # implement the metalasso algorithm
  fit = metalasso(C[,ind],Y,obs=sizes,lam1=lam,maxit=maxiter,tol=tol)

  # inverse wavelet transform to obtain the original coefficients beta
  eta = list()
  for(m in 1:ngroup){
    eta[[m]] = rep(0,prod(dims))
    eta[[m]][ind] = fit$coe[[m]]
  }
  if(length(dims)==2){
    beta = lapply(eta,reconstr.2D,dims=dims,fil.num=fil.num,family=family,
                  min.level=min.level)
  }else if(length(dims)==3){
    beta = lapply(eta,reconstr.3D,dims=dims,family=family,fil.num=fil.num,
                  min.level=min.level)
  }

  # calculate estimate of the original intercepts beta0
  beta0=rep(NA,ngroup)
  for(i in 1:ngroup){
    Ym = Y[(cumsize[i]+1):cumsize[(i+1)]]
    Xm = X[(cumsize[i]+1):cumsize[(i+1)]]
    Ym.hat = sapply(Xm,function(a) sum(a*beta[[i]]))
    beta0[i] = mean(Ym-Ym.hat)
  }

  # return results
  return(list(beta = beta,
              beta0 = beta0,
              family = family,
              fil.num = fil.num,
              min.level = min.level,
              lam = lam,
              iteration = fit$iteration,
              converge = fit$converge
  ))
}


#############################################
# cross validation on multiSOS
#############################################
cv.multiSOS <- function(X,Y,family,fil.num,sizes,maxiter,tol,levels,lams,k){
  if(missing(family)){family="DaubLeAsymm"}
  if(missing(fil.num)){fil.num=4}
  if(missing(sizes)){cat("Error: missing group sizes.","\n")}
  if(missing(maxiter)){maxiter=200}
  if(missing(tol)){tol=1e-8}

  if(missing(levels)){cat("missing DWT levels")}
  if(missing(lams)){cat("missing shrinkage lambdas")}
  n.level = length(levels)
  n.lam = length(lams)
  if(missing(k)){k=5}

  ngroup = length(sizes)
  cumsize = c(0,cumsum(sizes))

  # split the data into k lists proportional to group sizes
  split = vector("list",k)
  batch = round(sizes/k)
  cumbatch = c(0,cumsum(batch))
  for(i in 1:ngroup){
    order = sample(1:sizes[i])
    for(j in 1:(k-1)){
      split[[j]] = c(split[[j]],(order[((j-1)*batch[i]+1):(j*batch[i])]+cumsize[i]))
    }
    split[[k]] = c(split[[k]],(order[-(1:((k-1)*batch[i]))]+cumsize[i]))
  }

  MSE = matrix(0,n.level,n.lam)
  for(a in 1:n.level){
    for(b in 1:n.lam){
      min.level = levels[a]
      lam = lams[b]

      message(cat("start iteration for min.level=",min.level,", lam=",lam))

      for(i in 1:k){
        train = sort(unlist(split[-i]))
        test = sort(unlist(split[[i]]))
        if(i<k){
          fit = multiSOS(X[train], Y[train], family=family, fil.num=fil.num,
                        min.level=min.level, sizes = (sizes-batch), lam=lam,
                        maxiter=maxiter, tol=tol)
          Ytest.hat=rep(NA,sum(batch))
          for(j in 1:ngroup){
            Ym.hat = sapply(X[test][(cumbatch[j]+1):cumbatch[(j+1)]],
                            function(c) sum(c*fit$beta[[j]])) + fit$beta0[j]
            Ytest.hat[(cumbatch[j]+1):cumbatch[(j+1)]] = Ym.hat
          }
          mse = sum((Y[test]-Ytest.hat)^2)


        }else{
          fit = multiSOS(X[train], Y[train], family=family, fil.num=fil.num,
                        min.level=min.level, sizes = ((k-1)*batch), lam=lam,
                        maxiter=maxiter, tol=tol)
          Ytest.hat=rep(NA,sum(sizes-(k-1)*batch))
          for(j in 1:ngroup){
            Ym.hat = sapply(X[test][((cumsize[j]-(k-1)*cumbatch[j])+1):(cumsize[(j+1)]-(k-1)*cumbatch[(j+1)])],
                            function(c) sum(c*fit$beta[[j]])) + fit$beta0[j]
            Ytest.hat[((cumsize[j]-(k-1)*cumbatch[j])+1):(cumsize[(j+1)]-(k-1)*cumbatch[(j+1)])] = Ym.hat
          }
          mse = sum((Y[test]-Ytest.hat)^2)
        }

        message(cat("fold ",i,", mse=",mse))
        MSE[a,b] = MSE[a,b] + mse
      }
      message(cat("min.level=",min.level,", lam=",lam,", MSE=",MSE[a,b],"\n"))
    }
  }

  row.names(MSE) = paste0("min.level=",levels)
  colnames(MSE) = paste0("lambda=",lams)

  return(list(MSE = MSE,
              best.level = levels[which(MSE==min(MSE),arr.ind=T)[1,1]],
              best.lam = lams[which(MSE==min(MSE),arr.ind=T)[1,2]],
              family = family,
              fil.num = fil.num
  ))
}
