library(wavethresh)

###############################################################################################
# 2D discreet wavelet transform with selected min.level
# built based on wavethresh::imwd and wavethresh::imwr
# dimensions of input image should be power of 2 or padded with 0 to make up to power of 2
###############################################################################################
decomp.2D <- function(x,family="DaubLeAsymm",fil.num=4,min.level=2){
  if((class(x)!="list")||(class(x[[1]])!="matrix")){
    cat("input should be a list of square matrices","\n")
  }
  nsample = length(x)
  dims = dim(x[[1]])
  upper.level = min(log2(dims))
  if(upper.level<=min.level){
    cat("invalid min.level, should be less than min(log2(dim(x[[1]]))).","\n")
  }
  dim1d = prod(dims)
  wdt = matrix(NA,nsample,dim1d)

  for(i in 1:nsample){
    input = x[[i]]
    wd0 = imwd(input,filter.number=fil.num,family=family)

    for(j in (upper.level-1):min.level){
      rep.sub = c(eval(parse(text=paste0("wd0$w",j,"L4"))),eval(parse(text=paste0("wd0$w",j,"L1"))),
                  eval(parse(text=paste0("wd0$w",j,"L2"))),eval(parse(text=paste0("wd0$w",j,"L3"))))
      wdt[i,1:(4^(j+1))] = rep.sub
    }

  }

  return(wdt)
}

#######################################################################################
# inverse 2D discrete wavelet transform with selected min.level
# inverse procedure to decomp.2D but for single input vector
#######################################################################################
reconstr.2D <- function(eta,dims,family="DaubLeAsymm",fil.num=4,min.level=2){
  if(length(eta)!=prod(dims)){cat("length of eta does not match the given dims.","\n")}
  upper.level = min(log2(dims))
  wr0 = imwd(matrix(rep(0,prod(dims)),dims[1],dims[2]),filter.number=fil.num,family=family,RetFather=F)

  if(min.level>1){
    sub.wd = imwd(matrix(eta[1:(4^min.level)],2^min.level,2^min.level),family=family,filter.number=fil.num)
    wr0$w0Lconstant = sub.wd$w0Lconstant
    for(j in 0:(min.level-1)){
      value = eval(parse(text=paste0("sub.wd$w",j,"L1")))
      eval(parse(text=paste0("wr0$w",j,"L1 = value")))
      value = eval(parse(text=paste0("sub.wd$w",j,"L2")))
      eval(parse(text=paste0("wr0$w",j,"L2 = value")))
      value = eval(parse(text=paste0("sub.wd$w",j,"L3")))
      eval(parse(text=paste0("wr0$w",j,"L3 = value")))
    }
    for(j in min.level:(upper.level-1)){
      size = 4^j
      eval(parse(text=paste0("wr0$w",j,"L1 = eta[(size+1):(2*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L2 = eta[(2*size+1):(3*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L3 = eta[(3*size+1):(4*size)]")))
    }
    wr = imwr(wr0)

  }else if(min.level==1){
    a = eta[1]; b = eta[2]; c = eta[3]; d = eta[4]
    wr0$w0Lconstant = wr0$w0L4 = (a+b+c+d)/2
    wr0$w0L1 = (a+b-c-d)/2
    wr0$w0L2 = (a+c-b-d)/2
    wr0$w0L3 = (a+d-b-c)/2
    for(j in 1:(upper.level-1)){
      size = 4^j
      eval(parse(text=paste0("wr0$w",j,"L1 = eta[(size+1):(2*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L2 = eta[(2*size+1):(3*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L3 = eta[(3*size+1):(4*size)]")))
    }
    wr = imwr(wr0)

  }else{
    for(j in 0:(upper.level-1)){
      size = 4^j
      eval(parse(text=paste0("wr0$w",j,"L1 = eta[(size+1):(2*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L2 = eta[(2*size+1):(3*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L3 = eta[(3*size+1):(4*size)]")))
    }
    wr = imwr(wr0)
  }

  return(wr)
}



###############################################################################################
# 3D discreet wavelet transform with selected min.level
# built based on wavethresh::wd3D and wavethresh::wr3D
# dimensions of input image should be power of 2 or padded with 0 to make up to power of 2
###############################################################################################
decomp.3D <- function(x,family="DaubLeAsymm",fil.num=4,min.level=2){
  if((class(x)!="list")||(class(x[[1]])!="array")|(length(dim(x[[1]]))!=3)){
    cat("input should be a list of 3D arrays","\n")
  }
  nsample = length(x)
  dims = dim(x[[1]])
  if(min(log2(dims))<=min.level){
    cat("invalid min.level, should be less than min(log2(dim(x[[1]]))).","\n")
  }
  dim1d = prod(dims)
  wdt = matrix(NA,nsample,dim1d)

  for(i in 1:nsample){
    input = x[[i]]
    wd0 = wd3D(input,filter.number=fil.num,family=family)
    pos = 2^min.level
    subwd = wd3D(array(data=rep(0,pos^3),dim=c(pos,pos,pos)),filter.number=fil.num,family=family)
    subwd$a = wd0$a[1:pos,1:pos,1:pos]
    if(pos==1){
      inv.subwd = subwd$a
    }else{
      inv.subwd = wr3D(subwd)
    }
    wd = wd0$a
    wd[1:pos,1:pos,1:pos] = inv.subwd

    wdt[i,] = c(wd)
  }

  return(wdt)
}


#######################################################################################
# inverse 3D discrete wavelet transform with selected min.level
# inverse procedure to decomp.3D but for single input vector
#######################################################################################
reconstr.3D <- function(eta,dims,family="DaubLeAsymm",fil.num=4,min.level=2){
  if(length(eta)!=prod(dims)){cat("length of eta does not match the given dims.","\n")}
  wr0 = wd3D(array(data=rep(0,prod(dims)),dim=dims),filter.number=fil.num,family=family)
  wr0$a = array(data=eta,dim=dims)

  if(min.level>0){
    pos = 2^min.level
    subwr = wr0$a[1:pos,1:pos,1:pos]
    inv.subwr = wd3D(subwr,filter.number=fil.num,family=family)
    wr0$a[1:pos,1:pos,1:pos] = inv.subwr$a
    wr = wr3D(wr0)
  }else{
    wr = wr3D(wr0)
  }


  return(wr)
}


#################################################
# function to get wavelet transform matrix B
#################################################
getB <- function(dims,family,fil.num,min.level){
  p = prod(dims)
  X = list()
  if(length(dims)==2){
    for(k in 1:p){
      xk = rep(0,p)
      xk[k] = 1
      X[[k]] = matrix(xk,dims[1],dims[2])
    }
  }else if(length(dims)==3){
    for(k in 1:p){
      xk = rep(0,p)
      xk[k] = 1
      X[[k]] = array(xk,dim=dims)
    }
  }
  if(length(dims)==2){
    B = decomp.2D(X, family=family, fil.num=fil.num, min.level=min.level)
  }else if(length(dims)==3){
    B = decomp.3D(X, family=family, fil.num=fil.num, min.level=min.level)
  }

  return(B)
}
