gen_3d <- function (order, x, y, z, xi, xj, xk, yi, yj, yk, zi, zj, zk, p_array){
  if (order == 0){
    xx = x + (xi + yi + zi)/3
    yy = y + (xj + yj + zj)/3
    zz = z + (xk + yk + zk)/3
    if (length(p_array[[1]]) == 0){
      p_array = list(list(xx, yy, zz))
      return(p_array)
    } else {
      p_array = append(p_array, list(list(xx, yy, zz)))    }
    return(p_array)
    } else {
      p_array = gen_3d(order-1, x, y, z, yi/2, yj/2, yk/2, zi/2, zj/2, zk/2, xi/2, xj/2, xk/2, p_array)

      p_array = gen_3d(order-1, x + xi/2, y + xj/2, z + xk/2,  zi/2, zj/2, zk/2, xi/2, xj/2, xk/2,
             yi/2, yj/2, yk/2, p_array)
      p_array = gen_3d(order-1, x + xi/2 + yi/2, y + xj/2 + yj/2, z + xk/2 + yk/2, zi/2, zj/2, zk/2,
             xi/2, xj/2, xk/2, yi/2, yj/2, yk/2, p_array)
      p_array = gen_3d(order-1, x + xi/2 + yi, y + xj/2+ yj, z + xk/2 + yk, -xi/2, -xj/2, -xk/2, -yi/2,
             -yj/2, -yk/2, zi/2, zj/2, zk/2, p_array)
      p_array = gen_3d(order-1, x + xi/2 + yi + zi/2, y + xj/2 + yj + zj/2, z + xk/2 + yk +zk/2, -xi/2,
             -xj/2, -xk/2, -yi/2, -yj/2, -yk/2, zi/2, zj/2, zk/2, p_array)
      p_array = gen_3d(order-1, x + xi/2 + yi + zi, y + xj/2 + yj + zj, z + xk/2 + yk + zk, -zi/2, -zj/2,
             -zk/2, xi/2, xj/2, xk/2, -yi/2, -yj/2, -yk/2, p_array)
      p_array = gen_3d(order-1, x + xi/2 + yi/2 + zi, y + xj/2 + yj/2 + zj , z + xk/2 + yk/2 + zk, -zi/2,
             -zj/2, -zk/2, xi/2, xj/2, xk/2, -yi/2, -yj/2, -yk/2, p_array)
      p_array = gen_3d(order-1, x + xi/2 + zi, y + xj/2 + zj, z + xk/2 + zk, yi/2, yj/2, yk/2, -zi/2, -zj/2,
             -zk/2, -xi/2, -xj/2, -xk/2, p_array)
      return(p_array)
    }
}

hilbert_3d <- function(order){
  n = 2^(order)
  hilbert_curve = list(list())
  hilbert_curve = gen_3d(order, 0, 0, 0, n, 0, 0, 0, n, 0, 0, 0, n, hilbert_curve)

  res <- rapply(hilbert_curve[], as.integer,  how = "replace")

  return(res)
}


sf_objects <-function (sf_obj){
  coords_mat = st_coordinates(sf_obj)
  M1 = coords_mat[,1:3]
  if (sum(colnames(M1) != c("X", "Y", "Z"))!=0){
    M1[,3] = 0
  }
  colnames(M1) <- c("x", "y", "z")
  return (M1)
}

sort <- function(geodf, origin, radius, bins){
  M1 = sf_objects(geodf)
  M2 = sweep(M1, 2, origin)
  bin_interval = ((radius*2) / bins)
  offset = as.integer(radius/bin_interval)
  binned = array(,c(bins, bins, bins, 3))
  for (i in 1:nrow(M2)){
    x = as.integer(M2[i,][1]/bin_interval) + offset
    y = as.integer(M2[i,][2]/bin_interval) + offset
    z = as.integer(M2[i,][3]/bin_interval) + offset
    if (x > bins-1 || x < 0){
      next
    }
    if (y > bins-1 || y < 0){
      next
    }
    if (z > bins-1 || x < 0){
      next
    }

    binned[x,y,z,] = M2[i,]
  }

  order = np.log2(32)
  curve  = hilbert_3d(order)

  for (k in 1:length(curve)) {
    curve[[k]][1]
    #x1 = binned
  }

}

sf_obj = st_crs(nc)



########################
##### NEW VERSION ######
########################

partition <- function(A, st, en, od, ax, di){
  i = st
  j = en
  while (TRUE){
    while(i<j & as.logical(BitTest(A[[i+1]][ax+1], od)) == di){
      i = i + 1
    }
    while(i<j & as.logical(BitTest(A[[j+1]][ax+1], od))!= di){
      j = j - 1
    }
    if (j <= i){
      return (list(i, A))
    }
    temporal = A[[i+1]]
    A[[i+1]] = A[[j+1]]
    A[[j+1]] = temporal
  }
}


BitTest <-function (x, od){
  result = bitwAnd (x, (1 *(2**od) ) )
  return (as.integer(as.logical(result)))
}


BitFlip<-function(b, pos){
  b = bitwXor(b, 1 * (2**pos))
  return (b)
}



HSort <- function (input, st, en, od, c, e, d, di, cnt){
  if (en<=st){
    return (input)
  }
  vec = partition(input, st, en, od, (d + c) %% n, BitTest(e, (d + c) %% n))
  p = as.integer(vec[1])
  input = vec[[2]]

  if (c == n - 1){
    if (od ==0){
      return (input)
    }

    d2 = (d + n + n - (if (di){2} else{cnt + 2} )) %% n
    e = BitFlip(e, d2)
    e = BitFlip(e, (d + c) %% n)
    input = HSort(input, st, p - 1, od - 1, 0, e, d2, FALSE, 0)

    e = BitFlip(e, (d + c) %% n)
    e = BitFlip(e, d2)
    d2 = (d + n + n - (if (di){cnt + 2} else {2})) %% n
    input = HSort(input, p, en, od - 1, 0, e, d2, FALSE, 0)
    return (input)

  } else{
    input = HSort(input, st, p - 1, od, c + 1, e, d, FALSE, (if (di){1} else{cnt + 1}) )
    e = BitFlip(e, (d + c) %% n)
    e = BitFlip(e, (d + c + 1) %% n)
    input = HSort(input, p, en, od, c + 1, e, d, TRUE, (if (di){1} else{cnt + 1}) )
    e = BitFlip(e, (d + c + 1) %% n)
    e = BitFlip(e, (d + c) %% n)
    return (input)
  }
}


N = 4  # 9 points
n = 3  # 2 dimension
m = 3  # order of Hilbert curve


st=0
en=N - 1
od=m - 1
c=0
e=0
d=0
di=FALSE
cnt=0

#input =list(c(-72, 2, 1), c(-2, 4, 1) ,c(-3, 4,1), c(-2, 5, 1), c(3, 5, 1),
#                    c(1, 6, 1), c(-32, 6,1), c(5, -60,1), c(3, -80,1) )


input = list(c(2,2,1),c(2,4,1),c(30,4,1), c(2,5,1) )

outcome = HSort(input, st=0,  en=N - 1, od=m - 1, c=0, e=0, d=0, di=FALSE, cnt=0)
outcome
