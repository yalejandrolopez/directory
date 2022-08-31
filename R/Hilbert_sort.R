
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
