# GCD
gcd = function(a, b){
  while (b != 0) {
    r = a %% b
    a = b
    b = r
  }
  return (a)
}


# Bezout's Coefficients
bezoutCoef = function(m, n){
  g = gcd(m, n)
  m = m / g
  n = n / g
  v1 = c(1, 0)
  v2 = c(0, 1)
  while (n > 1){
    r = m %% n
    q = m %/% n
    v3 = v1 - q * v2
    v1 = v2
    v2 = v3
    m = n
    n = r
  }
  return (v2)
}
#bezoutCoef(173, 167)
#bezoutCoef(603, 626)

inverse = function(a, n) {
  if (gcd(a, n) != 1) {
    print("No inverse")
  }
  coeff = bezoutCoef(a, n)
  inverse = coeff[1]
  if (inverse < 0) {
    inverse = n + inverse
  }
  return (inverse)
}


gcdSteps = function(a, b) {
  step = 0
  while (b != 0) {
    r = a %% b
    a = b
    b = r
    step = step + 1
  }
  return (step)
}


# n must be positive integer
highestPower = function(n, b) {
  degree = 0
  while (n >= (b ^ degree)) {
    degree = degree + 1
  }
  return (degree-1)
}


# comupte the number of digits of an integer
# digitNum = function(b) {
#  return (floor(log10(b) + 1))
# }

digitNum = function(a) {
  return (highestPower(a, 10) + 1)
}


highestPowerDivisible = function(n, b) {
  degree = 0
  while (n %% (b ^ degree) == 0) {
    degree = degree + 1
  }
  return (degree-1)
}



limitOfStep = function(a, lob) {
  limitList = NULL
  for (i in 1:length(lob)) {
    b = lob[i]
    step = gcdSteps(a, b)
    digits = digitNum(b)
    limit = step / digits
    limitList = c(limitList, limit)
  }
  return (limitList)
}


#blist = NULL
#for (i in 10:1000010) {
#  blist = c(blist, i)
#}

# blist = seq(from = 10, to = 10010, by = 1)
# ycoord = limitOfStep(53751, blist)
# plot(x = blist, y = ycoord, type = "p")
#plot(x = blist, y = limitOfStep(2437377, blist), type = "p")
#hist(x = limitOfStep(5437, blist))



# n must be postive integer
binaryPowerRep = function(n) {
  powerList = NULL
  while (n > 0) {
    power = highestPower(n, 2)
    powerList= c(powerList, power)
    n = n - 2^(power)
  }
  return (powerList)
}


modExp = function(a, e, n) {
  result = 1
  if (e == 0) {
    return (result)
  }
  else {
    powerList = binaryPowerRep(e)
    powerList = sort(powerList, decreasing = FALSE)
    
    startPower = 0
    temp = 1
    for (i in (1 : length(powerList))) {
      power = powerList[i]
      for (j in (startPower : power)) {
        if (j == 0) {
          temp = a
        } else {
          temp = (temp ^ 2) %% n
        }
      }
      startPower = power + 1
      result = (result * temp) %% n
    }
    return (result)
  }
}

# modExp(504, 1109, 2881)
# modExp(1874, 1109, 2881)
# modExp(347, 1109, 2881)
# modExp(515, 1109, 2881)
# modExp(2088, 1109, 2881)
# modExp(2356, 1109, 2881)
# modExp(736, 1109, 2881)
# modExp(468, 1109, 2881)

discreteLog = function(b, r, n) {
  logs = NULL
  for (i in (1 : 50)) {
    if ((modExp(b, i, n) - r) %% n == 0) {
      logs = c(logs, i)
    }
  }
  return (logs)
}



# reduced residue system
coprimeClass = function(n) {
  classes = NULL
  for (i in (0 : (n-1))) {
    if (gcd(i, n) == 1) {
      classes = c(classes, i)
    }
  }
  return (classes)
}

eulerPhi = function(n) {
  return (length(coprimeClass(n)))
}

# ord=0 => ord=phi(n)
ord = function(a, n) {
  i = 1
  while (modExp(a, i, n) != 1) {
    i = i + 1
  }
  return (i)
}


primitiveRoot = function(n) {
  classes = coprimeClass(n)
  roots = NULL
  phi = eulerPhi(n)
  for (i in (1 : length(classes))) {
    a = classes[i]
    d = ord(a, n)
    if (d == phi) {
      roots = c(roots, a)
    }
  }
  return (roots)
}

# compute the index i such that r^i = a (mod n) where a is from integers with gcd(a, n) = 1
ind = function(r, a, n) {
  i = 0
  while(modExp(r, i, n) != a) {
    i = i + 1
  }
  return (i)
}

# display the index table for every a in {a | gcd(a, n) = 1}
indexTable = function(r, n) {
  indices = NULL
  reducedResidues = coprimeClass(n)
  for (a in reducedResidues) {
    i = ind(r, a, n)
    indices = c(indices, i)
  }
  return (indices)
}
