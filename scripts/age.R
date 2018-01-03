
#P(t | l) = P(l | t) P(t)
dt <- function(t, l, n, b){
    f <- function(tc)dtl(tc, l) * d_mutant_age2(t-tc, n, b)
    integrate(f, 0, t-1e-14)$value
}

dtl <- function(t, l) pt(t) * plt(t, l) / pl(l)
pt <- function(t) dexp(t, 1)
plt <- function(t, l) dexp(l * 1e-4, t)
pl <- function(l){
    f <- function(t)pt(t) * plt(t, l)
    integrate(f, 0, Inf)$value
}

#eq 5.1 in griffiths tavare 2003
#n : sample size
#b : subsample size
#t : time
#d_mutant_age: exact eqn
#d_mu
d_mutant_age2 <- function(t, n, b){
    sapply(t, d_mutant_age_detapprox_st, n, b)
}

d_mutant_age <- function(t, n, b){
    sapply(t, d_mutant_age_single_t, n, b)
}
d_mutant_age_single_t <- function(t, n, b){
    k <- 2:n
    a <- sum(exp(log(k) + p_nk(b, n, k, log=T) + log(PA_n(k, t, n))))
    b <- sum(exp(log(k) + p_nk(b, n, k, log=T) + lET_k(k)))
    
    a / b
}

d_mutant_age_detapprox_st<- function(t, n, b){
    ka <- EA_n(t, n)
    a <- exp(log(ka) + p_nk_cont(b, n, ka, log=T) )
    k <- 2:n
    b <- sum(exp(log(k) + p_nk(b, n, k, log=T) + lET_k(k)))
    
    a / b
}


p_nk <- function(b, n, k, log=F){
    if (log){
        lchoose(n-b-1, k-2) - lchoose(n-1, k-1)
    } else {
        choose(n-b-1, k-2) / choose(n-1, k-1)
    }
}

p_nk_cont <- function(b, n, k, log=T){
    a = lgamma(n-b) - lgamma(k-1) - lgamma(n-b - k + 2)
    b = lgamma(n) - lgamma(k) - lgamma(n - k +1)
    a - b
}

ET_k <- function(k){
    1 / choose(k, 2)
}
lET_k <- function(k){
    - lchoose(k, 2)
}

PA_n <- function(k, t, n){
    sapply(k, PA_n_single_k, t, n)
}
PA_n_single_k <- function(k, t, n){
    j <- k:n
    sign <- (-1)^(j-k)

    log_rho_jt <- -choose(j, 2) * t
    a <- log(2 * j - 1) + lrfact(k, j-1) + lffact(n, j)
    b <- lfactorial(k) + lfactorial(j-k) + lrfact(n, j)

    log_terms <- log_rho_jt + a - b

    sum(sign * exp(log_terms))
}

EA_n <- function(t, n){
    n / (n + (1-n)*exp(-t/2) )
}

rfact <- function(x,y){
    if(y==0) return (1)
    mapply(function(x, y)prod(x:(x+y-1)),
           x, y)
}
lrfact <- function(x,y){
    res <- mapply(function(x, y)sum(log(x:(x+y-1))),
           x, y)
    res[y==0] <- 0
    res
}
lffact <- function(x,y){
    res <- mapply(function(x, y)sum(log(x:(x-y+1))),
           x, y)
    res[y==0] <- 0
    res
}
