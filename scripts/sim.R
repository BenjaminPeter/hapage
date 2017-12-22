require(dplyr)
init <- function(L){
    chrom0 <- tibble(x0=0, x1=L, l=L, n=1, parent=0, age=0, gen=0, id=0)
}

#this version keeps all the generated fragments
sim1 <- function(last_gen, r=1){
    last_gen_id <- tail(last_gen$gen, 1)
    new_gen_id <- last_gen_id + 1
    n_haps <- nrow(last_gen)
    max_hap_id <- max(last_gen$id)

    #rec step
    recs <- rpois(n_haps, last_gen$n * last_gen$l * r)
    recs <- pmin(recs, last_gen$n) #at most 1 rec event per chunk, while this is an approximation it only matters in the first few generations
    n_rec <- sum(recs)

    if(n_rec > 0){
        #set up data for all the new chunks by recombination
        frag_length <- rep(last_gen$x1-last_gen$x0, recs) 
        rec_loc <- runif(n_rec, 0, frag_length) + rep(last_gen$x0, recs)

        x0 <- c(rec_loc,                rep(last_gen$x0, recs))
        x1 <- c(rep(last_gen$x1, recs), rec_loc)
        l <- x1 - x0 
        n <- 1
        parent_id <- rep(rep(last_gen$id, recs), 2)
        age <- 1
        new_ids <- 1:(2 * n_rec) + max_hap_id
        new_chroms <- tibble(x0, x1, l, n, age=age, id=new_ids,
                                 parent=parent_id, gen=new_gen_id)


        chroms <- last_gen %>% mutate(n=n-recs, age=age+1) %>% 
            bind_rows(new_chroms) %>%
            filter(n > 0)
    } else {
        chroms <- last_gen #%>% filter(n>0)
        chroms$age <- chroms$age + 1

    }
    #propagation step
    #this may need to be changed to account for demography, dilution, etc.
    chroms$n <- rpois(nrow(chroms), chroms$n)
    chroms$gen <- new_gen_id 
    chroms[chroms$n > 0,]
}

simn <- function(L=1e6, n_gen=1500, r=1e-8){
    data <- list()
    data[[1]] <- init(L)
    for(i in 2:(n_gen+1)){
        data[[i]] <- sim1(data[[i-1]], r)
        if(nrow(data[[i]]) == 0) break #extinct
    }
    do.call(bind_rows, data)
}

prun <- function(n=10000, r=1e-8, L=1e6, n_gen=2000){
    library(parallel)
    res <- mclapply(1:n, function(i){if(i%%1000 ==0)print(i); simn(r=r, L=L, n_gen=n_gen)},
                    mc.cores=detectCores())
    res 
}

bigrun <- function(){
    x <- prun(1e6)
    saveRDS(x, "data/test1M.rds")
}

