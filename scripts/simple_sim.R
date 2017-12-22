require(dplyr)

prun <- function(n=1000){
    library(parallel)
    cl <- makeCluster(detectCores()-1)  
    #cl <- makeCluster(2)  
    #get library support needed to run the code
    #clusterEvalQ(cl,source("sims.R"))
    #put objects in place that might be needed for the code
    #clusterExport(cl,c("myData"))
    #... then parallel replicate...
    res <- mclapply(1:n, function(i)sim_single_frag_ngen(r=1e-8, L=1e6, n_gen=1600),
                    mc.cores=8)
    #res <- mclapply(1:n, function(i){x <- runif(i)})
    #stop the cluster
    stopCluster(cl)

    res 
}
sim_single_frag_ngen <- function(L=1e8, n_gen=1500, r=1){
    data <- list()
    data[[1]] <- data.frame(x0=0, x1=L, l=L, n=1, parent=0, gen=1)
    for(i in 2:n_gen){
        data[[i]] <- sim_single_frag_1gen(data[[i-1]], r)
        if(nrow(data[[i]]) == 0) break #extinct
    }
    do.call(rbind, data)
}
sim_single_frag_1gen <- function(last_gen, r=1){
    last_gen_id <- tail(last_gen$gen, 1)
    new_gen_id <- last_gen_id + 1
    n_haps <- nrow(last_gen)

    #rec step
    recs <- rpois(n_haps, last_gen$n * last_gen$l * r)
    recs <- pmin(recs, last_gen$n) #at most 1 rec event per chunk, while this is an approximation it only matters in the first few generations
    n_rec <- sum(recs)

    if(n_rec > 0){
        #set up data for all the new chunks by recombination
        frag_length <- rep(last_gen$x1-last_gen$x0, recs) 
        rec_loc <- runif(n_rec, 0, frag_length) + rep(last_gen$x0, recs)
        keep_left <- runif(n_rec) < 0.5


        x0 <- ifelse(!keep_left, rec_loc, rep(last_gen$x0, recs))
        x1 <- ifelse( keep_left, rec_loc, rep(last_gen$x1, recs))
        l <- x1 - x0 
        n <- 1
        parent <- rep(1:n_haps, recs) 
        new_chroms <- data.frame(x0, x1, l, n, parent, gen=new_gen_id)

        chroms <- last_gen %>% mutate(n=n-recs) %>% 
            bind_rows(new_chroms) %>%
            filter(n > 0)
    } else {
        chroms <- last_gen #%>% filter(n>0)

    }
    #propagation step
    #this may need to be changed to account for demography, dilution, etc.
    pop_size <- 0
    n_new <- rpois(nrow(chroms), chroms$n)

    #from here on it's just WF
    #n_new <- c(rmultinom(1, pop_size, chroms$n))
    chroms$n <- n_new
    chroms$gen <- new_gen_id 
    chroms[chroms$n > 0,]
}
