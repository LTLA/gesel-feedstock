library(Matrix)
dir.create("embeddings", showWarnings=FALSE)

compute_chunked_neighbors <- function(mat, l2, num_neighbors, chunk_size, chunk_id) {
    chunk_start <- (chunk_id - 1L) * chunk_size + 1L
    chunk_end <- min(chunk_start + chunk_size - 1L, ncol(mat))
    chunk <- chunk_start : chunk_end

    d2 <- l2 - 2 * crossprod(mat, mat[,chunk,drop=FALSE])
    best <- matrix(0L, num_neighbors, length(chunk))
    dist <- matrix(0, num_neighbors, length(chunk))

    # Can't escape the quadratic runtime at this point. Oh well.
    start <- 1L
    for (i in chunk) {
        current <- d2[,start]
        o <- order(current)
        neighbors <- head(o, num_neighbors + 1)
        neighbors <- head(setdiff(neighbors, i), num_neighbors)
        best[,start] <- neighbors
        dist[,start] <- current[neighbors]
        start <- start + 1
    }

    list(index = t(best), distance2 = t(dist) + l2[chunk])
}

all.files <- list.files("assets", pattern="_set2gene.tsv.gz")
for (x in all.files) {
    lines <- readLines(file.path("assets", x))
    segments <- strsplit(lines, "\t")
    segments <- lapply(segments, function(x) cumsum(as.integer(x)) + 1L)

    i <- unlist(segments, use.names=FALSE)
    j <- rep(seq_along(segments), lengths(segments))
    mat <- sparseMatrix(i = i, j = j, x = 1, dims=c(max(i), length(segments)))

    # Converting set occupancy to proportions. This ensures that the distances
    # between large sets are comparable to those distances between small sets.
    mat <- t(t(mat) / colSums(mat))

    # We do an exact neighbor search as PCA is too lossy. Each axis ends up 
    # defining a family of gene sets, and if you just take the top 50 axes,
    # you basically just end up with the top 50 most different gene sets.
    l2 <- colSums(mat^2)

    perplexity <- 30
    num_neighbors <- round(perplexity * 3)

    # Operating in chunks for speed. We also use the sparse matrix
    # multiplication to compute distances more efficiently; the risk
    # of catastrophic cancellation is acceptable, given that we're
    # going to just do a visualization anyway.
    chunk_size <- 1000L
    chunk_ids <- seq_len(ceiling(ncol(mat) / chunk_size))

    cl <- parallel::makeCluster(4, type="FORK")
    out <- parallel::parLapply(cl, 
        chunk_ids, 
        compute_chunked_neighbors, 
        mat=mat, 
        l2=l2, 
        chunk_size=chunk_size, 
        num_neighbors=num_neighbors
    )
    parallel::stopCluster(cl)

    best <- do.call(rbind, lapply(out, function(x) x$index))
    dist <- do.call(rbind, lapply(out, function(x) x$distance))

    best <- t(best)
    dist <- t(dist)

    # TODO: figure out why scran.chan doesn't like tied distances,
    # which is why this jittering bit is needed.
    set.seed(99)
    for (i in seq_len(ncol(dist))) {
        dist[,i] <- sqrt(pmax(0, sort(jitter(dist[,i]))) + 0.01)
    }

    # Forgive me, father, for I have sinned. But I just couldn't wait
    # for Rtsne to get its shit together, hence this ':::' import.
    tsne_out <- scran.chan:::run_tsne(
        nnidx = best - 1L,
        nndist = dist,
        perplexity = perplexity,
        interpolate = FALSE,
        max_depth = 7L,
        max_iter = 500L,
        seed = 42,
        nthreads = 1 
    )

    species <- sub("_.*", "", x)
    png(paste0("embeddings/", species, "_tsne.png"), width=6, height=6, units="in", res=150)
    plot(tsne_out[1,], tsne_out[2,], pch=16, cex=0.1, xlab="t-SNE 1", ylab="t-SNE 2", main = species)
    dev.off()

    handle <- gzfile(paste0('embeddings/', species, '_tsne.tsv.gz'), open="wb")
    write.table(t(tsne_out), handle, sep="\t")
    close(handle)
}
