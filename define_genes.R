# This defines equivalence classes for each gene, namely all genes that an
# Ensembl ID and Entrez ID. The idea here is to be less tied to a single
# primary identifier, but treat both Ensembl and Entrez as equal partners.
#
# Here, we're generally using Ensembl 108 for our EnsDb objects, and
# whatever happens to be the current state of the OrgDb objects.

library(AnnotationHub)
ahub <- AnnotationHub(cache="ahub_cache", ask=FALSE)

lists <- list(
    `10090` = list(ensdb = "AH109367", orgdb = "AH107060"),
    `9606` = list(ensdb = "AH109336", orgdb = "AH107059"),
    `6239` = list(ensdb = "AH109275", orgdb = "AH107064"),
    `10116` = list(ensdb = "AH109438", orgdb = "AH107062"),
    `7227` = list(ensdb = "AH109306", orgdb = "AH107058"),
    `7955` = list(ensdb = "AH109309", orgdb = "AH107067"),
    `9598` = list(ensdb = "AH109433", orgdb = "AH107055")
)

library(BiocParallel)
dir.create("genes", showWarnings=FALSE)
dump <- function(x, out) {
    if (is.list(x)) {
        dump <- unlist(bplapply(x, paste, collapse="\t", BPPARAM=MulticoreParam()))
    } else {
        dump <- x
    }
    handle <- gzfile(file.path("genes", out), open="wb")
    writeLines(dump, con=handle)
    close(handle)
}

library(igraph)
for (species in names(lists)) {
    ensdb <- ahub[[lists[[species]]$ensdb]]
    orgdb <- ahub[[lists[[species]]$orgdb]]

    from.o <- select(orgdb, keys=keys(orgdb), columns="ENSEMBL")
    from.e <- select(ensdb, keys=keys(ensdb), columns="ENTREZID")

    edges <- rbind(
        DataFrame(ensembl = from.o$ENSEMBL, entrez = from.o$ENTREZID),
        DataFrame(ensembl = from.e$GENEID, entrez = from.e$ENTREZID)
    )
    edges <- edges[!is.na(edges$ensembl) & !is.na(edges$entrez),]
    edges <- unique(edges)

    g <- make_graph(rbind(edges$ensembl, edges$entrez), directed=FALSE)
    classes <- components(g)
    pooled <- names(classes$membership)

    is.entrez <- grepl("^[0-9]", pooled)
    f <- factor(classes$membership)
    stopifnot(identical(levels(f), as.character(seq_along(levels(f)))))

    names.entrez <- names(classes$membership)[is.entrez]
    f.entrez <- f[is.entrez]
    by.entrez <- split(names.entrez, f.entrez, drop=FALSE)
    dump(by.entrez, paste0(species, "_entrez.tsv.gz"))

    names.ensembl <- names(classes$membership)[!is.entrez]
    f.ensembl <- f[!is.entrez]
    by.ensembl <- split(names.ensembl, f.ensembl, drop=FALSE)
    dump(by.ensembl, paste0(species, "_ensembl.tsv.gz"))

    stopifnot(identical(names(by.entrez), names(by.ensembl)))
    by.entrez <- unname(by.entrez)
    by.ensembl <- unname(by.ensembl)

    entrez2sym <- select(orgdb, keys=unique(names.entrez), columns="SYMBOL")
    entrez2sym <- entrez2sym[!is.na(entrez2sym$SYMBOL),]
    sym.by.entrez <- split(entrez2sym$SYMBOL, entrez2sym$ENTREZID)

    ensembl2sym <- select(ensdb, keys=unique(names.ensembl), columns="SYMBOL")
    ensembl2sym <- ensembl2sym[!is.na(ensembl2sym$SYMBOL),]
    sym.by.ensembl <- split(ensembl2sym$SYMBOL, ensembl2sym$GENEID)

    # Can't figure out how to do this faster... whatever.
    by.sym <- bplapply(seq_along(by.entrez), FUN=function(i) {
        paste(
            union(
                unlist(sym.by.entrez[by.entrez[[i]]], use.names=FALSE), 
                unlist(sym.by.ensembl[by.ensembl[[i]]], use.names=FALSE)
            ),
            collapse="\t"
        )
    }, BPPARAM=MulticoreParam())
    dump(by.sym, paste0(species, "_symbols.tsv.gz"))
}
