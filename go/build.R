GO.from.org <- function(orgdb, prefix) {
    mappings <- select(orgdb, keytype="GO", keys=keys(orgdb, "GO"), columns="ENTREZID")

    library(BiocParallel)
    output <- split(mappings$ENTREZID, mappings$GO)
    output <- bplapply(output, function(x) {
        current <- unique(sort(x))
        paste(current, collapse="\t")
    }, BPPARAM=MulticoreParam())

    # Saving the names and descriptions.
    library(GO.db)
    info <- select(GO.db, keys=names(output), column="TERM")
    payload <- sprintf("%s\t%s\t%s", names(output), info$TERM[match(names(output), info$GOID)], output)

    con <- gzfile(paste0(prefix, ".gmt.gz"), open="wb")
    write(payload, file=con)
    close(con)

    # Emitting relevant session information.
    print(metadata(orgdb))
}

library(org.Hs.eg.db)
GO.from.org(org.Hs.eg.db, "human-GO")

library(org.Mm.eg.db)
GO.from.org(org.Mm.eg.db, "mouse-GO")

library(org.Dm.eg.db)
GO.from.org(org.Dm.eg.db, "fly-GO")

library(org.Ce.eg.db)
GO.from.org(org.Ce.eg.db, "worm-GO")

library(org.Rn.eg.db)
GO.from.org(org.Rn.eg.db, "rat-GO")

library(org.Pt.eg.db)
GO.from.org(org.Pt.eg.db, "chimp-GO")

library(org.Dr.eg.db)
GO.from.org(org.Dr.eg.db, "zebrafish-GO")
