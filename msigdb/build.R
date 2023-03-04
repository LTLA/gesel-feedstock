# Need to prune out gene sets with difficult licensing terms.

version <- "v2022.1"
ignore <- c(
    "c2.cp.kegg", # difficult license
    "c2.cp.biocarta", # difficult license
    "c2.all", # includes the above
    paste0("c2.cp.", version), # includes the above
    "c3.all", # too much stuff
    "c4.all", # too much stuff
    "c5.all", # too much stuff
    "c5.go", # duplicate of existing GO terms
    "c7.all", # too much stuff
    "msigdb" # too much stuff, and also includes the difficult licenses.
)

library(org.Hs.eg.db)
mapping <- mapIds(org.Hs.eg.db, keys=keys(org.Hs.eg.db), keytype="ENTREZID", column="ENSEMBL")

dir.create("output")
dir <- sprintf("msigdb_%s.Hs_files_to_download_locally/msigdb_%s.Hs_GMTs", version, version)

for (x in list.files(dir, pattern="\\.entrez\\.gmt$")) {
    if (any(startsWith(x, ignore))) {
        next 
    }

    contents <- strsplit(readLines(file.path(dir, x)), "\t")
    output <- character(length(contents))
    for (i in seq_along(contents)) {
        current <- contents[[i]]
        ids <- mapping[tail(current, -2)]
        output[i] <- paste(c(head(current, 2), ids[!is.na(ids)]), collapse="\t")
    }

    handle <- gzfile(file.path("output", paste0(x, ".gz")), open="wb")
    write(output, file=handle)
    close(handle)
}
