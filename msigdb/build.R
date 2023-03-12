# Need to prune out gene sets with difficult licensing terms
# or are redundant with each other.

version <- "v2023.1"
ignore <- c(
    "c2.cp.kegg", # difficult license
    "m2.cp.kegg", # difficult license
    "c2.cp.biocarta", # difficult license
    "m2.cp.biocarta", # difficult license
    "c2.all", # includes the above
    paste0("c2.cp.", version), # includes the above
    "c3.all", # too much stuff
    "c4.all", # too much stuff
    "c5.all", # too much stuff
    "c5.go", # duplicate of existing GO terms
    "c7.all", # too much stuff
    "msigdb" # too much stuff, and also includes the difficult licenses.
)

zip_and_dump <- function(indir, outdir) {
    dir.create(outdir, showWarnings=FALSE)
    for (x in list.files(indir, pattern="\\.entrez\\.gmt$")) {
        if (any(startsWith(x, ignore))) {
            next 
        }
        contents <- readLines(file.path(indir, x))
        handle <- gzfile(file.path(outdir, paste0(x, ".gz")), open="wb")
        writeLines(contents, con=handle)
        close(handle)
    }
}

dir <- sprintf("msigdb_%s.Hs_files_to_download_locally/msigdb_%s.Hs_GMTs", version, version)
zip_and_dump(dir, "human_output")

dir <- sprintf("msigdb_%s.Mm_files_to_download_locally/msigdb_%s.Mm_GMTs", version, version)
zip_and_dump(dir, "mouse_output")
