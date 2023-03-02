library(jsonlite)
library(S4Vectors)
library(BiocFileCache)
bfc <- BiocFileCache("gmt_cache", ask=FALSE)

manifest <- jsonlite::fromJSON("manifest.json", simplifyVector=FALSE)
descriptions <- list()
collections <- list()
all.genes <- all.ids <- integer(0)
species.genes <- list()
counter <- 0L

for (i in seq_along(manifest)) {
    current <- manifest[[i]]
    gc()

    # Checking validity.
    stopifnot(isSingleString(current$title), !grepl("[\n\t]", current$title))
    stopifnot(isSingleString(current$description), !grepl("[\n\t]", current$description))
    stopifnot(isSingleString(current$maintainer), !grepl("[\n\t]", current$maintainer))
    stopifnot(isSingleString(current$species), !grepl("[\n\t]", current$species))
    stopifnot(isSingleString(current$source), !grepl("[\n\t]", current$source))
    stopifnot(isSingleString(current$url), !grepl("[\n\t]", current$url))

    # Loading the GMT file.
    gmt.path <- bfcrpath(bfc, current$url)
    fragments <- strsplit(readLines(gmt.path), "\t")

    set.names <- vapply(fragments, FUN=function(x) x[1], FUN.VALUE="")
    set.descriptions <- vapply(fragments, FUN=function(x) x[2], FUN.VALUE="")
    cursets <- lapply(fragments, FUN=tail, n=-2)
    unlisted <- unlist(cursets)

    # Checking the species.
    regex <- c(
        `Homo sapiens` = "^ENSG[0-9]{11}$",
        `Mus musculus` = "^ENSMUSG[0-9]{11}$",
        `Rattus norvegicus` = "^ENSRNOG[0-9]{11}$",
        `Pan troglodytes` = "^ENSPTRG[0-9]{11}$",
        `Caenorhabditis elegans` = "^WBGene[0-9]{8}$",
        `Danio rerio` = "^ENSDARG[0-9]{11}$",
        `Drosophila melanogaster` = "^FBgn[0-9]{7}$"
    )
    cur.species <- match.arg(current$species, names(regex))
    stopifnot(all(grepl(regex[[cur.species]], unlisted)))
    species.genes[[cur.species]] <- union(species.genes[[cur.species]], unlisted)

    # Filling the values.
    descriptions[[i]] <- DataFrame(name=set.names, description=set.descriptions, size=lengths(cursets))
    collections[[i]] <- DataFrame(number = length(cursets), title=current$title, description=current$description, species=cur.species, maintainer=current$maintainer, source=current$source)

    all.genes <- c(all.genes, unlisted)
    all.ids <- c(all.ids, rep(seq_along(cursets) + counter, lengths(cursets)))
    counter <- counter + length(cursets)
}

descriptions <- do.call(rbind, descriptions)
collections <- do.call(rbind, collections)
stopifnot(anyDuplicated(collections$title) == 0)

u.genes <- unique(all.genes)
all.genes <- match(all.genes, u.genes)
by.gene <- split(all.ids, factor(all.genes, seq_along(u.genes)))
by.set <- split(all.genes, factor(all.ids, seq_len(nrow(descriptions))))

#################################
# Trying our damned best to map genes to one or more symbols and/or Entrez IDs.

gathered <- list(
    `Mus musculus` = "AH109367",
    `Homo sapiens` = "AH109336",
    `Caenorhabditis elegans` = "AH109275",
    `Rattus norvegicus` = "AH109438",
    `Drosophila melanogaster` = "AH109306",
    `Danio rerio` = "AH109309",
    `Pan troglodytes` = "AH109433"
)

library(AnnotationHub)
found.symbol <- found.symbol.ids <- found.entrez <- found.entrez.ids <- character(0)
ahub <- AnnotationHub(cache="ensdb_cache")

for (x in names(species.genes)) {
    gc()
    ensdb <- ahub[[gathered[[x]]]]
    
    mapped.symbol <- select(ensdb, keys=species.genes[[x]], keytype="GENEID", columns="SYMBOL")
    keep <- !is.na(mapped.symbol$SYMBOL)
    found.symbol.ids <- c(found.symbol.ids, mapped.symbol$GENEID[keep])
    found.symbol <- c(found.symbol, mapped.symbol$SYMBOL[keep])

    mapped.entrez <- select(ensdb, keys=species.genes[[x]], keytype="GENEID", columns="ENTREZID")
    keep <- !is.na(mapped.entrez$ENTREZID)
    found.entrez.ids <- c(found.entrez.ids, mapped.entrez$GENEID[keep])
    found.entrez <- c(found.entrez, mapped.entrez$ENTREZID[keep])
}

u.genes <- unlist(species.genes, use.names=FALSE)
found.symbol <- c(u.genes, found.symbol) # get the Ensembl ID at the front.
found.symbol.ids <- c(u.genes, found.symbol.ids)
symbol.mapping <- split(found.symbol, factor(found.symbol.ids, levels=u.genes))

found.entrez <- c(u.genes, found.entrez) # get the Ensembl ID at the front.
found.entrez.ids <- c(u.genes, found.entrez.ids)
entrez.mapping <- split(found.entrez, factor(found.entrez.ids, levels=u.genes))

#################################
# Build a search index for the descriptions and names.

tokenizer <- function(x) {
    out <- gsub("[^a-zA-Z0-9-]", " ", tolower(x))
    tokens <- strsplit(out, "\\s+")
    ids <- rep(seq_along(tokens), lengths(tokens))
    by.token <- split(ids, unlist(tokens))
    by.token[!(names(by.token) %in% c("", "-"))]
}

by.dtoken <- tokenizer(descriptions$description)
by.ntoken <- tokenizer(descriptions$name)

#################################

dir <- "assets"
dir.create(dir)

saveTabbedIndices <- function(y, path, include.names = FALSE) {
    x <- vapply(y, function(z) {
        z <- sort(z) # convert to diffs to reduce integer size
        z <- c(z[1] - 1L, diff(z)) # get to 0-based indexing.
        paste(z, collapse="\t")
    }, "")
    write(x, file=file.path(dir, path))

    strlen <- nchar(x, type="bytes") # deal with UTF-8 chars.
    handle <- gzfile(file.path(dir, paste0(path, ".ranges.gz")), open="wb")
    if (!include.names || is.null(names(y))) {
        write(strlen, file=handle, ncolumns=1)
    } else {
        write.table(data.frame(X=names(y), Y=strlen), col.names=FALSE, row.names=FALSE, quote=FALSE, file=handle, sep="\t")
    }
    close(handle)
}

saveTabbedIndices(by.gene, path="gene2set.tsv")
saveTabbedIndices(by.set, path="set2gene.tsv")
saveTabbedIndices(by.ntoken, path="tokens-names.tsv", include.names=TRUE)
saveTabbedIndices(by.dtoken, path="tokens-descriptions.tsv", include.names=TRUE)

saveLines <- function(lines, path) {
    write(lines, file=file.path(dir, path))
    handle <- gzfile(file.path(dir, paste0(path, ".ranges.gz")))
    write(nchar(lines, type="bytes"), file=handle, ncolumns=1) # deal with UTF-8 chars.
    close(handle)
}

{
    collections.start <- c(0L, head(cumsum(collections$number), -1))
    tstripped <- gsub("\t|\n", " ", collections$title)
    dstripped <- gsub("\t|\n", " ", collections$description)
    collected <- sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t$s", collections$id, collections.start, collections$number, tstripped, dstripped, collections$species, collections$maintainer, collections$source)
    saveLines(collected, "collections.tsv")

    collected <- sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s", collections$id, collections$number, tstripped, dstripped, collections$species, collections$maintainer, collections$source)
    handle <- gzfile(file.path(dir, "collections.tsv.gz"))
    write(collected, file=handle, ncolumns=1)
    close(handle)
}

{
    descriptions.collections <- rep(seq_along(collections$number) - 1L, collections$number)
    stopifnot(length(descriptions.collections) == nrow(descriptions))
    descriptions.internal <- sequence(collections$number) - 1L
    stopifnot(length(descriptions.internal) == nrow(descriptions))

    nstripped <- gsub("\t|\n", " ", descriptions$name) 
    dstripped <- gsub("\t|\n", " ", descriptions$description) 
    collected <- sprintf("%s\t%s\t%s\t%s\t%s", nstripped, dstripped, descriptions$size, descriptions.collections, descriptions.internal)
    saveLines(collected, "sets.tsv")

    collected <- sprintf("%s\t%s\t%s", nstripped, dstripped, descriptions$size)
    handle <- gzfile(file.path(dir, "sets.tsv.gz"))
    write(collected, file=handle, ncolumns=1)
    close(handle)
}

collected <- vapply(symbol.mapping, function(x) paste(gsub("\t|\n", " ", x), collapse="\t"), "")
handle <- gzfile(file.path(dir, "symbol2gene.tsv.gz"))
write(collected, file=handle)
close(handle)

collected <- vapply(entrez.mapping, function(x) paste(gsub("\t|\n", " ", x), collapse="\t"), "")
handle <- gzfile(file.path(dir, "entrez2gene.tsv.gz"))
write(collected, file=handle)
close(handle)
