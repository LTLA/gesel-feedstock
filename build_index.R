library(jsonlite)
library(S4Vectors)
library(BiocFileCache)
bfc <- BiocFileCache("gmt_cache", ask=FALSE)

manifest <- jsonlite::fromJSON("manifest.json", simplifyVector=FALSE)

species.descriptions <- list()
species.collections <- list()
species.counter <- list()
species.mapping <- list()
species.sets <- list()
species.genes <- list()

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
    stopifnot(current$id %in% c("entrez", "ensembl", "symbol"))

    # Setting up a species.
    cur.species <- current$species
    if (!(cur.species %in% names(species.genes))) {
        species.mapping[[cur.species]] <- list()
        species.counter[[cur.species]] <- 0L
        species.genes[[cur.species]] <- integer()
        species.sets[[cur.species]] <- integer()
        species.descriptions[[cur.species]] <- list()
        species.collections[[cur.species]] <- list()
    }

    # Loading the GMT file.
    gmt.path <- bfcrpath(bfc, current$url)
    fragments <- strsplit(readLines(gmt.path), "\t")

    set.names <- vapply(fragments, FUN=function(x) x[1], FUN.VALUE="")
    set.descriptions <- vapply(fragments, FUN=function(x) x[2], FUN.VALUE="")
    cursets <- lapply(fragments, FUN=tail, n=-2)
    unlisted <- unlist(cursets, use.names=FALSE)

    counter <- species.counter[[cur.species]]
    set.ids <- rep(seq_along(cursets) + counter, lengths(cursets))

    # Mapping to the gene IDs.
    gene.mappings <- species.mapping[[cur.species]]
    if (!(current$id %in% names(gene.mappings))) {
        gene.path <- bfcrpath(bfc, paste0("https://github.com/LTLA/gesel-feedstock/releases/download/genes-v1.0.0/", cur.species, "_", current$id, ".tsv.gz"))
        all.lines <- readLines(gene.path)
        collected <- vector("list", length(all.lines))
        keep <- all.lines != ""
        collected[keep] <- strsplit(all.lines[keep], "\t")
        gene.mappings[[current$id]] <- split(rep(seq_along(collected), lengths(collected)), unlist(collected))
        species.mapping[[cur.species]] <- gene.mappings
    }
    known.ids <- gene.mappings[[current$id]]

    m <- match(unlisted, names(known.ids))
    keep <- !is.na(m)
    stopifnot(mean(keep) >= 0.95)
    gene.ids <- known.ids[m[keep]]
    set.ids <- rep(set.ids[keep], lengths(gene.ids))
    gene.ids <- unlist(gene.ids, use.names=FALSE)
    stopifnot(length(set.ids) == length(gene.ids))

    # Filling the values.
    j <- length(species.descriptions[[cur.species]]) + 1L
    species.descriptions[[cur.species]][[j]] <- DataFrame(name=set.names, description=set.descriptions, size=lengths(cursets))
    species.collections[[cur.species]][[j]] <- DataFrame(number = length(cursets), title=current$title, description=current$description, species=cur.species, maintainer=current$maintainer, source=current$source)

    species.genes[[cur.species]] <- c(species.genes[[cur.species]], gene.ids)
    species.sets[[cur.species]] <- c(species.sets[[cur.species]], set.ids)
    species.counter[[cur.species]] <- counter + length(set.names)
}

#################################
# Defining useful functions.

tokenizer <- function(x) {
    out <- gsub("[^a-zA-Z0-9-]", " ", tolower(x))
    tokens <- strsplit(out, "\\s+")
    ids <- rep(seq_along(tokens), lengths(tokens))
    by.token <- split(ids, unlist(tokens))
    by.token[!(names(by.token) %in% c("", "-"))]
}

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

saveLines <- function(lines, path, ...) {
    write(lines, file=file.path(dir, path))

    nc <- nchar(lines, type="bytes") # deal with UTF-8 chars.

    handle <- gzfile(file.path(dir, paste0(path, ".ranges.gz")), open="wb")
    write.table(data.frame(X=nc, ...), file=handle, row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")
    close(handle)
}

#################################
# Looping across species.

for (species in names(species.descriptions)) {
    descriptions <- do.call(rbind, species.descriptions[[species]])
    collections <- do.call(rbind, species.collections[[species]])
    stopifnot(anyDuplicated(collections$title) == 0)

    all.genes <- species.genes[[species]]
    all.ids <- species.sets[[species]]
    by.gene <- split(all.ids, factor(all.genes, seq_along(species.mapping[[species]][[1]])))
    by.set <- split(all.genes, factor(all.ids, seq_len(nrow(descriptions))))

    # Build a search index for the descriptions and names.
    by.dtoken <- tokenizer(descriptions$description)
    by.ntoken <- tokenizer(descriptions$name)
    saveTabbedIndices(by.ntoken, path=paste0(species, "_tokens-names.tsv"), include.names=TRUE)
    saveTabbedIndices(by.dtoken, path=paste0(species, "_tokens-descriptions.tsv"), include.names=TRUE)

    # Creating the indices.
    saveTabbedIndices(by.gene, path=paste0(species, "_gene2set.tsv"))
    saveTabbedIndices(by.set, path=paste0(species, "_set2gene.tsv"))

    {
        tstripped <- gsub("\t|\n", " ", collections$title)
        dstripped <- gsub("\t|\n", " ", collections$description)
        common <- sprintf("%s\t%s\t%s\t%s\t%s", tstripped, dstripped, collections$species, collections$maintainer, collections$source)
        saveLines(common, paste0(species, "_collections.tsv"), size=collections$number)

        collected <- sprintf("%s\t%s", common, collections$number)
        handle <- gzfile(file.path(dir, paste0(species, "_collections.tsv.gz")))
        write(collected, file=handle, ncolumns=1)
        close(handle)
    }

    {
        nstripped <- gsub("\t|\n", " ", descriptions$name) 
        dstripped <- gsub("\t|\n", " ", descriptions$description) 
        common <- sprintf("%s\t%s", nstripped, dstripped)
        saveLines(common, paste0(species, "_sets.tsv"), size=descriptions$size)

        collected <- sprintf("%s\t%s", common, descriptions$size)
        handle <- gzfile(file.path(dir, paste0(species, "_sets.tsv.gz")))
        write(collected, file=handle, ncolumns=1)
        close(handle)
    }
}
