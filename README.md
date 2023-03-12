# Build gene sets to feed gesel

This repository contains code to build indices for **gesel**, the client-side gene set search interface.
The indices themselves are available on the [Releases page](https://github.com/LTLA/gesel-feedstock/releases);
this can be fetched by applications directly or via a CORS proxy.

## Overview of resources

Each species contains a separate copy of the files described in this section.
Files from a particular species will be prefixed with that species' NCBI taxonomy ID, e.g., `9606_ensembl.tsv.gz`.
For brevity, we will omit the prefix in the rest of this section.

### Gene mappings

Genes are defined as abstract "equivalence classes" that can be associated with one, zero or multiple identifiers or symbols.
Each equivalence class is defined as a component of the graph constructed from the relationships between Ensembl and Entrez identifiers.

- `ensembl.tsv.gz` is a Gzip-compressed tab-separated file where each line corresponds to a gene equivalence class, and the fields are Ensembl identifiers associated with that gene.
  An empty line indicates that the equivalence class contains no Ensembl IDs.
- `entrez.tsv.gz` is a Gzip-compressed tab-separated file where each line corresponds to a gene equivalence class, and the fields are Entrez identifiers associated with that gene.
  An empty line indicates that the equivalence class contains no Entrez IDs.
- `symbol.tsv.gz` is a Gzip-compressed tab-separated file where each line corresponds to a gene equivalence class, and the fields are symbols associated with that gene.
  An empty line indicates that the equivalence class contains no symbols.

All files have the same number of lines as they represent difference aspects of the same underlying array of equivalence classes.
Each gene's identity (i.e., the "gene ID") is defined as the 0-based index of the corresponding line in either file.

### Collection details

`collections.tsv.gz` is a Gzip-compressed tab-separated file where each line corresponds to a gene set collection and contains the following fields:

- `title`: the title of the collection.
- `description`: the description of the collection.
- `species`: the species involved in the collection.
- `maintainer`: the maintainer of the collection.
- `source`: the source URL for the collection.
- `number`: the number of gene sets in this collection.

Each collection's identity (i.e., the "collection ID") is defined as the 0-based index of the corresponding line in `collections.tsv.gz`.

`collections.tsv` is an uncompressed tab-separated file that has the same number of lines and order of collections as `collections.tsv.gz`.
It contains all fields in `collections.tsv.gz` except for `number`.

`collections.tsv.ranges.gz` is a Gzip-compressed file where each line corresponds to a collection in `collections.tsv`.
Each line contains the following fields:

- `bytes`: the number of bytes taken up by the corresponding line in `collections.tsv` (excluding the newline).
- `number`: the number of gene sets in this collection.

Applications can either download `collections.tsv.gz` to obtain information about all collections,
or they can download `collections.tsv.ranges.gz` and perform HTTP range requests on `collections.tsv` to obtain information about individual collections.
The former pays a higher up-front cost for easier batch processing.
To reduce the download size, we do not store `start` in `collections.tsv.gz`, as these can be computed easily on the client. 

### Set details

`sets.tsv.gz` is a Gzip-compressed tab-separated file where each line corresponds to a gene set and contains the following fields:

- `name`: the name of the set.
- `description`: the description of the set.
- `size`: the number of genes in the set.

Each set's identity (i.e., the "set ID") is defined as the 0-based index of the corresponding line in `sets.tsv.gz`.
Sets from the same collection are always stored in consecutive lines, ordered by their position within that collection.

`sets.tsv` is an uncompressed tab-separated file that has the same number of lines and order of sets as `sets.tsv.gz`.
It contains all fields in `sets.tsv.gz` except for `size`.

`sets.tsv.ranges.gz` is a Gzip-compressed file where each line corresponds to a set in `sets.tsv`.
Each line contains two tab-separated fields:

- `bytes`: the number of bytes taken up by the corresponding line in `sets.tsv` (excluding the newline).
- `size`: the number of genes in the set.

Applications can either download `sets.tsv.gz` to obtain information about all sets,
or they can download `sets.tsv.ranges.gz` and perform HTTP range requests on `sets.tsv` to obtain information about individual sets.
The former pays a higher up-front cost for easier batch processing.
To reduce the download size, we do not store `collection` and `position` in `sets.tsv.gz`, as these can be computed easily on the client. 

### Mappings between sets and genes

`set2gene.tsv` is a tab-separated file where each line corresponds to a gene set in the same order as `sets.tsv.gz`.
On each line, the first field contains the gene ID of the first gene in the set.
All subsequent fields contain increments from the preceding ID, i.e., computing the cumulative sum across all fields yields the array of gene IDs for this set.

`set2gene.tsv.ranges.gz` is a Gzip-compressed file where each line corresponds to a set in `set2gene.tsv`.
Each line contains an integer specifying the number of bytes taken up by the corresponding line in `set2gene.tsv` (excluding the newline).
This can be used for HTTP range requests to obtain the composition of each set.

`gene2set.tsv` is a tab-separated file where each line corresponds to a gene in the same order as `symbol2gene.tsv.gz`.
On each line, the first field contains the set ID of the first set containing that gene.
All subsequent fields contain increments from the preceding ID, i.e., computing the cumulative sum across all fields yields the array of IDs of sets containing this gene.

`gene2set.tsv.ranges.gz` is a Gzip-compressed file where each line corresponds to a set in `gene2set.tsv`.
Each line contains an integer specifying the number of bytes taken up by the corresponding line in `gene2set.tsv` (excluding the newline).
This can be used for HTTP range requests to obtain the identities of the sets containing a particular gene.

`set2gene.tsv.gz` is a Gzip-compressed version of `set2gene.tsv`.
Similarly, `gene2set.tsv.gz` is a Gzip-compressed version of `gene2set.tsv`.
Applications can either download these `*.tsv.gz` files to obtain all relationships up-front,
or they can download `*.ranges.gz` and perform HTTP range requests on the corresponding `*.tsv` to obtain each individual relationship.

### Text search tokens

`tokens-names.tsv` is a tab-separated file where each line corresponds to a token.
On each line, the first field contains the set ID of the first set where the name contains the corresponding token.
All subsequent fields contain increments from the preceding ID, i.e., computing the cumulative sum across all fields yields the array of set IDs that contain this token in its name.

`tokens-descriptions.tsv` is a tab-separated file where each line corresponds to a token.
On each line, the first field contains the set ID of the first set where the description contains the corresponding token.
All subsequent fields contain increments from the preceding ID, i.e., computing the cumulative sum across all fields yields the array of set IDs that contain this token in its description.

`tokens-names.tsv.ranges.gz` is a Gzip-compressed file where each line corresponds to a set in `tokens-names.tsv`.
Each line contains:

- `token`: a token string.
- `number`: an integer specifying the number of bytes taken up by the corresponding line in `tokens-names.tsv` (excluding the newline).

The same logic applies to `tokens-descriptions.tsv.ranges.gz` for `tokens-descriptions.tsv`.
This can be used for HTTP range requests to obtain the identities of the sets that match tokens in their names or descriptions.

The tokenization strategy is very simple - every contiguous stretch of ASCII alphanumeric characters or dashes (`-`) is treated as a token.
Query strings should be processed in the same manner to generate tokens for matching against the indices.
Handling of `?` or `*` wildcards is at the discretion of the client implementation.

`tokens-names.tsv.gz` is a Gzip-compressed version of `tokens-names.tsv`.
Similarly, `tokens-descriptions.tsv.gz` is a Gzip-compressed version of `tokens-descriptions.tsv`.
Applications can either download these `*.tsv.gz` files to obtain all relationships up-front,
or they can download `*.ranges.gz` and perform HTTP range requests on the corresponding `*.tsv` to obtain each individual relationship.

## Contributing gene sets

Make a [pull request](https://github.com/LTLA/gesel-feedstock/pulls) and add an entry to [`manifest.json`](manifest.json) to point to a GMT file of your choice.
Each entry in the array represents a collection with the following metadata:

- `title`: the title of the collection.
  This should not contain tabs or newlines.
- `description`: the description of the collection.
  This should not contain tabs or newlines.
- `species`: the NCBI taxonomy ID for the species.
- `maintainer`: the name of the maintainer of the collection.
- `source`: the source of the collection.
  This may reference an article or the code used to generate the collection, and is intended for human readers.
- `id`: the type of identifier.
  This should be one of `"entrez"`, `"ensembl"` or `"symbol"`; the former two are more reliable and preferred.
- `url`: the URL to the collection's GMT file.
  This should be downloadable.
  The GMT file should use Ensembl identifiers for all genes.
