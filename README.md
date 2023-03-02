# Build gene sets to feed gesel

This repository contains code to build indices for **gesel**, the client-side gene set search interface.

To contribute, make a [pull request](https://github.com/LTLA/gesel-feedstock/pulls) and add an entry to [`manifest.json`](manifest.json) to point to a GMT file of your choice.
Each entry in the array represents a collection with the following metadata:

- `title`: the title of the collection.
  This should not contain tabs or newlines.
- `description`: the description of the collection.
  This should not contain tabs or newlines.
- `species`: the full name of the species.
- `maintainer`: the name of the maintainer of the collection.
- `source`: the source of the collection.
  This may reference an article or the code used to generate the collection, and is intended for human readers.
- `url`: the URL to the collection's GMT file.
  This should be downloadable.
  The GMT file should use Ensembl identifiers for all genes.
