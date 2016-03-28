# Processing human Evolutionary Rate Covariation data

This resource processes the [evolutionary rate covaration](http://csb.pitt.edu/erc_analysis/Methods.php) (ERC) data from [Priedigkeit et al 2015](https://doi.org/10.1371/journal.pgen.1004967 "Evolutionary Signatures amongst Disease Genes Permit Novel Methods for Gene Prioritization and Construction of Informative Gene-Based Networks"). ERC assesses whether two genes have a similar evolutionary history. See the corresponding [Thinklab discussion](http://doi.org/10.15363/thinklab.d57 "Selecting informative ERC (evolutionary rate covariation) values between genes") for more information.

## Files

`erc.ipynb` converts the raw ERC data into a tidy TSV (`erc_mam33.tsv.gz`) containing each gene pair mapped to Entrez Gene. `entrez-group.R` converts to using Entrez Gene as the primary gene identifiers. `erc_mam33-entrez.tsv.gz` is not tracked, since it's 821 MB, but the subset of gene pairs with ERC values >= 0.6 is stored in `erc_mam33-entrez-gt-0.6.tsv.gz`. 

## License

All original content in this repository is licensed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/ "CC0 1.0 Universal: Public Domain Dedication"). ERC data is used with permission.
