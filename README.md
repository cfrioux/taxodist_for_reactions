# Distance calculation between taxa and metabolic elements

## Description

In [MetaCyc](http://metacyc.org/), pathways are associated to taxa through *expected taxonomic ranges*. Reactions are themselves associated to pathways. This way, we can link reactions and the *expected taxonomic ranges*. 

This program relies on a PADMet version of the Metacyc database (see [padmet-utils](https://github.com/AuReMe/padmet-utils/blob/master/connection/pgdb_to_padmet.py) to build a PADMet file from the PGDB dat files [downloadable](https://biocyc.org/download.shtml) on BioCyc).

This program retrieves reactions that are associated to pathways from a PADMet file, and the expected taxonomic ranges of these pathways. Each taxon is fetched using the NCBI API to obtain information such as its lineage. Alternatively, this data is loaded from a JSON file (*an example is provided in the repository*), that can also be generated (see usage) to prevent calling the API at each run. 

For each reaction, the distance between itself and the closest of the expected taxonomic ranges is retrieved. The distance calculation relies on the topology of the taxonomic tree according to the lineages provided by [NCBI Taxonomy Browser]
(https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi). More details on the distance algorithm will be provided soon in this documentation.

For a given organism, the program can calculate the distance between the organism and the reactions.

## Dependencies

* [PADMet Library](https://github.com/AuReMe/padmet)
* [Biopython](https://biopython.org/) ([Entrez](https://biopython.org/DIST/docs/api/Bio.Entrez-module.html) module)

## Usage

```
dist_pwy_rxn.py --email [--padmet] [--fromjson] [--tojson] [--orga] 

optional arguments:
  -h, --help           show this help message and exit
  --email EMAIL        email is required to call NCBI API
  --orga ORGA          organism of interest
  --fromjson FROMJSON  json input to retrieve data
  --tojson TOJSON      json output to store data
  --padmet PADMET      padmet file with raw data on taxonomy and pathways 
  ```


/!\ The email required in the program inputs are used for calls to the NCBI API, see [Biopython Entrez](https://biopython.org/DIST/docs/api/Bio.Entrez-module.html) documentation and [NCBI Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25497/) policies. The runtime of the program is limited by the calls to the API: 3 requests per second. 
