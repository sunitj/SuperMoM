# IMG/M ER Tar balls
This set of scripts load your IMG/M ER annotated data the database.

## What does it do exactly?
* Parses an IMG tar-ball.
* Consolidates data from all files in the tar-ball and creates a single tab-delimited file for all the data that needs to be loaded in the database.
* Creates the Neo4j database.

## Schema
![image](images/graphDB_schema.png)

## Screen Shot of the current database
![Sample](images/Sample.png)

## How to use it?
* Calculate the length and GC content by using the [length+GC.pl](https://github.com/Geo-omics/scripts/blob/master/SeqTools/length%2BGC.pl) script.


```
perl length+GC.pl -f img_scaffolds.fna -gc > file.length_gc
```

* Make sure you have all the Perl dependencies installed, see: [Dependencies](https://github.com/sunitj/SuperMoM#perl).
* Make sure the Neo4j server is running.
* Populate the database by running the `perl` script `createDB.pl`, like so:

```
perl createDB.pl -gff file.gff -phylodist file.phylodist -gene_prod file.geneproducts -map file.map -lgc file.length_gc
```

**WARNING:** The steps mentioned above have only been tested on a small sample set. Use these scripts at your own risk.

## ToDo:
* Streamline the script
* Test the database
* Figure out bottlenecks and resource usage