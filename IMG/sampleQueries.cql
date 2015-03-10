// Top Organisms in the Database
MATCH (sn:SciName)<-[r:BELONGS_TO]-() RETURN sn.sciName, count(r) ORDER BY count(r) DESC;

// Genes on each Contig
MATCH (ig:ImgGenes)-->(gc:GeneCoords)-[]->(ic:ImgContigs)-->(sn:SciName) RETURN ic.imgContigName,gc.geneStart,gc.geneStop,ig.geneProduct,sn.sciName ORDER BY ic.imgContigName;
<<<<<<< Updated upstream

// Number of Scaffolds in a Genus
MATCH (a)-[r:`IN_GENUS`]->(b) RETURN b.taxaGenus,count(r) ORDER BY count(r) DESC LIMIT 25
=======
>>>>>>> Stashed changes
