// Load IMG Data

CREATE CONSTRAINT ON (project:Project) ASSERT project.projectID IS UNIQUE;
CREATE CONSTRAINT ON (contig:ImgContigs) ASSERT contig.name IS UNIQUE;
CREATE CONSTRAINT ON (gene:ImgGenes) ASSERT gene.name IS UNIQUE;
CREATE CONSTRAINT ON (locus:GeneCoords) ASSERT locus.name IS UNIQUE;
CREATE CONSTRAINT ON (gsource:GeneSource) ASSERT gsource.name IS UNIQUE;
CREATE CONSTRAINT ON (domain:TaxaDomain) ASSERT domain.name IS UNIQUE;
CREATE CONSTRAINT ON (phylum:TaxaPhylum) ASSERT phylum.name IS UNIQUE;
CREATE CONSTRAINT ON (class:TaxaClass) ASSERT class.name IS UNIQUE;
CREATE CONSTRAINT ON (torder:TaxaOrder) ASSERT torder.name IS UNIQUE;
CREATE CONSTRAINT ON (family:TaxaFamily) ASSERT family.name IS UNIQUE;
CREATE CONSTRAINT ON (genus:TaxaGenus) ASSERT genus.name IS UNIQUE;
CREATE CONSTRAINT ON (sci:SciName) ASSERT sci.name IS UNIQUE;

//CREATE INDEX ON :ImgContigs(ImgContigName);
//CREATE INDEX ON :ImgContigs(originalContigName);
//CREATE INDEX ON :ImgGenes(geneImgName);

USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:/Users/sunitj/Neo4j/IMG/Staging/testRun2.tsv' AS line FIELDTERMINATOR '\t'
MERGE (contig:Contigs {name:LOWER(line.ImgContigID)}) 
SET contig.length = TOINT(line.ContigLength),
	contig.gc = TOFLOAT(line.ContigGC), 
	contig.origName = LOWER(line.OriginalContigID)

MERGE (gene: Genes {name: LOWER(line.ImgGeneID)})
SET gene.length = TOINT(line.GeneLength), 
	gene.product = LOWER(line.GeneProduct)

MERGE (locus: Locus {name:LOWER(line.ImgGeneID)})
SET locus.geneStart = TOINT(line.GeneStart), 
	locus.geneStop = TOINT(line.GeneStop),
	locus.geneStrand = TOSTRING(line.GeneStrand),
	locus.geneTaxaPercID = TOFLOAT(line.TaxaPercID)
							
MERGE (gsource: GeneSource {name:LOWER(line.GeneSource)})
MERGE (domain: TaxaDomain {name:LOWER(line.TaxaDomain)})
MERGE (phylum: TaxaPhylum {name:LOWER(line.TaxaPhylum)})
MERGE (class: TaxaClass {name:LOWER(line.TaxaClass)})
MERGE (torder: TaxaOrder {name:LOWER(line.TaxaOrder)})
MERGE (family: TaxaFamily {name:LOWER(line.TaxaFamily)})
MERGE (genus: TaxaGenus {name: LOWER(line.TaxaGenus)})
MERGE (sci: TaxaSciName {name: LOWER(line.Sci_Name)})

// Relationships
MERGE (sci)-[:IN_GENUS]->(genus)
MERGE (genus)-[:IN_FAMILY]->(family)
MERGE (family)-[:IN_ORDER]->(torder)
MERGE (torder)-[:IN_CLASS]->(class)
MERGE (class)-[:IN_PHYLUM]->(phylum)
MERGE (phylum)-[:IN_DOMAIN]->(domain)
MERGE (gene)-[:HAS_SOURCE]->(gsource)
MERGE (gene)-[:PRESENT_ON]->(locus)
MERGE (locus)-[:LOCATED_ON]->(contig)
MERGE (locus)-[:IN_ORGANISM]->(sci)
;
