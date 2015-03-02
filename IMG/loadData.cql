// Load IMG Data

CREATE CONSTRAINT ON (project:Project) ASSERT project.projectID IS UNIQUE;
CREATE CONSTRAINT ON (contig:ImgContigs) ASSERT contig.name IS UNIQUE;
CREATE CONSTRAINT ON (gene:ImgGenes) ASSERT gene.name IS UNIQUE;
CREATE CONSTRAINT ON (coord:GeneCoords) ASSERT coord.locus IS UNIQUE;
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
SET gene.length = TOINT(line.GeneLength) WHERE gene.length IS NOT NULL, 
	gene.product = LOWER(line.GeneProduct) WHERE gene.product IS NOT NULL

MERGE (coord: Coords {locus:LOWER(line.ImgGeneID)})
SET coord.geneStart = TOINT(line.GeneStart) WHERE coord.geneStart IS NOT NULL, 
	coord.geneStop = TOINT(line.GeneStop)WHERE coord.geneStop IS NOT NULL,
	coord.geneStrand = TOSTRING(line.GeneStrand) WHERE coord.geneStrand IS NOT NULL,
	coord.geneTaxaPercID = TOFLOAT(line.TaxaPercID) WHERE coord.TaxaPercID IS NOT NULL
							
MERGE (gsource: GeneSource {name:LOWER(line.GeneSource)} WHERE gsource IS NOT NULL)
MERGE (domain: TaxaDomain {name:LOWER(line.TaxaDomain)} WHERE domain IS NOT NULL)
MERGE (phylum: TaxaPhylum {name:LOWER(line.TaxaPhylum)} WHERE phylum IS NOT NULL)
MERGE (class: TaxaClass {name:LOWER(line.TaxaClass)} WHERE class IS NOT NULL)
MERGE (torder: TaxaOrder {name:LOWER(line.TaxaOrder)} WHERE torder IS NOT NULL)
MERGE (family: TaxaFamily {name:LOWER(line.TaxaFamily)} WHERE family IS NOT NULL)
MERGE (genus: TaxaGenus {name: LOWER(line.TaxaGenus)} WHERE genus IS NOT NULL)
MERGE (sci: TaxaSciName {name: LOWER(line.Sci_Name)} WHERE sci IS NOT NULL)

// Relationships
MERGE (sci)-[:IN_GENUS]->(genus) WHERE sci OR genus IS NOT NULL
MERGE (genus)-[:IN_FAMILY]->(family) WHERE family OR genus IS NOT NULL
MERGE (family)-[:IN_ORDER]->(torder) WHERE family OR torder IS NOT NULL
MERGE (torder)-[:IN_CLASS]->(class) WHERE class OR torder IS NOT NULL
MERGE (class)-[:IN_PHYLUM]->(phylum) WHERE class OR phylum IS NOT NULL
MERGE (phylum)-[:IN_DOMAIN]->(domain) WHERE domain OR phylum IS NOT NULL
MERGE (gene)-[:HAS_SOURCE]->(gsource) WHERE gene OR gsource IS NOT NULL
MERGE (coord)-[:ANNOTATED_AS]->(gene) WHERE gene OR coord IS NOT NULL
MERGE (coord)-[:LOCATED_ON]->(contig) WHERE coord IS NOT NULL
MERGE (coord)-[:IN_ORGANISM]->(sci) WHERE sci OR coord IS NOT NULL
;
