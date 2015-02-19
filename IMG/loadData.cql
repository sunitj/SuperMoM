// Load IMG Data

CREATE CONSTRAINT ON (p:ProjectName) ASSERT p.projectName IS UNIQUE;
CREATE CONSTRAINT ON (ic:ImgContigs) ASSERT ic.ImgContigName IS UNIQUE;
CREATE CONSTRAINT ON (coord:GeneCoords) ASSERT coord.geneCoordName IS UNIQUE;
CREATE CONSTRAINT ON (ig:ImgGenes) ASSERT ig.geneImgName IS UNIQUE;
CREATE CONSTRAINT ON (td:TaxaDomain) ASSERT td.taxaDomain IS UNIQUE;
CREATE CONSTRAINT ON (tp:TaxaPhylum) ASSERT tp.taxaPhylum IS UNIQUE;
CREATE CONSTRAINT ON (tc:TaxaClass) ASSERT tc.taxaClass IS UNIQUE;
CREATE CONSTRAINT ON (torder:TaxaOrder) ASSERT torder.taxaOrder IS UNIQUE;
CREATE CONSTRAINT ON (tf:TaxaFamily) ASSERT tf.taxaFamily IS UNIQUE;
CREATE CONSTRAINT ON (tg:TaxaGenus) ASSERT tg.taxaGenus IS UNIQUE;

//CREATE INDEX ON :ImgContigs(ImgContigName);
//CREATE INDEX ON :ImgContigs(originalContigName);
//CREATE INDEX ON :ImgGenes(geneImgName);
//CREATE INDEX ON :TaxaDomain(taxaDomain);
//CREATE INDEX ON :TaxaPhylum(taxaPhylum);
//CREATE INDEX ON :TaxaClass(taxaClass);
//CREATE INDEX ON :TaxaOrder(taxaOrder);
//CREATE INDEX ON :TaxaFamily(taxaFamily);
//CREATE INDEX ON :TaxaGenus(taxaGenus);

LOAD CSV WITH HEADERS FROM 'file:/Users/sunitj/Neo4j/IMG/Staging/random300_testRun2.tsv' AS line FIELDTERMINATOR '\t'
MERGE (ic:ImgContigs {imgContigName:LOWER(line.ImgContigID), contigLength:TOINT(line.ContigLength), contigGC: TOFLOAT(line.ContigGC), originalContigName: LOWER(line.OriginalContigID)})
MERGE (ig: ImgGenes {geneImgName: LOWER(line.ImgGeneID), geneLength: TOINT(line.GeneLength), geneProduct: LOWER(line.GeneProduct), geneTaxaPercID:TOFLOAT(line.TaxaPercID)})
MERGE (coord: GeneCoords {geneCoordName:LOWER(line.ImgGeneID),geneStart:TOINT(line.GeneStart), geneStop:TOINT(line.GeneStop), geneStrand:TOSTRING(line.GeneStrand)})
MERGE (td: TaxaDomain {taxaDomain:LOWER(line.TaxaDomain)})
MERGE (tp: TaxaPhylum {taxaPhylum:LOWER(line.TaxaPhylum)})
MERGE (tc: TaxaClass {taxaClass:LOWER(line.TaxaClass)})
MERGE (tor: TaxaOrder {taxaOrder:LOWER(line.TaxaOrder)})
MERGE (tf: TaxaFamily {taxaFamily:LOWER(line.TaxaFamily)})
MERGE (tg: TaxaGenus {taxaGenus: LOWER(line.TaxaGenus)})
MERGE (sn: SciName {sciName: LOWER(line.Sci_Name)})

// Relationships
MERGE (ig)-[loc:LOCATED_AT]->(coord)-[ann:ANNOTATED_AS]->(ic)
MERGE (td)<-[dp:IN_DOMAIN]-(tp)<-[pc:IN_PHYLUM]-(tc)<-[co:IN_CLASS]-(tor)<-[of:IN_ORDER]-(tf)<-[fg:IN_FAMILY]-(tg)<-[fs:IN_GENUS]-(sn)
//MERGE (coord)-[ctd:IN_DOMAIN]->(td)
//MERGE (coord)-[ctp:IN_PHYLUM]->(tp)
//MERGE (coord)-[ctc:IN_CLASS]->(tc)
//MERGE (coord)-[ctor:IN_ORDER]->(tor)
//MERGE (coord)-[ctf:IN_FAMILY]->(tf)
//MERGE (coord)-[ctg:IN_GENUS]->(tg)
//MERGE (coord)-[csn:IN_SPECIES]->(sn)
MERGE (ic)-[csp:BELONGS_TO]->(sn)
;
