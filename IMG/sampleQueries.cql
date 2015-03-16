// Find Housekeeping Genes by species
MATCH (s:Source)<-[]-(l:Locus)-[]->(t:Taxa) WHERE s.type="housekeeping" RETURN t.id as Species,collect(DISTINCT s.id) as List, count(DISTINCT s.id) as HK_Count ORDER BY HK_Count DESC;
