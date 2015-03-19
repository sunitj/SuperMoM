#!/usr/bin/perl

=head1 DESCRIPTION

createDB.pl -- Do this.

=head1 USAGE

perl createDB.pl

=head2 Neo4j specific options

    -port   <INT>   port where neo4j server is running; (default: 7474)

=head2 Script specific options

    -DIR	<CHAR>	Path to uncompressed IMG tarball; (default: "./")
    -project	<CHAR>	File name with multiple key/value pairs OR just an ID number.
    
    -version -v	<BOOLEAN>	version of the current script
    -help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Wed Dec  3 12:54:05 EST 2014)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

## Core Modules ##
use strict;
# use warnings;
use Getopt::Long;
use File::Basename;
use FileHandle;

## External Modules ##
use REST::Neo4p;
use REST::Neo4p::Batch;
use REST::Neo4p::Schema;
use Bio::SeqIO;
use Term::ProgressBar;

my $help;
my $version=fileparse($0).".pl\tv0.9.96";

my $port=7474;
my $project;
my %HKG; # 36 Housekeeping Genes; Cicarelli et al Science 2006
my $DIR="./";
GetOptions(
    'port:i'=>\$port,
    'DIR:s'=>\$DIR,
    'project=s'=>\$project,
    'v|version'=>sub{print $version."\n"; exit;},
    'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

# Connect to the Neo4j Server
my $server='http://127.0.0.1:'.$port;
eval {
    REST::Neo4p->connect($server);
};
ref $@ ? $@->rethrow : die $@ if $@;

# Get all File names in the given directory
unless ($DIR=~ m/\/$/){$DIR=$DIR."/";}
my @FILES=<$DIR*>;

print "[Update] Getting required files from Directory:\t$DIR\n";
my ($cog, $ec, $faa, $fna, $geneProd, $gff, $ko, $mapTxt, $pfam, $phylodist, $config);
foreach my $f(@FILES){
    $cog=$f if ($f=~ /.*.cog.txt/);
    $ec=$f if ($f=~ /.*.ec.txt/);
    $faa=$f if ($f=~ /.*.faa/);
    $fna=$f if ($f=~ /.*.fna/);
    $geneProd=$f if ($f=~ /.*.gene_product.txt/);
    $gff=$f if ($f=~ /.*.gff/);
    $ko=$f if ($f=~ /.*.ko.txt/);
    $mapTxt=$f if ($f=~ /.*.map.txt/);
    $pfam=$f if ($f=~ /.*.pfam.txt/);
    $phylodist=$f if ($f=~ /.*.phylodist.txt/);
}

die "[FATAL] Need project id."if (! $project);


# Read all the files

## Read GFF3 files
print "[Update] Reading GFF file:\t$gff\n";
my (%gff_attributes, %annotation);
my $GFF3=FileHandle->new();
open($GFF3, "<",$gff)||die $!;
while(my $line=<$GFF3>){
    chomp $line;
    $line=lc($line);
    parseGFF3($line);
}
close $GFF3;

## Read the Phylodist File
print "[Update] Reading PhyloDist file:\t$phylodist\n";
my %PhyloDist;
my $PD=FileHandle->new();
open($PD, "<",$phylodist) || die $!;
while(my $line=<$PD>){
    chomp $line;
    $line=lc($line);
    my @stuff=split(/\t/, $line);
    my $locusID=$stuff[0];
    my $perc_id=$stuff[3];
    my(@lineage)=split(/;/, $stuff[4]);
    # lineage = domain;phylum;class;order;family;genus;species;taxon_name
    my $sciName=pop(@lineage);
    my $species=pop(@lineage); # discarded Species Name
    $PhyloDist{$locusID}{"DOMAIN"}=$lineage[0];
    $PhyloDist{$locusID}{"PHYLUM"}=$lineage[1];
    $PhyloDist{$locusID}{"CLASS"}=$lineage[2];
    $PhyloDist{$locusID}{"ORDER"}=$lineage[3];
    $PhyloDist{$locusID}{"FAMILY"}=$lineage[4];
    $PhyloDist{$locusID}{"GENUS"}=$lineage[5];
    $PhyloDist{$locusID}{"SPECIES"}=$sciName;
    $PhyloDist{$locusID}{"PERCENT"}=$perc_id;
    #join("\t", @lineage, $sciName, $perc_id);
}
close $PD;

## Read the Gene_Product file
print "[Update] Reading Gene Product file:\t$geneProd\n";
&loadHouseKeepingGenes;
my %GeneProd;
my $hypoProt=0;
my $GP=FileHandle->new();
open($GP, "<",$geneProd) || die $!;
while(my $line=<$GP>){
    chomp $line;
    $line=lc($line);
    my ($locusID, $product, $source)=split(/\t/, $line);
    $GeneProd{$locusID}{"PRODUCT"}=$product;
    $GeneProd{$locusID}{"SOURCE"}=$source;
    if ($HKG{uc($source)}) {
	$GeneProd{$locusID}{"TYPE"}="housekeeping";
    }
    
}
close $GP;

## Read the map.txt
print "[Update] Reading Names Map file:\t$mapTxt\n";
my %nameMap;
my $MAP=FileHandle->new();
open($MAP, "<",$mapTxt) || die $!;
while(my $line=<$MAP>){
    chomp $line;
    $line=lc($line);
    my ($scaffoldID, $imgScaffoldID)=split(/\t/, $line);
    $nameMap{$imgScaffoldID}=$scaffoldID;
}
close $MAP;

## Calculate the length and GC from scaffold file
print "[Update] Reading scaffold fasta file:\t$fna\n";
print "[Update] Calculating Length and GC content.\n";
my %lenGC;
my $FASTA=FileHandle->new();
open($FASTA, "<",$fna) || die $!;
$/=">";
while(my $line=<$FASTA>){
    chomp $line;
    next unless $line;
    $line=lc($line);
    my ($header, @sequence)=split(/\n/, $line);
    my $seq=join("",@sequence);
    my $seqLen=length($seq);
    
    die $header." Length is 0\n" if ($seqLen==0);
    
    my ($g, $c);
    while ( $seq =~ /g/ig ) { $g++ }
    while ( $seq =~ /c/ig ) { $c++ }
    
    my $GC = (($g+$c)/$seqLen)*100;
    $lenGC{$header}{"LEN"}=$seqLen;
    $lenGC{$header}{"GC"}=sprintf( "%.4f", $GC);

}
close $FASTA;
$/="\n";

print "[Update] All files read!\n";
### NODES
my %project_nodes;
my %contig_nodes;
my %locus_nodes;
my %source_nodes;
my %taxa_nodes;

### RELATIONS
my %contig2project;
my %locus2contig;
my %locus2source;
my %locus2taxa;

### SCHEMA ###


### INDEXES ###
my $idx = REST::Neo4p->get_index_by_name('node', 'my_nodes') ||
	  REST::Neo4p::Index->new('node', 'my_nodes');
my $rel_idx = REST::Neo4p->get_index_by_name('relationship', 'my_relationships') ||
	  REST::Neo4p::Index->new('relationship', 'my_relationships');

my %projectMetadata;
if (-s $project) {
    my $META=FileHandle->new();
    open($META, "<",$project) || die $!;
    while (my $line=<$META>) {
	chomp $line;
	next unless $line;
	next if ($line=~/^#/);
	$line=lc($line);
	my($tag, $value)=split(/\t/, $line);
	$projectMetadata{$tag}=$value;
	$project=$value if ($tag eq "id");
    }
    close $META;
    die "[FATAL] No project id found!\n" if(! $projectMetadata{"id"});
}

($project_nodes{$project})= $idx->find_entries(id=>$project);
unless($project_nodes{$project}){
    $project_nodes{$project}=REST::Neo4p::Node->new({id=>$project});
    $project_nodes{$project}->set_labels("Project");
    $idx->add_entry($project_nodes{$project}, id=>$project);
    if (scalar(keys(%projectMetadata)>0)) {
	foreach (keys %projectMetadata){
	    next if ($_ eq "id");
	    my $value=lc($projectMetadata{$_});
	    my $key=lc($_);
	    $project_nodes{$project}->set_property({$key=>$value});
	}
    }
    else{
	$project_nodes{$project}->set_property({name=>"Another Nameless Project",
						type=>"Metagenome"});
    }
}
###

# When did the uploading start
my $BEGIN=FileHandle->new();
open($BEGIN, ">","BEGIN") || die $!;
close $BEGIN;

# Create a Progress bar
my $numContigs=keys(%annotation);
my $progress = Term::ProgressBar->new ({count => $numContigs ,name => 'Populating Database'});
my $currentContig=0;

# Create a Log file
my $logFile=$project.".neo4j.loading.log";
my $LOG=FileHandle->new();
open($LOG, ">",$logFile) || die $!;
print $LOG "\# $version\n";

## Populate the database and create a tabulr file at the same time.
foreach my $contig(keys %annotation){
    next unless ($lenGC{$contig}{"LEN"});
    print $LOG "Contig:\t$contig\n";
    ($contig_nodes{$contig})= $idx->find_entries(id=>$contig);
    unless($contig_nodes{$contig}){
        $contig_nodes{$contig}=REST::Neo4p::Node->new({id=>$contig,
						       length=>$lenGC{$contig}{"LEN"},
						       gc=>$lenGC{$contig}{"GC"},
						       name=>$nameMap{$contig}
						       });
        $contig_nodes{$contig}->set_labels("Scaffolds");
#	$contig_nodes{$contig}->set_property({blah=>$blah}) if ($blah);
        $idx->add_entry($contig_nodes{$contig}, id=>$contig);
    }
    
    # RELATION: (Scaffolds)-[:BELONGS_TO]->(Project)
    $contig2project{$project}{$contig}=$rel_idx->find_entries(id=>$contig."_".$project);
    unless($contig2project{$project}{$contig}){
        $contig2project{$project}{$contig}=$contig_nodes{$contig}->relate_to($project_nodes{$project}, "BELONGS_TO");
	$rel_idx->add_entry($contig2project{$project}{$contig}, id=>$contig."_".$project)
    }
    print $LOG "\t\t\tREL:\t($contig)-[:BELONGS_TO]->($project)\n";
    
    foreach my $gene(keys %{$annotation{$contig}}){
	print $LOG "\tGene:\t$gene\n";
	($locus_nodes{$gene})= $idx->find_entries(id=>"Locus_".$gene);
	unless($locus_nodes{$gene}){
	    $locus_nodes{$gene}=REST::Neo4p::Node->new({id=>$gene});
	    $locus_nodes{$gene}->set_labels("Locus");
	    if ($annotation{$contig}{$gene}{"START"}){
		$locus_nodes{$gene}->set_property({length=>$annotation{$contig}{$gene}{"LEN"}});
	    }
	    $idx->add_entry($locus_nodes{$gene}, name=>"Locus_".$gene);
	}
	print $LOG "\t\t\tLocus:\t$gene\n";

	if ($annotation{$contig}{$gene}{"START"}) {
	    # RELATION: (Locus)-[:LOCATED_ON]->(Contig)
	    $locus2contig{$contig}{$gene}=$rel_idx->find_entries(id=>$gene."_".$contig);
	    unless($locus2contig{$contig}{$gene}){
		$locus2contig{$contig}{$gene}=$locus_nodes{$gene}->relate_to($contig_nodes{$contig}, "LOCATED_ON");
		foreach (keys %{$annotation{$contig}{$gene}}){
		    next if $_ eq "LEN";
		    my $value=lc($annotation{$contig}{$gene}{$_});
		    my $key=lc($_);
		    $locus2contig{$contig}{$gene}->set_property({$key=>$value});
		}
		$rel_idx->add_entry($locus2contig{$contig}{$gene}, id=>$gene."_".$contig);
	    }
	}
	else{
	    $annotation{$contig}{$gene}{"START"}=0;
	    $annotation{$contig}{$gene}{"STOP"}=0;
	    $annotation{$contig}{$gene}{"TYPE"}=0;
	    $annotation{$contig}{$gene}{"LEN"}=0;
	    $annotation{$contig}{$gene}{"STRAND"}=0;
	}
	print $LOG "\t\t\tREL:\t($gene)-[:LOCATED_ON]->($contig)\n";
	
	my $source=$GeneProd{$gene}{"SOURCE"};
	if ($source) {
	    ($source_nodes{$gene})= $idx->find_entries(id=>"Source_".$source);
	    unless($source_nodes{$gene}){
		$source_nodes{$gene}=REST::Neo4p::Node->new({id=>$source});
		$source_nodes{$gene}->set_labels("Source");
		foreach (keys %{$GeneProd{$gene}}){
		    next if $_ eq "SOURCE";
		    my $value=lc($GeneProd{$gene}{$_});
		    my $key=lc($_);
		    $source_nodes{$gene}->set_property({$key=>$value});
		}
		$idx->add_entry($source_nodes{$gene}, id=>"Source_".$source);
	    }
	    print $LOG "\t\t\tSource:\t$source\n";
	    
	    # RELATION: (Locus)-[:HAS_SOURCE]->(Source)
	    $locus2source{$source}{$gene}=$rel_idx->find_entries(id=>$source."_".$gene);
	    unless($locus2source{$source}{$gene}){
		$locus2source{$source}{$gene}=$locus_nodes{$gene}->relate_to($source_nodes{$gene}, "HAS_SOURCE");
		$rel_idx->add_entry($locus2source{$source}{$gene}, id=>$source."_".$gene);
	    }
	    print $LOG "\t\t\tREL:\t($gene)-[:HAS_SOURCE]->($source)\n";
	}
	else{
	    $source="Unknown"
	}
	
	my $product;
	unless($GeneProd{$gene}{"PRODUCT"}){
	    $product="Unknown";
	}
	
	if ($PhyloDist{$gene}{"DOMAIN"}){
	    my $species=$PhyloDist{$gene}{"SPECIES"};
	    ($taxa_nodes{$gene})= $idx->find_entries(id=>$species);
	    unless($taxa_nodes{$gene}){
		$taxa_nodes{$gene}=REST::Neo4p::Node->new({id=>$species});
		$taxa_nodes{$gene}->set_labels("Taxa");
		foreach (keys %{$PhyloDist{$gene}}){
		    next if $_ eq "SPECIES";
		    next if $_ eq "PERCENT";
		    my $value=lc($PhyloDist{$gene}{$_});
		    my $key=lc($_);
		    $taxa_nodes{$gene}->set_property({$key=>$value});
		}
		$idx->add_entry($taxa_nodes{$gene}, id=>$species);
	    }
	    print $LOG "\t\t\tTaxa:\t$species\n";
	    
	    # RELATION (Locus)-[:IN_ORGANISM]->(Taxa)
	    $locus2taxa{$species}{$gene}=$rel_idx->find_entries(id=>$species."_".$gene);
	    unless($locus2taxa{$species}{$gene}){
		$locus2taxa{$species}{$gene}=$locus_nodes{$gene}->relate_to($taxa_nodes{$gene}, "IN_ORGANISM");
		$locus2taxa{$species}{$gene}->set_property({identity=>$PhyloDist{$gene}{"PERCENT"}});
		$rel_idx->add_entry($locus2taxa{$species}{$gene}, id=>$species."_".$gene);
	    }
	    print $LOG "\t\t\tREL:\t($gene)-[:IN_ORGANISM]->($species)\n";
	    
	}
	else{
	    foreach (keys %{$PhyloDist{$gene}}){
		next if $_ eq "PERCENT";
		$PhyloDist{$gene}{$_}="Unknown";
	    }
	    $PhyloDist{$gene}{"PERCENT"}=0;
	}
	
#	print $OUT $printLine."\n";
    }


    # Update Progress bar
    $currentContig++;
    $progress->update($currentContig);

}
close $LOG;
# close $OUT;

# When did the loading end
my $END=FileHandle->new();
open($END, ">","END") || die $!;
close $END;


##############################
######   Sub-Routines   ######
##############################
sub trim{
    my $line=shift;
    chomp $line;
    $line=~ s/^\s+//;
    $line=~ s/\s+$//;
    return $line;
}

sub parseGFF3{
#http://gmod.org/wiki/GFF
# contig, source, type, start,stop,score,strand, phase,attributes
    my $line=shift;
    my ($contig, $source, $type, $start,$stop,$score,$strand, $phase,$attribs)=split(/\t/, $line);
#    $contig=$nameMap{$contig};
    my(@attributes)=split(/\;/, $attribs);
     my ($locusID, $ID, $Name,$Alias, $Parent, $Target, $Gap, $Derives_from, $Note, $Dbxref, $Onto, $repeat_type, $repeat_unit, $repeat_fam, $product);
    foreach my $att(@attributes){
	if ($att=~/^id\=(.*)/){
	    $ID= $1;
	}
        if ($att=~/locus_tag\=(.*)/){
            $locusID=$1;
            next;
    	}
    }
    if ($locusID){
        foreach my $att(@attributes){
            next if ($att=~/^id/);
            next if ($att=~/^locus_tag/);
            my($tag, $data)=split(/\=/,$att);
	    $gff_attributes{$tag}++;
	    $annotation{$contig}{$locusID}{$tag}=$data;
        }
    }
    else{
	foreach my $att(@attributes){
	    if ($Parent){
		$locusID=$Parent."__exon"
	    }
	    elsif($type=~/repeat/){
		$locusID=$ID."__".
		($repeat_type ? $repeat_type : "Unknown")."__".
		($repeat_unit ? $repeat_unit : "Unknown")."__"; # rpt_type=CRISPR;rpt_unit=13023..13055;rpt_family=blah
            }
            elsif($ID){
		$locusID=$ID."__".$type;
            }
	    else{
		$locusID="WTF__".$type;
	    }
        }
    }
    
    $strand=0 if (($strand != 1) && ($strand != -1));
    my $geneLen=abs($stop-$start)+1;
    $annotation{$contig}{$locusID}{"START"}=$start;
    $annotation{$contig}{$locusID}{"STOP"}=$stop;
    $annotation{$contig}{$locusID}{"TYPE"}=$type;
    $annotation{$contig}{$locusID}{"LEN"}=$geneLen;
    $annotation{$contig}{$locusID}{"STRAND"}=$strand;
    #$annotation{$contig}{$locusID}=join("\t", $start, $stop, $geneLen,$strand);
    
    return;
}

sub loadHouseKeepingGenes{
%HKG=(
"COG0080"=>"L11",
"COG0081"=>"L1",
"COG0087"=>"L3",
"COG0091"=>"L22",
"COG0093"=>"L14",
"COG0094"=>"L5",
"COG0097"=>"L6P/L9E",
"COG0102"=>"L13",
"COG0197"=>"L16/L10E",
"COG0200"=>"L15",
"COG0256"=>"L18",
"COG0048"=>"S12",
"COG0049"=>"S7",
"COG0052"=>"S2",
"COG0092"=>"S3",
"COG0096"=>"S8",
"COG0098"=>"S5",
"COG0099"=>"S13",
"COG0100"=>"S11",
"COG0103"=>"S9",
"COG0184"=>"S15P/S13E",
"COG0186"=>"S17",
"COG0522"=>"S4",
"COG0016"=>"Phenylalanyl-tRNA synthethase alpha subunit",
"COG0018"=>"Arginyl-tRNA synthetase",
"COG0060"=>"Isoleucyl-tRNA synthetase",
"COG0124"=>"Histidyl-tRNA synthetase",
"COG0143"=>"Methionyl-tRNA synthetase",
"COG0172"=>"Seryl-tRNA synthetase",
"COG0201"=>"Preprotein translocase subunit SecY",
"COG0495"=>"Leucyl-tRNA synthetase",
"COG0525"=>"Valyl-tRNA synthetase",
"COG0202"=>"DNA-directed RNA polymerase, alpha subunit/40 kD subunit",
"COG0085"=>"DNA-directed RNA polymerase, beta subunit/140 kD subunit",
"COG0012"=>"Predicted GTPase",
"COG0533"=>"Metal-dependent protease"
);
}
