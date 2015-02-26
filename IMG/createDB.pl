#!/usr/bin/perl

=head1 DESCRIPTION

createDB.pl -- Do this.

=head1 USAGE

perl createDB.pl

=head2 Options
    -gff    <CHAR>  GFF file
    -phylodist  <CHAR>  Phylodist file  
    -port   <INT>   port where neo4j server is running.
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
use warnings;
use Getopt::Long;
use File::Basename;
use FileHandle;

## External Modules ##
use REST::Neo4p;
use REST::Neo4p::Batch;
use Bio::SeqIO;

my $help;
my $version=fileparse($0).".pl\tv0.0.1b";

my $port=7474;
my ($gff, $phylodist, $geneProd, $mapTxt, $crisprTxt, $out, $lgc);
GetOptions(
    'port:i'=>\$port,
    'gff=s'=>\$gff,
    'gene_prod=s'=>\$geneProd,
    'phylodist=s'=>\$phylodist,
    'map=s'=>\$mapTxt,
    'crispr=s'=>\$crisprTxt,
    'out=s'=>\$out,
    'lgc=s'=>\$lgc,
    'v|version'=>sub{print $version."\n"; exit;},
    'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

# Connect to the Neo4j Server
#my $server='http://127.0.0.1:'.$port;
#eval {
#    REST::Neo4p->connect($server);
#};
#ref $@ ? $@->rethrow : die $@ if $@;

# Read all the files
## Read GFF3 files
my (%gff_attributes, %annotation);
my $GFF3=FileHandle->new();
open($GFF3, "<",$gff)||die $!;
while(my $line=<$GFF3>){
    chomp $line;
    parseGFF3($line);
}
close $GFF3;

## Read the Phylodist File
my %PhyloDist;
my $PD=FileHandle->new();
open($PD, "<",$phylodist) || die $!;
while(my $line=<$PD>){
    chomp $line;
    my @stuff=split(/\t/, $line);
    my $locusID=$stuff[0];
    my $perc_id=$stuff[3];
    my(@lineage)=split(/;/, $stuff[4]);
    # lineage = domain;phylum;class;order;family;genus;species;taxon_name
    my $sciName=pop(@lineage);
    my $species=pop(@lineage); # discarded Species Name
    $PhyloDist{$locusID}=join("\t", @lineage, $sciName, $perc_id);
}
close $PD;

## Read the Gene_Product file
my %GeneProd;
my $GP=FileHandle->new();
open($GP, "<",$geneProd) || die $!;
while(my $line=<$GP>){
    chomp $line;
    my ($locusID, $product, $source)=split(/\t/, $line);
    $GeneProd{$locusID}=join("\t", $product, $source);
}
close $GP;

## Read the map.txt
my %nameMap;
my $MAP=FileHandle->new();
open($MAP, "<",$mapTxt) || die $!;
while(my $line=<$MAP>){
    chomp $line;
    my ($scaffoldID, $imgScaffoldID)=split(/\t/, $line);
    $nameMap{$imgScaffoldID}=$scaffoldID;
}
close $MAP;

### Read the crispr.txt
#my %crispr;
#my $CRISP=FileHandle->new();
#open($CRISP, "<",$crisprTxt) || die $!;
#while(my $line=<$CRISP>){
#    chomp $line;
#    my ($imgScaffoldID, $crisprNum, $pos, $repeat, $spacer, $tool)=split(/\t/, $line);
#    $crispr{$imgScaffoldID}=join("\t", $repeat, $spacer);
#}
#close $CRISP;

## Read the length+GC output
my %lenGC;
my $LGC=FileHandle->new();
open($LGC, "<",$lgc) || die $!;
while(my $line=<$LGC>){
    chomp $line;
    my ($imgScaffoldID, $gc, $length)=split(/\t/, $line);
    $lenGC{$imgScaffoldID}=join("\t", $length, $gc);
}
close $LGC;

# Write data to a TSV for batch load to database using `LOAD CSV` from the neo4j shell
# Until I can figure out how to do it from here.
my $OUT=FileHandle->new();
open($OUT, ">",$out) || die $!;

# Header
print $OUT "ImgContigID\tOriginalContigID\tContigLength\tContigGC\t";
print $OUT "ImgGeneID\tGeneStart\tGeneStop\tGeneLength\tGeneStrand\t";
print $OUT "GeneProduct\tGeneSource\t";
print $OUT "TaxaDomain\tTaxaPhylum\tTaxaClass\tTaxaOrder\tTaxaFamily\tTaxaGenus\tSci_Name\tTaxaPercID\t";
# print $OUT "Repeat_Seq\tSpacer_Seq";
print $OUT "\n";

foreach my $contig(keys %annotation){
    foreach my $gene(keys %{$annotation{$contig}}){
	
	unless($lenGC{$contig}){
	    $lenGC{$contig}=join("\t", qw(0 0));
	}
	
	unless($annotation{$contig}{$gene}){
	    $annotation{$contig}{$gene}=join("\t", qw(0 0 0 .));
	}
	
	unless ($GeneProd{$gene}){
	    $GeneProd{$gene}=join("\t", qw(Unknown Unknown));
	}
	
	unless ($PhyloDist{$gene}){
	    $PhyloDist{$gene}=join("\t", (map 'Unknown',(1..7)),"0");
	}
	
	# Header: "IMG_contig_ID\tOriginal_contig_ID\tContig_Length\tContig_GC\t
	# IMG_Gene_ID\tGene_Start\tGene_Stop\tGene_Length\tGene_Strand\t
	# Gene_Product\tGene_Source\t
	# Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tTaxa_Perc_ID\t
	## Repeat_Seq\tSpacer_Seq\n";
	my $printLine = $contig."\t".$nameMap{$contig}."\t".$lenGC{$contig}."\t".$gene."\t";
	$printLine .= $annotation{$contig}{$gene}."\t";
	$printLine.=$GeneProd{$gene}."\t";
	$printLine.=$PhyloDist{$gene}."\t";
#	print $OUT ($crispr{$contig} ? $crispr{$contig} : "\t\t" );
#	$printLine=~ s/\s+$//g;
	print $OUT $printLine."\n";
    }
}
close $OUT;

## Loading the Data sequentially
#my $contig_names=REST::Neo4p::Index->new('node','contig_names');
#my $gene_names=REST::Neo4p::Index->new('node','gene_names');


# Load the data in a batch
# Ref: http://www.slideshare.net/majensen1/dcpm-meetup (Slide 12)
#batch{   
#}

##############################
######   Sub-Routines   ######
##############################
sub parseGFF3{
#http://gmod.org/wiki/GFF
# contig, source, type, start,stop,score,strand, phase,attributes
    my $line=shift;
    my ($contig, $source, $type, $start,$stop,$score,$strand, $phase,$attribs)=split(/\t/, $line);
#    $contig=$nameMap{$contig};
    my(@attributes)=split(/\;/, $attribs);

     my ($locusID, $ID, $Name,$Alias, $Parent, $Target, $Gap, $Derives_from, $Note, $Dbxref, $Onto, $repeat_type, $repeat_unit, $repeat_fam, $product);
    foreach my $att(@attributes){
	if ($att=~/^ID\=(.*)/){
	    $ID= $1;
	    next;
	}
        if ($att=~/locus_tag\=(.*)/){
            $locusID=$1;
            next;
    	}
    }
    if ($locusID){
        foreach my $att(@attributes){
            next if ($att=~/^ID/);
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
            else{
		$locusID=$ID."__".$type;
            }
        }
    }
    #$annotation{$contig}{$locusID}{"START"}=$start;
    #$annotation{$contig}{$locusID}{"STOP"}=$stop;
    #$annotation{$contig}{$locusID}{"TYPE"}=$type;
    #$annotation{$contig}{$locusID}{"LEN"}=(abs($stop-$start));
    #$annotation{$contig}{$locusID}{"STRAND"}=$strand;
    my $geneLen=abs($stop-$start);
    $annotation{$contig}{$locusID}=join("\t", $start, $stop, $geneLen,$strand);
    
    return;
}
