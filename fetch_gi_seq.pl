#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use threads;
use threads::shared;
use WWW::Mechanize;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $gi_file          = '';
my $length_threshold = 10;
my $outfile          = undef;
my $failed_file      = undef;
my $parallel         = 5;
my $batch_number     = 100;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'length=i'   => \$length_threshold,
    'infile=s'   => \$gi_file,
    'outfile=s'  => \$outfile,
    'failed=s'   => \$failed_file,
    'parallel=i' => \$parallel,
    'batch=i'    => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

unless ( defined $outfile ) {
    $outfile = $gi_file;
    $outfile =~ s/\.\w+?$/\.fasta/;
    print "Write to $outfile\n";
}
unless ( defined $failed_file ) {
    $failed_file = $gi_file;
    $failed_file =~ s/\.\w+?$/\.failed/;
}

my $http_proxy = 'http://127.0.0.1:8088/';

#----------------------------------------------------------#
# Download
#----------------------------------------------------------#
$|++;

# get gi file
my @gids;
{
    open my $infh, "<", $gi_file;
    my $file_content = do { local $/; <$infh> };
    close $infh;
    @gids = $file_content =~ /gi\:(\d+)/g;
}

while ( @gids > $batch_number ) {
    my @working = splice @gids, 0, $batch_number;

    # download $batch_number gi and write then to $outfile
    # $batch_number * 5 is retrying number
    manager( \@working, $batch_number * 5 );
}
manager( \@gids );

exit;

#----------------------------------------------------------#
# Subroutine
#----------------------------------------------------------#
sub manager {
    my $gids  = shift;
    my $retry = shift;

    my @working_gids;
    share(@working_gids);
    @working_gids = @$gids;
    my $content : shared;

    my $download_fasta = sub {
        my $thread_id   = shift;
        my $thread_name = "thread $thread_id: ";

        my $mech = WWW::Mechanize->new;
        $mech->proxy( [ 'http', 'ftp' ], $http_proxy ) if $http_proxy;
        $mech->stack_depth(0);    # no history to save memory

        my $url_part_1 = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?";
        my $url_part_2 = "dopt=fasta&dispmax=5&sendto=t&db=nucleotide";
        my $url_part_3 = "&qty=1&c_start=1&list_uids=";
        my $url_part_4 = "&from=begin&to=end";

        my $error = 0;
    GET: while (1) {
            my $gi;
            {
                lock @working_gids;
                last GET if @working_gids == 0;
                $gi = shift @working_gids;
            }

            my $url = join '',
                ( $url_part_1, $url_part_2, $url_part_3, $gi, $url_part_4 );

            print $thread_name, "getting gi: $gi\n";
            $mech->get($url);

            my $html_content = $mech->content;
            if ( $html_content =~ /^>.+\n\n$/s ) {
                $error = 0;
                if ( length $html_content > $length_threshold ) {
                    lock $content;
                    $content .= $html_content;
                    print $thread_name, "Done.\n";
                }
                else {
                    print $thread_name, "Seq too short. Jump to next\n";
                }
                redo GET;
            }
            else {
                $error++;
                lock @working_gids;
                unshift @working_gids, $gi;
                if ( $error < $retry ) {
                    print $thread_name, "Broken file downloaded. Retry...\n";
                    redo GET;
                }
                else {
                    print $thread_name, "Connection lost.\n";
                    last GET;
                }
            }
        }

        return;
    };

    my @workers;
    for ( 1 .. $parallel ) {
        $workers[$_] = threads->new( $download_fasta, $_ );
    }
    for ( 1 .. $parallel ) {
        $workers[$_]->join;
    }

    open my $outfh, ">>", $outfile
        or die "Cannot append to file: $!";
    print {$outfh} $content;
    close $outfh;

    if ( @working_gids > 0 ) {
        open my $failed_fh, ">>", $failed_file
            or die "Cannot append to file: $!";
        print {$failed_fh} join '', map {"gi:$_\n"} @working_gids;
        close $failed_fh;
    }

    return;
}

__END__

=head1 NAME

    fetch_gi_seq.pl - Fetch fasta file of a gi from NCBI

=head1 SYNOPSIS

    fetch_gi_seq.pl [options]
        Options:
            --help              brief help message
            --man               full documentation
            --length            length threshold
            --infile            gi filename
            --outfile            output dir of fasta files
            --length            length threshold
            --parallel          run in parallel mode
            --retry             retry number

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
