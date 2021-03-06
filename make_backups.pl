#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Class;
use Archive::Zip qw( :ERROR_CODES :CONSTANTS );

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

my $backup_dir = 'F:\Software\AppData';

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'd|backup_dir' => \$backup_dir,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# directories to be backuped
#----------------------------------------------------------#
my %backup_of = (

    strawberry => {
        dir    => [ dir('C:\strawberry'), ],
        action => sub {
            print "Clean useless files in Perl directory...\n";
        },
        delete => [
            dir('C:\strawberry\cpan\build'),
            dir('C:\strawberry\cpan\sources\authors\id'),
            dir('C:\strawberry\perl\html'),
            dir('C:\strawberry\perl\man'),
        ],
    },

    Mozilla     => { dir => [ dir( $ENV{APPDATA}, 'Mozilla' ) ], },
    ActiveState => { dir => [ dir( $ENV{APPDATA}, 'ActiveState' ) ], },

    #VanDyke     => { dir => [ dir( $ENV{APPDATA}, 'VanDyke' ) ], },
    Launchy       => { dir => [ dir( $ENV{APPDATA}, 'Launchy' ) ], },
    #BeyondCompare => { dir => [ dir( $ENV{APPDATA}, 'Scooter Software' ) ], },
    #XShell        => { dir => [ dir( $ENV{APPDATA}, 'NetSarang' ) ], },

    Perforce => { dir => [ dir('D:\Tools\Perforce') ], },
    FlashFXP => { dir => [ dir('D:\Tools\FlashFXP') ], },

    #freeime  => { dir => [ dir('D:\Tools\freeime') ], },
    Scripts => { dir => [ dir('E:\Scripts') ], },
    lab => { dir => [ dir('E:\lab') ], },
    #EverNote => {
    #    dir =>
    #        [ dir( $ENV{USERPROFILE}, 'Documents', 'My EverNote Files' ), ],
    #},
    #SlickEdit => {
    #    dir =>
    #        [ dir( $ENV{USERPROFILE}, 'Documents', 'My SlickEdit Config' ), ],
    #},
    StartMenu => {
        dir => [
            dir( $ENV{ProgramData}, 'Microsoft', 'Windows', 'Start Menu' ),
        ],
    },
    ultraedit => {
        dir => [
            dir( $ENV{APPDATA}, 'IDMComp' ),
            dir('D:\Program Files\IDM Computer Solutions\UltraEdit'),
        ],
    },
    #bioinfo => {
    #    dir => [
    #        dir('D:\Tools\clustalw1.83.XP'), dir('D:\Tools\clustalx1.81'),
    #        dir('D:\Tools\HYPHY'),           dir('D:\Tools\muscle'),
    #        dir('D:\Tools\paml'),            dir('D:\Tools\PAUP'),
    #        dir('D:\Tools\phylip'),          dir('D:\Tools\Primer3'),
    #        dir('D:\Tools\ProSeq'),          dir('D:\Tools\readseq'),
    #        dir('D:\Tools\t_coffee'),
    #    ],
    #},
    NetProxy => {
        dir => [
            dir('D:\Tools\CCProxy'),  dir('D:\Tools\FreeCap'),
            dir('D:\Tools\mproxy12'), dir('D:\Tools\QQWry'),
        ],
        file => [ file('D:\proxy.pac') ],
        reg  => {
            freecap      => q{"HKEY_CURRENT_USER\Software\Bert's Software"},
            fc_uninstall => q{"HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\Windows}
                . q{\CurrentVersion\Uninstall\FreeCap_is1"},
        },
    },
);

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
BACKUP: for my $key ( keys %backup_of ) {
    my $zipped_file = file( $backup_dir, "$key.zip" );
    if ( -e $zipped_file ) {
        print "=== $zipped_file already exists, skip it\n\n";
        next BACKUP;
    }

    print "=== Start backuping [$key] as $zipped_file\n";

    my $zip = Archive::Zip->new;
    my @temp_files;

    # execute actions before backup
    if ( exists $backup_of{$key}->{action} ) {
        print "Execute actions of [$key]\n";
        $backup_of{$key}->{action}->();
    }

    # export registry entry to .reg files
    if ( exists $backup_of{$key}->{reg} ) {
        print "Export registry of [$key]\n";
        my $reg = $backup_of{$key}->{reg};
        for my $regkey ( keys %$reg ) {
            my $reg_file    = "$regkey.reg";
            my $regkey_path = $reg->{$regkey};
            my $cmd         = "regedit /e $reg_file $regkey_path";
            print $cmd, "\n";
            system $cmd;

            push @temp_files, $reg_file;
            $zip->addFile($reg_file);
        }
    }

    # delete directories
    if ( exists $backup_of{$key}->{delete} ) {
        for my $del ( @{ $backup_of{$key}->{delete} } ) {
            next if !-e $del;
            if ( !-d $del ) {
                warn "$del is not a directory!\n";
                warn "Stop backuping [$key]\n";
                next BACKUP;
            }
            print "Delete directory tree $del\n";
            $del->rmtree;
        }
    }

    # backup directories
    if ( exists $backup_of{$key}->{dir} ) {
        for my $dir ( @{ $backup_of{$key}->{dir} } ) {
            if ( !-d $dir ) {
                warn "$dir is not a directory!\n";
                warn "Stop backuping [$key]\n";
                next BACKUP;
            }
            print "Add directory tree $dir to archive\n";
            $zip->addTree( "$dir", $dir->dir_list(-1) );    # stringify
        }
    }

    # backup files
    if ( exists $backup_of{$key}->{file} ) {
        for my $file ( @{ $backup_of{$key}->{file} } ) {
            if ( !-e $file ) {
                warn "$file does not exist!\n";
                warn "Stop backuping [$key]\n";
                next BACKUP;
            }
            print "Add file $file to archive\n";
            $zip->addFile( "$file", $file->basename );    # stringify
        }
    }

    # comments attached to the archive
    $backup_of{$key}->{time} = scalar localtime;
    $zip->zipfileComment( "Backup comments:\n" . Dump $backup_of{$key} );

    # write backup file to disks
    unless ( $zip->writeToFileNamed("$zipped_file") == AZ_OK ) {
        die "write error\n";
    }

    # clear temporary files
    for (@temp_files) {
        unlink $_;
    }

    print "=== Finish backuping [$key]\n\n";
}

exit;

__END__

=head1 NAME

make_backups.pl - Backup directories and files in your machine

=head1 SYNOPSIS

perl make_backups.pl [options] [file ...]
  Options:
    --help              brief help message
    --man               full documentation
    -d, --backup_dir    where to store backup files

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
