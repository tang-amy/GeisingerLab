# Data Visualization Using Circos
## Installation
This user guide is based on the documentation provided by circos: [http://circos.ca/documentation/tutorials/configuration/installation/].
However, the official documentation is a bit complicated and hard to follow, which is why I wrote this user guide.
### Installing Circos
Circos is a command line tool (e.g. no GUI) and is perl based. It would be easier to run it on a mac computer (or linux). If you are using Mac OS, perl is most likely pre-installed. 

To check if perl is installed, type the following command in the terminal:
```
which perl
```
If the terminal returns something like: 
```
/usr/bin/perl
```
It means perl is already installed.

To check the version of perl:
```
perl --version
```
Example output:
```
This is perl 5, version 18, subversion 4 (v5.18.4) built for darwin-thread-multi-2level
(with 2 registered patches, see perl -V for more detail)

Copyright 1987-2013, Larry Wall

Perl may be copied only under the terms of either the Artistic License or the
GNU General Public License, which may be found in the Perl 5 source kit.

Complete documentation for Perl, including FAQ lists, should be found on
this system using "man perl" or "perldoc perl".  If you have access to the
Internet, point your browser at http://www.perl.org/, the Perl Home Page.
```
To install Circos, first make a directory where you want the software to be installed.
Download the latest version of Circos on: http://circos.ca/software/download/.
Install circos following instructions on:
Circos Tutorials: Configuration and Installation - Distribution and Installation

http://circos.ca/documentation/tutorials/configuration/perl_and_modules/index.mhtml#circos-libgd-gd

Configuration (.conf files): specifies the features of the graph (cytogenetic bands, labels, ticks and, of course, data)

Karyotype:  defines the names, sizes and colors of chromosomes that you will use in the image.




