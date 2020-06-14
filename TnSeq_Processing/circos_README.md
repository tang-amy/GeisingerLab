# Data Visualization Using Circos
## Installation
This user guide is based on the documentation provided by circos: http://circos.ca/documentation/tutorials/configuration/installation/.
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
1. To install Circos, first make a directory where you want the software to be installed. For example, if you want to creat a folder in the home directory, execute the following in terminal: 
```
cd ~
mkdir circos
```
2. Download the latest version of Circos on: http://circos.ca/software/download/. Move the .tgz file to the folder created in step 1.

3. De

Configuration (.conf files): specifies the features of the graph (cytogenetic bands, labels, ticks and, of course, data)

Karyotype:  defines the names, sizes and colors of chromosomes that you will use in the image.




