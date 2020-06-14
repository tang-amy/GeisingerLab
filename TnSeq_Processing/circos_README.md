# Data Visualization Using Circos
## Installation
This user guide is based on the documentation provided by circos: http://circos.ca/documentation/tutorials/configuration/installation/.
However, the official documentation is a bit complicated and hard to follow, which is why I wrote this user guide.
### Installing Circos
Circos is a command line tool (e.g. no GUI) and is perl based. It would be easier to run it on a mac computer (or linux). If you are using Mac OS, perl is most likely pre-installed. 

To check if perl is installed, type the following command in the terminal:
```bash
which perl
```
If the terminal returns something like: 
```bash
/usr/bin/perl
```
It means perl is already installed.

To check the version of perl:
```bash
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
```bash
# "~" stands for home directory, you can install Circos to anywhere you want if you know what you are doing.
cd ~
mkdir circos
```
2. Download the latest version of Circos on: http://circos.ca/software/download/. Move the .tgz file to the folder created in step 1.

3. Unpack the .tgz file (named circos-0.67-pre4.tgz for example):
```bash
# first change working directory to circos
cd ~/circos
# unpack file
tar -xzvf circos-0.67-pre4.tgz
# now list files in the current directory, you should be able to see a folder named circos-0.67-pre4, in addition to circos-0.67-pre4.tgz
ls
# now you can detele the .tgz file
rm -f circos-0.67-pre4.tgz
```

4. Circos is ready to use if you finish installing the required perl packages. You can run circos by executing `./bin/circos`, but you can also add it to the shell path so you can run it anywhere by simply typing `circos`.
```
# check what's in the current PATH
echo $PATH
# add circos to the PATH (suppose circos is installed in /home/circos/circos-0.67-pre4)
export PATH=/home/circos/circos-0.67-pre4/bin:$PATH
# now you should be able to see /home/circos/circos-0.67-pre4/bin in the PATH
echo $PATH
# try execute circos
circos
# you will likely see error messages, this is normal as the required perl packages are not installed yet.
```
Alternatively, use the alias command so you don't have to type in './bin/circos every time':
```
sudo alias circos=/home/circos/circos-0.67-pre4/bin/circos
```
### Installing Perl Dependencies

Configuration (.conf files): specifies the features of the graph (cytogenetic bands, labels, ticks and, of course, data)

Karyotype:  defines the names, sizes and colors of chromosomes that you will use in the image.




