# Data Visualization Using Circos
## Installation
This user guide is based on the [documentation](http://circos.ca/documentation/tutorials/configuration/installation/) provided by circos.
However, the official documentation is a bit complicated and hard to follow, which is why I wrote this user guide.
### Install Circos
Circos is a command line tool (e.g. no GUI) and is perl based. It would be easier to run it on a mac computer (or linux). If you are using Mac OS, perl is most likely pre-installed (see Installing Perl Dependencies for detail). 

1. To install Circos, first make a directory where you want the software to be installed. For example, if you want to creat a folder in the home directory, execute the following in terminal: 
```bash
# "~" stands for home directory, you can install Circos to anywhere you want if you know what you are doing.
cd ~
mkdir circos
```
2. Download the [latest version of Circos](http://circos.ca/software/download/). Move the .tgz file to the folder created in the previous step.

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
```bash
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
Alternatively, use the alias command so you don't have to type in `./bin/circos` every time:
```bash
sudo alias circos=/home/circos/circos-0.67-pre4/bin/circos
```
### Install Perl Dependencies
1. First, check if perl is installed by typing the following command in the terminal:
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

2. There are many perl modules (packages) required for circos. The default perl does not come with these modules, so we need to install them manually. To see which modules are missing:
```bash
circos -module
# If returns "command not found: circos", it means you didn't link circos with its path correctly, see step 4 in Installing Circos.
```
This will return something like (with a lot more missing items):
```
missing  1.29 Carp
missing  0.36 Clone
ok       2.63 Config::General
ok       3.40 Cwd
ok      2.145 Data::Dumper
ok       2.52 Digest::MD5
ok       2.84 File::Basename
ok       3.40 File::Spec::Functions
ok       0.23 File::Temp
ok       1.51 FindBin
ok       0.39 Font::TTF::Font
missing  2.53 GD
missing  0.2 GD::Polyline
ok       2.39 Getopt::Long
```
3. Now let's install these packages using the cpan install command, except for GD and GD::Polyline which are a bit tricky to install. 
```bash
# add the names of the missing modules listed above (except for GD and GD::Polyline)
sudo cpan install Carp Clone
# enter your password when prompted
# check if the modules are installed correctly after done:
circos -module
```
4. The last thing is to install GD and GD::Polyline. Mostly likely `cpan install` will not install them correctly.
You will have to install the modules that GD is dependent on, including `libpng`, `freetype`, `libgd` and `jpegsrc`, then install the GD module. Follow the [instructions](http://circos.ca/documentation/tutorials/configuration/perl_and_modules/) on the circos website (section INSTALLING libpng, freetype, libgd AND gd) about how to install these packagess. 
Basically, download the tar.gz files for these modules, unzip the files, go to each of the folders and run `./configure -prefix=/usr/local`, `make`, and `make install` as stated on the link. If you encounter permission issue, add sudo in front of your command, e.g. `sudo make install`. (This installation process is called "build from source", which is sort of a standard way for installing softwares in Unix environments.)

5. Check again if all the perl dependencies are satisfied.
```bash
circos -module
```
If you see no "missing" module, congratulations the installation is done!

6. To see if circos runs correctly, let's test it by running the example dataset. Suppose your working directory is where circos is installed, run the following command:
```bash
circos -conf example/etc/circos.conf
# this by default should generate two output files `circos.png` and `circos.svg`
# open the image:
open circos.png
# you could also double click the image file in Finder
```

## Plot Data with Circos
Circos is highly versatile and can be used to generate various types of plots. See detailed tutorials on the [circos website](http://circos.ca/documentation/tutorials/).
### Command Line Usage
The basic usage of Circos is by typing the following command:
```
circos -conf circos.conf
```
An [example .conf file](https://github.com/tang-amy/GeisingerLab/blob/master/TnSeq_Processing/circos_plot/sample_with_3_barcode_tracks_mariner_1000_tn10_6000.conf) used for plotting *A. baumannii* Tn insertions and essential genes is available in this repository. **Do not run this .conf file directly!** You will have to modify the details in this file (file paths, tracks to plot, etc.).

### Key File Formats Used

**Configuration** (.conf): specifies the features of the graph (cytogenetic bands, labels, ticks and, of course, data)

**Karyotype** (.txt): defines the names, sizes and colors of chromosomes that you will use in the image. For *A. baumannii* there is only one chromosome and I have uploaded the [*A. baumannii* NZ_CP012004 karyotype file](https://github.com/tang-amy/GeisingerLab/blob/master/TnSeq_Processing/circos_plot/karyotype.Ab17978.txt) in this repository.

**Data File** (.txt): data that you wish to plot as different tracks. For scatter / line / histogram plots, the data files contains four requried columns and one optional column. See [Data Files](http://circos.ca/documentation/tutorials/configuration/data_files/) section on circos tutorial.

### Basic Configuration
Inside the configuration file, there are different sections (enclosed in `<section>` and `</section>`) that defines different features of the graph (e.g. ticks, labels, sub plots). You do not have to write the configuration file from scratch, as there are existing configurations in the `circos/etc/` folder.
#### Karyotype
This should be the first line of the configuration file. It specifies the location of the karyotype file. **Change the file path accordingly.**
```perl
karyotype=/home/circos/karyotype.Ab17978.txt
```
#### Image
Specifies the setting for image output, and the output name. In this case, the default setting is used by importing the `image.conf` from `etc/`. Output name can be specified by `file* =`, followed by the prefix you wish to use. The output is by default generated in the working directory. 
```perl
<image>
<<include etc/image.conf>>
file* = thickness_10_black_genome_with_3_barcode_tracks_mariner1000_tn10_6000
</image>
```
#### Ideogram
The ideogram is a circular visualization of the genome. For polyploid, the ideogram can plot different chromosomes in different colors. In the case of *A.baumannii*, only one chromosome will be plotted.
```perl
<ideogram>
<spacing>
default = 0.005r
</spacing>

radius  = 0.90r
thickness = 10p
fill = yes
stroke_color = dgrey
stroke_thickness = 2p
show_label = no
</ideogram>
```
#### Ticks
The genome positions can be indicated by adding ticks to the ideogram.
```perl
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p

# the tick label is derived by multiplying the tick position
# by 'multiplier' and casting it in 'format':
#
# sprintf(format,position*multiplier)

multiplier       = 1e-6

# %d   - integer
# %f   - float
# %.1f - float with one decimal
# %.2f - float with two decimals
#
# for other formats, see http://perldoc.perl.org/functions/sprintf.html

format           = %.2f

<tick>
spacing        = 0.04u
size           = 10p
</tick>

<tick>
spacing        = 0.2u
size           = 18p
show_label     = yes
label_size     = 40p
label_offset   = 18p
format         = %.2f
</tick>

</ticks>
```
#### Plots
This section specifies the tracks you wish to plot. Inside `<plots>` and `</plots>` you can add multiple tracks by writing multiple `<plot>` and `</plot>` blocks. The settings can be applied globally if it is outside individual `<plot>` and `</plot>` blocks, or can be applied specifically to one plot if written inside the block. You could also show data selectively by defining the `<rule>` and `</rule>` block (e.g. to plot only values greater than certain values). The example configuration contains 5 tracks. **Change the file paths accordingly.**

```perl
<plots>

<plot>
# Track 1: combined essentiality barcodes
type = histogram 
file = /Users/yunfei/circos-0.69-9/crispr_paper_circos/data_files_for_plotting/circos_mariner_tn10_essentiality.txt
r1 = 0.25r
r0 = 0.2r
# r1 is the position of the outer circle of the track
# r0 is the position of the inner circle of the track

stroke_type = outline
thickness = 1
color = red
extend_bin = no
</plot>

<plot>
# Track 2: Mariner essentiality barcodes
type = histogram 
file = /Users/yunfei/circos-0.69-9/crispr_paper_circos/data_files_for_plotting/circos_mariner_essentiality.txt 
r1 = 0.35r
r0 = 0.3r

stroke_type = outline
thickness = 1
color = grey
extend_bin = no

</plot>

<plot>
# Track 3: Tn10 essentiality barcodes
type = histogram 
file = /Users/yunfei/circos-0.69-9/crispr_paper_circos/data_files_for_plotting/circos_tn10_essentiality.txt
r1 = 0.45r
r0 = 0.4r

stroke_type = outline
thickness = 1
color = grey
extend_bin = no
</plot>

<plot>
# Track 4: Mariner reads
type = histogram 
file = /Users/yunfei/circos-0.69-9/crispr_paper_circos/data_files_for_plotting/linear_transitcombinedSS83ATAJB83_div200.txt
r1 = 0.72r
r0 = 0.5r

stroke_type = outline
thickness = 1
color = green
extend_bin = no

<rules>
<rule>
condition = var(value) > 1000
condition = 1
show = no
</rule>
</rules>

</plot>

<plot>
# Track 5: Tn10 reads
type = histogram 
file = /Users/yunfei/circos-0.69-9/crispr_paper_circos/data_files_for_plotting/linear_denomenator_700_TN10-TTR-combined_wig_output-20200609_copy.txt
r1 = 0.94r
r0 = 0.72r

stroke_type = outline
thickness = 1
color = blue
extend_bin = no

<rules>
<rule>
condition = var(value) > 6000
condition = 1
show = no
</rule>
</rules>
</plot>

</plots>

```
#### Others
The example configuration uses the default setting for colors and fonts, as well as the default housekeeping confuration.

```perl
<<include etc/colors_fonts_patterns.conf>> 
<<include etc/housekeeping.conf>> 
```
#### Output
By default, the output is `.png` and `.svg` image files. The `.png` file is smaller with limited resolution, and is good for preview purposes. The `.svg` file is a vector graph, meaning that it remains smooth no matter how large you zoom in, and it is better for publication purposes.
