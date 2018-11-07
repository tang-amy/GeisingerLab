Transit Installation and Troubleshooting
=============================
_On: MacOS, Ubuntu Linux_
</br>

MacOS/Linux Installation
-------------------------------
Make sure Python2 is installed and the python variable path goes to Python2.
```bash
python --version
```
If the version is Python 2.7.x, then the pathing variable is correct. 
</br>

Install pip with the following command:
```bash
sudo easy install pip 
``` 
</br>

Use pip to install Transit (this should install of the necessary dependencies):
```bash
pip install tnseq-transit
``` 

Running Transit
-------------------------------
<i> The following for running Transit is taken from their website. If the above instructions were followed for installation, Transit should be installed as Python package, so Transit can be called from the command line/terminal. If the installation method was different, please see [Transit's official documentation](https://transit.readthedocs.io/en/latest/transit_running.html). </i>

### Running in GUI Mode: ###
Use command line/terminal and type the following command:
```bash
transit
```

### Running in Command Line/Headless Mode: ###
The following are the required arguments:
```bash
transit <method> <comma-separated .wig files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]
``` 
__IMPORTANT:__ When running Transit in Command Line/Headless Mode, you need to have the absolute path names (unless Transit is in the same folder as all of the input and output files). 
To avoid typing in or copy-and-pasting in all of the files, create a directory of only the wig files which are being analyzed. 

Run the following code in the command line/terminal:
```bash
for f in PATH/*.wig; do  full_path=`pwd`/$f; all_inputs="${all_inputs}${all_inputs:+,}$full_path"; done
# check that the variable is set correctly
echo $all_inputs 
```
When running this code, it is important to note, the _PATH_ for the file should not include the leading '/' character if you are entering a directory name.

Run transit with this string of inputs (comma-separated wig files) that was created:
```bash
transit <method> $all_inputs <annotation .prot_table or GFF3> <output file> [Optional Arguments]
```

After running transit, clear the variable name:
```bash 
unset all_inputs
``` 

General Troubleshooting
------------------------------
For general troubleshooting, visit [Transit's website](https://transit.readthedocs.io/en/latest/transit_install.html#troubleshooting). Below contains fixes to specific errors/issues that we had run into when installing/running Transit.

Troubleshooting Installation
-------------------------------
<h3> Troubleshooting on MacOS: </h3>

__Error message:__ 
> Cannot install [pkg]. It is a distutils installed project and thus we cannot accurately determine which files belong to it which would lead to only a partial uninstall.

Keep in track the name of the package that pip cannot install. In the following code, where the word [pkg] is used, type the name and version of the package that pip cannot install. 
```bash
sudo -H pip install --ignore-installed [pkg]
```
The entire tnseq-transit package can be reinstalled with this command (replace [pkg] with _tnseq-transit_). Any existing installations/dependencies should be left alone.

<h3> Troubleshooting on Linux: </h3>

__Error:__ contains a line similar to _def getDefaultPublisher() -> Publisher:_
If installed through Linux's package manager (_apt-get install_), there is a versioning error within the package for PyPubSub. 
You can either:
* Delete the Transit folder that was created upon installation and follow the above instructions for installation.
* Uninstall and reinstall the version of PyPubSub (3.3.0) that Transit is dependent on.
    ```bash
	sudo -H pip uninstall PyPubSub
	sudo -H pip install PyPubSub==3.3.0
    ```z

Troubleshooting Running Transit
-------------------------------
<h3> Specific errors not mentioned on Transit's official documentation: </h3>

__Error Message:__ _ValueError: invalid literal for int() with base 10: '>129734'_ 

This is caused by incorrect parsing of the genome bank (.gbk) file for the prot table (.prot_table file) which leaves ">" or "<" signs in which cannot later on be converted to numbers. As long as all of the "<" and ">" are removed, this error should no longer occur. 

The prot table (.prot_table file) can be remade using this [make_prot_table.py script](https://github.com/tang-amy/GeisingerLab/blob/wip/TnSeqProcessing/src/make_prot_table.py).   

__Error Message:__ _IOError: [Errno 2] No such file or directory:..._

Transit needs to have the absolute file pathing (cannot use "~" in place of /Users/...) in order to read/write the files unless it is in the same directory as the input and output files. For easier way to input a long list of files, see _Running Transit: Running in Command Line/Headless Mode_.
