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
<i> The following for running Transit is taken from their website. If the above instructions were followed for installation, because Transit was installed as Python package, it can be called from the command line/terminal. If the installation method was different, please see [Transit's official documentation](https://transit.readthedocs.io/en/latest/transit_running.html). </i>

### Running in GUI Mode: ###
Use command line/terminal and type the following commands:
```bash
transit
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

__Error:__

<p> If installed through Linux's package manager (_apt-get install_), there may be a dependency problem within the package. Delete the Transit folder that was created upon installation and follow the above instructions for installation. </p>

Troubleshooting Running Transit
-------------------------------
<h3> Specific errors not mentioned on Transit's official documentation: </h3>

__Error Message:__ _ValueError: invalid literal for int() with base 10: '>129734'_ 

This is caused by incorrect parsing of the genome bank (.gbk) file for the prot table (.prot_table file) which leaves ">" or "<" signs in which cannot later on be converted to numbers. As long as all of the "<" and ">" are removed, this error should no longer occur. 

The prot table (.prot_table file) can be remade using the make_prot_table.py on this Github Repository.   
