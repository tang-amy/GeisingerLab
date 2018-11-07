Transit Installation Set Up and Troubleshooting
=============================
_On: MacOS, Ubuntu Linux_
</br>

MacOS/Linux Installation/Running GUI Mode
-------------------------------
Make sure Python2 is installed and the python variable path goes to Python2.
'''bash
python --version
'''
If the version is Python 2.7.x, then the pathing variable is correct. 
</br>

Install pip with the following command:
'''bash
sudo easy install pip 
''' 
</br>

Use pip to install Transit (this should install of the necessary dependencies):\
'''bash
pip install tnseq-transit
''' 
</br>

Run Transit in GUI Mode by typing _transit_ into the command line/terminal.

MacOS Troubleshooting
-------------------------------
__Error message:__ 
> Cannot install [pkg]. It is a distutils installed project and thus we cannot accurately determine which files belong to it which would lead to only a partial uninstall.

Keep in track the name of the package that pip cannot install. In the following code, where the word [pkg] is used, type the name and version of the package that pip cannot install. 
'''bash
sudo -H pip install --ignore-installed [pkg]
'''
The entire tnseq-transit package can be reinstalled with this command (replace [pkg] with _tnseq-transit_). Any existing installations/dependencies should be left alone.

Linux Troubleshooting
-------------------------------
__Error:__
<p> If installed through Linux's package manager (_apt-get install_), there may be a dependency problem within the package. Delete the transit folder that was created upon installation and follow the above instructions for installation. </p>

General Troubleshooting
-------------------------------
For more general troubleshooting, visit [Transit's website](https://transit.readthedocs.io/en/latest/transit_install.html#troubleshooting).