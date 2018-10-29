VirtualBox Set Up and Troubleshooting
=============================
_On: Host MacOS, Guest Linux (Ubuntu)_
</br>

VirtualMachine SetUp
-----------------------------
<h3> Specs </h3>
<ul>
<li> Ram: Half of available ram (16GB, 16384 MB) </li>
<li> CPUs: Half of # of total CPUs (8) </li>
<li> Video memory: highest possible to reduce lag/pixelation (128 MB) </li>
<li> Storage/memory: doesn't need to be over 128GB but can be set to 256GB (use a shared folder to store data files) </li>
</ul>
<h3> VBoxGuestAdditions: </h3>
<p> Necessary for full-screen and folder-sharing between Guest and Host machines.
<ul>
<li> Install the correct version of VBoxGuestAdditions (iso) according to version of VirtualBox being used. </li>
<li> Mount the iso using the menu bar (install VboxGuestAdditions iso) or by terminal. </li>
</ul>

VirtualBox/VirtualMachine Troubleshooting
--------------------------------------------------------
<h3> Gatekeeping On MacOS </h3>
<p> Mac OSX security has gatekeeping to block kernels from unknown developers and downloaded/outside kernels. An error message saying that the .kext file is corrupted, .kext can't be run may be displayed </p>

__Follow steps below in host (Mac) terminal:__
```bash
sudo su # goes into root user account
sudo spctl --master-disable # disables gatekeeping, allows everything
exit
```
<p> After allowing apps to be downloaded from anywhere, install VirtualBox. 

<b><i>Go to Settings and Preferences, in Security and Privacy, change allowed apps to be apps downloaded from AppStore/Apple and from known developers.</i></b> </p>

<h3> Folder sharing </h3>
<p> Once the folder being mounted is created, it may be able to be mounted via VirtualBox app in the settings of whatever VM the folder is being mounted onto. For the options, check: <b> auto-mount </b> and <b>make permanent</b>. The folder should appear on the desktop as a drive and in the <i>/media</i> folder. </p>

__If the folder doesn't appear but is mounted in VM's settings:__
```bash 
id g # check groups that the user is in
sudo usermod -a -G vboxsf $USER # to add user to necessary group for the VBox Shared Folder
```
__If folder is not appearing even though user is in vboxsf group:__
```bash
#turn on rc.local in ubuntu
printf '%s\n' '#!/bin/bash' 'exit 0' | sudo tee -a /etc/rc.local
sudo chmod +x /etc/rc.local
sudo reboot

sudo systemctl edit --full rc-local
# add to rc.local: 
# sudo mount -t vboxsf <Shared Folder specified in VB setting> <mount point path in guest system>
sudo reboot
```





