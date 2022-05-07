# WGS analysis by Breseq
## Download file from TUFC to RC cluster
To download large files on the cluster, use the -xfer node.
```
ssh -Y username@xfer.discovery.neu.edu
```
Enter your NU account password when prompted. Next, redirect to the destination directory for storing raw files. Make one by `mkdir newfolder` when needed.
```
cd /work/geisingerlab
mkdir newfolder
cd newfolder
```
In `newfolder`, execute the following command to download from the ftp server. The username and password can be found in the shared folder on google drive.
```
wget -r -np -R "index.html*" ftp://130.64.74.72/RESULTS/ --user edward.geisinger --password YEmTZuz --auth-no-challenge 2>download_log.err &
```
Replace `RESULTS` with the actual folder name of the sequencing result. This can be found by right clicking the folder in the ftp server opened in Finder, and clicking `Get Info`. Remember to include the `/` at the end, otherwise the command may not recognize it as a folder and as a result may attempt to download the whole server. 

Check the status of the download by checking the `download_log.err` file. Press `G` on the keyboard to move the the last line of the file where it shows the latest status of the download. 
```
less download_log.err
```
Exit the -xfer node.
```
exit
```

## Preprocess raw files for breseq analysis
Connect to the cluster (login node).
```
ssh -Y username@login-00.discovery.neu.edu
```
Enter your NU account password when prompted. Next, redirect to the directory containing the raw files.
```
cd /work/geisinger/newfolder
```
