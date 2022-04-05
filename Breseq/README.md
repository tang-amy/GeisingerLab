# Download file from TUFC to RC cluster
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
wget -r -np -R "index.html*" ftp://130.64.74.72/RESULTFOLDER/ --user edward.geisinger --password YEmTZuz --auth-no-challenge 2>download_log.err &
```
Change `RESULTFOLDER` to the actual folder name of the sequencing result. This can be found by right click the folder in the ftp server opened in Finder, and click `Get Info`. Remember to include the `/` at the end, otherwise the command may not recognize it as a folder and may attempt to download the whole server. 
