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
