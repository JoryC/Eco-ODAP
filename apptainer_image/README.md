# apptainer_image

[//]: # (To combine the parts of the apptainer image together, execute the command in the apptainer_image directory: cat part_* > r-odaf_default.sif)

[//]: # (Verify the checksums using: `md5sum r-odaf_default.sif` and, `sha256sum r-odaf_default.sif`.)

[//]: # (The checksums should be: `16aa94fb9cd9a5d890ff699113cb31b8  r-odaf_default.sif` and, `725239ad404aa200d7f293cea67dc179a06f40e20e8213dde4e62ba5f627729e  r-odaf_default.sif`, respectively, if combined properly.)

The apptainer image can be stored once and re-used over and over. If storing on DRAC servers, I suggest saving it to a group-sharable directory in the projects directory so all users in your group can access the file. Currently **Method 3 is preferred**.

Note: Methods 1 and 2 will require you to install additional conda environments interactively in the apptainer shell... which is not recommended.

Tip: to avoid using docker as root, add your user to the docker group `sudo usermod -aG docker $USER`.
Tip: Use a file transfer tool like FileZilla to copy the SIF image to a long-term storage location on DRAC servers like your userhome/projects/ location...

## Method 1 - Build from the original Main R-ODAF Dockerfile

You can build the Main R-ODAF docker image following the [documentation](https://github.com/R-ODAF/Main). Once built, the docker image can be [converted to an apptainer image](https://docs.alliancecan.ca/wiki/Apptainer#Building_an_Apptainer_image).

Specifically, download the original Main R-ODAF [DOCKER](https://github.com/R-ODAF/Main) file locally to this directory on a machine with `sudo` privledges and docker installed. 
List your images using the docker images command, then create a tar archive file using docker save and the image ID:

```
user@my_pc:Eco-ODAP/apptainer_image$ docker build -t r-odaf_default .

user@my_pc:Eco-ODAP/apptainer_image$ docker images
REPOSITORY                        TAG               IMAGE ID       CREATED          SIZE
r-odaf_default                  sometag            5a08b442ac65   1 hour ago      2.5GB

user@my_pc:Eco-ODAP/apptainer_image$ docker save 5a08b442ac65 -o r-odaf_default.tar
```

Next, convert the tar file into an Apptainer container using the `docker-archive` bootstrap agent (apptainer must be installed):

```
user@my_pc:Eco-ODAP/apptainer_image$ apptainer build r-odaf_default.sif docker-archive:r-odaf_default.tar
INFO:    Starting build...
Getting image source signatures
Copying blob sha256:7b940d9c2334bf5f2c547da9c9bc62a48fd2127f3461e2a08f113da7d1f152a1
  102.44 MiB / 102.44 MiB [==================================================] 5s
Copying blob sha256:caf639e1a8776734d5adbf271b2a5f3823b8d4f236d62c8f1a51d38e2f07e3db
  23.89 KiB / 23.89 KiB [====================================================] 0s
Copying config sha256:8c12b8d4e752fb3c5620f208f2136b8f61fdcb91b2a0e6f9cefa1d3db5dfe104
  4.27 KiB / 4.27 KiB [======================================================] 0s
Writing manifest to image destination
Storing signatures
INFO:    Creating SIF file...
INFO:    Build complete: r-odaf_default.sif
```

## Method 2 - Download the Main R-ODAF apptainer image

Send me an e-mail at jory.curry@ec.gc.ca and ask for me to share with you a copy of the original Main R-ODAF apptainer image file through [Google Drive](https://drive.google.com/file/d/1d6QhN1aYcHrZV3WBsWMePOhsnJKHl1hT/view?usp=drive_link)!

## Method 3 - Create the apptainer image from the provided Dockerfile

Install docker and apptainer on your local machine with sudo priveledges.
Ensure the build context is correct; the Directory structure should look like:
```
/path/to/cloned/repo/
├── apptainer_image/
│   └── Dockerfile
└── conda_environments/
    ├── r-odaf_default_updated_my_r_pkgs_environment.yml
    ├── r-odaf_default_mymultiqc_environment.yml
    └── r-odaf_default_DESeq2_report_env_environment.yml
```

Build the Dockerfile:
```
user@my_pc:~$ cd /path/to/cloned/Eco-ODAP/ #Change directories to the top of the cloned directory
user@my_pc:/path/to/cloned/Eco-ODAP/$ docker build -t r-odaf_default -f apptainer_image/Dockerfile .
```

Save the docker container:
```
user@my_pc:/Eco-ODAP/$ docker save -o r-odaf_default.tar r-odaf_default
```

Convert the docker image to an apptainer image:
```
user@my_pc:~$ cd /path/to/cloned/Eco-ODAP/apptainer_image/ #Change directories to the apptainer_image directory
user@my_pc:Eco-ODAP/apptainer_image$ apptainer build r-odaf_default.sif docker-archive:r-odaf_default.tar
INFO:    Starting build...
Copying blob 89169d87dbe2 done   |
Copying blob 29b75e9b3d5e done   |
Writing manifest to image destination
2024/10/28 11:35:19  info unpack layer: sha256:17ea5ece76fbd585e33798aceb5695859529dce90c725b3631f14b10418174b4
2024/10/28 11:35:21  warn rootless{usr/bin/ping} ignoring (usually) harmless EPERM on setxattr "security.capability"
2024/10/28 11:35:24  warn rootless{usr/sbin/arping} ignoring (usually) harmless EPERM on setxattr "security.capability"
2024/10/28 11:35:24  warn rootless{usr/sbin/clockdiff} ignoring (usually) harmless EPERM on setxattr "security.capability"
2024/10/28 11:35:26  info unpack layer: sha256:38c96f6ab21f39325a06e92164129850242fb75deb582b4142d6ca0df6ecd16e
INFO:    Creating SIF file...
INFO:    Build complete: r-odaf_default.sif
```

Use the apptainer:
```
user@my_pc:Eco-ODAP/apptainer_image$ apptainer exec --fakeroot --bind .:/mnt r-odaf_default.sif \
bash -c "find /mnt -type d -exec chmod 750 {} \; && \
find /mnt -type f -exec chmod 750 {} \; && \
source /opt/miniconda3/etc/profile.d/conda.sh && \
conda activate && \
echo "Hello World!" > helloworld.txt"
```