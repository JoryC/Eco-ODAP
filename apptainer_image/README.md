# apptainer_image

[//]: # (To combine the parts of the apptainer image together, execute the command in the apptainer_image directory: cat part_* > r-odaf_default.sif)

[//]: # (Verify the checksums using: `md5sum r-odaf_default.sif` and, `sha256sum r-odaf_default.sif`.)

[//]: # (The checksums should be: `16aa94fb9cd9a5d890ff699113cb31b8  r-odaf_default.sif` and, `725239ad404aa200d7f293cea67dc179a06f40e20e8213dde4e62ba5f627729e  r-odaf_default.sif`, respectively, if combined properly.)

The apptainer image can be stored once and re-used over and over. If storing on DRAC servers, I suggest saving it to a group-sharable directory in the projects directory so all users in your group can access the file. 

## Method 1

You can build the Main R-ODAF docker image following the [documentation](https://github.com/R-ODAF/Main). Once built, the docker image can be [converted to an apptainer image](https://docs.alliancecan.ca/wiki/Apptainer#Building_an_Apptainer_image).

Name the image `r-odaf_default.sif`

## Method 2

Send me an e-mail at jory.curry@ec.gc.ca and ask for me to share with you a copy of the apptainer image file through [Google Drive](https://drive.google.com/file/d/1d6QhN1aYcHrZV3WBsWMePOhsnJKHl1hT/view?usp=drive_link)!