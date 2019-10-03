##################################
#####   Use Git Easily       #####
#####   Author: Wei Si       #####
#####   Date: Jan. 7. 2019   #####
##################################

##### Give permission by this command: chmod u+x easyGit.sh;
##### ./easyGit.sh ( using "Update" as the annotation );
##### ./easyGit.sh Annotation ( using other annotation defined by yourself );


#!/bin/bash
## commit easily;

if [ $# -eq 0 ]
then
    Annotation="Update";
else
    Annotation=$1;
fi

#git pull
git add .
git commit -m "${Annotation}"
git push

