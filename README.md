# danuser/common
[![build status](https://git.biohpc.swmed.edu/danuser/common/badges/master/build.svg)](https://git.biohpc.swmed.edu/danuser/common/commits/master)

## Purpose
 Core functions and MovieData APIs shared across applications.

For Users
---------
+ Basic user tutorial is provided here: 
    + https://git.biohpc.swmed.edu/danuser/common/tree/master/toolfun/movieManagement/tutorial.
+ User documentation is provided in subdirectories listed as `doc`:
    + https://git.biohpc.swmed.edu/danuser/common/tree/master/toolfun/movieManagement/doc.
    + https://git.biohpc.swmed.edu/danuser/common/tree/master/toolfun/movieManagement/interface/doc.

For Developers
--------------
+ Changes to the `master` branch are resticted.  
+ Please commit changes to the `testing` branch. 
+ The `master` branch will be merged on a monthly basis with `testing`.

+ New commits will trigger the gitlab CI pipeline [testing: [![build status](https://git.biohpc.swmed.edu/danuser/common/badges/testing/build.svg)](https://git.biohpc.swmed.edu/danuser/common/commits/testing)]  [master: [![build status](https://git.biohpc.swmed.edu/danuser/common/badges/master/build.svg)](https://git.biohpc.swmed.edu/danuser/common/commits/master)]
   + Test datasets for the automated CI pipeline are maintained in the BioHPC lamella cloud via the "danuserweb" account.
   + See more details here https://git.biohpc.swmed.edu/danuser/ci-scripts.
