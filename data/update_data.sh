#!/bin/bash

cd /Volumes/BOOTCAMP/Users/tsz/Desktop/opfo_images/

rsync -R -P */*/*/*/dor[1-9]_Multifocus*.tif /Users/tsz/Documents/unil/opfo_main/4_traits/data/img/

rsync -R -P */*/*/[a-z][a-z][a-z][1-4].xlsx /Users/tsz/Documents/unil/opfo_main/4_traits/data/img/
