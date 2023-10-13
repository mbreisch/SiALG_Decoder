#!/bin/bash

#Application path location of applicaiton

Dependencies=`pwd`/Dependencies

source ${Dependencies}/root/install/bin/thisroot.sh

export LD_LIBRARY_PATH=`pwd`/lib:${Dependencies}/zeromq-4.0.7/lib:${Dependencies}/boost_1_66_0/install/lib:${Dependencies}/pqxx/install/lib:$LD_LIBRARY_PATH/lib:${Dependencies}/root/install/lib

export SEGFAULT_SIGNALS="all"
