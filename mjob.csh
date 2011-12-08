#!/bin/csh
# Submit batch jobs
# This must be used with MkRunBat.sh

cd _pwd_
eval `scramv1 runtime -csh`
cd -
&
tar zcvf _output_file_.tgz *.log *.pdf *.txt *.root
rfcp _output_file_.tgz _storage_
