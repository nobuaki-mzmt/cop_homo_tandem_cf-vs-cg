#!/bin/bash

# define variable
date=`date '+%Y%m%d'`

from_filename="draft_CoptoHomoTandem.docx"
from_exe="docx"
to_filename="draft_CoptoHomoTandem.md"
to_exe="gfm"

current_branch=`git rev-parse --abbrev-ref HEAD`


# command
pandoc -s ${from_filename} -f ${from_exe} -t ${to_exe} ${options} -o ${to_filename}


gitfunc() {
  git add -A
  git commit -m "Auto commit by docx2md"
  git push --set-upstream origin ${currentbranch}
}

if [[ $# -gt 0 ]]; then
  echo "catch commandline args"
  if [[ $1 == "push" ]]; then
    echo "git add / commit / push"
    gitfunc
  fi
fi