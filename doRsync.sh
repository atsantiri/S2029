#!/bin/bash

# Identity file
if [ "$1" == "ibc" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf_ibc"
elif [ "$1" == "mlg" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf"
elif [ "$1" == "ld" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf_ld"
elif [ "$1" == "at" ]; then
  SSH_KEY="$HOME/.ssh/rcmp_at_triumf_at"
else
  echo "Error: first argument must be either ibc, mlg, ld or at"
  exit 1
fi

# From where
remote="rcmp@142.90.96.93"

# Skip first argument
shift

# And now read directiories if passed
if [ "$#" -eq 0 ]; then
  dirs=("Raw" "Cluster" "Data" "Filter" "Merger") # default case
else
  dirs=("$@")
fi

#Act on dirs
for d in "${dirs[@]}"; do
  echo "Rsyncing ${d}"
  rsync -e "ssh -i ${SSH_KEY}" -avz --progress "${remote}:/home/rcmp/S2029/RootFiles/${d}/" "./RootFiles/${d}/"
done
