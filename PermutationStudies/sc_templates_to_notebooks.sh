#!/bin/bash

# Convert template notebooks to executable notebooks

# Number of permutations to perform
n_perm=10000
# Thresholds for the geospatial clustering
thresholds=("0.05" "0.1" "0.2" "0.5" "1" "2" "4")

mkdir -p notebook_empty

# households
echo "Households"
sed -e "s/%%data_key%%/${data_key}/g" \
  -e "s/%%n_perm%%/${n_perm}/g" \
  -e "s/%%cutoff%%/${cutoff}/g" \
  template/Perm-hh.ipynb > notebook_empty/Perm-hh-${n_perm}.ipynb

# buildings
echo "Buildings"
sed -e "s/%%n_perm%%/${n_perm}/g" \
  template/Perm-bd.ipynb > notebook_empty/Perm-bd-${n_perm}.ipynb

# buildings2
echo "Buildings2"
sed -e "s/%%n_perm%%/${n_perm}/g" \
  template/Perm-bd2.ipynb > notebook_empty/Perm-bd2-${n_perm}.ipynb

# locations
echo "Locations"
for threshold in "${thresholds[@]}"; do
  sed -e "s/%%threshold%%/${threshold}/g" \
    -e "s/%%n_perm%%/${n_perm}/g" \
    template/Perm-lc.ipynb > notebook_empty/Perm-lc-${threshold}-${n_perm}.ipynb
done

# locations2
echo "Locations2"
for threshold in "${thresholds[@]}"; do
  sed -e "s/%%threshold%%/${threshold}/g" \
    -e "s/%%n_perm%%/${n_perm}/g" \
    template/Perm-lc2.ipynb > notebook_empty/Perm-lc2-${threshold}-${n_perm}.ipynb
done
