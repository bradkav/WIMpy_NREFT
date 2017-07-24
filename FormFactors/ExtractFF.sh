#!/bin/sh
sed -e 's/\\right.*$//' $1 > tmp.txt
sed -i '' 's/.*\\left(//' tmp.txt
sed -i '' '/&&/d' tmp.txt
sed -i '' 's/y^7.*$//' tmp.txt
sed -i '' 's/y^7//' tmp.txt
sed -i '' 's/y^6//' tmp.txt
sed -i '' 's/y^5//' tmp.txt
sed -i '' 's/y^4//' tmp.txt
sed -i '' 's/y^3//' tmp.txt
sed -i '' 's/y^2//' tmp.txt
sed -i '' 's/y//' tmp.txt
sed -i '' 's/\+/\ \+/' tmp.txt
sed -i '' 's/\-/\ \-/' tmp.txt
sed -i '' 's/^[ \t]*//' tmp.txt
mv tmp.txt FormFactors.txt