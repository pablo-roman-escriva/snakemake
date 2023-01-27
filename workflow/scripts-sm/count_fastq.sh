#!/bin/bash 
lines=$(cat $1 | wc -l)
count=$(expr $lines / 4) 
echo -e "$count" > $2