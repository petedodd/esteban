#!/bin/bash
for i in $(ls *pdf); do
	convert $i ${i/pdf}png
done
