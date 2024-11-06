#!/bin/sh
grep -o -E "^>\w+.*$" promoters.fa | tr -d ">" > seq_id.txt

