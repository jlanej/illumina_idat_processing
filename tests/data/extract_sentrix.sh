#!/usr/bin/env bash
# usage: ./extract_sentrix.sh SampleSheet.csv > out.tsv

awk -F',' '
BEGIN { indata=0; }
$0=="[Data]" { indata=1; next }
indata==1 {
  if (hdr==0) {
    for (i=1;i<=NF;i++) {
      gsub(/^[ \t]+|[ \t]+$/,"",$i)
      if ($i=="SentrixBarcode_A") sb=i
      else if ($i=="SentrixPosition_A") sp=i
      else if ($i=="Sample_Group") sg=i   # change to Sample_Name if desired
    }
    hdr=1
    next
  }

  if (sb && sp && sg && $sb!="" && $sp!="" && $sg!="")
    print $sb "_" $sp "\t" $sg
}
' "$1"
