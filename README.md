# filterBamClip #


### Installation
```
pip install git+http://github.com:wckdouglas/filterBamClip.git
```

```
usage: filterSoftClip.py [-h] -i INBAM [-o OUTBAM] [-s SINGLE_END]
                         [-b BOTH_END]

Filter alignments from sam file by softclip ratio

optional arguments:
  -h, --help            show this help message and exit
  -i INBAM, --inbam INBAM
                        input bam file
  -o OUTBAM, --outbam OUTBAM
                        output bam file (defulat: - )
  -s SINGLE_END, --single_end SINGLE_END
                        Maximum ratio of the whole alignment being clipped in
                        one end (default: 0.2)
  -b BOTH_END, --both_end BOTH_END
                        Maximum ratio of the whole alignment being clipped in
                        sum(each end) (default : 0.5)
```
