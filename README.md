# polymc

When you are an engineer and want to work on haplotyping algorithm in matlab, you need a .mat file instead of a fragment as a text file.
The matlab code convert_frag_mat.m is for converting a fragment file, the output of extractHAIRS from HAPCUT package to a sparse matrix.
This version is dedicated to biallelic variants which can be used for polypoids as well as diploids.
The input is a fragment file like 

```
Header
1 NC_001133-318 1 0000101001 AIGGGGDIII
1 NC_001133-288 1 0000010000 ?HHHEGFGHI
1 NC_001133-348 2 000101001 HHGGEEHII
```

A sample .mat file is a standard matrix file.
```
[-1,-1,1,0,0,0,0,0;
-1,-1,-1,-1,0,0,0,0;
-1,-1,-1,-1,0,0,0,0]
```


By Sina Majidian Dec 2018

Iran University of Science and Technology
