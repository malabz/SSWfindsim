# Find similar segments using SIMD/OpenMP Striped Smith-Waterman algorithm

This project is based on https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.

## Arugments
```text
--center   FILE      center file name (Required)
--seq      FILE      sequence file name (Required)
--out      FILE      output file name (Required)
--matrix   FILE      matrix file name
--thread   N         threads in this program (default 1)
--dna                this file is DNA sequence file
--protein            this file is protein seqeuence file
--help               print help message
--version            show program version
```

## Input and output

- Input: two files, center (target) file is the base sequence, must have only 1 sequence; sequence (query) file has many ($x$) sequences.
- Output: one file with $x$ lines, every line contains:
- * `sequence_id`, which means the number of sequence;
- * `query_begin`, `query_end`, which means the place of center sequence;
- * `target_begin`, `target_end`, which means the place of query sequence;
- * `score`, SW score;
- * `alignment_list`, which print the sequence alignment ways: `I` insert, `D` delete and `M` match.
- The order of line is: `sequence_id query_begin query_end target_begin target_end score alignment list`.
- * Sample: find the similarity of target sequence `AAAATC` and query sequence `AAAAATC` with treated them as `DNA` sequences:
- * Will print `0 1 6 0 5 12 M6`.
