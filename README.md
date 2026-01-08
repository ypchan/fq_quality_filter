# fq_quality_filter
To extract sanger-quality reads from long-reads fastq files.
## Compile
```bash
gcc -O3 -march=native -pipe fq_quality_filter.c -o fq_quality_filter -lz -pthread
```

## Usage
```bash
./fq_quality_filter 
```
```markdown
fq_quality_filter

Usage:
  ./fq_quality_filter -i in.fastq[.gz] -o out.fastq[.gz|-] [--threads 1] \
     -q INT --fraction FLOAT   (repeatable)

Examples:
  ./fq_quality_filter -i in.fq.gz -o out.fq.gz --threads 8 -q 20 --fraction 0.90 -q 30 --fraction 0.80
  ./fq_quality_filter -i in.fq.gz -o - --threads 8 -q 20 --fraction 0.90 -q 30 --fraction 0.80 | pigz -p 16 > out.fq.gz

Notes:
  - Auto-detect phred 33/64 by sampling quality lines.
  - Assumes 4-line FASTQ records.

```
