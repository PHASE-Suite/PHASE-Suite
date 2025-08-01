#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
set -o errtrace

###############################################################################
# Script Name: ecNGS Data Analysis Pipeline for PacBio Mutation Identification Analysis
# Author: Jaime Miranda
# Date: July 30, 2025
# Version: 2.0
#
# Description:
# This script automates a multi-step error corrected Next-Generation Sequencing
# (ecNGS) data analysis pipeline. It processes raw BAM files to identify somatic mutations
# providing a comprehensive overview of sequence variants within different sample subsets.
#
# The pipeline integrates various bioinformatics tools to perform:
# 1. Alignment and initial filtering of sequencing reads.
# 2. Conversion of aligned reads to BED format and calculation of base counts.
# 3. Splitting BAM files into chunks and extracting mismatches (mutations).
# 4. Joining forward and reverse strand mismatch data.
# 5. Merging results from all samples and applying final quality filters.
#
# Input Requirements:
# - Original BAM files: Located in the 'original-bams/' directory.
# - Reference Genome: A BWA-indexed FASTA file (e.g., hg38, mm10)
# - Sam2tsv.jar: A Java tool from JVarkit for converting SAM to TSV.
#
# Output Generated:
# The script creates several subdirectories under the main 'workdir' to store
# intermediate and final results:
# - filtered-bams/: Intermediate filtered BAM files.
# - output-bams/: Intermediate unfiltered BAM files.
# - fwd-rev/: Intermediate forward/reverse strand mismatch data.
# - bed/: BED files per sample and per barcode, containing genomic regions.
# - bases/: Total base counts per sample and per barcode.
# - filtered-results/: Per-sample filtered mutation results.
# - v1-joined/: Combined and globally filtered mutation results.
# - sample_times.tsv: A summary of processing times for each sample.
#
# Key Dependencies (External Tools):
# - samtools: For BAM/SAM manipulation.
# - bwa: For read alignment.
# - bedtools: For converting BAM to BED and genomic interval operations.
# - Java: To run sam2tsv.jar.
# - GNU Parallel: For parallel execution of tasks across multiple cores/samples.
# - awk, grep, sort, uniq, cut, tr: Standard Unix text processing utilities.
#
# Parallelization Strategy:
# - The script extensively uses GNU Parallel to process multiple samples concurrently
#   and to parallelize certain steps within individual sample processing (e.g.,
#   bwa alignment).
#
# Error Handling:
# - `set -euo pipefail`: Ensures the script exits immediately on any error,
#   uninitialized variable, or if a command in a pipeline fails.
# - `IFS=$'\n\t'`: Prevents word splitting issues.
# - `set -o errtrace`: Ensures traps are inherited by shell functions.
# - `error_handler()`: A trap function that provides detailed error messages
#   including the function name and line number where the error occurred.
#
# Configuration:
# - Important parameters like thread counts, number of splits, and file paths
#   are defined at the beginning for easy modification.
#
# Workflow Overview:
# The pipeline first processes each sample independently (Steps 1-4) in parallel
# up to the 'join_strands' stage, after all samples are processed and their
# results are merged and globally filtered (Step 5).
#
###############################################################################

########## ERROR HANDLER ##########
# Capture any non‐zero exit and report which function (and line) failed.
error_handler() {
  local exit_code=$?
  local line_no=${BASH_LINENO[0]}
  local func="${FUNCNAME[1]:-MAIN}"   # If not in a function, show MAIN
  echo >&2 "✗ ERROR in function '${func}' at or near line ${line_no} (exit code ${exit_code})."
  exit "$exit_code"
}
trap error_handler ERR

#### CONFIGURATION ####
threads_per_job=8
jobs=3
splits=4

ref_fa="/path/to/reference/genome.fa"

workdir="$(pwd)"
bam_dir="$workdir/original-bams"
out_dirs=( bases bed output-bams filtered-bams fwd-rev filtered-results v1-joined )

# File to accumulate per‐sample runtimes
summary_file="$workdir/sample_times.tsv"
: > "$summary_file"    # truncate/create

#### PREPARE DIRECTORIES ####
for d in "${out_dirs[@]}"; do
  mkdir -p "$workdir/$d"
done

#########################
# 1) Align & filter     #
#########################
align_and_filter() {
  local bamfile=$1
  local sample=${bamfile%.bam}

  # 1) Convert BAM to FASTQ, align to reference with BWA, and write output and filtered BAM
  samtools fastq -@ 2 "$bam_dir/$bamfile" \
    | bwa mem -t"$((threads_per_job-4))" -B21 -O11 "$ref_fa" - \
    | tee \
        >(samtools view -b -@ 1 -F 4079 -q 60 -o "$workdir/filtered-bams/${sample}.filtered.bam") \
        | samtools view -b -@ 1 -o "$workdir/output-bams/${sample}.bam"
}

####################################
# 2) BED-making & base counting    #
####################################
make_bed_and_count() {
  local sample=$1
  local filt="$workdir/filtered-bams/${sample}.filtered.bam"

  # 1) Convert filtered BAM to BED and use awk to match fwd/rev read pairs
  bedtools bamtobed -i "$filt" \
    | awk '
        /fwd/ { n=$4; sub(/fwd/,"",n); f[n]=$1"\t"$2"\t"$3 }
        /rev/ { n=$4; sub(/rev/,"",n); r[n]=$1"\t"$2"\t"$3 }
        END {
          for(n in f) if(n in r){
            split(f[n],F,"\t"); split(r[n],R,"\t")
            if(F[1]==R[1] && F[2]<R[3] && R[2]<F[3]){
              s=(F[2]>R[2]?F[2]:R[2])
              e=(F[3]<R[3]?F[3]:R[3])
              if (e > s)
                print n, F[1],s,e
            }
          }
        }' OFS=$'\t' \
    | tee "$workdir/bed/${sample}.bed" \
    | awk '{ sum += ($4 - $3) } END { print sum }' \
      > "$workdir/bases/${sample}.bases"
   
  # 3) Remove the intermediate filtered BAM  
  rm -f "$filt"
}

#########################################
# 3) Split & extract mismatches         #
#########################################
split_and_extract() {
  local sample=$1
  local bam="$workdir/output-bams/${sample}.bam"
  local hdr="$workdir/output-bams/${sample}.header.sam"
  local prefix="$workdir/output-bams/${sample}.chunk_"

  # 1) Extract the BAM header
  samtools view -H "$bam" > "$hdr"
  
  # 2) Count total alignments in the BAM and compute lines per chunk based on $splits
  local total chunk
  total=$(samtools view -c "$bam")
  chunk=$(( (total + splits - 1) / splits ))
  
  # 3) Split the alignment lines (no header) into "splits" chunks
  samtools view "$bam" | split -l "$chunk" - "$prefix"

  # 4) For each chunk, run sam2tsv then AWK to extract “/fwd” and “/rev” mismatches
  for f in ${prefix}*; do
    for st in fwd rev; do
      (
        cat "$hdr" "$f" \
          | sam2tsv -R "$ref_fa" \
          | awk -v st="$st" '
              $1 ~ ("/" st "$") && $10 == "M" {
                alt = toupper($6)
                ref = toupper($9)

                if (alt != ref && ref != "N") {
                  readname = $1
                  sub("/" st "$", "", readname)
                  key = readname ":" $4 ":" $8 ":" ref ":" alt
                  print key, $2, $3
                }
              }' OFS=$'\t'
      ) > "$workdir/fwd-rev/${sample}_$(basename "$f").${st}" &
    done
  done
  wait

  # 5) Remove the temporary header, chunk files and output-bam
  rm -f "$hdr" ${prefix}*
  rm -f "$bam"
}

####################################
# 4) Join strands & finalize       #
####################################
join_strands() {
  local sample=$1
  local out="$workdir/filtered-results/${sample}.FR.joined.v1"

  # 1) Concatenate all per‐chunk “.fwd” and “.rev” into .all.fwd / .all.rev
  cat "$workdir/fwd-rev/${sample}_"*.fwd > "$workdir/fwd-rev/${sample}.all.fwd"
  cat "$workdir/fwd-rev/${sample}_"*.rev > "$workdir/fwd-rev/${sample}.all.rev"

  # 2) Sort each by the five‐part key (field 1)
  sort -k1,1 "$workdir/fwd-rev/${sample}.all.fwd" > "$workdir/fwd-rev/${sample}.all.fwd.sorted"
  sort -k1,1 "$workdir/fwd-rev/${sample}.all.rev" > "$workdir/fwd-rev/${sample}.all.rev.sorted"

  # 3) Join on field 1 (the five‐part key)
  join -t $'\t' -j1 \
    "$workdir/fwd-rev/${sample}.all.fwd.sorted" \
    "$workdir/fwd-rev/${sample}.all.rev.sorted" \
    > "$workdir/fwd-rev/${sample}.FR.joined"
  
  # 4) Split that key (tr ":" "\t"), then reorder columns  
  tr ":" "\t" < "$workdir/fwd-rev/${sample}.FR.joined" \
    | awk '{ print $1,$6,$7,$8,$9,$2,$3,$4,$5 }' OFS=$'\t' \
    | sort -S80% --parallel=4 -u \
    > "$out"
    
  # 5) Clean up all intermediate files  
  rm -f \
    "$workdir/fwd-rev/${sample}.FR.joined" \
    "$workdir/fwd-rev/${sample}_"* \
    "$workdir/fwd-rev/${sample}.all."*
}

####################################
# Wrapper for one sample           #
####################################
# Description: Orchestrates the execution of alignment, bed-making, mismatch
#              extraction, and strand joining for a single sample.
# Arguments: $1 - BAM filename (e.g., sample.bam)
# Outputs: Updates various output directories and logs timing.
process_sample() {
  local bamfile=$1
  local sample=${bamfile%.bam}
  echo "===== [${sample}] Starting ====="
  local start_ts=$(date +%s)

  align_and_filter    "$bamfile"
  make_bed_and_count  "$sample"
  split_and_extract   "$sample"
  join_strands        "$sample"

  local end_ts=$(date +%s)
  local elapsed=$((end_ts - start_ts))
  # Format elapsed into HH:MM:SS
  local formatted
  printf -v formatted "%02d:%02d:%02d" \
    $((elapsed/3600)) \
    $(( (elapsed%3600) / 60 )) \
    $(( elapsed%60 ))

  printf ">>> [%-10s] Done in %s <<<\n\n" "$sample" "$formatted"

  # Append to the shared summary file (TSV: sample\tHH:MM:SS)
  echo -e "${sample}\t${formatted}" >> "$summary_file"
}

# Export functions and variables for GNU parallel
export -f align_and_filter make_bed_and_count split_and_extract join_strands process_sample
export threads_per_job splits ref_fa sam2tsv_jar workdir bam_dir summary_file

######################################
# Main entry for samples (Steps 1-4) #
######################################
cd "$bam_dir"
parallel --jobs "$jobs" process_sample ::: *.bam

echo "All samples have been submitted and finished."
echo
echo "===== Per‐sample timing summary ====="
printf "Sample\tRuntime (HH:MM:SS)\n"
sort "$summary_file"       # print sorted by sample name; change if you prefer original order

###############################
# 5) Combine & final filter  #
###############################
echo "===== Merging all samples & applying final filters ====="

# Merge all per-sample .FR.joined.v1
cat "$workdir/filtered-results/"*.FR.joined.v1 > "$workdir/v1-joined/all.joined"

# Strand-/MAPQ-/duplication filters
sort -k6 -S80% "$workdir/v1-joined/all.joined" \
  | uniq -u -f5 \
  | awk '$2==0 || $2==16' \
  | awk '$4==0 || $4==16' \
  | awk '$3==60 && $5==60' \
  > "$workdir/v1-joined/all.joined.1"

# Reorder fields & drop any remaining duplicate keys
cut -f1,6-9 "$workdir/v1-joined/all.joined.1" \
  | awk '{ print $2, $3, $4, $5, $1 }' OFS=$'\t' \
  | sort -k5 -S80% \
  | uniq -u -f4 \
  > "$workdir/v1-joined/all.joined.2"

echo "===== All done. Final file: v1-joined/all.joined.2 ====="

