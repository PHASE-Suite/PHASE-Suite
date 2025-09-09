#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
set -o errtrace

###############################################################################
# Script Name: PHASE Suite - ecNGS Data Processing Pipeline
# Author: Jaime Miranda
# Date: September 9, 2025
# Version: 3.0
#
# Description:
# This script is the central engine of the PHASE Suite. It automates a multi-step,
# error-corrected Next-Generation Sequencing (ecNGS) data analysis pipeline
# that is controlled by external configuration files. The script processes raw
# PacBio HiFi BAM files to identify high-confidence somatic mutations and
# demultiplexes them by barcode, preparing the data for the PHASE Explorer
# visualization app.
#
# The pipeline integrates various bioinformatics tools to perform:
# 1. Alignment and initial filtering of sequencing reads.
# 2. Consensus BED region generation and per-sample base counting.
# 3. Mismatch extraction from forward and reverse strands.
# 4. Error correction via strand-joining logic.
# 5. Global merging and filtering of results from all samples.
# 6. Barcode-based demultiplexing of the final mutation and BED files.
#
# Execution Requirements:
# 1. A main configuration file (e.g., pipeline.conf).
# 2. A barcode configuration file (e.g., my_barcodes.cfg).
# 3. Original BAM files located in an 'original-bams/' directory.
# 4. A BWA-indexed reference genome, with the path specified in the config.
#
# Output Generated:
# The script creates several subdirectories. The final, user-facing outputs
# are located in the 'indexing/' directory.
#
# Final Outputs (in indexing/):
# - mutations/*.muts: Demultiplexed, high-confidence mutation lists.
# - bases/*.bases: Demultiplexed total base counts for normalization.
# - bed/*.bed: Demultiplexed BED files defining consensus regions.
#
# Intermediate Outputs:
# - filtered-bams/, output-bams/, fwd-rev/, bed/, bases/, filtered-results/
# - v1-joined/: Contains the globally merged and filtered master mutation list.
# - sample_times.tsv: A summary of processing times for each sample.
#
# Key Dependencies (External Tools):
# - samtools: For BAM/SAM manipulation.
# - bwa: For read alignment.
# - bedtools: For genomic interval operations.
# - sam2tsv (from jvarkit): For converting SAM to TSV to find mismatches.
# - seqkit: For efficient FASTA searching during demultiplexing.
# - GNU Parallel: For parallel execution of tasks across multiple cores.
#
# Workflow Overview:
# The pipeline first processes each raw BAM file in parallel to identify mutations
# on a per-sample basis. Once all samples are complete, their results are merged and
# subjected to a filtering stage to create a final list of mutations.
# In the final stage, this list is demultiplexed, assigning each mutation back to
# its specific barcode of origin for every sample.
#
# Usage:
# bash phase_pipeline.sh <main_config_file> <barcode_config_file>
#
# Example:
# bash phase_pipeline.sh pipeline.conf my_barcodes.cfg
#
###############################################################################

########## ERROR HANDLER ##########
error_handler() {
  local exit_code=$?
  local line_no=${BASH_LINENO[0]}
  local func="${FUNCNAME[1]:-MAIN}"
  echo >&2 "✗ ERROR in function '${func}' at or near line ${line_no} (exit code ${exit_code})."
  exit "$exit_code"
}
trap error_handler ERR

#### CONFIGURATION LOADING ####

# 1) Check for correct number of arguments
if [[ $# -ne 2 ]]; then
    echo >&2 "Usage: $0 <main_config_file> <barcode_config_file>"
    echo >&2 "Example: $0 pipeline.conf my_barcodes.cfg"
    exit 1
fi

main_config_file=$1
barcode_config_file=$2

# 2) Load Main Configuration
if [[ ! -f "$main_config_file" ]]; then
    echo >&2 "✗ ERROR: Main configuration file not found at '$main_config_file'"
    exit 1
fi
source "$main_config_file"
echo "===== Loaded main settings from '$main_config_file' ====="

# 3) Load Barcode Configuration
if [[ ! -f "$barcode_config_file" ]]; then
    echo >&2 "✗ ERROR: Barcode file not found at '$barcode_config_file'"
    exit 1
fi

declare -A barcodes=()
while read -r bc_name bc_seq; do
    # Skip empty lines or lines that are comments
    [[ -z "$bc_name" || "$bc_name" =~ ^# ]] && continue
    barcodes["$bc_name"]="$bc_seq"
done < "$barcode_config_file"

if [[ ${#barcodes[@]} -eq 0 ]]; then
    echo >&2 "✗ ERROR: No barcodes were loaded from '$barcode_config_file'."
    echo >&2 "Ensure the file is not empty and is formatted correctly (name<space>sequence)."
    exit 1
fi
echo "===== Loaded ${#barcodes[@]} barcodes from '$barcode_config_file' ====="


#### SETUP & PREPARATION ####
workdir="$(pwd)"
bam_dir="$workdir/original-bams"
out_dirs=(
    bases bed output-bams filtered-bams fwd-rev filtered-results v1-joined
    indexing indexing/bases indexing/bed indexing/mutations indexing/temp
)

# File to accumulate per‐sample runtimes
summary_file="$workdir/sample_times.tsv"
: > "$summary_file" # truncate/create

echo "===== Preparing work directories... ====="
for d in "${out_dirs[@]}"; do
  mkdir -p "$workdir/$d"
done

# Create per-barcode temp directories for list files
for bc_name in "${!barcodes[@]}"; do
    mkdir -p "$workdir/indexing/temp/${bc_name}"
done

##############################################################################
#                                                                            #
#                      PART 1: CORE MUTATION PIPELINE                        #
#                                                                            #
##############################################################################

#########################
# 1) Align & filter     #
#########################
align_and_filter() {
  local bamfile=$1
  local sample=${bamfile%.bam}

  # 1) Convert BAM to FASTQ, align to reference, and write output and filtered BAM
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
    | awk '{ sum += ($4 - $3) } END { print sum+0 }' \
      > "$workdir/bases/${sample}.bases"
      
  # 2) Remove the intermediate filtered BAM 
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
  rm -f "$hdr" ${prefix}* "$bam"
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
  rm -f "$workdir/fwd-rev/${sample}.FR.joined" "$workdir/fwd-rev/${sample}_"* "$workdir/fwd-rev/${sample}.all."*
}

####################################
# Wrapper for one sample           #
####################################
process_sample() {
  local bamfile=$1
  local sample=${bamfile%.bam}
  echo "===== [${sample}] Starting Core Pipeline ====="
  local start_ts=$(date +%s)

  align_and_filter    "$bamfile"
  make_bed_and_count  "$sample"
  split_and_extract   "$sample"
  join_strands        "$sample"

  local end_ts=$(date +%s)
  local elapsed=$((end_ts - start_ts))
  local formatted
  printf -v formatted "%02d:%02d:%02d" \
    $((elapsed/3600)) \
    $(( (elapsed%3600) / 60 )) \
    $(( elapsed%60 ))

  printf ">>> [${sample}] Core pipeline done in %s <<<\n\n" "$formatted"
  echo -e "${sample}\t${formatted}" >> "$summary_file"
}

##############################################################################
#                                                                            #
#                        PART 2: DEMULTIPLEXING                              #
#                                                                            #
##############################################################################

demultiplex_sample() {
  local sample=$1
  local fasta_tmp="$workdir/indexing/temp/${sample}.fasta"
  echo "===== [${sample}] Starting Demultiplexing ====="

  # Create FASTA once per sample for efficiency
  samtools fasta "$bam_dir/${sample}.bam" > "$fasta_tmp"

  for bc_name in "${!barcodes[@]}"; do
    local bc_seq="${barcodes[$bc_name]}"
    local list_file="$workdir/indexing/temp/${bc_name}/${sample}.list"
    
    seqkit grep -s -p "$bc_seq" -R 1:20 "$fasta_tmp" | \
      grep ">" | tr -d ">fwdrev" | sort -u > "$list_file"

    if [[ ! -s "$list_file" ]]; then
      echo "  - No reads found for barcode ${bc_name}."
      > "$workdir/indexing/mutations/${sample}.${bc_name}.muts"
      > "$workdir/indexing/bed/${sample}.${bc_name}.bed"
      echo 0 > "$workdir/indexing/bases/${sample}.${bc_name}.bases"
      continue
    fi
    
    echo "  - Filtering results for barcode ${bc_name}..."
    
    awk 'FNR==NR{a[$1];next} ($5 in a)' "$list_file" "$workdir/v1-joined/all.joined.2" > "$workdir/indexing/mutations/${sample}.${bc_name}.muts"
    awk 'FNR==NR{a[$1];next} ($1 in a)' "$list_file" "$workdir/bed/${sample}.bed" > "$workdir/indexing/bed/${sample}.${bc_name}.bed"
    awk '{ sum += ($4 - $3) } END { print sum+0 }' "$workdir/indexing/bed/${sample}.${bc_name}.bed" > "$workdir/indexing/bases/${sample}.${bc_name}.bases"
  done

  rm -f "$fasta_tmp"
  echo ">>> [${sample}] Demultiplexing done <<<"
}

##############################################################################
#                                                                            #
#                             MAIN EXECUTION                                 #
#                                                                            #
##############################################################################

# --- Part 1: Run core pipeline for all samples ---
echo "#####################################################"
echo "# Starting Core Mutation Pipeline for All Samples   #"
echo "#####################################################"
export -f align_and_filter make_bed_and_count split_and_extract join_strands process_sample
export threads_per_job splits ref_fa sam2tsv_jar workdir bam_dir summary_file
cd "$bam_dir"
sample_files=(*.bam)
cd "$workdir"
parallel --jobs "$jobs" process_sample ::: "${sample_files[@]}"

echo "All samples have completed the core pipeline."
echo
echo "===== Per‐sample timing summary ====="
printf "Sample\tRuntime (HH:MM:SS)\n"
sort "$summary_file"

# --- Part 2: Combine, merge, and apply final filters globally ---
echo
echo "#####################################################"
echo "# Merging All Samples & Applying Final Filters      #"
echo "#####################################################"
cat "$workdir/filtered-results/"*.FR.joined.v1 > "$workdir/v1-joined/all.joined"

sort -k6 -S80% "$workdir/v1-joined/all.joined" \
  | uniq -u -f5 \
  | awk '$2==0 || $2==16' \
  | awk '$4==0 || $4==16' \
  | awk '$3==60 && $5==60' \
  > "$workdir/v1-joined/all.joined.1"

cut -f1,6-9 "$workdir/v1-joined/all.joined.1" \
  | awk '{ print $2, $3, $4, $5, $1 }' OFS=$'\t' \
  | sort -k5 -S80% \
  | uniq -u -f4 \
  > "$workdir/v1-joined/all.joined.2"

echo "Global filtering complete. Final mutation file: v1-joined/all.joined.2"

# --- Part 3: Demultiplex all samples by barcode ---
echo
echo "#####################################################"
echo "# Starting Demultiplexing Process for All Samples   #"
echo "#####################################################"
export -f demultiplex_sample
export -v barcodes workdir bam_dir
sample_names=("${sample_files[@]%.bam}")
parallel --jobs "$jobs" demultiplex_sample ::: "${sample_names[@]}"

echo
echo "===== Pipeline finished successfully! ====="
