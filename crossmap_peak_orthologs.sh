#!/bin/bash
# crossmap_peak_orthologs.sh: Map peaks between genome versions using CrossMap and HALPER.
# This script uses CrossMap to map peaks from one genome assembly to another,
# followed by HALPER to identify 1-1 orthologous regions.

#SBATCH --job-name=crossmap
#SBATCH --ntasks-per-core=1
#SBATCH --error=logs/crossmap_%A_%a.err.txt
#SBATCH --output=logs/crossmap_%A_%a.out.txt

# Default argument values (can be overridden by args)
OVERWRITE='FALSE'; NAME=''; SNP='FALSE'; PARALLEL='16';
# HALPER defaults
MIN_LEN=50; PROTECT_DIST=10; MAX_FRAC=1.5

# Constants (should not be overridden by args)
TMP_LABEL="tmp$(date +%s)"
# Make a temporary dir, portable on Linux and Mac
TMP_DIR=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')

function usage()
{
    echo "crossmap_peak_orthologs.sh takes in a bed/narrowpeak file"
    echo "and uses CrossMap + HALPER to map to the orthologous regions in a target genome."
    echo "Summits of narrowpeak files will be used to anchor peaks."
    echo "Example run call:"
    echo ""
    echo "sbatch crossmap_peak_orthologs.sh -b myPeaks.narrowPeak -o /dir/to/output/ -c rheMac10ToRheMac8.over.chain"
    echo ""
    echo "Required parameters:"   
    echo "--input-bed-file  -b FILENAME     bed or narrowPeak file, can be gzipped, {.bed/.bed.gz/.narrowpeak/.narrowpeak.gz}"
    echo "--output-dir      -o /OUT/DIR/    path to main output directory"
    echo "--chain-file      -c CHAINFILE    chain file for liftover e.g. {rheMac10ToRheMac8.over.chain}"
    echo ""
    echo "Optional parameters:"   
    echo "--name-base       -n NAME         Mapped file name prefix"
    echo "--force-overwrite -f              Whether to overwrite intermediate files"
    echo "--snp                             Map SNPs rather than peaks (change min peak length=1/protect dist=0)"
    echo "--parallel        -p THREADS      Number of threads for sorting (default: 16)"
    echo "-max_frac                         As in orthologFind.py"
    echo "-min_len                          As in orthologFind.py"
    echo "-protect_dist                     As in orthologFind.py"
    echo ""
}

# read in command line arguments
while [[ $1 != "" ]]; do
    case $1 in
        # the required parameters
        -b | --input-bed-file ) shift
                                BEDFILE=$1
                                echo "Bed file is $BEDFILE."
                                ;;
        -c | --chain-file )     shift
                                CHAINFILE=$1
                                echo "Chain file is $CHAINFILE."
                                ;;
        -o | --output-dir )     shift
                                OUTDIR=$1
                                echo "Output dir is $OUTDIR."
                                ;;
        # some optional parameters
        -n | --name-base )      shift
                                NAME=$1
                                echo "Run name is $NAME."
                                ;;
        -p | --parallel )       shift
                                PARALLEL=$1
                                echo "Using $PARALLEL threads for sorting."
                                ;;
        --snp )                 SNP='TRUE';
                                echo "Mapping SNPs rather than peaks. Setting min length=1 and protect_dist=0."
                                MIN_LEN=1
                                PROTECT_DIST=0
                                ;;
        -f | --force-overwrite ) OVERWRITE='TRUE'
                                echo "Overwrite intermediate files."
                                ;;
        -min_len )              shift
                                MIN_LEN=$1
                                echo "-min_len is $MIN_LEN."
                                ;;
        -protect_dist )         shift
                                PROTECT_DIST=$1
                                echo "-protect_dist is $PROTECT_DIST."
                                ;;
        -max_frac )             shift
                                MAX_FRAC=$1
                                echo "-max_frac is $MAX_FRAC."
                                ;;                                
        -preserve )             shift
                                PRESERVE=$1
                                echo "-preserve is $PRESERVE."
                                ;;
        -h | --help )           usage
                                exit 1
                                ;;
        *)
        usage
        exit 1
    esac
    shift
done

function detect_genomes_from_chain() {
    # Extract source and target genome info from chain filename
    CHAIN_BASE=$(basename $CHAINFILE)
    if [[ $CHAIN_BASE =~ ([A-Za-z0-9]+)To([A-Za-z0-9]+) ]]; then
        SRC_GENOME=${BASH_REMATCH[1]}
        TGT_GENOME=${BASH_REMATCH[2]}
        echo "Detected mapping from $SRC_GENOME to $TGT_GENOME from chain filename"
    else
        echo "WARNING: Could not detect source and target from chain filename."
        echo "Using generic names in output files."
        SRC_GENOME="source"
        TGT_GENOME="target"
    fi
}

function check_params()
{
    # check required parameters
    if [[ -z "$CHAINFILE" ]]; then
        echo 'Please set chain file for liftover.'; exit 1
    elif [[ -z "$OUTDIR" ]]; then
        echo 'Please set output directory.'; exit 1
    elif [[ -z "$BEDFILE" ]]; then
        echo 'Please set input bed file.'; exit 1
    elif [[ ${NAME} == '' ]]; then
        NAME=$(basename $BEDFILE | sed 's/.gz$//g;s/.narrowPeak$//g;s/.bed$//g')
    fi
    
    if [ ! -f "$CHAINFILE" ]; then
        echo "Chain file $CHAINFILE not found."; exit 1
    elif [ ! -f "$BEDFILE" ]; then
        echo "Input bed/narrowpeak file does not exist: $BEDFILE"; exit 1
    fi
}

function check_bed()
{
    # unzip gzipped bed file
    INPUTBED=${TMP_DIR}/${NAME}.input.bed
    if file --mime-type "${BEDFILE}" | grep -q gzip$ ;
        then zcat $BEDFILE > $INPUTBED
        else cat $BEDFILE > $INPUTBED
    fi
}

function num_null_values(){
    # $1: path to file
    # $2: column number (integer)
    # echo the number of "-1" values in this column of this file
    res=$(cut "-f$2-$2" $1 | grep "\-1" | wc -l)
    echo $res
}

function num_columns(){
    # $1: path to file
    # echo number of columns
    awk '{print NF; exit}' $1
}

function format_bed()
{
    ##########################################
    UNIQUEBED=${TMP_DIR}/${NAME}.unique.bed
    # check if basic 3 column bed file w/o name column
    if [[ $(awk '{print NF; exit}' ${INPUTBED}) -lt 4 ]]; 
        then echo "Bed file doesn't have name column. Adding"
        echo "USING CHR:START-END in the NAME column."
        awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $1 ":" $2 "-" $3 ":" int(($3 - $2) / 2), "0", "."}' $INPUTBED > $UNIQUEBED    
    # name column found, check if duplicate names
    elif [[ $(awk '++A[$4] > 1 { print "true"; exit 1 }' $INPUTBED) ]]; then
        echo "Non-unique bed peak names detected. Giving unique names now."
        if [[ ($(num_columns $INPUTBED) == 10) && ($(num_null_values $INPUTBED 10) == 0) ]]; 
            # If there are 10 columns, and the 10th column has no null values "-1", then use the summit values.
            then echo "Appending CHR:START-END:SUMMIT to NAME column."
            awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $1 ":" $2 "-" $3 ":" $10, $5, $6, $7, $8, $9, $10}' $INPUTBED > $UNIQUEBED
        else 
            # Otherwise, don't use the summit column.
            echo "Appending CHR:START-END to NAME column."
            awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $1 ":" $2 "-" $3, $5, $6, $7, $8, $9, $10}' $INPUTBED > $UNIQUEBED
        fi
    else 
        echo "Bed peak names are unique; moving onwards."
        cat $INPUTBED > $UNIQUEBED
    fi

    # make sure the scores column is numeric, strand column is +|-|.
    echo "Formatting bed score and strand columns."

    # Ensure that score column is numeric
    awk 'BEGIN {FS="\t"; OFS="\t"} { if($5 !~ /^[0-9]+\.{0,1}[0-9]*$/){$5 = 0} print;}' $UNIQUEBED > ${TMP_DIR}/${NAME}.unique2.bed

    # Ensure that strand column is "+", "-", or "."
    awk 'BEGIN {FS="\t"; OFS="\t"} { if($6 !~ /\+|\-|\./){$6 = "."} print;}' ${TMP_DIR}/${NAME}.unique2.bed > ${TMP_DIR}/${NAME}.unique3.bed
    mv ${TMP_DIR}/${NAME}.unique3.bed $UNIQUEBED; rm ${TMP_DIR}/${NAME}.unique2.bed

    # get the 6-column bed for broad peak
    SIMPLEBED=${TMP_DIR}/${NAME}.6col.bed
    cut -f1-6 $UNIQUEBED > $SIMPLEBED
}

function get_summits()
{
    # get summits from narrowpeak file
    SUMMITFILE=${TMP_DIR}/${NAME}.summits.bed
    echo "Extracting bed/narrowpeak summits for mapping."
    if [[ $(awk '{print NF; exit}' ${UNIQUEBED}) -lt 3 ]]; 
        then echo "Too few columns to be a bed file."; cleanup; exit 1
    elif [[ ($(num_columns $INPUTBED) == 10) && ($(num_null_values $INPUTBED 10) == 0) ]]; 
        # If there are 10 columns, and the 10th column has no null values "-1", then use the summit values.
        then echo "Summits detected. Using the summits."
        awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2+$10, $2+$10+1, $4, $5, $6}' $UNIQUEBED > $SUMMITFILE
    else [[ $(awk '{print NF; exit}' ${UNIQUEBED}) -gt 3 ]]
        # Otherwise, take the mean value between start and end of each peak
        echo "No summits found, taking the mean between start and end as the summits."
        awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, int(($2+$3)/2), int(($2+$3)/2)+1, $4, $5, $6}' $UNIQUEBED > $SUMMITFILE
    fi
}

function prepare_dirs()
{
    # Make output dir if it doesn't exist
    mkdir -p ${OUTDIR}
    mkdir -p ${OUTDIR}/logs
}

function map_summits()
{
    echo "Mapping summits using chain file."
    CROSSMAP_SUMMIT_FILE=${TMP_DIR}/${NAME}.summits.crossmap.unsorted.bed
    CROSSMAP_SUMMIT_SORTED=${TMP_DIR}/${NAME}.summits.crossmap.sorted.bed
    
    if [[ ! -f ${OUTDIR}/${NAME}.summits.crossmap.sorted.bed.gz || $OVERWRITE == 'TRUE' ]]; then
        # Map summits using CrossMap
        CrossMap bed $CHAINFILE $SUMMITFILE $CROSSMAP_SUMMIT_FILE
        sort --parallel $PARALLEL -k 1,1 -k2,2n $CROSSMAP_SUMMIT_FILE > $CROSSMAP_SUMMIT_SORTED
    else
        # Using previously mapped summits
        echo "Using previously mapped summits."
        gunzip -c ${OUTDIR}/${NAME}.summits.crossmap.sorted.bed.gz > $CROSSMAP_SUMMIT_SORTED
    fi
}

function map_peaks()
{
    echo "Mapping peaks using chain file."
    CROSSMAP_PEAK_FILE=${TMP_DIR}/${NAME}.peaks.crossmap.unsorted.bed
    CROSSMAP_PEAK_SORTED=${TMP_DIR}/${NAME}.peaks.crossmap.sorted.bed
    
    if [[ ! -f ${OUTDIR}/${NAME}.peaks.crossmap.sorted.bed.gz || $OVERWRITE == 'TRUE' ]]; then
        # Map peaks using CrossMap
        CrossMap bed $CHAINFILE $SIMPLEBED $CROSSMAP_PEAK_FILE
        sort --parallel $PARALLEL -k 1,1 -k2,2n $CROSSMAP_PEAK_FILE > $CROSSMAP_PEAK_SORTED
    else
        # Using previously mapped peaks
        echo "Using previously mapped peaks."
        gunzip -c ${OUTDIR}/${NAME}.peaks.crossmap.sorted.bed.gz > $CROSSMAP_PEAK_SORTED
    fi
}

function run_halper()
{    
    echo "Applying HALPER to identify 1-1 orthologous peaks in the target genome."
    OUTFILE=${TMP_DIR}/${NAME}.orthologs.tmp.narrowPeak
    FINAL_OUTFILE=${TMP_DIR}/${NAME}.${TGT_GENOME}.HALPER.narrowPeak
    
    args="-max_frac $MAX_FRAC -min_len $MIN_LEN -protect_dist $PROTECT_DIST \
    -tFile $CROSSMAP_PEAK_SORTED -sFile $CROSSMAP_SUMMIT_SORTED \
    -narrowPeak -qFile $UNIQUEBED -oFile $OUTFILE"
    
    # Add optional arguments if they are set
    if ! [[ -z ${PRESERVE+x} ]]; then
        # change PRESERVE from comma-delimited to space-delimited for args
        PRESERVE=$(echo $PRESERVE | tr "," " ")
        args="${args} -preserve $PRESERVE"
    fi
    
    python -m orthologFind $args
    echo "Mapped $(line_count $OUTFILE) peaks out of $(line_count $INPUTBED) total peaks."

    # Sort the output and copy to destination
    sort -k 1,1 -k2,2n $OUTFILE | uniq -u > $FINAL_OUTFILE
    
    # Archive and copy the results
    gzip --force $CROSSMAP_PEAK_SORTED $CROSSMAP_SUMMIT_SORTED $FINAL_OUTFILE
    
    # Copy the final files to the output directory
    rsync -aq ${CROSSMAP_PEAK_SORTED}.gz ${OUTDIR}/${NAME}.peaks.crossmap.sorted.bed.gz
    rsync -aq ${CROSSMAP_SUMMIT_SORTED}.gz ${OUTDIR}/${NAME}.summits.crossmap.sorted.bed.gz
    rsync -aq ${FINAL_OUTFILE}.gz ${OUTDIR}/${NAME}.${TGT_GENOME}.HALPER.narrowPeak.gz
    
    echo "Successfully saved outputs to ${OUTDIR}"
}

function line_count()
{
    res=$(wc -l < $1)
    echo $res
}

function map_snps()
{
    echo "Mapping SNP locations using chain file."
    CROSSMAP_PEAK_FILE=${TMP_DIR}/${NAME}.snps.crossmap.unsorted.bed
    CROSSMAP_PEAK_SORTED=${TMP_DIR}/${NAME}.snps.crossmap.sorted.bed
    
    if [[ ! -f ${OUTDIR}/${NAME}.snps.crossmap.sorted.bed.gz || $OVERWRITE == 'TRUE' ]]; then
        # Map SNPs using CrossMap
        CrossMap bed $CHAINFILE $SIMPLEBED $CROSSMAP_PEAK_FILE
        sort --parallel $PARALLEL -k 1,1 -k2,2n $CROSSMAP_PEAK_FILE > $CROSSMAP_PEAK_SORTED
    else
        # Using previously mapped SNPs
        echo "Using previously mapped SNPs."
        gunzip -c ${OUTDIR}/${NAME}.snps.crossmap.sorted.bed.gz > $CROSSMAP_PEAK_SORTED
    fi
    
    echo "Mapped $(line_count ${CROSSMAP_PEAK_SORTED}) SNPs out of $(line_count $INPUTBED) total SNPs."
    gzip --force ${CROSSMAP_PEAK_SORTED}
    rsync -aq ${CROSSMAP_PEAK_SORTED}.gz ${OUTDIR}/${NAME}.snps.crossmap.sorted.bed.gz
    echo "Successfully saved SNP mapping to ${OUTDIR}/${NAME}.snps.crossmap.sorted.bed.gz"
}

function cleanup()
{
    rm -rf $TMP_DIR
}

function main()
{
    # Prepare files and directories
    check_params
    detect_genomes_from_chain
    check_bed
    format_bed
    prepare_dirs
    
    # Perform mapping based on whether we're handling SNPs or peaks
    if [[ "${SNP}" == 'TRUE' ]]; then 
        map_snps
    else 
        get_summits
        map_summits
        map_peaks
        run_halper
    fi
    
    # Clean up temporary files
    cleanup
    
    echo "Process completed successfully."
}

# Execute the main function
main
