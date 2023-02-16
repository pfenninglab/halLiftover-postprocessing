#!/bin/bash
# halper_map_peak_orthologs.sh: Map orthologs with halLiftover and postprocess with HALPER.
# Adapted from a script by BaDoi Phan:
# /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh

# ### NOTE slurm flags ###
#
# On the Lane cluster:
# sbatch -p pool1 --mem 10000 -w compute-1-11 halper_map_peak_orthologs.sh [flags]
# We use -w compute-1-11 because this node has a Cactus alignment in /scratch already.
#
# On the PSC:
# sbatch -p RM-shared --mem 2000 halper_map_peak_orthologs.sh [flags]
#
# To run on multiple target species in parallel using a slurm array job, please use:
# --array=1-[number of target species]
# Please note this may fail on the Lane cluster due to a limited length slurm queue.

#SBATCH --job-name=halliftover
#SBATCH --ntasks-per-core=1
#SBATCH --error=logs/halliftover_%A_%a.out.txt
#SBATCH --output=logs/halliftover_%A_%a.out.txt


# CACTUSFILE=/projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/10plusway-master.hal
source activate hal

# Default argument values (can be overridden by args)
CACTUSFILE=/scratch/cactus_alignments/241-mammalian-2020v2.hal
OVERWRITE='FALSE'; NAME=''; SNP=''; 
# HALPER defaults
MIN_LEN=50; PROTECT_DIST=10; MAX_FRAC=1.5

# Constants (should not be overridden by args)
TMP_LABEL="tmp$(date +%s)"
# Make a temporary dir, portable on Linux and Mac
# from https://unix.stackexchange.com/a/84980
TMP_HAL_DIR=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')
CACTUS_DIR=/scratch/cactus_alignments

function usage()
{
    echo "halper_map_peak_orthologs.sh takes in a bed/narrowpeak file in one genome coordinate"
    echo "and uses halLiftover + HALPER to map to the orthologous regions in a target genome."
    echo "Summits of narrowpeak files will be used to anchor peaks."
    echo "Example run call:"
    echo ""
    echo "sbatch halper_map_peak_orthologs.sh -b myPeaks.narrowPeak -o /dir/to/output/ -s sourceSpecies -t targetSpecies "
    echo ""
    echo "On the Lane cluster, please use the following sbatch flags to use a node that already has the cactus alignment:"
    echo "    sbatch -p pfen1 -w compute-1-11 --mem 10000 halper_map_peak_orthologs.sh [other args]"
    echo ""
    echo "Required parameters:"   
    echo "--input-bed-file  -b FILENAME     bed or narrowPeak file, can be gzipped, {.bed/.bed.gz/.narrowpeak/.narrowpeak.gz}"
    echo "--output-dir      -o /OUT/DIR/    path to main output directory"
    echo "--source-species  -s SPECIES      sourceSpecies in the Cactus file e.g. {Homo_sapiens, Mus_musculus, Macaca_mulatta}"
    echo "--target-species  -t SPECIES1,SPECIES2"
    echo "                                  comma separated list of targetSpecies in the Cactus file e.g. {Human, Mouse, Rhesus, etc}"
    echo "                                  OR path to a file with one target species per line"
    echo ""
    echo ""
    echo "Optional parameters:"   
    echo "--name-base       -n NAME         Mapped file name prefix after halliftover and Halper"
    echo "--cactus-file     -c HALFILE      Cactus multi-species alignment file {msaFile.hal}"
    echo "--force-overwrite -f              Whether to overwrite intermediate halliftover files."
    echo "--snp                             Map SNPs rather than peaks (change min peak length=1/protect dist=0)."
    echo "--keepChrPrefix                   Only keep orthologs where the chromosome name starts with this prefix."
    echo "--halPath                         Path to a halLiftover binary, in case you don't want to build halLiftover."
    echo "-max_frac                         As in orthologFind.py"
    echo "-min_len                          As in orthologFind.py"
    echo "-protect_dist                     As in orthologFind.py"
    echo ""
    echo ""
    echo "Note: genomes are whatever you get from running:"
    echo "  halStats --genomes file.hal"
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
        -s | --source-species )      shift
                                SOURCE=$1
                                echo "Source species $SOURCE."
                                ;;
        -t | --target-species )      shift
                                TARGETS=$1
                                echo "Target species $TARGETS."
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
        -c | --cactus-file )    shift
                                CACTUSFILE=$1
                                echo "Multi-species cactus alignment file is $CACTUSFILE."
                                ;;
        --snp )                 shift
                                SNP='TRUE';
                                echo "Mapping SNPs rather than peaks. Setting min length=1."
                                ;;
        -f | --force-overwrite ) shift
                                OVERWRITE='TRUE'
                                echo "Overwrite intermediate halLiftover files."
                                ;;
        --keepChrPrefix )       shift
                                KEEPCHRPREFIX=$1
                                echo "Keep chr prefix is $KEEPCHRPREFIX."
                                ;;
        --halPath )             shift
                                HALLIFTOVER_PATH=$1
                                echo "halLiftover binary path is $HALLIFTOVER_PATH."
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

function check_params()
{
    ############################################################
    # check required parameters, and some sensibilites. ########
    if [[ -z "$SOURCE" ]];
        then echo 'Please set source species e.g. {Human, Mouse, Rhesus, etc}.'; exit 1
    elif [[ -z "$TARGETS" ]];
        then echo 'Please set target species e.g. {Human, Mouse, Rhesus, etc}.'; exit 1
    elif [[ -z "$OUTDIR" ]];
        then echo 'Please set output directory.'; exit 1
    elif [[ "$SOURCE" == "$TARGETS" ]];
        then echo "Source and target species cannot be the same."; exit 1
    elif [[ ${NAME} == '' ]]
        then NAME=$(basename $BEDFILE | sed 's/.gz$//g;s/.narrowPeak$//g;s/.bed$//g')
    fi
    if [ ! -f "$CACTUSFILE" ]; then
        echo "Cactus file $CACTUSFILE not found."; exit 1
    elif [ ! -f "$BEDFILE" ]; then
        echo "Input bed/narrowpeak file does not exist: $BEDFILE"; exit 1
    fi
}

function check_bed()
{
    # unzip gzipped bed file
    INPUTBED=${TMP_HAL_DIR}/${NAME}.input.${SOURCE}.${TMP_LABEL}.bed
    if file --mime-type "${BEDFILE}" | grep -q gzip$ ;
        then zcat $BEDFILE > $INPUTBED
        else cat $BEDFILE > $INPUTBED
    fi
}

# if a path to a halLiftover binary was passed, then use that binary
# otherwise, use halLiftover directly from the command
function run_halLiftover()
{
    if [[ $HALLIFTOVER_PATH ]]; then
        $HALLIFTOVER_PATH "$@"
    else
        halLiftover "$@"
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
    UNIQUEBED=${TMP_HAL_DIR}/${NAME}.unique.${SOURCE}.${TMP_LABEL}.bed
    # check if basic 3 column bed file w/o name column
    if [[ $(awk '{print NF; exit}' ${INPUTBED}) -lt 4 ]]; 
        then echo "Bed file doesn't have name column. Adding"
        echo "USING CHR:START-END in the NAME column."
        awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $1 ":" $2 "-" $3 ":" int(($3 - $2) / 2), "0", "."}' $INPUTBED > $UNIQUEBED    
    # name column found, check if duplicate names
    elif [[ $(awk '++A[$4] > 1 { print "true"; exit 1 }' $INPUTBED) ]]; then
        echo "Non-unique bed peak names detected. Giving unique names now."
        echo "Bed file doesn't have name column. Adding"
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

    #################################################################
    # make sure the scores column is numeric, strand column is +|-|.
    echo "Formatting bed score and strand columns."

    # Ensure that score column is numeric
    # Adapted from https://stackoverflow.com/a/33706695
    # Regex: 1 or more digits, 0 or 1 period, 0 or more digits
    # If score is not numeric, then set it to 0, otherwise keep it
    awk 'BEGIN {FS="\t"; OFS="\t"} { if($5 !~ /^[0-9]+\.{0,1}[0-9]*$/){$5 = 0} print;}' $UNIQUEBED > ${TMP_HAL_DIR}/${NAME}.unique2.${SOURCE}.${TMP_LABEL}.bed

    # Ensure that strand column is "+", "-", or "."
    awk 'BEGIN {FS="\t"; OFS="\t"} { if($6 !~ /\+|\-|\./){$6 = "."} print;}' ${TMP_HAL_DIR}/${NAME}.unique2.${SOURCE}.${TMP_LABEL}.bed > ${TMP_HAL_DIR}/${NAME}.unique3.${SOURCE}.${TMP_LABEL}.bed
    mv ${TMP_HAL_DIR}/${NAME}.unique3.${SOURCE}.${TMP_LABEL}.bed $UNIQUEBED; rm ${TMP_HAL_DIR}/${NAME}.unique2.${SOURCE}.${TMP_LABEL}.bed

    ######################################
    # get the 6-column bed for broad peak
    SIMPLEBED=${TMP_HAL_DIR}/${NAME}.6col.${SOURCE}.${TMP_LABEL}.bed; cut -f1-6 $UNIQUEBED > $SIMPLEBED
}

function get_summits()
{
    ##################################################################
    # get summits from narrowpeak file or max cactus depth from cactus
    SUMMITFILE=${TMP_HAL_DIR}/${NAME}.summit.${SOURCE}.${TMP_LABEL}.bed
    echo "Extracting bed/narrowpeak summits to halLiftover."
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
    ########################################################################
    # Make output dir if it doesn't exist
    mkdir -p ${OUTDIR}

    # If the /scratch dir exists, copy the cactus file there, then use the cactus file in /scratch
    # Otherwise, use the cactus file in the place the user provided.
    if [ -d "/scratch" ]; then
        echo "/scratch dir found; using cactus file in /scratch dir"
        if [ ! -f "$CACTUS_DIR/$(basename $CACTUSFILE)" ]; then
            echo "Copying cactus file to /scratch dir. This will take a long time."
            rsync -Paq $CACTUSFILE $CACTUS_DIR
        fi
        CACTUSFILE=$CACTUS_DIR/$(basename $CACTUSFILE)
    else
        echo "no /scratch dir, using cactus file $CACTUSFILE"
    fi
}

function lift_summits()
{
    ######################################################
    echo "Lifting ${SOURCE} summits to ${TARGET} genome."
    HALLIFTEDSFILE=${NAME}.${SOURCE}To${TARGET}.halLiftover.sFile.bed
    if [[ ! -f ${OUTDIR}/${HALLIFTEDSFILE}.gz || $OVERWRITE == 'TRUE' ]]
        # hal-liftover the summits de novo
        then run_halLiftover $CACTUSFILE $SOURCE $SUMMITFILE $TARGET ${TMP_HAL_DIR}/$HALLIFTEDSFILE
    else
        # using previously hal-liftover summits
        echo "The file ${HALLIFTEDSFILE}.gz exist without permission to overwrite."
        echo "Unzipping and using the previous halliftover sFile."
        gunzip -c ${OUTDIR}/${HALLIFTEDSFILE}.gz > ${TMP_HAL_DIR}/${HALLIFTEDSFILE}
    fi
}

function lift_peaks()
{
    ###################################################
    echo "Lifting ${SOURCE} peaks to ${TARGET} genome."
    HALLIFTEDTFILE=${NAME}.${SOURCE}To${TARGET}.halLiftover.tFile.bed
    if [[ ! -f ${OUTDIR}/${HALLIFTEDTFILE}.gz || $OVERWRITE == 'TRUE' ]]
        # hal-liftover the peaks de novo
        then run_halLiftover $CACTUSFILE $SOURCE $SIMPLEBED $TARGET ${TMP_HAL_DIR}/$HALLIFTEDTFILE
    else
        # using previously hal-liftover peaks
        echo "The file ${HALLIFTEDTFILE}.gz exists without permission to overwrite."
        echo "Unzipping and using the previous halliftover tFile."
        gunzip -c ${OUTDIR}/${HALLIFTEDTFILE}.gz > ${TMP_HAL_DIR}/${HALLIFTEDTFILE}
    fi
}

function run_halper()
{    
    ########################################################################
    echo "Apply HALPER to identify 1-1 orthologous in the ${TARGET} genome."
    OUTFILE=${TMP_HAL_DIR}/${NAME}.${SOURCE}To${TARGET}.orthologs.${SOURCE}To${TARGET}.${TMP_LABEL}.narrowPeak
    args="-max_frac $MAX_FRAC -min_len $MIN_LEN -protect_dist $PROTECT_DIST \
    -tFile ${TMP_HAL_DIR}/$HALLIFTEDTFILE -sFile ${TMP_HAL_DIR}/$HALLIFTEDSFILE \
    -narrowPeak -qFile $UNIQUEBED -oFile $OUTFILE"
    # Add optional arguments if they are set
    if ! [[ -z ${KEEPCHRPREFIX+x} ]]; then
        args="${args} -keepChrPrefix $KEEPCHRPREFIX"
    fi
    if ! [[ -z ${PRESERVE+x} ]]; then
        # change PRESERVE from comma-delimited to space-delimited for args
        PRESERVE=$(echo $PRESERVE | tr "," " ")
        args="${args} -preserve $PRESERVE"
    fi
    python -m orthologFind $args
    echo "Mapped $(line_count $OUTFILE) peaks out of $(line_count $INPUTBED) total peaks."

    #################################################################
    # sort hal-liftOvered mouse peaks to hg38 and copy to destination
    # if [[ ${SNP} == "TRUE" ]]; then rm ${TMP_HAL_DIR}/$HALLIFTEDTFILE; fi
    OUTFILE2=${TMP_HAL_DIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak
    sort -k 1,1 -k2,2n $OUTFILE | uniq -u > $OUTFILE2
    gzip --force ${TMP_HAL_DIR}/$HALLIFTEDTFILE ${TMP_HAL_DIR}/$HALLIFTEDSFILE $OUTFILE2
    rsync --remove-source-files -Paq ${TMP_HAL_DIR}/${NAME}.${SOURCE}To${TARGET}*.gz ${OUTDIR}/
    rm ${TMP_HAL_DIR}/${NAME}*${SOURCE}To${TARGET}.${TMP_LABEL}*
    echo "Done."; 
}

function line_count()
{
    res=$(wc -l < $1)
    echo $res
}

function map_snps()
{
    ####################################################################
    echo "Only the ${SOURCE} bed files of the SNPs to ${TARGET} genome."
    HALLIFTEDTFILE=${NAME}.${SOURCE}To${TARGET}.halLiftover.tFile.bed
    OUTFILE=${TMP_HAL_DIR}/${NAME}.${SOURCE}To${TARGET}.snp.bed
    if [[ ! -f ${OUTDIR}/${HALLIFTEDTFILE}.gz || $OVERWRITE == 'TRUE' ]]
        # hal-liftover the peaks de novo
        then run_halLiftover $CACTUSFILE $SOURCE $SIMPLEBED $TARGET $OUTFILE
        sort -k 1,1 -k2,2n $OUTFILE | uniq -u > ${TMP_HAL_DIR}/$HALLIFTEDTFILE
        rm $OUTFILE
    else
        # using previously hal-liftover peaks
        echo "The file ${HALLIFTEDTFILE}.gz exists without permission to overwrite."
        echo "Unzipping and using the previous halliftover tFile."
        gunzip -c ${OUTDIR}/${HALLIFTEDTFILE}.gz > ${TMP_HAL_DIR}/${HALLIFTEDTFILE}
    fi
    echo "Mapped $(line_count ${TMP_HAL_DIR}/$HALLIFTEDTFILE) SNPs out of $(line_count $INPUTBED) total SNPs."
    ######################################
    gzip --force ${TMP_HAL_DIR}/$HALLIFTEDTFILE
    rsync --remove-source-files -Paq ${TMP_HAL_DIR}/$HALLIFTEDTFILE.gz ${OUTDIR}/
}

function cleanup()
{
    rm -r $TMP_HAL_DIR
}

function liftover_target()
{
    if [[ "${SNP}" == 'TRUE' ]]; then map_snps
    else get_summits; lift_summits; lift_peaks; run_halper; fi
}

##############################
# prepare bed files for input
check_params
check_bed
format_bed
prepare_dirs

##########################################
# perform liftover for each target species

# Get array of targets, either from comma-delimited string, or by reading from a file.
if [ ! -f "$TARGETS" ]; then
    # interpret TARGETS as a comma-delimited string and read elements into array
    TARGETS=( $(echo $TARGETS | tr "," "\n") )
else
    # TARGETS is a file, read lines from it and store as array
    echo "Reading in targets line-by-line from $TARGETS"
    unset -v lines
    while IFS= read -r; do
        lines+=("$REPLY")
    done <$TARGETS
    [[ $REPLY ]] && lines+=("$REPLY")
    TARGETS=("${lines[@]}")
fi

# If this is not a slurm array job, then loop over targets sequentially.
# Otherwise, process the single target that corresponds to this array task ID.
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    for TARGET in "${TARGETS[@]}"; do
        liftover_target
    done
else
    # access the target that corresponds to this array ID (0-indexed)
    TARGET=${TARGETS[$(($SLURM_ARRAY_TASK_ID - 1))]}
    liftover_target
fi

##################
# final clean up
cleanup
