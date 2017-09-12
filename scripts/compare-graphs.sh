#!/usr/bin/env bash
# compare-graphs.sh: compare a set of graph against each other using toil-vg mapeval on AWS

set -e

# What toil-vg should we install?
TOIL_VG_PACKAGE="git+https://github.com/vgteam/toil-vg.git@f737ddd68313fc9bb7f84046d2adb8e81b42a8e4"

# What Toil should we use?
TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:3.10.1

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] OUTPUT_PATH KEYPAIR_NAME REGION_NAME GRAPH [GRAPH [GRAPH ...]] \n"
    printf "Options:\n\n"
    printf "\t-p PACKAGE\tUse the given Python package specifier to install toil-vg.\n"
    printf "\t-t CONTAINER\tUse the given Toil container in the cluster (default: ${TOIL_APPLIANCE_SELF}).\n"
    exit 1
}

while getopts "hp:t:" o; do
    case "${o}" in
        p)
            TOIL_VG_PACKAGE="${OPTARG}"
            ;;
        t)
            TOIL_APPLIANCE_SELF="${OPTARG}"
            ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "4" ]]; then
    # Too few arguments
    usage
fi

OUTPUT_PATH="${1}"
shift
KEYPAIR_NAME="${1}"
shift
REGION_NAME="${1}"
shift

GRAPH_NAMES=( )
while [[ "$#" -gt "0" ]]; do
    # Put all the args as graph names
    GRAPH_NAMES+=("$1")
    shift
done

# What's our unique run ID? Should be lower-case and start with a letter for maximum compatibility.
# See <https://gist.github.com/earthgecko/3089509>
RUN_ID="run$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 32 | head -n 1)"

# What cluster should we use?
CLUSTER_NAME="${RUN_ID}"

# Where do we keep our input files
INPUT_STORE="https://cgl-pipeline-inputs.s3.amazonaws.com/vg_cgl/bakeoff"

# Where do we save our results from the various jobs responsible for writing them?
OUTPUT_STORE="aws:us-west-2:cgl-pipeline-inputs/vg_cgl/comparison-script/runs/${RUN_ID}"
OUTPUT_STORE_URL="s3://cgl-pipeline-inputs/vg_cgl/comparison-script/runs/${RUN_ID}"

# Where do we store our jobs?
JOB_TREE="aws:us-west-2:${RUN_ID}"

echo "Running run ${RUN_ID} as ${KEYPAIR_NAME} to compare ${GRAPH_NAMES[*]} on ${REGION_NAME} into ${OUTPUT_PATH}"

function get_input_url() {
    # Prints the input URL to download for the given file name
    local BASE_FILENAME="${1}"
    shift
    echo "${INPUT_STORE}/${BASE_FILENAME}"
}

function get_graph_url() {
    # Prints the base URL for the given graph
    local BASE_GRAPHNAME="${1}"
    shift
    get_input_url "${BASE_GRAPHNAME}-${REGION_NAME}"
}

# Make sure we don't leave the cluster running on exit.
function kill_cluster() {
    set +e
    aws s3 rm --recursive "${OUTPUT_STORE_URL}"
    toil clean "${JOB_TREE}"
    toil destroy-cluster "${CLUSTER_NAME}" -z us-west-2a
}
trap kill_cluster EXIT

# Convert just graph stems to full base urls
GRAPH_URLS=()
for GRAPH_STEM in "${GRAPH_NAMES[@]}"; do
    GRAPH_URLS+=(`get_graph_url "${GRAPH_STEM}"`)
done

TOIL_APPLIANCE_SELF="${TOIL_APPLIANCE_SELF}" toil launch-cluster "${CLUSTER_NAME}" --nodeType=t2.medium -z us-west-2a "--keyPairName=${KEYPAIR_NAME}"

# We need to manually install git to make pip + git work...
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt update
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt install git -y

# For hot deployment to work, toil-vg needs to be in a virtualenv that can see the system Toil
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" virtualenv --system-site-packages venv

toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install numpy
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install scipy
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install scikit-learn
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install "${TOIL_VG_PACKAGE}"

# We need the master's IP to make Mesos go
MASTER_IP="$(toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" hostname -i)"

toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg mapeval \
    --fasta `get_input_url "${REGION_NAME}.fa"` \
    --index-bases "${GRAPH_URLS[@]}" \
    --gam-names "${GRAPH_NAMES[@]}" \
    --gam_input_reads `get_input_url "comparison-${REGION_NAME}.gam"` \
    --bwa --bwa-paired --vg-paired \
    --mapeval-threshold 200 \
    --realTimeLogging --logInfo \
    "${JOB_TREE}" \
    "${OUTPUT_STORE}" \
    `get_input_url "comparison-${REGION_NAME}.pos"` \
    --batchSystem mesos --provisioner=aws "--mesosMaster=${MASTER_IP}:5050" --nodeType=t2.large
    
mkdir -p ./out
aws s3 sync "${OUTPUT_STORE_URL}" "${OUTPUT_PATH}"

# Cluster, tree, and output will get cleaned up by the exit trap