import argparse, sys, os, stat
from path import path
from util import merge_items

SCRIPT = "/software/R-3.1.0/bin/Rscript /nfs/users/nfs_s/sd11/repo/dirichlet/dp_combined/RunDP_pipeline.R"

def _write_script(filename, lines):
    outf = open(filename, "w")
    for line in lines: outf.write(line+"\n")
    outf.close()
    # Make executable
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)
    
def _get_r_wrapper_lines(phase):
    assert phase=="cons" or phase=="trees", "The phase parameter can only be trees or cons"
    
    # Create the run command, which differs per phase
    if phase == "cons":
        run_cmd = merge_items(["${CMD} ${DATASET} ${ITERS} ${BURNIN} ${DATPATH} ${PURITYFILE} \"",phase,"\" ${PARALLEL} ${THREADS} ${BINSIZE} \"1\" ${NBLOCKS}"], sep="")
    else:
        run_cmd = merge_items(["${CMD} ${DATASET} ${ITERS} ${BURNIN} ${DATPATH} ${PURITYFILE} \"",phase,"\" ${PARALLEL} ${THREADS} ${BINSIZE} $LSB_JOBINDEX"], sep="")
    
    # Create the other lines
    lines = ["#$LSB_JOBINDEX",
             "DATASET=$1",
             "ITERS=$2",
             "BURNIN=$3",
             "DATPATH=$4",
             "PURITYFILE=$5",
             "ANALYSIS=$6",
             "PARALLEL=$7",
             "THREADS=$8",
             "BINSIZE=$9",
             "NBLOCKS=${10}",
             merge_items(["CMD=",SCRIPT], sep=""),
             run_cmd]
    return(lines)

def generateDPrunScript(run_dir, dp_in_dir, dp_master_file, projectname):
    
    '''
    This function will create a series of shell scripts that make life easier when running dirichlet clustering pipelines.
    
    within the directory that is supplied by the run_dir parameter this will be created:
        - submit.block.sh : for the block parallel tree based method, when run this creates two LSF jobs for the trees and cons phase respectively
        - submit.nd.sh : for running the nD clustering, when run this creates a single LSF job
        - resubmit.block.sh : for resubmitting the cons step of the tree based method, when run this creates a single LSF job for the cons phase
        These would ideally never be called from the command line and only through the other scripts:
        - R wrapper script RunBlockTreeDP_trees.sh 
        - R wrapper script RunBlockTreeDP_cons.sh
    '''
    
    
    ''' Write block parallel run scripts '''
    # Write R wrappers first
    lines = _get_r_wrapper_lines("trees")
    _write_script(path.joinpath(run_dir, "RunBlockTreeDP_trees.sh"), lines)
    lines = _get_r_wrapper_lines("cons")
    _write_script(path.joinpath(run_dir, "RunBlockTreeDP_cons.sh"), lines)
    
    # Write wrapper around R wrappers
    '''
    Example
        SAMPLE=$1
        QUEUE="basement"
        JOBNAME="tree_pros"
        ACCEPTEDHOSTS="-m vr-2-3-10 vr-2-3-02 vr-2-3-05 vr-2-3-08 vr-2-3-15 vr-2-3-13"
        #CMD="Rscript ~/repo/dirichlet/dp_combined/RunDP_simulated.R"
        
        PARAMS="${SAMPLE} 200 30 /lustre/scratch110/sanger/sd11/dirichlet/prostate_mets/Data/ /lustre/scratch110/sanger/sd11/dirichlet/prostate_mets/prostate_mets.txt tree_dp true 5 0.05"
        NBLOCKS=10
        
        MEMTREE=75000
        MEMCONS=75000
        
        bsub -M ${MEMTREE} -R "select[mem>${MEMTREE}] rusage[mem=${MEMTREE}] span[hosts=1]" -n 5 -J "${JOBNAME}_t[1-${NBLOCKS}]" -q "${QUEUE}" -o $PWD/logs/${JOBNAME}_t.%J.%I.out -e $PWD/logs/${JOBNAME}_t.%J.%I.err "./RunBlockTreeDP_trees.sh ${PARAMS}"
        bsub -w"ended(${JOBNAME}_t[1-${NBLOCKS}])" -M ${MEMCONS} -R "select[mem>${MEMCONS}] rusage[mem=${MEMCONS}]" -J "${JOBNAME}_c" -q "${QUEUE}" -o $PWD/logs/${JOBNAME}_c.%J.out -e $PWD/logs/${JOBNAME}_c.%J.err "./RunBlockTreeDP_cons.sh ${PARAMS} ${NBLOCKS}"
    '''
    lines = ["SAMPLE=$1", 
             "QUEUE=\"basement\"",
             merge_items(["JOBNAME=\"tree_",projectname,"\""], sep=""),
             merge_items(["PARAMS=\"${SAMPLE} 200 30", dp_in_dir, dp_master_file, "tree_dp true 5 0.05\""]),
             "NBLOCKS=\"10\"",
             "MEMTREE=\"15000\"",
             "MEMCONS=\"15000\"",
             merge_items(["bsub -M ${MEMTREE} -R \"select[mem>${MEMTREE}] rusage[mem=${MEMTREE}] span[hosts=1]\" -n 5 -J \"${JOBNAME}_t[1-${NBLOCKS}]\" -q \"${QUEUE}\" -o $PWD/logs/${JOBNAME}_t.%J.%I.out -e $PWD/logs/${JOBNAME}_t.%J.%I.err \"./RunBlockTreeDP_trees.sh ${PARAMS}\""]),
             merge_items(["bsub -w\"ended(${JOBNAME}_t[1-${NBLOCKS}])\" -M ${MEMCONS} -R \"select[mem>${MEMCONS}] rusage[mem=${MEMCONS}]\" -J \"${JOBNAME}_c\" -q \"${QUEUE}\" -o $PWD/logs/${JOBNAME}_c.%J.out -e $PWD/logs/${JOBNAME}_c.%J.err \"./RunBlockTreeDP_cons.sh ${PARAMS} ${NBLOCKS}\""])]
    _write_script(path.joinpath(run_dir, "submit.block.sh"), lines)


    ''' Write resubmit block parallel run script 
    
    Example:
        resubmit.prostate_mets.block.sh
        
        SAMPLE=$1
        QUEUE="long"
        JOBNAME="tree_pros"
        ACCEPTEDHOSTS="-m vr-2-3-10 vr-2-3-02 vr-2-3-05 vr-2-3-08 vr-2-3-15 vr-2-3-13"
        #CMD="Rscript ~/repo/dirichlet/dp_combined/RunDP_simulated.R"
        
        PARAMS="${SAMPLE} 200 30 /lustre/scratch110/sanger/sd11/dirichlet/prostate_mets/Data/ /lustre/scratch110/sanger/sd11/dirichlet/prostate_mets/prostate_mets.txt tree_dp true 5 0.05"
        NBLOCKS=10
        
        MEMCONS=75000
        
        bsub -M ${MEMCONS} -R "select[mem>${MEMCONS}] rusage[mem=${MEMCONS}]" -J "${JOBNAME}_c" -q "${QUEUE}" -o $PWD/logs/${JOBNAME}_c.%J.out -e $PWD/logs/${JOBNAME}_c.%J.err "./RunBlockTreeDP_cons.sh ${PARAMS} ${NBLOCKS}"
    '''
    lines = ["SAMPLE=$1", 
             "QUEUE=\"basement\"",
             merge_items(["JOBNAME=\"tree_",projectname,"\""], sep=""),
             merge_items(["PARAMS=\"${SAMPLE} 200 30", dp_in_dir, dp_master_file, "tree_dp true 5 0.05\""]),
             "NBLOCKS=\"10\"",
             "MEMCONS=\"15000\"",
             merge_items(["bsub -w\"ended(${JOBNAME}_t[1-${NBLOCKS}])\" -M ${MEMCONS} -R \"select[mem>${MEMCONS}] rusage[mem=${MEMCONS}]\" -J \"${JOBNAME}_c\" -q \"${QUEUE}\" -o $PWD/logs/${JOBNAME}_c.%J.out -e $PWD/logs/${JOBNAME}_c.%J.err \"./RunBlockTreeDP_cons.sh ${PARAMS} ${NBLOCKS}\""])]
    _write_script(path.joinpath(run_dir, "resubmit.block.sh"), lines)
    
    
    ''' Write the nD run script
    
    Example:
        submit.prostate_mets.nd.sh
    
        QUEUE="normal"
        JOBNAME="nd_pros"
        # -m "vr-2-3-10 vr-2-3-02 vr-2-3-05 vr-2-3-08 vr-2-3-15 vr-2-3-13"
        CMD="Rscript ~/repo/dirichlet/dp_combined/RunDP_pipeline.R"
        PARAMS="1 1250 250 /lustre/scratch110/sanger/sd11/dirichlet/prostate_mets/Data/ /lustre/scratch110/sanger/sd11/dirichlet/prostate_mets/prostate_mets.txt nd_dp false 1 NA 1 1"
        MEMORY="17000"
        bsub -J ${JOBNAME} -q ${QUEUE} -M ${MEMORY} -R 'span[hosts=1] select[mem>'${MEMORY}'] rusage[mem='${MEMORY}']' -o $PWD/logs/${JOBNAME}.%J.out -e $PWD/logs/${JOBNAME}.%J.err "${CMD} ${PARAMS}"
    
    '''
    lines = ["SAMPLE=$1",
             "QUEUE=\"normal\"",
             merge_items(["JOBNAME=\"nd_",projectname,"\""], sep=""),
             merge_items(["CMD=\"", SCRIPT, "\""], sep=""),
             merge_items(["PARAMS=\"${SAMPLE} 1250 250", dp_in_dir, dp_master_file, "nd_dp false 1 NA 1 1\""]),
             "MEMORY=\"15000\"",
             "bsub -J ${JOBNAME} -q ${QUEUE} -M ${MEMORY} -R 'span[hosts=1] select[mem>'${MEMORY}'] rusage[mem='${MEMORY}']' -o $PWD/logs/${JOBNAME}.%J.out -e $PWD/logs/${JOBNAME}.%J.err \"${CMD} ${PARAMS}\""]
    _write_script(path.joinpath(run_dir, "submit.nd.sh"), lines)


def main(argv):
    parser = argparse.ArgumentParser(prog='generateDPrunScript',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Full path to dirichlet master file")
    parser.add_argument("-d", required=True, type=str, help="Full path to directory where dirichlet input files are stored")
    parser.add_argument("-r", required=True, type=str, help="Full path to directory where pipelines will be run")
    parser.add_argument("-p", required=True, type=str, help="Short project name to be used in the bsub command names for each tracking through bjobs+grep")
    args = parser.parse_args()
    generateDPrunScript(args.r, args.d, args.i, args.p)
    
if __name__ == '__main__':
    main(sys.argv[0:])