#!/bin/bash

# Use this script to submit a set of dreadds simulation jobs to PBS.
# It displays a tty form (using dialog) that allows the number and
# shape of jobs to be specified, as well as the parameter space to explore.
# The parameter space is specified by supplying lists of values for
# --species.max-dispersal-radius
# --species.niche-breadth
# --gene-flow-max-distance
#
# While the dreadds config file allow for the specification of
# multiple species, the simulations performed here are limited to one
# species. This is mainly to keep the user interface from being
# unwieldy.

dreadds=dreadds # the executable to run (in $PATH or specify absolute pathname)
utils=/g/data/ka2/util/bin # dialog lives here, if not available in $PATH.
                           # See https://invisible-island.net/dialog/
rc_file=$HOME/.dreaddsrc

here=$(cd $(dirname "$0") && pwd)
this=$here/$(basename "$0")

# Put $executable in the same directory as this script or in $HOME/bin
export PATH=$PATH:$here:$HOME/bin

usage() {
    [ "$1" ] && echo "$1" >&2
    echo "Usage: $0 config-file path_to_output_dir" >&2
    exit 1
}

readvars() {
    read n_jobs
    read qsub_args
    read n_runs
    read time_steps
    read mdr_vec
    read snb_vec
    read gfmd_vec
    read post_cmd    
}

readrc() {
    read _comment
    readvars
}

submit_jobs() {
    [ $# != 2 ] && usage
    config=${1:?}
    out_dir=${2:?}

    mkdir -p "$out_dir"
    cd "$out_dir" || exit 1
    out_dir=$PWD

    # Check we can read config file now we are in out_dir as that's the pwd of the job
    [ -r "$config" ] || usage "Cannot read $config"
    executable=$(which "$dreadds") || exit 1

    # Default job parameters
    n_jobs=10
    qsub_args="-l ncpus=16,walltime=1:00:00"
    n_runs=50 # per CPU 
    time_steps=100
    mdr_vec="3 4 5"
    snb_vec="5,40 6,50 7,60"
    gfmd_vec="5 6 7"
    post_cmd=""
    # Override from rc file
    [ -s "$rc_file" ] && readrc < "$rc_file"
    
    exec 3>&1
    args=$(PATH="$utils:$PATH" dialog \
               --form "Job configuration. Max 1 species. Use tab to select OK|cancel" \
	       0 0 0 \
               "Number of jobs to submit:" \
               1 2 "$n_jobs" 1 28 5 5 \
               "qsub arguments"  \
               2 2 "$qsub_args" 2 28 40 100 \
               "Runs per CPU:" \
               3 2 "$n_runs"  3 28  5 5 \
               "Time steps per run:" \
               4 2 "$time_steps" 4 28 5 5 \
               "species.max-dispersal-radius (space separated list)" \
               5 2 "$mdr_vec" 6 2 60 100 \
               "species.niche-breadth (space separated list of CSV)" \
               7 2 "$snb_vec" 8 2 60 100 \
               "gene-flow-max-distance  (space separated list)" \
               9 2 "$gfmd_vec" 10 2 60 100 \
               "Command to run in output directory after each run" \
               11 2 "$post_cmd" 12 2 60 100 \
	       2>&1 1>&3)
    x=$?
    echo
    [ "$x" = 0 ] || exit
    exec 3>&-

    # Update job args from dialog results and save for next time
    readvars <<< "$args"
    (echo "# Generated and used by $this"; echo "$args") > "$rc_file"
    
    for i in $(seq 1 $n_jobs); do
	# note "--" to run executable with args
        job_name="$(basename "$this")--$(basename "$config")--$i"
        echo "qsub $qsub_args -N $job_name ..."
        # Note that pbsdsh doesn't handle empty args, so we have to
        # pass non-empty strings for each arg.
	qsub \
            $qsub_args \
	    -N "$job_name" \
	    -- "$this" \
            "$executable" \
	    "$config" \
	    "${n_runs:-1}" \
	    "${time_steps:--}" \
	    "${mdr_vec:--}" \
	    "${snb_vec:--}" \
	    "${gfmd_vec:--}" \
	    "${out_dir:-.}" \
            "${post_cmd:-:}"
    done
}


spawn_tasks() {
    # This runs when we're invoked inside a PBS job Note that $this
    # refers to the original pathname of this script, not the
    # temporary PBS job script due to qsub … -- $executable [$arg … ]
    node_count=${PBS_NCPUS:?}
    for node_n in $(seq 1 $node_count); do
	pbsdsh -n $node_n --  "$this" --task "${PBS_JOBID}-$node_n" "$@" &
    done
    wait
}


run_simulation() {
    task_name=${1:?}; shift
    executable=${1:?}; shift
    config=${1:?}; shift
    n_runs=${1:?}; shift
    time_steps=${1:?}; shift
    mdr_vec=${1:?}; shift
    snb_vec=${1:?}; shift  # each value is comma separated tuple. v1,v2 v1,v2, ...
    gfmd_vec=${1:?}; shift
    out_dir=${1:?}; shift
    post_cmd=${1:?}; shift

    scratch_dir=${PBS_JOBFS:?}/$task_name
    mkdir -p "$scratch_dir"
    
    cd "$scratch_dir" || exit 1

    for mdr in $mdr_vec; do
	for snb in $snb_vec; do
	    for gfmd in $gfmd_vec; do
		for i in $(seq 1 $n_runs); do
		    dir="${task_name}-mdr$mdr-nb$snb-gf$gfmd-$i"
		    mkdir -p "$dir" || exit 1
		    args=()
                    [ "$time_steps" != - ] && args+=("--iterations=$time_steps")
		    [ "$mdr" != - ] && args+=("--species.max-dispersal-radius=$mdr")
		    [ "$snb" != - ] &&  args+=("--species.niche-breadth=$snb")
		    [ "$gfmd" != - ] && args+=("--gene-flow-max-distance=$gfmd")
                    echo Running  $executable \
			-c "$config" \
			-o "$dir" \
			"${args[@]}"
		    $executable \
			-c "$config" \
			-o "$dir" \
			"${args[@]}" || \
                        exit 1 # bail if any one fails as the rest will probably fail
                    (cd "$out_dir" && eval "$post_cmd")
		    tar cfz "$out_dir/$dir.tar.gz" "$dir" && \
			rm -r "$dir"
		done
	    done
	done
    done
}


if [ "$1" = "--task" ]; then
    # Worker process in job started by pbsdsh above.
    shift
    export PATH=$PBS_O_PATH
    run_simulation "$@"
elif [ "$PBS_JOBID" ]; then
    # Master process in job.
    export PATH=$PBS_O_PATH
    spawn_tasks "$@" # args from stuff after "$this" in qsub
else
    # Running from command line    
    submit_jobs "$@"
fi