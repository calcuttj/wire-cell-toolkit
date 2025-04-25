#!/usr/bin/env bats

# Run tests related to applying imaging to a rich uboone event.

bats_load_library "wct-bats.sh"

# Define the archive format file extension for temporary files.  Can
# be "tar" or "tar.gz" to exercise WCT streaming compression.  Use
# just "tar" to make running faster.
arrext="tar.gz"

dag_file="dag.json"
log_file="wire-cell.log"
img_numpy_file="clusters-numpy.$arrext"
img_json_file="clusters-json.$arrext"

setup_file () {
    skip_if_no_input

    cd_tmp file

    # was: uboone/celltree/celltreeOVERLAY-event6501.tar.bz2
    local infile="$(input_file frames/celltreeOVERLAY-event6501.tar.bz2)"
    [[ -s "$infile" ]]

    local cfg_file="$(relative_path test-wct-uboone-img.jsonnet)"
    [[ -s "$cfg_file" ]]

    run_idempotently -s "$cfg_file" -t "$dag_file" -- \
                     compile_jsonnet "$cfg_file" "$dag_file" \
                     -A infile="$infile" -A outpat="clusters-%s.$arrext" -A formats="json,numpy"
    [[ -s "$dag_file" ]]

    run_idempotently -s "$dag_file" -s "$infile" \
                     -t "$img_numpy_file" -t "$img_json_file" -t "$log_file" -- \
                     wire-cell -l "$log_file" -L debug -c "$dag_file"
    [[ -s "$log_file" ]]
}


# bats test_tags=history
@test "make history" {
    cd_tmp file

    [[ -s "$img_numpy_file" ]]
    [[ -s "$img_json_file" ]]
    if [ "$arrext" = "tar" ] ; then
        gzip -k "$img_numpy_file"
        gzip -k "$img_json_file"
        saveout -c history "${img_numpy_file}.gz"
        saveout -c history "${img_json_file}.gz"
    else
        saveout -c history "$img_numpy_file"
        saveout -c history "$img_json_file"
    fi
}

# bats test_tags=dotify
@test "dotify dag" {
    cd_tmp file

    local dag_viz="dag.pdf"
    run_idempotently -s "$dag_file" -t "$dag_viz" -- \
                     dotify_graph "$dag_file" "$dag_viz"
    saveout -c plots "$dag_viz"
}

@test "check wire-cell log file" {
    cd_tmp file

    local errors="$(egrep ' W | E ' $log_file)"
    info "$errors"
    [[ -z "$errors" ]]
}

function do_blobs () {
    local name="$1"; shift
    local what="$1"; shift
    local args=( $@ )

    cd_tmp file

    declare -A logs

    for fmt in json numpy
    do
        local log="${what}-${name}-${fmt}.log"
        logs[$fmt]="$log"

        local dat="clusters-${fmt}.$arrext"
        run_idempotently -s "$dat" -t "$log" -- wcpy img $what "${args[@]}" -o "$log" "$dat"
        [[ -s "$log" ]]
    done

    local delta="$(diff -u ${logs[*]})"
    echo "$delta"
    [[ -z "$delta" ]]
}

@test "inspect blobs quietly" {
    do_blobs "quietly" inspect
}
@test "inspect blobs verbosely" {
    do_blobs "verbosely" inspect --verbose

    local lfile=inspect-verbosely-numpy.log
    [[ -s "$lfile" ]]
    local got=$(grep 'nn for m' "$lfile" | grep -v 'b=1\b')
    info "$got"
    [[ -n "$got" ]]

}
@test "dump blobs" {
    do_blobs "nominal" dump-blobs
}

# bats test_tags=plots
@test "plot blobs" {
    cd_tmp file
    run_idempotently -t blob-views.pdf -s clusters-numpy.$arrext -- \
                     wcpy img plot-blobs -o blob-views.pdf --plot views clusters-numpy.$arrext 
}

@test "valid cluster graph schema" {
    cd_tmp file

    local moo="$(which moo)"
    if [ -z "$moo" ] ; then
        skip "No moo available to check schema (see https://brettviren.github.io/moo/moo.html#install)"
    fi

    local sfile=$(relative_path ../../aux/docs/cluster-graph-schema.jsonnet)
    [[ -s "$sfile" ]]

    local tfile=clusters-json.$arrext
    [[ -s "$tfile" ]]
    local dfile="cluster_6501_graph.json"

    run_idempotently -s "$tfile" -t "$dfile" -- \
                     tar -xf $tfile
    [[ -s "$dfile" ]]
    touch -r "$tfile" "$dfile"  # override in-tar attributes

    local lfile="cluster_6501_graph.log"

    run_idempotently -s "$sfile" -s "$dfile" -t "$lfile" -- \
                     $moo validate -o $lfile -t wirecell.cluster.Cluster -s "$sfile" "$dfile"
}

@test "plot view projections" {
    run_idempotently -s clusters-json.$arrext -t found-projection.pdf -- \
                     wcpy img blob-activity-mask \
                     -o found-projection.pdf --found clusters-json.$arrext --vmin=0.001 \
                     --channel-lines 296,298-315,317,319-327,336,337,343-345,348-351,376-400,410-445,447-484,501-503,505-520,522-524,536-559,561-592,595-598,600-632,634-652,654,656-671,864-911,3936-3983,7136-7199,7201-7214,7216-7263



}

@test "check view projections" {
    run_idempotently -s clusters-json.$arrext -t found-projection.json -- \
                     wcpy img blob-activity-stats \
                     -f json -o found-projection.json clusters-json.$arrext
    skip_if_missing jq
    info "$(cat found-projection.json)"
    [[ "$(cat found-projection.json | jq '.pqtot > 0.75')" = "true" ]]
    [[ "$(cat found-projection.json | jq '.pqfound > 0.90')" = "true" ]]
    [[ "$(cat found-projection.json | jq '.nchan')" = "8256" ]] 
    [[ "$(cat found-projection.json | jq '.nslice')" = "2398" ]] 
    
}

function roundtrip () {
    local ifmt=$1; shift
    local ofmt=$1; shift
    local ifile=clusters-${ifmt}.$arrext
    local ofile=clusters-${ifmt}-${ofmt}.$arrext
    local dag=roundtrip-${ifmt}-${ofmt}
    local lfile=$dag.log
    local cfg=$(relative_path test-wct-uboone-img-roundtrip.jsonnet)

    cd_tmp file

    run_idempotently -s "$cfg" -t "${dag}.json" -- \
                     compile_jsonnet $cfg ${dag}.json \
                     -A infile=$ifile \
                     -A outfile=$ofile \
                     -A infmt=$ifmt \
                     -A outfmt=$ofmt
    
    run_idempotently -s "${dag}.json" -t "${dag}.png" -- \
                     dotify_graph "${dag}.json" "${dag}.png"

    run_idempotently -s "${dag}.json" -s "$ifile" -t "$ofile" -- \
                     wire-cell -l $lfile -L debug -c "${dag}.json"

    local errors="$(egrep ' W | E ' $lfile)"
    echo "$errors"
    [[ -z "$errors" ]]
    
    ilog=roundtrip-inspect-${ifmt}2${ofmt}-ilog-${ifmt}.log
    olog=roundtrip-inspect-${ifmt}2${ofmt}-olog-${ofmt}.log
    run_idempotently -s "$ifile" -t "$ilog" -- \
                     wcpy img inspect -o "$ilog" "$ifile"
    run_idempotently -s "$ofile" -t "$olog" -- \
                     wcpy img inspect -o $olog $ofile

    local delta="$(diff -u $ilog $olog)"
    echo "$delta"
    [[ -z "$delta" ]]
}

#@test "round trip json to json" {
#    roundtrip json json
#}
#@test "round trip numpy to numpy" {
#    roundtrip numpy numpy
#}
#@test "round trip numpy to json" {
#    roundtrip numpy json
#}
#@test "round trip json to numpy" {
#    roundtrip json numpy
#}

function roundtrip2 () {
    local ifmt=$1; shift
    local ofmt=$1; shift
    local ifile=clusters-${ifmt}.$arrext
    local ofile1=clusters2-${ifmt}-${ofmt}.$arrext
    local ofile2=clusters2-${ifmt}-${ofmt}-${ifmt}.$arrext


    local dag1=roundtrip2-${ifmt}-${ofmt}
    local dag2=roundtrip2-${ifmt}-${ofmt}-${ifmt}

    local lfile1=${dag1}.log
    local lfile2=${dag2}.log
    local cfg=$(relative_path test-wct-uboone-img-roundtrip.jsonnet)

    cd_tmp file

    run_idempotently -s "$cfg" -t "${dag1}.json" -- \
                     compile_jsonnet ${cfg} ${dag1}.json \
                     -A infile=$ifile \
                     -A outfile=$ofile1 \
                     -A infmt=$ifmt \
                     -A outfmt=$ofmt

    run_idempotently -s "${dag1}.json" -t "${dag1}.png" -- \
                     dotify_graph "${dag1}.json" "${dag1}.png"

    run_idempotently -s "${dag1}.json" -s "$ifile" -t "$ofile1" -- \
                     wire-cell -l $lfile1 -L debug -c ${dag1}.json
    local errors="$(egrep ' W | E ' $lfile1)"
    echo "$errors"
    [[ -z "$errors" ]]


    run_idempotently -s "$cfg" -t "${dag2}.json" -- \
                     compile_jsonnet ${cfg} ${dag2}.json \
                     -A infile=$ofile1 \
                     -A outfile=$ofile2 \
                     -A infmt=$ofmt \
                     -A outfmt=$ifmt

    run_idempotently -s "${dag2}.json" -t "${dag2}.png" -- \
                     dotify_graph "${dag2}.json" "${dag2}.png"

    run_idempotently -s "${dag2}.json" -s "$ofile1" -t "$ofile2" -- \
                     wire-cell -l $lfile2 -L debug -c ${dag2}.json
    local errors="$(egrep ' W | E ' $lfile2)"
    echo "$errors"
    [[ -z "$errors" ]]
    

    ilog=${dag1}-inspect.log
    olog=${dag2}-inspect.log
    run_idempotently -s "$ifile" -t "$ilog" -- \
                     wcpy img inspect -o $ilog $ifile
    run_idempotently -s "$ofile2" -t "$olog" -- \
                     wcpy img inspect -o $olog $ofile2

    local delta="$(diff -u $ilog $olog)"
    echo "$delta"
    [[ -z "$delta" ]]
}

@test "round trip json to tensor" {
    roundtrip2 json tensor
}

@test "round trip numpy to tensor" {
    roundtrip2 numpy tensor
}