HIBRIDON_ROOT_PATH="$1"
TEMP_PATH="$2"

RETURNCODE_SUCCESS=0
RETURNCODE_FAILURE=1

function error()
{
    local error_message="$1"
    echo "ERROR: $error_message"
}

function command_is_available()
{
    local command="$1"  # eg 'csh'
    which $command > /dev/null
    if [ "$?" = "$RETURNCODE_SUCCESS" ]
    then
        echo 'true'
    else
        echo 'false'
    fi
}


function show_compare_results()
{
    local hibridon_root_path="$1"
    local temp_path="$2"

    mkdir -p "$temp_path"
    local comp_exe_path="$temp_path/comp_tests.x"
    gfortran $hibridon_root_path/ci/comp_tests.f90 -o "$comp_exe_path"
    echo ''
    echo '************************'
    echo 'Comparing outputs from tests...'
    echo '************************'

    files=("Csbtest1.ics" "Cctest1.ics" "Ccbrstest1.ics" "Ccbtest1.ics" "Ccrstest1.ics" "Cstest1.ics") 
    #files=`cd "$hibridon_root_path"/testnew/ && ls *.ics`
    # Pour l'instant test seulement sur les fichiers ics qui ne posent pas de problème
    # de dépassement du nombre de colones.
    for file in "${files[@]}"
    do
        echo "$file"...
        $comp_exe_path $hibridon_root_path/tests/$file $hibridon_root_path/testnew/$file
        if [ $? != 0 ]
        then
            error "hibtest failed"
            return $RETURNCODE_FAILURE
        fi
    done
    echo '************************'
}

function test_build()
{
    local hibridon_root_path="$1"
    local temp_path="$2"

    local old_path=$PATH
    export PATH=$PATH:$hibridon_root_path/bin

    if [ "$(command_is_available tcsh)" != 'true' ]
    then
        error "failed to find tcsh, which is the language makefirst script is written in"
        return $RETURNCODE_FAILURE
    fi

    hibtest
    if [ $? != 0 ]
    then
        error "hibtest failed"
        return $RETURNCODE_FAILURE
    fi
    export PATH=$old_path

    mkdir -p "$temp_path"
    local compare_result_file_path="$temp_path/compare_results.txt"
    show_compare_results "$hibridon_root_path" "$temp_path" | tee "$compare_result_file_path"
    if [ $? != 0 ]
    then
        error "show_compare_results failed"
        return $RETURNCODE_FAILURE
    fi

    grep -q 'FAILED' $compare_result_file_path
    if [ $? = 0 ]
    then
        error "some tests gave results that are too different from the expected results (see $compare_result_file_path for details)"
        return $RETURNCODE_FAILURE
    fi
    return $RETURNCODE_FAILURE
}

test_build "$HIBRIDON_ROOT_PATH" "$TEMP_PATH"
