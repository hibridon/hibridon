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

function build()
{
    local hibridon_root_path="$1"
    local old_path=$PATH
    export PATH=$PATH:$hibridon_root_path/bin

    if [ "$(command_is_available tcsh)" != 'true' ]
    then
        error "failed to find tcsh, which is the language makefirst script is written in"
        return $RETURNCODE_FAILURE
    fi

    hib_makeobj
    if [ $? != 0 ]
    then
        error "hib_makeobj failed"
        return $RETURNCODE_FAILURE
    fi

    export PATH=$old_path
}


build "$HIBRIDON_ROOT_PATH"
