HIBRIDON_ROOT_PATH="$1"
TEMP_PATH="$2"
CONFIG_ID="$3"  # 'fr.univ-rennes1.ipr.physix.gfortran' or 'fr.univ-rennes1.ipr.physix.ifort'

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
    local config_id="$2" # 'fr.univ-rennes1.ipr.physix.gfortran' or 'fr.univ-rennes1.ipr.physix.ifort'
    local old_path=$PATH
    export PATH=$PATH:$hibridon_root_path/bin

    case $config_id in 
        'fr.univ-rennes1.ipr.physix.ifort')
            local ifort_module='compilers/ifort/19.1.1'
            module load "$ifort_module"
            if [ $? != 0 ]
            then
                error "failed to load module $ifort_module"
            fi
            ;;
    esac

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


build "$HIBRIDON_ROOT_PATH" "$CONFIG_ID"
