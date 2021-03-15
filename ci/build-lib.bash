HIBRIDON_ROOT_PATH="$1"

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

function unattended_makefirst()
{
    local hibridon_root_path="$1"
    local old_path=$PATH
    export PATH=$PATH:$hibridon_root_path/bin

    if [ "$(command_is_available tcsh)" != 'true' ]
    then
        error "failed to find tcsh, which is the language makefirst script is written in"
        return $RETURNCODE_FAILURE
    fi

    if [ "$(command_is_available expect)" != 'true' ]
    then
        error "failed to find expect command, which is needed"
        return $RETURNCODE_FAILURE
    fi

    local la_lib_path='/usr/lib/x86_64-linux-gnu/'
    local blas_lib_name='blas'
    local lapack_lib_name='lapack'
    local lib_name=''

    # check that the library dev files are here
    for lib_name in $blas_lib_name $lapack_lib_name
    do
        local lib_path="${la_lib_path}/lib${lib_name}.a"
        if [ ! -f "$lib_path" ]
        then
            error "failed to find $lib_path"
            return $RETURNCODE_FAILURE
        fi
    done

	expect <<-EOF
	spawn $hibridon_root_path/bin/makefirst
	expect -ex "Have you read the license agreement (/home/graffy/work/hibridon/LICENSE) \[y/n\]?  "
	send "y\r"
	expect -ex "Does your PATH include /home/graffy/work/hibridon/bin \[y/n\]?  "
	send "y\r"
	expect -ex "Your fortran compiler is ifort -O3 -save, is this correct? \[y/n\]  "
	send "n\r"
	expect -ex "Enter the name of your fortran compiler: (g95, gfortran, pgf95, pgf90 ...  "
	send "gfortran\r"
	expect -ex "Are you using Intel's mklib for lapack and blas? \[y/n\]  "
	send "n\r"
	expect -ex "Are you using Intel x86 (686 386) Architecture with PGI Lapack? \[y/n\]  "
	send "n\r"
	expect -ex "Enter the path for your BLAS library  "
	send "$la_lib_path\r"
	expect -ex "Enter the name  s of your libraries (-lxxx, -lyyy, etc, including the hypen)  "
	send -- "-l$lapack_lib_name -l$blas_lib_name\r"
	expect eof;
	catch wait result
	exit [lindex \$result 3]
	EOF
	if [ $? != 0 ]
	then
		error "makefirst failed"
		return $RETURNCODE_FAILURE
	fi

    export PATH=$old_path
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

function test_build()
{
    local hibridon_root_path="$1"
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
        error "hib_makeobj failed"
        return $RETURNCODE_FAILURE
    fi

    export PATH=$old_path
}

function build_lib()
{
    local hibridon_root_path="$1"

    unattended_makefirst $hibridon_root_path
    if [ $? != 0 ]
    then
        error "unattended_makefirst failed"
        return $RETURNCODE_FAILURE
    fi
    build $hibridon_root_path
    if [ $? != 0 ]
    then
        error "build failed"
        return $RETURNCODE_FAILURE
    fi
    test_build $hibridon_root_path
    if [ $? != 0 ]
    then
        error "test_build failed"
        return $RETURNCODE_FAILURE
    fi

}

# unattended_makefirst $HIBRIDON_ROOT_PATH
# build $HIBRIDON_ROOT_PATH
# test_build $HIBRIDON_ROOT_PATH
build_lib $HIBRIDON_ROOT_PATH
