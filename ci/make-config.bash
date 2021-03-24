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

function unattended_makefirst()
{
    local hibridon_root_path="$1"
    local config_id="$2" # 'fr.univ-rennes1.ipr.physix.gfortran' or 'fr.univ-rennes1.ipr.physix.ifort'
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

    local la_lib_path=''
    local blas_lib_name=''
    local lapack_lib_name=''
    local fcomp_exe_name=''

    case "$config_id" in
        'fr.univ-rennes1.ipr.physix.gfortran')
            fcomp_exe_name='gfortran'
            local check_la_lib_location='true'
            la_lib_path='/usr/lib'
            local os_name=$(lsb_release --short --id)
            case $os_name in
                'Debian'|'Ubuntu')
                    # On debian or ubuntu, the default libblas.a (which could point on blas, atlas or openblas, depending on update-alternatives) is sometimes in '/usr/lib', sometimes in '/usr/lib/x86_64-linux-gnu/' but the linker knows its location, as explained in https://askubuntu.com/questions/52617/what-is-usr-lib-i386-linux-gnu-for:
                    # > This change was made to enable installing versions of the same library compiled for different architectures (e.g. on an AMD64 system, one version might go in /usr/lib/x86_64-linux-gnu while the other goes in i386-linux-gnu).
                    # > Both the standard linker and dynamic linker know about these directories, so the change should be invisible for most applications. If the application is searching for actual library files manually, then it will need modification.
                    # > Details of the changes to Ubuntu can be found here:
                    # > https://wiki.ubuntu.com/MultiarchSpec
                    check_la_lib_location='false'
                    ;;
                *)
            esac
            blas_lib_name='blas'
            lapack_lib_name='lapack'

            if [ "$check_la_lib_location" = 'true' ]
            then
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
            fi

			expect <<-EOF
			spawn $hibridon_root_path/bin/makefirst
			expect -ex "Have you read the license agreement ($(hibriddir)/LICENSE) \[y/n\]?  "
			send "y\r"
			expect -ex "Does your PATH include $(hibriddir)/bin \[y/n\]?  "
			send "y\r"
			expect -ex "Your fortran compiler is ifort -O3 -save, is this correct? \[y/n\]  "
			send "n\r"
			expect -ex "Enter the name of your fortran compiler: (g95, gfortran, pgf95, pgf90 ...  "
			send "$fcomp_exe_name\r"
			expect -ex "Are you using Intel's mklib for lapack and blas? \[y/n\]  "
			send "n\r"
			expect -ex "Are you using Intel x86 (686 386) Architecture with PGI Lapack? \[y/n\]  "
			send "n\r"
			expect -ex "Enter the path for your BLAS library  "
			send "$la_lib_path/\r"
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

            ;;
        'fr.univ-rennes1.ipr.physix.ifort')
            local ifort_module='compilers/ifort/19.1.1'
            module load "$ifort_module"
            if [ $? != 0 ]
            then
                error "failed to load module $ifort_module"
            fi
            fcomp_exe_name='ifort'
            la_lib_path='/usr/lib'
            blas_lib_name='blas'
            lapack_lib_name='lapack'

			expect <<-EOF
			spawn $hibridon_root_path/bin/makefirst
			expect -ex "Have you read the license agreement ($(hibriddir)/LICENSE) \[y/n\]?  "
			send "y\r"
			expect -ex "Does your PATH include $(hibriddir)/bin \[y/n\]?  "
			send "y\r"
			expect -ex "Your fortran compiler is ifort -O3 -save, is this correct? \[y/n\]  "
			send "y\r"
			expect -ex "Are you using Intel's mklib for lapack and blas? \[y/n\]  "
			send "y\r"
			expect -ex "Are you using Intel IA64 Architecture? \[y/n\]  "
			send "n\r"
			expect eof;
			catch wait result
			exit [lindex \$result 3]
			EOF
			if [ $? != 0 ]
			then
						error "makefirst failed"
						return $RETURNCODE_FAILURE
			fi

            ;;
        *)
            error "unhandled config id : $config_id"
            return $RETURNCODE_FAILURE
    esac


    export PATH=$old_path
}

unattended_makefirst "$HIBRIDON_ROOT_PATH" "$CONFIG_ID"
