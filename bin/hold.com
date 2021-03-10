 lslpp -h -qcOu xlfcmp | awk -F: '/COMMIT/ {print $3}' | head -1 | awk -F. '{printf "%d.%d\n", $1, $2}'
 lslpp -h -qcOu xlfcmp.obj | awk -F: '/COMPLETE/ {print $4}' | head -1 | awk -F. '{printf "%d.%d\n", $1, $2}'
