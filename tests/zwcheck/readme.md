samples for unit testing `zwcheck.py`

These samples were obtained with the following command:

on `debian`:
```sh
unbuffer make > >(tee make.stdout) 2> >(tee make.stderr >&2)
```

on `osx` (where the comand `unbuffer` was not available):
```sh
make > >(tee make.stdout) 2> >(tee make.stderr >&2)
```


