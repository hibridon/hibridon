# BATCH

To run the Hibridon package in batch mode, set up a run stream with BATCH as the first line.
This supresses the input line query Hibridon.

All other commands will be read from the run stream and executed. An example would be:

```
BATCH
PRINTS,job1
INP={input}
OUT={outfile}
JOB={jobname}
RSTART=3.5
RUN
EXIT
```

Note that on many systems it may be necessary to specify a pathname as well as a file name for  `{input}`, `{outfile}`, `{jobname}`  unless the corresponding files are located in your current directory.

Alternatively, running scattering calculations in the background can be accomplished by first defining a command file which contains the input commands. Following the above example, we would define a command file `runner.com` , say, which would contain the following lines:

```
INP={input}
PRINTS,job1
OUT={outfile}
JOB={jobname}
RSTART=3.5
RUN
EXIT
```

You would execute the command
```bash
hib_xx -k <kmax> < runner.com > runner.log &
# or
hib_xx -k <kmax> -c runner.com > runner.log &
```

Where `hib_xx` is the name of your executable code, and the file  `runner.log`  will contain a log of the execution.

⚠️ Note that using the `-c runner.com` option automatically sets Hibridon in BATCH mode. 

