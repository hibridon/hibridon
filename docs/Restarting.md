If you have a long-running calculation which is interrupted by a power failure or other system crash, you may be able to restart by
- Setting the [flag](Flags) RSFLAG = .TRUE.
- and then re-running the job with the original input data, unchanged except for RSFLAG

Similarly, if you finish a calculation, and wish to extend it to include addition partial waves, you can
- Set the [flag](Flags) RSFLAG = .TRUE.
- Increase the value of [JTOT2](JTOT)
- then re-run the job with the original input data, otherwise unchanged

Restarting will work only if either (or both) of the flags [WRPART](WRPART) or [WRXSEC](PRXSEC-and-WRXSEC) were .TRUE. in the input to the first job.