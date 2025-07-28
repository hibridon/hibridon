
# CHLIST

If `CHLIST` = .TRUE., the channels included are listed at every value of $\mathbf{J}_{\rm tot}$ (or $\bar{L}$).

If `CHLIST` = .FALSE., then the channels are listed in the output file only at the first value of $\mathbf{J}_{\rm tot}$.


Note that the flag `CHLIST` controls printing the channel list only in the output file.
The only way to print out the list of channels on your terminal screen is to set the flag [`BASTST`](./bastst.md) = .TRUE., in which case the job will cease after the list of channels is printed, and await a new command.

When the number of channels is large, you will certainly want to suppress printing the list of channels.