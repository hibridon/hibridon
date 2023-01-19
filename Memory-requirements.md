The major data memory requirements associated with the Hibridon code are:

-  kmax x kmax scratch matrices used in the propagation, where is the maximum number of channels set when the code is executed, with the -k or --kmax option.

- 3 kmax x kmax matrices, which hold the R-independent potential matrix elements and associated indices

- 25 vectors of length kmax

To reduce to 5 the total number of scratch matrices required, you can [link](link) the code with the -b option.