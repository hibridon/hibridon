# POT

In cases in which you need to input data from a file to calculate the [potential](../potlist.md), you can do so by first writing a subroutine `LOAPOT(iunit,filname)`, in which a file named `filname` is opened on logical unit `iunit` and subsequently read. This subroutine must be then linked with the hibridon code.

When you subsequently execute the hibridon code, you need only enter the command

```
POT=<filname>
```

to cause the program to transfer control to `LOAPOT` so that the file <filname> can be read.  Control is then returned to the Hibridon driver.
