# asard test suite
This folder contains some tests to ensure processing works as intended.
The tests need some test data, which is not delivered with the software.
You can request the test data from the package maintainer [John Truckenbrodt](mailto:john.truckenbrodt@dlr.de).
The tests look up the test folder using an environment variable `ASARD_TESTDATA`.
This folder needs to contain all data and the variable needs to be set to it.  

The script `run_slurm.sh` can be used to run all tests in parallel on an HPC platform using the SLURM scheduler.
It needs to be modified to set the `ASARD_TESTDATA` variable and a valid e-mail address for receiving SLURM notifications.
Some more configuration may be necessary on your system, e.g. the cluster and partition.