To execute a run case:

1 - choose the batch file you want to use:
    EX: "storage5.bat" will generate all folders needed to run that dataset
    -Will fire up all processes: Genie host, slaves, pest++, p-sensan

2 - sit back and wait

3a- Kill all the leftover genie instances with KillRun.bat
	-also kills other processes associated with the model run so you dont have to manually do it

3b - Clean up folders and files using Clean.bat



NOTE:  sup_storage5 is identical to storage5 data except there is a super 
parameter specified in the storage5.pst file in the sup_storage5 template