cd ..
$env:Path += ";C:\GNAT\2021\bin"; gprbuild lte.gpr
cd run
copy lt.exe lt.exe.bak
copy ..\obj\enso_opt.exe lt.exe
