del lt.exe.%1.dat.par 
copy ..\..\pukpr.github.io\examples\MLR\%1\ts.dat %1.dat
copy .\lt.exe.par.warne_backwards1920pt5 lt.exe.par
python3 .\fill_month_gaps.py %1.dat
python3 .\sub24.py %1.dat
set CLIMATE_INDEX=%1.dat
set METRIC=DTW
.\lt
copy .\lt.exe.par .\lt.exe.par.%1


