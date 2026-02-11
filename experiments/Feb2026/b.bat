set DATA=111
mkdir %DATA%
copy lt.exe.par.%DATA%_cc   %DATA%\lt.exe.par
copy %DATA%.dat             %DATA%
copy lt.exe                 %DATA%
copy lt.exe.resp            %DATA%
