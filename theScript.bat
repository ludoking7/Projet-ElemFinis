@echo off
set /a b=%1


if %b%==1 (cd .\ProjectPreProcessor\build
    mingw32-make
    .\myFem.exe 
    cd ..\..\)

if %b%==2 (cd .\Project\build
    mingw32-make
    .\myFem.exe
    cd ..\..\)

if %b%==3 (cd .\ProjectPostProcessor\build
    mingw32-make
    .\myFem.exe
    cd ..\..\)

if %b%==0 (cd .\ProjectPreProcessor\build
    mingw32-make
    .\myFem.exe
    cd ..\..\

    cd .\Project\build
    mingw32-make
    .\myFem.exe
    cd ..\..\

    cd .\ProjectPostProcessor\build
    mingw32-make
    .\myFem.exe
    cd ..\..\)

