@echo off

g++ src/*.cpp src/glad.c -Iinclude -Iglfw3 -Iopengl32 -o MagneticFields_C++.exe

if %errorlevel% equ 0 (
    echo.
    echo [Succes] Compilation complete
    echo.
    MagneticFields_C++.exe
) else (
    echo.
    echo [Error] Compilation Failed
    echo Ensure your msys64/minGW compiler is properly configured in your system PATH variables
    echo.
    pause
)