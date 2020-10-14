# Microsoft Developer Studio Generated NMAKE File, Format Version 4.20
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101
# TARGTYPE "Win32 (x86) QuickWin Application" 0x0107
# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102
# TARGTYPE "Win32 (x86) Static Library" 0x0104

!IF "$(CFG)" == ""
CFG=expresslib - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to expresslib - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "express - Win32 Release" && "$(CFG)" !=\
 "express - Win32 Debug" && "$(CFG)" != "scifit - Win32 Release" && "$(CFG)" !=\
 "scifit - Win32 Debug" && "$(CFG)" != "wexpress - Win32 Release" && "$(CFG)" !=\
 "wexpress - Win32 Debug" && "$(CFG)" != "expresslib - Win32 Release" &&\
 "$(CFG)" != "expresslib - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "Express.mak" CFG="expresslib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "express - Win32 Release" (based on\
 "Win32 (x86) QuickWin Application")
!MESSAGE "express - Win32 Debug" (based on "Win32 (x86) QuickWin Application")
!MESSAGE "scifit - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "scifit - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE "wexpress - Win32 Release" (based on\
 "Win32 (x86) Dynamic-Link Library")
!MESSAGE "wexpress - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "expresslib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "expresslib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "wexpress - Win32 Debug"

!IF  "$(CFG)" == "express - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\Express.exe"

CLEAN : 
	-@erase "$(INTDIR)\Datutil.obj"
	-@erase "$(INTDIR)\DATUTIL_MODULE.mod"
	-@erase "$(INTDIR)\Express.obj"
	-@erase "$(INTDIR)\Express1.obj"
	-@erase "$(INTDIR)\Express_Ext_Module.mod"
	-@erase "$(INTDIR)\Express_LOGICAL_Module.mod"
	-@erase "$(INTDIR)\Express_Module.mod"
	-@erase "$(INTDIR)\Expressl.obj"
	-@erase "$(INTDIR)\INTER_IO_module.mod"
	-@erase "$(INTDIR)\interIO.obj"
	-@erase "$(INTDIR)\STRING_OP_module.mod"
	-@erase "$(INTDIR)\strop.obj"
	-@erase "$(INTDIR)\test.obj"
	-@erase "$(OUTDIR)\Express.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=fl32.exe
# ADD BASE F90 /Ox /I "Release/" /c /nologo /MW
# ADD F90 /Ox /I "Release/" /c /nologo /MW
F90_PROJ=/Ox /I "Release/" /c /nologo /MW /Fo"Release/" 
F90_OBJS=.\Release/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/Express.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS=.\Release/
CPP_SBRS=.\.

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

MTL=mktyplib.exe
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/Express.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386 /nodefaultlib:"console.lib"
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386 /nodefaultlib:"console.lib"
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:windows /incremental:no\
 /pdb:"$(OUTDIR)/Express.pdb" /machine:I386 /nodefaultlib:"console.lib"\
 /out:"$(OUTDIR)/Express.exe" 
LINK32_OBJS= \
	"$(INTDIR)\Datutil.obj" \
	"$(INTDIR)\Express.obj" \
	"$(INTDIR)\Express1.obj" \
	"$(INTDIR)\Expressl.obj" \
	"$(INTDIR)\interIO.obj" \
	"$(INTDIR)\strop.obj" \
	"$(INTDIR)\test.obj"

"$(OUTDIR)\Express.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "express - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\Express.exe"

CLEAN : 
	-@erase "$(INTDIR)\Datutil.obj"
	-@erase "$(INTDIR)\datutil_module.mod"
	-@erase "$(INTDIR)\Express.obj"
	-@erase "$(INTDIR)\Express1.obj"
	-@erase "$(INTDIR)\express_ext_module.mod"
	-@erase "$(INTDIR)\express_logical_module.mod"
	-@erase "$(INTDIR)\express_module.mod"
	-@erase "$(INTDIR)\Expressl.obj"
	-@erase "$(INTDIR)\inter_io_module.mod"
	-@erase "$(INTDIR)\interIO.obj"
	-@erase "$(INTDIR)\string_op_module.mod"
	-@erase "$(INTDIR)\strop.obj"
	-@erase "$(INTDIR)\test.obj"
	-@erase "$(OUTDIR)\Express.exe"
	-@erase "$(OUTDIR)\Express.ilk"
	-@erase "$(OUTDIR)\Express.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=fl32.exe
# ADD BASE F90 /Zi /I "Debug/" /c /nologo /MW
# ADD F90 /Zi /I "Debug/" /c /nologo /MW
F90_PROJ=/Zi /I "Debug/" /c /nologo /MW /Fo"Debug/" /Fd"Debug/Express.pdb" 
F90_OBJS=.\Debug/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/Express.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=.\Debug/
CPP_SBRS=.\.

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

MTL=mktyplib.exe
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/Express.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386 /nodefaultlib:"console.lib"
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386 /nodefaultlib:"console.lib"
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:windows /incremental:yes\
 /pdb:"$(OUTDIR)/Express.pdb" /debug /machine:I386 /nodefaultlib:"console.lib"\
 /out:"$(OUTDIR)/Express.exe" 
LINK32_OBJS= \
	"$(INTDIR)\Datutil.obj" \
	"$(INTDIR)\Express.obj" \
	"$(INTDIR)\Express1.obj" \
	"$(INTDIR)\Expressl.obj" \
	"$(INTDIR)\interIO.obj" \
	"$(INTDIR)\strop.obj" \
	"$(INTDIR)\test.obj"

"$(OUTDIR)\Express.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "scifit - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "scifit\Release"
# PROP BASE Intermediate_Dir "scifit\Release"
# PROP BASE Target_Dir "scifit"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "scifit\Release"
# PROP Intermediate_Dir "scifit\Release"
# PROP Target_Dir "scifit"
OUTDIR=.\scifit\Release
INTDIR=.\scifit\Release

ALL : "wexpress - Win32 Release" "$(OUTDIR)\scifit.exe"

CLEAN : 
	-@erase "$(INTDIR)\CALLBAC.obj"
	-@erase "$(INTDIR)\CALLBACK_MOD.MOD"
	-@erase "$(INTDIR)\FIGWND.OBJ"
	-@erase "$(INTDIR)\wexpress.res"
	-@erase "$(INTDIR)\wmain.obj"
	-@erase "$(INTDIR)\WndFig.mod"
	-@erase "$(OUTDIR)\scifit.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=fl32.exe
# ADD BASE F90 /Ox /I "scifit\Release/" /c /nologo
# ADD F90 /Ox /I "scifit\Release/" /c /nologo
F90_PROJ=/Ox /I "scifit\Release/" /c /nologo /Fo"scifit\Release/" 
F90_OBJS=.\scifit\Release/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/scifit.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS=.\scifit\Release/
CPP_SBRS=.\.

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

MTL=mktyplib.exe
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/wexpress.res" /d "NDEBUG" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/scifit.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:windows /incremental:no\
 /pdb:"$(OUTDIR)/scifit.pdb" /machine:I386 /out:"$(OUTDIR)/scifit.exe" 
LINK32_OBJS= \
	"$(INTDIR)\CALLBAC.obj" \
	"$(INTDIR)\FIGWND.OBJ" \
	"$(INTDIR)\wexpress.res" \
	"$(INTDIR)\wmain.obj" \
	"..\..\..\pandau\figure\Release\figure.lib" \
	"..\..\..\pandau\math95c\Release\math95c.lib" \
	"..\..\..\pandau\others\Release\others.lib" \
	".\wexpress\Release\wexpress.lib"

"$(OUTDIR)\scifit.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "scifit\Debug"
# PROP BASE Intermediate_Dir "scifit\Debug"
# PROP BASE Target_Dir "scifit"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "scifit\Debug"
# PROP Intermediate_Dir "scifit\Debug"
# PROP Target_Dir "scifit"
OUTDIR=.\scifit\Debug
INTDIR=.\scifit\Debug

ALL : "wexpress - Win32 Debug" "$(OUTDIR)\scifit.exe"

CLEAN : 
	-@erase "$(INTDIR)\CALLBAC.obj"
	-@erase "$(INTDIR)\CALLBACK_MOD.MOD"
	-@erase "$(INTDIR)\FIGWND.OBJ"
	-@erase "$(INTDIR)\wexpress.res"
	-@erase "$(INTDIR)\wmain.obj"
	-@erase "$(INTDIR)\Wndfig.mod"
	-@erase "$(OUTDIR)\scifit.exe"
	-@erase "$(OUTDIR)\scifit.ilk"
	-@erase "$(OUTDIR)\scifit.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=fl32.exe
# ADD BASE F90 /Zi /I "scifit\Debug/" /c /nologo
# ADD F90 /Zi /I "scifit\Debug/" /c /nologo
F90_PROJ=/Zi /I "scifit\Debug/" /c /nologo /Fo"scifit\Debug/"\
 /Fd"scifit\Debug/Express.pdb" 
F90_OBJS=.\scifit\Debug/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/scifit.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=.\scifit\Debug/
CPP_SBRS=.\.

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

MTL=mktyplib.exe
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
RSC_PROJ=/l 0x409 /fo"$(INTDIR)/wexpress.res" /d "_DEBUG" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/scifit.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:windows /incremental:yes\
 /pdb:"$(OUTDIR)/scifit.pdb" /debug /machine:I386 /out:"$(OUTDIR)/scifit.exe" 
LINK32_OBJS= \
	"$(INTDIR)\CALLBAC.obj" \
	"$(INTDIR)\FIGWND.OBJ" \
	"$(INTDIR)\wexpress.res" \
	"$(INTDIR)\wmain.obj" \
	"..\..\..\pandau\figure\Release\figure.lib" \
	"..\..\..\pandau\math95c\Release\math95c.lib" \
	"..\..\..\pandau\others\Release\others.lib" \
	".\wexpress\Debug\wexpress.lib"

"$(OUTDIR)\scifit.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "wexpress - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "wexpress\Release"
# PROP BASE Intermediate_Dir "wexpress\Release"
# PROP BASE Target_Dir "wexpress"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "wexpress\Release"
# PROP Intermediate_Dir "wexpress\Release"
# PROP Target_Dir "wexpress"
OUTDIR=.\wexpress\Release
INTDIR=.\wexpress\Release

ALL : "$(OUTDIR)\wexpress.dll" "$(OUTDIR)\Express_Module.mod"

CLEAN : 
	-@erase "$(INTDIR)\Express.obj"
	-@erase "$(INTDIR)\Express_Module.mod"
	-@erase "$(OUTDIR)\wexpress.dll"
	-@erase "$(OUTDIR)\wexpress.exp"
	-@erase "$(OUTDIR)\wexpress.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=fl32.exe
# ADD BASE F90 /Ox /I "wexpress\Release/" /c /nologo /MT
# ADD F90 /Ox /I "wexpress\Release/" /c /nologo /MT
F90_PROJ=/Ox /I "wexpress\Release/" /c /nologo /MT /Fo"wexpress\Release/" 
F90_OBJS=.\wexpress\Release/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/wexpress.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS=.\wexpress\Release/
CPP_SBRS=.\.

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

MTL=mktyplib.exe
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/wexpress.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /machine:I386
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:windows /dll /incremental:no\
 /pdb:"$(OUTDIR)/wexpress.pdb" /machine:I386 /out:"$(OUTDIR)/wexpress.dll"\
 /implib:"$(OUTDIR)/wexpress.lib" 
LINK32_OBJS= \
	"$(INTDIR)\Express.obj" \
	"..\..\..\pandau\others\Release\others.lib"

"$(OUTDIR)\wexpress.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "wexpress - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "wexpress\Debug"
# PROP BASE Intermediate_Dir "wexpress\Debug"
# PROP BASE Target_Dir "wexpress"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "wexpress\Debug"
# PROP Intermediate_Dir "wexpress\Debug"
# PROP Target_Dir "wexpress"
OUTDIR=.\wexpress\Debug
INTDIR=.\wexpress\Debug

ALL : "$(OUTDIR)\wexpress.dll" "$(OUTDIR)\Express_Module.mod"

CLEAN : 
	-@erase "$(INTDIR)\Express.obj"
	-@erase "$(INTDIR)\Express_Module.mod"
	-@erase "$(OUTDIR)\wexpress.dll"
	-@erase "$(OUTDIR)\wexpress.exp"
	-@erase "$(OUTDIR)\wexpress.ilk"
	-@erase "$(OUTDIR)\wexpress.lib"
	-@erase "$(OUTDIR)\wexpress.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=fl32.exe
# ADD BASE F90 /Zi /I "wexpress\Debug/" /c /nologo /MT
# ADD F90 /Zi /I "wexpress\Debug/" /c /nologo /MT
F90_PROJ=/Zi /I "wexpress\Debug/" /c /nologo /MT /Fo"wexpress\Debug/"\
 /Fd"wexpress\Debug/Express.pdb" 
F90_OBJS=.\wexpress\Debug/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /MTd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /MTd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/wexpress.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 
CPP_OBJS=.\wexpress\Debug/
CPP_SBRS=.\.

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

MTL=mktyplib.exe
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/wexpress.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /debug /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /debug /machine:I386
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\
 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\
 odbccp32.lib /nologo /subsystem:windows /dll /incremental:yes\
 /pdb:"$(OUTDIR)/wexpress.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)/wexpress.dll" /implib:"$(OUTDIR)/wexpress.lib" 
LINK32_OBJS= \
	"$(INTDIR)\Express.obj" \
	"..\..\..\pandau\others\Release\others.lib"

"$(OUTDIR)\wexpress.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "expresslib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "expresslib\Release"
# PROP BASE Intermediate_Dir "expresslib\Release"
# PROP BASE Target_Dir "expresslib"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "expresslib\Release"
# PROP Intermediate_Dir "expresslib\Release"
# PROP Target_Dir "expresslib"
OUTDIR=.\expresslib\Release
INTDIR=.\expresslib\Release

ALL : "$(OUTDIR)\expresslib.lib" "$(OUTDIR)\Express_Ext_Module.mod"

CLEAN : 
	-@erase "$(INTDIR)\Express.obj"
	-@erase "$(INTDIR)\Express1.obj"
	-@erase "$(INTDIR)\Express_Ext_Module.mod"
	-@erase "$(INTDIR)\Express_LOGICAL_Module.mod"
	-@erase "$(INTDIR)\Express_Module.mod"
	-@erase "$(INTDIR)\Expressl.obj"
	-@erase "$(OUTDIR)\expresslib.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=fl32.exe
# ADD BASE F90 /Ox /I "expresslib\Release/" /c /nologo
# ADD F90 /Ox /I "expresslib\Release/" /c /nologo
F90_PROJ=/Ox /I "expresslib\Release/" /c /nologo /Fo"expresslib\Release/" 
F90_OBJS=.\expresslib\Release/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/expresslib.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS=.\expresslib\Release/
CPP_SBRS=.\.

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/expresslib.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/expresslib.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Express.obj" \
	"$(INTDIR)\Express1.obj" \
	"$(INTDIR)\Expressl.obj"

"$(OUTDIR)\expresslib.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "expresslib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "expresslib\Debug"
# PROP BASE Intermediate_Dir "expresslib\Debug"
# PROP BASE Target_Dir "expresslib"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "c:\hou qing\lib"
# PROP Intermediate_Dir "c:\hou qing\include"
# PROP Target_Dir "expresslib"
OUTDIR=c:\hou qing\lib
INTDIR=c:\hou qing\include

ALL : "$(OUTDIR)\expresslib.lib" "..\..\..\include\Express_Ext_Module.mod"

CLEAN : 
	-@erase "$(INTDIR)\Express.obj"
	-@erase "$(INTDIR)\Express1.obj"
	-@erase "$(INTDIR)\Express_Ext_Module.mod"
	-@erase "$(INTDIR)\express_logical_module.mod"
	-@erase "$(INTDIR)\express_module.mod"
	-@erase "$(INTDIR)\Expressl.obj"
	-@erase "$(OUTDIR)\expresslib.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

F90=fl32.exe
# ADD BASE F90 /Z7 /I "expresslib\Debug/" /c /nologo
# ADD F90 /Z7 /I "c:\hou qing\include/" /c /nologo
F90_PROJ=/Z7 /I "c:\hou qing\include/" /c /nologo /Fo"c:\hou qing\include/" 
F90_OBJS="c:\hou qing\include/"

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /MLd /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/expresslib.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS="c:\hou qing\include/"
CPP_SBRS=.\.

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/expresslib.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/expresslib.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Express.obj" \
	"$(INTDIR)\Express1.obj" \
	"$(INTDIR)\Expressl.obj"

"$(OUTDIR)\expresslib.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

################################################################################
# Begin Target

# Name "express - Win32 Release"
# Name "express - Win32 Debug"

!IF  "$(CFG)" == "express - Win32 Release"

!ELSEIF  "$(CFG)" == "express - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\Express.f90

!IF  "$(CFG)" == "express - Win32 Release"

NODEP_F90_EXPRE=\
	".\Release\STRING_OP_module.mod"\
	
F90_MODOUT=\
	"Express_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\STRING_OP_module.mod"
   $(BuildCmds)

"$(INTDIR)\Express_Module.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\STRING_OP_module.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "express - Win32 Debug"

DEP_F90_EXPRE=\
	".\Debug/STRING_OP_module.mod"\
	
F90_MODOUT=\
	"Express_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express.obj" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"\
 "$(INTDIR)\string_op_module.mod"
   $(BuildCmds)

"$(INTDIR)\express_module.mod" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"\
 "$(INTDIR)\string_op_module.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\test.f90

!IF  "$(CFG)" == "express - Win32 Release"

NODEP_F90_TEST_=\
	".\Release\Express_Ext_Module.mod"\
	".\Release\Express_LOGICAL_Module.mod"\
	".\Release\Express_Module.mod"\
	".\Release\INTER_IO_module.mod"\
	".\Release\STRING_OP_module.mod"\
	

"$(INTDIR)\test.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\Express_Module.mod"\
 "$(INTDIR)\STRING_OP_module.mod" "$(INTDIR)\Express_Ext_Module.mod"\
 "$(INTDIR)\INTER_IO_module.mod" "$(INTDIR)\Express_LOGICAL_Module.mod"


!ELSEIF  "$(CFG)" == "express - Win32 Debug"

DEP_F90_TEST_=\
	".\Debug/Express_LOGICAL_Module.mod"\
	".\Debug/Express_Module.mod"\
	".\Debug/STRING_OP_module.mod"\
	".\Debug\express_ext_module.mod"\
	".\Debug\inter_io_module.mod"\
	

"$(INTDIR)\test.obj" : $(SOURCE) $(DEP_F90_TEST_) "$(INTDIR)"\
 "$(INTDIR)\express_module.mod" "$(INTDIR)\string_op_module.mod"\
 "$(INTDIR)\express_ext_module.mod" "$(INTDIR)\inter_io_module.mod"\
 "$(INTDIR)\express_logical_module.mod"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE="\hou qing\Houqing\Mathlib\strop\strop.f90"

!IF  "$(CFG)" == "express - Win32 Release"

F90_MODOUT=\
	"STRING_OP_module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\strop.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\STRING_OP_module.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "express - Win32 Debug"

F90_MODOUT=\
	"STRING_OP_module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\strop.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\string_op_module.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Express1.f90

!IF  "$(CFG)" == "express - Win32 Release"

NODEP_F90_EXPRES=\
	".\Release\DATUTIL_MODULE.mod"\
	".\Release\Express_LOGICAL_Module.mod"\
	".\Release\Express_Module.mod"\
	".\Release\STRING_OP_module.mod"\
	
F90_MODOUT=\
	"Express_Ext_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express1.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\Express_Module.mod"\
 "$(INTDIR)\STRING_OP_module.mod" "$(INTDIR)\DATUTIL_MODULE.mod"\
 "$(INTDIR)\Express_LOGICAL_Module.mod"
   $(BuildCmds)

"$(INTDIR)\Express_Ext_Module.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\Express_Module.mod" "$(INTDIR)\STRING_OP_module.mod"\
 "$(INTDIR)\DATUTIL_MODULE.mod" "$(INTDIR)\Express_LOGICAL_Module.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "express - Win32 Debug"

DEP_F90_EXPRES=\
	".\Debug/DATUTIL_MODULE.mod"\
	".\Debug/Express_LOGICAL_Module.mod"\
	".\Debug/Express_Module.mod"\
	".\Debug/STRING_OP_module.mod"\
	
F90_MODOUT=\
	"Express_Ext_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express1.obj" : $(SOURCE) $(DEP_F90_EXPRES) "$(INTDIR)"\
 "$(INTDIR)\express_module.mod" "$(INTDIR)\string_op_module.mod"\
 "$(INTDIR)\datutil_module.mod" "$(INTDIR)\express_logical_module.mod"
   $(BuildCmds)

"$(INTDIR)\express_ext_module.mod" : $(SOURCE) $(DEP_F90_EXPRES) "$(INTDIR)"\
 "$(INTDIR)\express_module.mod" "$(INTDIR)\string_op_module.mod"\
 "$(INTDIR)\datutil_module.mod" "$(INTDIR)\express_logical_module.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE="\hou qing\Houqing\Mathlib\Datautil\Datutil.f90"

!IF  "$(CFG)" == "express - Win32 Release"

F90_MODOUT=\
	"DATUTIL_MODULE"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Datutil.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\DATUTIL_MODULE.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "express - Win32 Debug"

F90_MODOUT=\
	"DATUTIL_MODULE"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Datutil.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\datutil_module.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE="\hou qing\Houqing\Mathlib\interface\interIO.f90"

!IF  "$(CFG)" == "express - Win32 Release"

NODEP_F90_INTER=\
	".\Release\STRING_OP_module.mod"\
	
F90_MODOUT=\
	"INTER_IO_module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\interIO.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\STRING_OP_module.mod"
   $(BuildCmds)

"$(INTDIR)\INTER_IO_module.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\STRING_OP_module.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "express - Win32 Debug"

DEP_F90_INTER=\
	".\Debug/STRING_OP_module.mod"\
	
F90_MODOUT=\
	"INTER_IO_module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\interIO.obj" : $(SOURCE) $(DEP_F90_INTER) "$(INTDIR)"\
 "$(INTDIR)\string_op_module.mod"
   $(BuildCmds)

"$(INTDIR)\inter_io_module.mod" : $(SOURCE) $(DEP_F90_INTER) "$(INTDIR)"\
 "$(INTDIR)\string_op_module.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Expressl.f90

!IF  "$(CFG)" == "express - Win32 Release"

NODEP_F90_EXPRESS=\
	".\Release\Express_Module.mod"\
	".\Release\STRING_OP_module.mod"\
	
F90_MODOUT=\
	"Express_LOGICAL_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Expressl.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\Express_Module.mod"\
 "$(INTDIR)\STRING_OP_module.mod"
   $(BuildCmds)

"$(INTDIR)\Express_LOGICAL_Module.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\Express_Module.mod" "$(INTDIR)\STRING_OP_module.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "express - Win32 Debug"

DEP_F90_EXPRESS=\
	".\Debug/Express_Module.mod"\
	".\Debug/STRING_OP_module.mod"\
	
F90_MODOUT=\
	"Express_LOGICAL_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Expressl.obj" : $(SOURCE) $(DEP_F90_EXPRESS) "$(INTDIR)"\
 "$(INTDIR)\express_module.mod" "$(INTDIR)\string_op_module.mod"
   $(BuildCmds)

"$(INTDIR)\express_logical_module.mod" : $(SOURCE) $(DEP_F90_EXPRESS)\
 "$(INTDIR)" "$(INTDIR)\express_module.mod" "$(INTDIR)\string_op_module.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
# End Target
################################################################################
# Begin Target

# Name "scifit - Win32 Release"
# Name "scifit - Win32 Debug"

!IF  "$(CFG)" == "scifit - Win32 Release"

!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE="\Hou Qing\pandau\figure\Release\figure.lib"

!IF  "$(CFG)" == "scifit - Win32 Release"

!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE="\Hou Qing\pandau\math95c\Release\math95c.lib"

!IF  "$(CFG)" == "scifit - Win32 Release"

!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE="\Hou Qing\pandau\others\Release\others.lib"

!IF  "$(CFG)" == "scifit - Win32 Release"

!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Project Dependency

# Project_Dep_Name "wexpress"

!IF  "$(CFG)" == "scifit - Win32 Release"

"wexpress - Win32 Release" : 
   $(MAKE) /$(MAKEFLAGS) /F ".\Express.mak" CFG="wexpress - Win32 Release" 

!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

"wexpress - Win32 Debug" : 
   $(MAKE) /$(MAKEFLAGS) /F ".\Express.mak" CFG="wexpress - Win32 Debug" 

!ENDIF 

# End Project Dependency
################################################################################
# Begin Source File

SOURCE=.\scifit\CALLBAC.f90

!IF  "$(CFG)" == "scifit - Win32 Release"

NODEP_F90_CALLB=\
	".\scifit\Release\WndFig.mod"\
	
F90_MODOUT=\
	"CALLBACK_MOD"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\CALLBAC.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\WndFig.mod"
   $(BuildCmds)

"$(INTDIR)\CALLBACK_MOD.MOD" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\WndFig.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

DEP_F90_CALLB=\
	".\Scifit\Debug\Wndfig.mod"\
	
F90_MODOUT=\
	"CALLBACK_MOD"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\CALLBAC.obj" : $(SOURCE) $(DEP_F90_CALLB) "$(INTDIR)"\
 "$(INTDIR)\Wndfig.mod"
   $(BuildCmds)

"$(INTDIR)\CALLBACK_MOD.MOD" : $(SOURCE) $(DEP_F90_CALLB) "$(INTDIR)"\
 "$(INTDIR)\Wndfig.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\scifit\FIGWND.F90
DEP_F90_FIGWN=\
	{$(INCLUDE)}"\fig_modules.mod"\
	{$(INCLUDE)}"\wghdef.mod"\
	

!IF  "$(CFG)" == "scifit - Win32 Release"

F90_MODOUT=\
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\FIGWND.OBJ" : $(SOURCE) $(DEP_F90_FIGWN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\WndFig.mod" : $(SOURCE) $(DEP_F90_FIGWN) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

F90_MODOUT=\
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig" \
	"WndFig"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\FIGWND.OBJ" : $(SOURCE) $(DEP_F90_FIGWN) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Wndfig.mod" : $(SOURCE) $(DEP_F90_FIGWN) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\scifit\wmain.f90

!IF  "$(CFG)" == "scifit - Win32 Release"

DEP_F90_WMAIN=\
	{$(INCLUDE)}"\wghdef.mod"\
	
NODEP_F90_WMAIN=\
	".\scifit\Release\CALLBACK_MOD.MOD"\
	".\scifit\Release\WndFig.mod"\
	

"$(INTDIR)\wmain.obj" : $(SOURCE) $(DEP_F90_WMAIN) "$(INTDIR)"\
 "$(INTDIR)\CALLBACK_MOD.MOD" "$(INTDIR)\WndFig.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"

DEP_F90_WMAIN=\
	".\scifit\Debug\CALLBACK_MOD.MOD"\
	".\Scifit\Debug\Wndfig.mod"\
	{$(INCLUDE)}"\wghdef.mod"\
	

"$(INTDIR)\wmain.obj" : $(SOURCE) $(DEP_F90_WMAIN) "$(INTDIR)"\
 "$(INTDIR)\CALLBACK_MOD.MOD" "$(INTDIR)\Wndfig.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\scifit\wexpress.rc
DEP_RSC_WEXPR=\
	".\Scifit\Bitmap1.bmp"\
	".\Scifit\Bmp00008.bmp"\
	".\Scifit\Bmp00009.bmp"\
	".\Scifit\Ok3.bmp"\
	

!IF  "$(CFG)" == "scifit - Win32 Release"


"$(INTDIR)\wexpress.res" : $(SOURCE) $(DEP_RSC_WEXPR) "$(INTDIR)"
   $(RSC) /l 0x409 /fo"$(INTDIR)/wexpress.res" /i "scifit" /d "NDEBUG"\
 $(SOURCE)


!ELSEIF  "$(CFG)" == "scifit - Win32 Debug"


"$(INTDIR)\wexpress.res" : $(SOURCE) $(DEP_RSC_WEXPR) "$(INTDIR)"
   $(RSC) /l 0x409 /fo"$(INTDIR)/wexpress.res" /i "scifit" /d "_DEBUG"\
 $(SOURCE)


!ENDIF 

# End Source File
# End Target
################################################################################
# Begin Target

# Name "wexpress - Win32 Release"
# Name "wexpress - Win32 Debug"

!IF  "$(CFG)" == "wexpress - Win32 Release"

!ELSEIF  "$(CFG)" == "wexpress - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\Express.f90
DEP_F90_EXPRE=\
	"\hou qing\include\STRING_OP_module.mod"\
	

!IF  "$(CFG)" == "wexpress - Win32 Release"

F90_MODOUT=\
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express.obj" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Express_Module.mod" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "wexpress - Win32 Debug"

F90_MODOUT=\
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module" \
	"Express_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express.obj" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Express_Module.mod" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE="\Hou Qing\pandau\others\Release\others.lib"

!IF  "$(CFG)" == "wexpress - Win32 Release"

!ELSEIF  "$(CFG)" == "wexpress - Win32 Debug"

!ENDIF 

# End Source File
# End Target
################################################################################
# Begin Target

# Name "expresslib - Win32 Release"
# Name "expresslib - Win32 Debug"

!IF  "$(CFG)" == "expresslib - Win32 Release"

!ELSEIF  "$(CFG)" == "expresslib - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\Expressl.f90
DEP_F90_EXPRESS=\
	"\hou qing\include\Express_Module.mod"\
	"\hou qing\include\STRING_OP_module.mod"\
	

!IF  "$(CFG)" == "expresslib - Win32 Release"

F90_MODOUT=\
	"Express_LOGICAL_Module" \
	"Express_LOGICAL_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Expressl.obj" : $(SOURCE) $(DEP_F90_EXPRESS) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Express_LOGICAL_Module.mod" : $(SOURCE) $(DEP_F90_EXPRESS)\
 "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "expresslib - Win32 Debug"

F90_MODOUT=\
	"Express_LOGICAL_Module" \
	"Express_LOGICAL_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Expressl.obj" : $(SOURCE) $(DEP_F90_EXPRESS) "$(INTDIR)"\
 "$(INTDIR)\express_module.mod"
   $(BuildCmds)

"$(INTDIR)\express_logical_module.mod" : $(SOURCE) $(DEP_F90_EXPRESS)\
 "$(INTDIR)" "$(INTDIR)\express_module.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Express1.f90

!IF  "$(CFG)" == "expresslib - Win32 Release"

DEP_F90_EXPRES=\
	"\hou qing\include\DATUTIL_MODULE.mod"\
	"\hou qing\include\STRING_OP_module.mod"\
	
NODEP_F90_EXPRES=\
	".\expresslib\Release\Express_LOGICAL_Module.mod"\
	".\expresslib\Release\Express_Module.mod"\
	
F90_MODOUT=\
	"Express_Ext_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express1.obj" : $(SOURCE) $(DEP_F90_EXPRES) "$(INTDIR)"\
 "$(INTDIR)\Express_Module.mod" "$(INTDIR)\Express_LOGICAL_Module.mod"
   $(BuildCmds)

"$(INTDIR)\Express_Ext_Module.mod" : $(SOURCE) $(DEP_F90_EXPRES) "$(INTDIR)"\
 "$(INTDIR)\Express_Module.mod" "$(INTDIR)\Express_LOGICAL_Module.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "expresslib - Win32 Debug"

DEP_F90_EXPRES=\
	"\hou qing\include\DATUTIL_MODULE.mod"\
	"\hou qing\include\Express_LOGICAL_Module.mod"\
	"\hou qing\include\Express_Module.mod"\
	"\hou qing\include\STRING_OP_module.mod"\
	
F90_MODOUT=\
	"Express_Ext_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express1.obj" : $(SOURCE) $(DEP_F90_EXPRES) "$(INTDIR)"\
 "$(INTDIR)\express_module.mod" "$(INTDIR)\express_logical_module.mod"
   $(BuildCmds)

"$(INTDIR)\Express_Ext_Module.mod" : $(SOURCE) $(DEP_F90_EXPRES) "$(INTDIR)"\
 "$(INTDIR)\express_module.mod" "$(INTDIR)\express_logical_module.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Express.f90
DEP_F90_EXPRE=\
	"\hou qing\include\STRING_OP_module.mod"\
	

!IF  "$(CFG)" == "expresslib - Win32 Release"

F90_MODOUT=\
	"Express_Module" \
	"Express_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express.obj" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\Express_Module.mod" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "expresslib - Win32 Debug"

F90_MODOUT=\
	"Express_Module" \
	"Express_Module"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\Express.obj" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\express_module.mod" : $(SOURCE) $(DEP_F90_EXPRE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
