; The name of the installer
Name "SDMCheckInstaller"

; The file to write
OutFile "SDMCheckInstaller.exe"

; Request application privileges for Windows Vista
RequestExecutionLevel admin

; Build Unicode installer
Unicode True

; The default installation directory
InstallDir $PROGRAMFILES\SDMCheck

;--------------------------------

; Pages

Page directory
Page instfiles

;--------------------------------

; The stuff to install
Section "Install" ;No components page, name is not important

  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  File SDMCheck.exe

  SetOutPath $INSTDIR\resources
  File /nonfatal /a /r "resources\"

SectionEnd

Section "Start Menu Shortcut"
  SetOutPath $INSTDIR
  CreateDirectory "$SMPROGRAMS\SDMCheck"
  CreateShortcut "$SMPROGRAMS\SDMCheck\SDMCheck.lnk" "$INSTDIR\SDMCheck.exe" "" "$INSTDIR\resources\Icon.ico" 0
SectionEnd
