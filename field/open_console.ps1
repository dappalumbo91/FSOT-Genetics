# Build and open the field console in the default browser.
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)
python scripts/build_field_console.py --open
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
