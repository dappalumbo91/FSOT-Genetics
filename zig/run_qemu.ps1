# Build freestanding genetics kernel and run under QEMU (serial log).
# Pattern: fsot-neuron-zig/run_qemu.ps1
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

$zigCmd = Get-Command zig -ErrorAction SilentlyContinue
$zig = $null
if ($zigCmd) { $zig = $zigCmd.Source }
if (-not $zig) {
    $cand = Get-ChildItem "$env:LOCALAPPDATA\Microsoft\WinGet\Packages" -Filter zig.exe -Recurse -ErrorAction SilentlyContinue |
        Select-Object -First 1 -ExpandProperty FullName
    if ($cand) { $zig = $cand }
}
if (-not $zig) { throw "zig not found on PATH" }

Write-Host "=== zig build kernel ==="
& $zig build kernel
if ($LASTEXITCODE -ne 0) { throw "zig build kernel failed" }

$kernelSrc = Join-Path $PSScriptRoot "zig-out\bin\fsot_genetics_kernel"
if (-not (Test-Path $kernelSrc)) {
    $alt = Get-ChildItem (Join-Path $PSScriptRoot "zig-out\bin") -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -like "fsot_genetics_kernel*" } | Select-Object -First 1
    if ($alt) { $kernelSrc = $alt.FullName }
}
if (-not (Test-Path $kernelSrc)) {
    Write-Host "FAIL: kernel binary not found under zig-out/bin"
    exit 2
}

$qemuCmd = Get-Command qemu-system-x86_64 -ErrorAction SilentlyContinue
$qemu = $null
if ($qemuCmd) { $qemu = $qemuCmd.Source }
if (-not $qemu -and (Test-Path "C:\Program Files\qemu\qemu-system-x86_64.exe")) {
    $qemu = "C:\Program Files\qemu\qemu-system-x86_64.exe"
}
if (-not $qemu) {
    Write-Host "WARN: qemu-system-x86_64 not found - kernel built at $kernelSrc"
    Write-Host "Install QEMU then re-run. Host gate: zig build host"
    exit 0
}

$kernel = Join-Path $env:TEMP "fsot_genetics_kernel"
$serialLog = Join-Path $env:TEMP "fsot_genetics_qemu_serial.log"
$errLog = Join-Path $env:TEMP "fsot_genetics_qemu_err.log"
Copy-Item -Force $kernelSrc $kernel
Remove-Item $serialLog, $errLog -ErrorAction SilentlyContinue

Write-Host "=== QEMU genetics product cell ==="
$argList = @(
    "-display", "none",
    "-serial", "file:$serialLog",
    "-no-reboot",
    "-m", "64M",
    "-kernel", $kernel
)
$p = Start-Process -FilePath $qemu -ArgumentList $argList -PassThru -WindowStyle Hidden -RedirectStandardError $errLog

$maxWaitSec = 60
$waited = 0
while (-not $p.HasExited -and $waited -lt $maxWaitSec) {
    Start-Sleep -Seconds 2
    $waited += 2
    if (Test-Path $serialLog) {
        $partial = Get-Content $serialLog -Raw -ErrorAction SilentlyContinue
        if ($partial -match "FSOT_STAGE_GENETICS_") {
            Start-Sleep -Seconds 1
            break
        }
    }
}
if (-not $p.HasExited) { Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue }

Write-Host "--- serial output ---"
if (Test-Path $serialLog) {
    Get-Content $serialLog
    $txt = Get-Content $serialLog -Raw
    if ($txt -match "FSOT_STAGE_GENETICS_FAIL") {
        Write-Host "=== QEMU GATE FAIL ==="
        exit 1
    }
    if ($txt -match "FSOT_STAGE_GENETICS_OK") {
        Write-Host "=== QEMU GATE PASS ==="
        exit 0
    }
    Write-Host "=== QEMU GATE INCOMPLETE (no stage marker) ==="
    exit 1
}
Write-Host "=== no serial log ==="
exit 1
