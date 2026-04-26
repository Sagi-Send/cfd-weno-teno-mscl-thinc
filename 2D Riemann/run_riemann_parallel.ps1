param(
    [ValidateSet("preview", "production")]
    [string]$Mode = "preview",
    [int]$Threads = 4,
    [double]$TFinal = 0.5,
    [switch]$SkipAnimation
)

$ErrorActionPreference = "Stop"

$caseDir = $PSScriptRoot
$buildDir = Join-Path $caseDir "build"
$exe = Join-Path $buildDir "ri-six-omp.exe"
$relativeExe = ".\build\ri-six-omp.exe"

if ($Mode -eq "preview") {
    $gridPhysical = 100
    $snapshots = 10
    $baseFramesName = "frames_preview"
    $baseOutputName = "riemann_density_3schemes_grid_106_preview"
} else {
    $gridPhysical = 800
    $snapshots = 100
    $baseFramesName = "frames"
    $baseOutputName = "riemann_density_3schemes_grid_806"
}

$timeTag = ("t{0:0.###}" -f $TFinal).Replace(".", "p")
if ([Math]::Abs($TFinal - 0.5) -lt 1.0e-12) {
    $framesName = $baseFramesName
    $outputName = "$baseOutputName.mp4"
} else {
    $framesName = "${baseFramesName}_${timeTag}"
    $outputName = "${baseOutputName}_${timeTag}.mp4"
}

$framesDir = Join-Path $caseDir $framesName
$outputPath = Join-Path $caseDir $outputName

New-Item -ItemType Directory -Force $buildDir | Out-Null
New-Item -ItemType Directory -Force $framesDir | Out-Null

Write-Host "Compiling OpenMP solver..."
Push-Location $caseDir
try {
    & gfortran -O3 -fopenmp -J "build" -o $relativeExe "ri-six.f90"
    if ($LASTEXITCODE -ne 0) {
        throw "gfortran failed with exit code $LASTEXITCODE"
    }
} finally {
    Pop-Location
}

$env:OMP_NUM_THREADS = "$Threads"

$schemes = @("TENO", "WENO", "MUSC")
$processes = @()

Write-Host "Launching $($schemes.Count) scheme processes with $Threads OpenMP threads each..."
foreach ($scheme in $schemes) {
    $stdoutLog = Join-Path $framesDir "$scheme.stdout.log"
    $stderrLog = Join-Path $framesDir "$scheme.stderr.log"
    $args = @($scheme, "$gridPhysical", "$snapshots", $framesName, "$TFinal")
    $processes += Start-Process -FilePath $exe `
        -ArgumentList $args `
        -WorkingDirectory $caseDir `
        -RedirectStandardOutput $stdoutLog `
        -RedirectStandardError $stderrLog `
        -WindowStyle Hidden `
        -PassThru
}

foreach ($process in $processes) {
    $process.WaitForExit()
    $process.Refresh()
}

$failed = $false
for ($idx = 0; $idx -lt $processes.Count; $idx++) {
    $exitCode = $processes[$idx].ExitCode
    if ($null -ne $exitCode -and $exitCode -ne 0) {
        Write-Host "$($schemes[$idx]) failed with exit code $($processes[$idx].ExitCode)"
        $failed = $true
    }
}
if ($failed) {
    throw "At least one scheme process failed. Check logs in $framesDir"
}

$expectedBytes = [int64]$snapshots * [int64]$gridPhysical * [int64]$gridPhysical * 8
foreach ($scheme in $schemes) {
    $dataPath = Join-Path $framesDir "rho_${scheme}_HLLC_grid_$($gridPhysical + 2 * 3)_frames_${snapshots}.bin"
    $metaPath = Join-Path $framesDir "rho_${scheme}_HLLC_grid_$($gridPhysical + 2 * 3)_frames_${snapshots}.meta"
    if (-not (Test-Path $dataPath)) {
        throw "Missing snapshot data file: $dataPath"
    }
    if (-not (Test-Path $metaPath)) {
        throw "Missing snapshot metadata file: $metaPath"
    }
    $actualBytes = (Get-Item $dataPath).Length
    if ($actualBytes -ne $expectedBytes) {
        throw "Unexpected size for ${dataPath}: $actualBytes bytes, expected $expectedBytes"
    }
}

Write-Host "All scheme processes completed."

if (-not $SkipAnimation) {
    $localFfmpeg = Join-Path $caseDir "tools\ffmpeg\bin\ffmpeg.exe"
    $ffmpegCommand = Get-Command ffmpeg -ErrorAction SilentlyContinue
    $ffmpegPath = $null
    if (Test-Path $localFfmpeg) {
        $ffmpegPath = $localFfmpeg
    } elseif ($ffmpegCommand) {
        $ffmpegPath = $ffmpegCommand.Source
    } else {
        try {
            $pythonFfmpeg = & python -c "import imageio_ffmpeg; print(imageio_ffmpeg.get_ffmpeg_exe())" 2>$null
            if ($LASTEXITCODE -eq 0 -and $pythonFfmpeg -and (Test-Path $pythonFfmpeg.Trim())) {
                $ffmpegPath = $pythonFfmpeg.Trim()
            }
        } catch {
            $pythonFfmpeg = $null
        }
    }
    if (-not $ffmpegPath) {
        throw "ffmpeg was not found. Install ffmpeg, install Python package imageio-ffmpeg, or place ffmpeg.exe at '$localFfmpeg', then rerun this script."
    }

    Write-Host "Rendering MP4 animation..."
    & python (Join-Path $caseDir "ri-animate.py") `
        --frames-dir $framesDir `
        --grid-total ($gridPhysical + 2 * 3) `
        --nsnap $snapshots `
        --tfinal $TFinal `
        --ffmpeg-path $ffmpegPath `
        --output $outputPath
    if ($LASTEXITCODE -ne 0) {
        throw "Animation rendering failed with exit code $LASTEXITCODE"
    }
    Write-Host "Animation written to $outputPath"
}
