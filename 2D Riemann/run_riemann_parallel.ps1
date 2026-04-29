param(
    [ValidateSet("preview", "production")]
    [string]$Mode = "preview",
    [ValidateSet("REFL", "OUTF", "COMPARE")]
    [string]$BoundaryMode = "REFL",
    [int]$Threads = 4,
    [double]$TFinal = 0.5,
    [switch]$SkipAnimation
)

$ErrorActionPreference = "Stop"

# Paths stay relative to the case folder so the script can be run from elsewhere.
$caseDir = $PSScriptRoot
$buildDir = Join-Path $caseDir "build"
$exe = Join-Path $buildDir "ri-six-omp.exe"
$relativeExe = ".\build\ri-six-omp.exe"

if ($Mode -eq "preview") {
    # Preview mode keeps the grid and frame count low for quick iteration.
    $gridPhysical = 100
    $snapshots = 10
    $baseFramesName = "frames_preview"
    $baseOutputName = "riemann_density_3schemes_grid_106_preview"
} else {
    # Production mode uses the report-resolution grid and longer animation stack.
    $gridPhysical = 800
    $snapshots = 100
    $baseFramesName = "frames"
    $baseOutputName = "riemann_density_3schemes_grid_806"
}

$timeTag = ("t{0:0.###}" -f $TFinal).Replace(".", "p")
$bcLabels = @{
    OUTF = "outflow"
    REFL = "reflection"
}
if ($BoundaryMode -eq "COMPARE") {
    # Comparison mode runs both boundary treatments and renders them together.
    $boundaryModes = @("OUTF", "REFL")
    $outputName = "${baseOutputName}_outflow_vs_reflection_${timeTag}.mp4"
} else {
    $boundaryModes = @($BoundaryMode)
    $outputName = "${baseOutputName}_$($bcLabels[$BoundaryMode])_${timeTag}.mp4"
}

$outputPath = Join-Path $caseDir $outputName

New-Item -ItemType Directory -Force $buildDir | Out-Null

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

function Get-FramesName {
    param([string]$BcMode)
    # Include boundary mode and final time so repeated runs do not overwrite each other.
    return "${baseFramesName}_$($bcLabels[$BcMode])_${timeTag}"
}

function Invoke-SchemeBatch {
    param([string]$BcMode)

    # Each scheme writes its own logs and density snapshots into the same run folder.
    $framesName = Get-FramesName $BcMode
    $framesDir = Join-Path $caseDir $framesName
    New-Item -ItemType Directory -Force $framesDir | Out-Null

    $processes = @()
    Write-Host "Launching $($schemes.Count) $BcMode scheme processes with $Threads OpenMP threads each..."
    foreach ($scheme in $schemes) {
        $stdoutLog = Join-Path $framesDir "$scheme.stdout.log"
        $stderrLog = Join-Path $framesDir "$scheme.stderr.log"
        $args = @($scheme, "$gridPhysical", "$snapshots", $framesName, "$TFinal", $BcMode)
        $processes += Start-Process -FilePath $exe `
            -ArgumentList $args `
            -WorkingDirectory $caseDir `
            -RedirectStandardOutput $stdoutLog `
            -RedirectStandardError $stderrLog `
            -WindowStyle Hidden `
            -PassThru
    }

    foreach ($process in $processes) {
        # Wait after launching all schemes so they can run in parallel.
        $process.WaitForExit()
        $process.Refresh()
    }

    $failed = $false
    for ($idx = 0; $idx -lt $processes.Count; $idx++) {
        $exitCode = $processes[$idx].ExitCode
        if ($null -ne $exitCode -and $exitCode -ne 0) {
            Write-Host "$BcMode $($schemes[$idx]) failed with exit code $($processes[$idx].ExitCode)"
            $failed = $true
        }
    }
    if ($failed) {
        throw "At least one $BcMode scheme process failed. Check logs in $framesDir"
    }

    $expectedBytes = [int64]$snapshots * [int64]$gridPhysical * [int64]$gridPhysical * 8
    foreach ($scheme in $schemes) {
        # Check both file presence and size before spending time on rendering.
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

    Write-Host "$BcMode scheme processes completed."
    return $framesDir
}

$completedFrameDirs = @{}
foreach ($bc in $boundaryModes) {
    $completedFrameDirs[$bc] = Invoke-SchemeBatch $bc
}

Write-Host "All requested scheme processes completed."

if (-not $SkipAnimation) {
    # Prefer a local ffmpeg, then PATH, then the imageio-ffmpeg Python package.
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
    if ($BoundaryMode -eq "COMPARE") {
        # Pass both frame directories so the renderer can build a two-row comparison.
        & python (Join-Path $caseDir "ri-animate.py") `
            --frames-dir-outflow $completedFrameDirs["OUTF"] `
            --frames-dir-reflection $completedFrameDirs["REFL"] `
            --grid-total ($gridPhysical + 2 * 3) `
            --nsnap $snapshots `
            --tfinal $TFinal `
            --ffmpeg-path $ffmpegPath `
            --output $outputPath
    } else {
        & python (Join-Path $caseDir "ri-animate.py") `
            --frames-dir $completedFrameDirs[$BoundaryMode] `
            --grid-total ($gridPhysical + 2 * 3) `
            --nsnap $snapshots `
            --tfinal $TFinal `
            --ffmpeg-path $ffmpegPath `
            --output $outputPath
    }
    if ($LASTEXITCODE -ne 0) {
        throw "Animation rendering failed with exit code $LASTEXITCODE"
    }
    Write-Host "Animation written to $outputPath"
}
