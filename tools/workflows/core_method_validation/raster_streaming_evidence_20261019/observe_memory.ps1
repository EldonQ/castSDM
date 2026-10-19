$ErrorActionPreference = 'Stop'
$env:R_LIBS_USER = 'E:/Package/.rlib_cast'
$env:LANG = 'C'
$env:LC_ALL = 'C'
$taskRoot = 'E:\Package\_toolaudit\castSDM_raster_20261019'
$taskScript = Join-Path $taskRoot 'process_memory_run.R'
& 'D:\R\R-4.4.3\bin\Rscript.exe' $taskScript prepare *> (Join-Path $taskRoot 'memory_prepare.log')
if ($LASTEXITCODE -ne 0) { throw 'Input preparation failed' }
$taskObservations = [System.Collections.Generic.List[object]]::new()
foreach ($taskMode in @('old', 'new')) {
    $taskPidPath = Join-Path $taskRoot "memory_$taskMode.pid"
    if (Test-Path -LiteralPath $taskPidPath) { Remove-Item -LiteralPath $taskPidPath }
    $taskProcess = Start-Process -FilePath 'D:\R\R-4.4.3\bin\Rscript.exe' -ArgumentList @($taskScript, $taskMode) -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskRoot "memory_$taskMode.log") -RedirectStandardError (Join-Path $taskRoot "memory_${taskMode}_error.log")
    $taskTimer = [System.Diagnostics.Stopwatch]::StartNew()
    while (-not $taskProcess.HasExited) {
        if (Test-Path -LiteralPath $taskPidPath) {
            $taskRId = [int](Get-Content -LiteralPath $taskPidPath -TotalCount 1)
            $taskRProcess = Get-Process -Id $taskRId -ErrorAction SilentlyContinue
            if ($null -ne $taskRProcess) {
                $taskObservations.Add([pscustomobject]@{ mode=$taskMode; elapsed_seconds=$taskTimer.Elapsed.TotalSeconds; working_set_bytes=$taskRProcess.WorkingSet64; private_bytes=$taskRProcess.PrivateMemorySize64 })
            }
        }
        Start-Sleep -Milliseconds 100
    }
    $taskProcess.WaitForExit()
    if ($taskProcess.ExitCode -ne 0) { throw "$taskMode memory run failed" }
}
$taskObservations | Export-Csv -LiteralPath (Join-Path $taskRoot 'process_memory_samples.csv') -NoTypeInformation
