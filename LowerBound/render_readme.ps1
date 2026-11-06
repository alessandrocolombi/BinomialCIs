# Rebuild the HTML documentation only; never starts a simulation.
$ErrorActionPreference = 'Stop'
$pandocCommand = Get-Command pandoc -ErrorAction SilentlyContinue
$pandocPath = if ($pandocCommand) { $pandocCommand.Source } else {
    $candidates = @(
        'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe',
        'C:/Program Files/RStudio/bin/pandoc/pandoc.exe'
    )
    $candidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
}
if (-not $pandocPath) { throw 'Pandoc not found. Add pandoc to PATH or install RStudio.' }
& $pandocPath (Join-Path $PSScriptRoot 'README.md') `
    '--from=markdown+tex_math_single_backslash' '--to=html5' '--standalone' `
    '--mathml' '--embed-resources' '--toc' '--toc-depth=2' `
    '--css' (Join-Path $PSScriptRoot 'readme.css') `
    '--metadata=title:LowerBound - Documentation' '--metadata=lang:en' `
    '--output' (Join-Path $PSScriptRoot 'README.html')
if ($LASTEXITCODE -ne 0) { throw "Pandoc failed with exit code $LASTEXITCODE" }
Write-Output (Join-Path $PSScriptRoot 'README.html')
