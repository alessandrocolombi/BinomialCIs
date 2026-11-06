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
$documents = @(
    @{ Source = 'README.md'; Output = 'README.html'; Title = 'LowerBound - Implementation'; From = 'markdown+tex_math_single_backslash'; Lang = 'en' },
    @{ Source = 'Notes/README.md'; Output = 'Notes/README.html'; Title = 'LowerBound - Theory'; From = 'markdown+tex_math_single_backslash'; Lang = 'en' },
    @{ Source = 'Notes/theory.tex'; Output = 'Notes/theory.html'; Title = 'LowerBound - Derivazioni formali'; From = 'latex'; Lang = 'it' }
)
foreach ($document in $documents) {
    & $pandocPath (Join-Path $PSScriptRoot $document.Source) `
        "--from=$($document.From)" '--to=html5' '--standalone' `
        '--mathml' '--embed-resources' '--toc' '--toc-depth=2' `
        '--css' (Join-Path $PSScriptRoot 'readme.css') `
        "--metadata=title:$($document.Title)" "--metadata=lang:$($document.Lang)" `
        '--output' (Join-Path $PSScriptRoot $document.Output)
    if ($LASTEXITCODE -ne 0) { throw "Pandoc failed for $($document.Source) with exit code $LASTEXITCODE" }
    Write-Output (Join-Path $PSScriptRoot $document.Output)
}
