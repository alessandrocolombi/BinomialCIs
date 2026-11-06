# User-facing documentation

Whenever creating or updating a README or other user-facing Markdown
documentation in this repository, also create or update a matching HTML file.
The user reads HTML rather than Markdown and needs rendered LaTeX formulas.
Prefer standalone HTML with embedded CSS and MathML so it works offline.
Keep the Markdown source and HTML synchronized and link the HTML in replies.
This requirement does not apply to internal instruction files such as AGENTS.md.

For LowerBound/README.md, use LowerBound/render_readme.ps1 (Pandoc) to rebuild
LowerBound/README.html. Do not run simulation scripts merely to rebuild docs.

# Simulation completion logs

For simulation or batch scripts created or modified, print a brief results
recap on successful completion (settings, elapsed time, key results and output
locations), followed by a final line exactly FINITO. Print this marker only
after all computation, output saves, and worker cleanup have succeeded.
Respect the user's request not to execute the VM simulation locally.
