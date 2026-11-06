from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "proofreading_report.md"
OUT_DIR = ROOT / "output" / "pdf"
TMP_DIR = ROOT / "tmp" / "pdfs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
TMP_DIR.mkdir(parents=True, exist_ok=True)


def protect_math_and_code(text: str):
    tokens = []

    def save(match):
        tokens.append(match.group(0))
        return f"@@TOKEN{len(tokens)-1}@@"

    text = re.sub(r"`[^`]*`", save, text)
    text = re.sub(r"\\\(.*?\\\)", save, text)
    return text, tokens


def escape_plain(text: str) -> str:
    text = (text.replace("—", "-").replace("–", "-")
                .replace("’", "'").replace("“", '"').replace("”", '"')
                .replace("…", "...").replace("→", r"$\to$")
                .replace("≥", r"$\ge$").replace("≤", r"$\le$"))
    replacements = {
        "&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_",
        "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def escape_code(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}\allowbreak{}",
        "_": r"\_\allowbreak{}",
        "{": r"\{\allowbreak{}",
        "}": r"\}\allowbreak{}",
        "%": r"\%\allowbreak{}",
        "#": r"\#\allowbreak{}",
        "&": r"\&\allowbreak{}",
        "^": r"\textasciicircum{}\allowbreak{}",
        "~": r"\textasciitilde{}\allowbreak{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def inline(text: str) -> str:
    text, tokens = protect_math_and_code(text)
    text = escape_plain(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", text)
    text = re.sub(r"\*(.+?)\*", r"\\emph{\1}", text)
    for i, token in enumerate(tokens):
        marker = f"@@TOKEN{i}@@"
        if token.startswith("`"):
            content = token[1:-1]
            rendered = r"\texttt{" + escape_code(content) + "}"
        else:
            rendered = token
        text = text.replace(marker, rendered)
    return text


lines = SOURCE.read_text(encoding="utf-8").splitlines()
body = []
in_display = False

for raw in lines:
    stripped = raw.strip()
    if stripped == r"\[":
        in_display = True
        body.append(r"\[")
        continue
    if stripped == r"\]":
        in_display = False
        body.append(r"\]")
        continue
    if in_display:
        body.append(raw)
        continue
    if not stripped:
        body.append("")
        continue
    if stripped.startswith("# "):
        body.append(r"\begin{center}\LARGE\bfseries " + inline(stripped[2:]) + r"\end{center}")
    elif stripped.startswith("## "):
        body.append(r"\section*{" + inline(stripped[3:]) + "}")
    elif stripped.startswith("### "):
        body.append(r"\subsection*{" + inline(stripped[4:]) + "}")
    elif re.match(r"^\d+\.\s", stripped):
        num, rest = stripped.split(".", 1)
        body.append(r"\par\medskip\noindent\textbf{" + num + ".} " + inline(rest.strip()))
    elif stripped.startswith("- "):
        body.append(r"\par\smallskip\noindent\hangindent=1.2em\hangafter=1 $\bullet$ " + inline(stripped[2:]))
    else:
        body.append(inline(stripped))

preamble = r"""\documentclass[10pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,bm}
\usepackage[a4paper,margin=2.15cm,headheight=14pt]{geometry}
\usepackage{xcolor}
\usepackage{microtype}
\usepackage{enumitem}
\usepackage{fancyhdr}
\usepackage{hyperref}
\definecolor{navy}{HTML}{17365D}
\definecolor{accent}{HTML}{B4432E}
\definecolor{muted}{HTML}{5B6673}
\hypersetup{colorlinks=true,linkcolor=navy,urlcolor=navy}
\setlength{\parindent}{0pt}
\setlength{\parskip}{4pt plus 1pt}
\setlength{\emergencystretch}{2em}
\pagestyle{fancy}
\fancyhf{}
\lhead{JASA draft - proofreading report}
\rhead{Technical and language review}
\cfoot{\thepage}
\renewcommand{\headrulewidth}{0.3pt}
\usepackage{titlesec}
\titleformat{\section}{\Large\bfseries\color{navy}}{}{0pt}{}
\titleformat{\subsection}{\large\bfseries\color{accent}}{}{0pt}{}
\begin{document}
"""

ending = r"""
\end{document}
"""

(TMP_DIR / "proofreading_report.tex").write_text(
    preamble + "\n".join(body) + ending,
    encoding="utf-8",
)
