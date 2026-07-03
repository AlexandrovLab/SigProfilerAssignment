import re


_OSF_TOC_RE = re.compile(r"^\s*@\[toc\]\([^)]+\)\s*$", re.IGNORECASE)


def on_page_markdown(markdown: str, /, *, page, config, files):
    # OSF wiki pages contain a non-standard marker like `@[toc](Sections)` which MkDocs
    # otherwise renders as a broken relative link. We keep the source files unchanged
    # (to match OSF) and strip the marker at build time.
    lines = markdown.splitlines()
    filtered = [line for line in lines if not _OSF_TOC_RE.match(line)]
    return "\n".join(filtered) + ("\n" if markdown.endswith("\n") else "")

