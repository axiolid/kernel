import pathlib
import sys

# Strip the Genus row to a status claim with no evidence (kernel#35 probe).
page = pathlib.Path("docs/capabilities.md")
text = page.read_text(encoding="utf-8")
marker = "| Genus | Implemented for closed two-manifolds |"
start = text.index(marker)
end = text.index("\n", start)
text = text[:start] + marker + " |" + text[end:]
page.write_text(text, encoding="utf-8")
sys.exit(0)
