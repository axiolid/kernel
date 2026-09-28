import { existsSync, readFileSync, readdirSync } from "node:fs";
import { join, resolve } from "node:path";

import { loadGlossary } from "../.vitepress/glossary-plugin.mjs";

const obsoleteRoute = ["m", "cs-axi", "oval-contract-map"].join("");

function fail(message) {
  throw new Error(`DOCS_UI_CHECK: ${message}`);
}

function requireIncludes(content, expected, subject) {
  if (!content.includes(expected)) fail(`${subject} is missing ${JSON.stringify(expected)}`);
}

function markdownFiles(root) {
  const files = [];
  const visit = (path) => {
    for (const entry of readdirSync(path, { withFileTypes: true })) {
      if ([".git", ".hermes", "node_modules", "target", ".vitepress"].includes(entry.name)) continue;
      const child = join(path, entry.name);
      if (entry.isDirectory()) visit(child);
      else if (entry.isFile() && child.endsWith(".md")) files.push(child);
    }
  };
  visit(root);
  return files;
}

function sourceCheck(root) {
  const docs = join(root, "docs");
  const config = readFileSync(join(docs, ".vitepress/config.ts"), "utf8");
  const css = readFileSync(join(docs, ".vitepress/theme/custom.css"), "utf8");
  const contract = join(docs, "architecture/openbim-geometry-boundary.md");
  const oldContract = join(docs, `architecture/${obsoleteRoute}.md`);
  const entries = loadGlossary(join(docs, "glossary.md"));

  if (entries.length < 8) fail(`glossary is unexpectedly small (${entries.length} entries)`);
  if (!existsSync(contract)) fail("neutral openbim.geometry boundary page is absent");
  const contractSource = readFileSync(contract, "utf8");
  requireIncludes(contractSource, "`openbim.geometry:tessellation.surface`", "neutral capability ID example");
  if (contractSource.includes("`geometry.tessellation.surface`")) fail("unqualified capability ID example remains");
  if (existsSync(oldContract)) fail("obsolete application-specific contract route still exists");

  requireIncludes(config, "createGlossaryPlugin(glossaryEntries)(md)", "VitePress config");
  requireIncludes(config, 'rel: "icon"', "VitePress config");
  requireIncludes(config, 'href: `${docsBase}mark.svg`', "VitePress config");
  requireIncludes(config, 'link: "/glossary"', "VitePress sidebar");
  requireIncludes(config, 'link: "/architecture/openbim-geometry-boundary"', "VitePress sidebar");
  requireIncludes(css, ".glossary-term:hover::after", "glossary tooltip CSS");
  requireIncludes(css, ".glossary-term:focus-visible::after", "glossary focus CSS");
  requireIncludes(css, "@media (hover: none)", "touch tooltip CSS");
  requireIncludes(css, ".glossary-term::after {\n    display: none;", "touch tooltip overflow guard");

  const forbidden = new RegExp(["M", "CS", "|Axi", "oval"].join(""), "i");
  const scanned = markdownFiles(docs);
  for (const path of scanned) {
    const match = readFileSync(path, "utf8").match(forbidden);
    if (match) fail(`application-specific boundary name remains in ${path}: ${match[0]}`);
  }

  console.log(`DOCS_UI_SOURCE_CHECK=PASS glossary=${entries.length} markdown=${scanned.length}`);
}

function builtPage(dist, route) {
  const candidates = route
    ? [join(dist, `${route}.html`), join(dist, route, "index.html")]
    : [join(dist, "index.html")];
  const path = candidates.find(existsSync);
  if (!path) fail(`built page missing for /${route}`);
  return readFileSync(path, "utf8");
}

function builtCheck(dist) {
  if (!existsSync(join(dist, "mark.svg"))) fail("favicon asset was not copied");
  const home = builtPage(dist, "");
  const architecture = builtPage(dist, "architecture");
  const glossary = builtPage(dist, "glossary");
  const contract = builtPage(dist, "architecture/openbim-geometry-boundary");

  requireIncludes(home, 'rel="icon"', "built home page");
  requireIncludes(home, 'href="/kernel/mark.svg"', "built favicon metadata");
  requireIncludes(architecture, 'class="glossary-term"', "built architecture page");
  requireIncludes(architecture, 'href="/kernel/glossary#dag"', "built DAG link");
  requireIncludes(architecture, "data-definition=", "built glossary tooltip");
  requireIncludes(architecture, "aria-label=", "built glossary accessible name");
  requireIncludes(glossary, "Directed acyclic graph", "built glossary page");
  requireIncludes(contract, "openbim.geometry:tessellation.surface", "built neutral capability ID");

  const obsoleteHtml = [
    join(dist, `architecture/${obsoleteRoute}.html`),
    join(dist, `architecture/${obsoleteRoute}/index.html`),
  ];
  if (obsoleteHtml.some(existsSync)) fail("obsolete contract route was built");

  const dagLinks = architecture.match(/href="\/kernel\/glossary#dag"/g) ?? [];
  if (dagLinks.length !== 1) fail(`architecture page has ${dagLinks.length} DAG glossary links; expected first use only`);

  console.log("DOCS_UI_BUILD_CHECK=PASS favicon=1 glossary=1 first_use_dag=1 neutral_boundary=1");
}

const args = process.argv.slice(2);
if (args[0] === "--dist") builtCheck(resolve(args[1] ?? ".vitepress/dist"));
else sourceCheck(resolve(args[0] ?? ".."));
