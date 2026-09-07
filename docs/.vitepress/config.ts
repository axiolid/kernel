import { fileURLToPath } from "node:url";

import { defineConfig } from "vitepress";

import { diagramPlugin } from "./diagram-plugin.ts";
import { createGlossaryPlugin, loadGlossary } from "./glossary-plugin.mjs";

const docsBase = "/kernel/";
const glossaryEntries = loadGlossary(fileURLToPath(new URL("../glossary.md", import.meta.url)));

const adrs = [
  ["0001-axiolid-ifc-split-and-kernel-contract", 1, "IFC split and kernel contract"],
  ["0002-hardware-abstraction-and-backend-selection", 2, "Hardware and backends"],
  ["0003-pure-rust-mesh-boolean", 3, "Pure-Rust mesh booleans"],
  ["0004-package-layout-and-backend-features", 4, "Package layout"],
  ["0009-layered-geometry-dag", 9, "Layered geometry DAG"],
  ["0011-native-accelerator-backends-out-of-tree", 11, "Native accelerators"],
  ["0012-scalar-reference-ownership", 12, "Scalar reference"],
  ["0013-deferred-performance-techniques", 13, "Deferred performance work"],
  ["0014-adopt-boolmesh-mesh-boolean", 14, "Boolmesh provider"],
  ["0015-adopt-earcut-polygon-triangulation", 15, "Earcut triangulation"],
  ["0016-predicate-ownership-and-adopted-implementations", 16, "Predicate ownership"],
  ["0017-solid-boolean-contract-before-implementation", 17, "Solid boolean contract"],
  ["0018-curve-evaluation-in-the-scalar-reference", 18, "Curve evaluation"],
  ["0019-validate-and-refine-nurbs-on-the-scalar-read-path", 19, "NURBS read path"],
  ["0020-exact-brep-kernel-model", 20, "Exact B-rep kernel model"],
  ["0021-capability-seams-live-in-the-kernel", 21, "Capability seams (superseded)"],
  ["0022-general-nurbs-kernel-capability", 22, "General NURBS kernel capability"],
  ["0023-solid-generation-is-an-l2-crate", 23, "Solid generation at L2"],
  ["0024-exact-brep-result-contracts", 24, "Exact B-rep result contracts"],
  ["0025-certified-nurbs-subdivision-oracle", 25, "Certified NURBS subdivision"],
  ["0026-certified-planar-nurbs-root-isolation", 26, "Planar NURBS roots"],
  ["0027-certified-nurbs-curve-surface-root-isolation", 27, "Curve/surface roots"],
  ["0028-certified-affine-surface-surface-tracing", 28, "Affine surface tracing"],
  ["0029-certified-trace-topology-integration", 29, "Trace topology integration"],
  ["0030-globally-certified-surface-projection", 30, "Certified surface projection"],
  ["0031-verified-periodic-curve-views", 31, "Verified periodic curve views"],
  ["0032-explicit-periodic-bspline-surfaces", 32, "Explicit periodic surfaces"],
  ["0033-mesh-plane-section-contract", 33, "Mesh plane sections"],
  ["0034-authored-open-profile-contract", 34, "Authored open profiles"],
  ["0035-nested-ownership-and-capability-contracts", 35, "Nested ownership layout"],
  ["0036-use-case-specific-compilation-closures", 36, "Use-case compilation closures"],
  ["0037-mapped-3d-verification-oracle", 37, "Mapped-3D verification oracle"],
  ["0038-constructed-intersection-curves", 38, "Constructed intersection curves"],
  ["0039-downstream-integration-contract", 39, "Downstream integration contract"],
  ["0040-versioned-c-abi-boundary", 40, "Versioned C ABI"],
  ["0041-cmake-and-native-release-bundles", 41, "CMake and native bundles"],
  ["0042-black-box-downstream-compatibility-gate", 42, "Black-box compatibility gate"],
  ["0043-release-version-and-changelog-automation", 43, "Release version/changelog automation"],
  ["0044-pointcloud-representation-and-reconstruction", 44, "Pointcloud representation and reconstruction"],
] as const;

const adrIcon = `<svg class="adr-icon" aria-hidden="true" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/><path d="M14 2v6h6"/><path d="M8 13h8M8 17h6"/></svg>`;

function adrSidebarItem([file, number, title]: (typeof adrs)[number]) {
  return {
    text: `<span class="adr-sidebar-label">${adrIcon}<span class="adr-sidebar-title">${title}</span><span class="adr-sidebar-number">${number}</span></span>`,
    link: `/adr/${file}`,
  };
}

export default defineConfig({
  title: "Axiolid",
  description: "A pure-Rust, format-agnostic geometry kernel.",
  base: docsBase,
  srcExclude: ["adr/_template.md"],
  markdown: {
    html: false,
    math: true,
    config(md) {
      diagramPlugin(md);
      createGlossaryPlugin(glossaryEntries)(md);
    },
  },
  cleanUrls: true,
  lastUpdated: true,
  head: [
    ["link", { rel: "icon", type: "image/svg+xml", href: `${docsBase}mark.svg` }],
    ["meta", { name: "theme-color", content: "#111827" }],
    ["meta", { property: "og:title", content: "Axiolid" }],
    ["meta", { property: "og:description", content: "A pure-Rust, format-agnostic geometry kernel." }],
  ],
  themeConfig: {
    logo: "/mark.svg",
    nav: [
      { text: "Guide", link: "/guide/getting-started" },
      { text: "Architecture", link: "/architecture" },
      { text: "Capabilities", link: "/capabilities" },
      { text: "Support map", link: "/support-map" },
      { text: "API", link: "https://docs.rs/axiolid" },
      { text: "GitHub", link: "https://github.com/axiolid/kernel" },
    ],
    sidebar: [
      {
        text: "Start here",
        items: [
          { text: "Overview", link: "/" },
          { text: "Getting started", link: "/guide/getting-started" },
          { text: "Downstream integration", link: "/guide/downstream-integration" },
          { text: "Geometry concepts", link: "/guide/geometry-concepts" },
          { text: "Glossary", link: "/glossary" },
          { text: "Capabilities", link: "/capabilities" },
          { text: "Support map", link: "/support-map" },
          { text: "Architecture", link: "/architecture" },
          { text: "Crate map", link: "/architecture/crate-map" },
          { text: "Closure profiles", link: "/architecture/closure-profiles" },
          { text: "2D-only consumers", link: "/architecture/2d-only-consumers" },
          { text: "C ABI v0.4", link: "/architecture/c-abi-v0.4" },
          { text: "Native distribution", link: "/architecture/native-distribution" },
          { text: "Dependency graph", link: "/architecture/dependency-graph" },
          { text: "Threading", link: "/architecture/threading" },
          { text: "openbim.geometry boundary", link: "/architecture/openbim-geometry-boundary" },
          { text: "Public crate reference", link: "/reference/crates" },
        ],
      },
      {
        text: "About",
        items: [
          { text: "Roadmap", link: "/ROADMAP" },
          { text: "Changelog", link: "/CHANGELOG" },
          { text: "Research", link: "/research/geometry-kernel-capability-comparison" },
          { text: "Licensing", link: "/guide/licensing" },
          { text: "Contributing", link: "/guide/contributing" },
          { text: "Where things go", link: "/contributing/where-things-go" },
          { text: "Crate migration", link: "/contributing/crate-migration" },
        ],
      },
      {
        text: "Architecture decisions",
        items: adrs.map(adrSidebarItem),
      },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/axiolid/kernel" }],
    footer: {
      message: "Released under the Mozilla Public License 2.0.",
      copyright: "Copyright © 2026 Axiolid contributors",
    },
    search: { provider: "local" },
  },
});
